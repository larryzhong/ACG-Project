#pragma once

#include <array>
#include <algorithm>
#include <cfloat>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

#include "core/ray.h"
#include "core/rng.h"
#include "core/vec2.h"
#include "core/vec3.h"
#include "core/color.h"
#include "math/aabb.h"
#include "math/transform.h"
#include "scene/hittable.h"
#include "scene/material.h"

class MeshData {
public:
    MeshData(std::vector<Vec3> positions,
             std::vector<Vec3> normals,
             std::vector<Vec2> uvs,
             std::vector<std::uint32_t> indices)
        : positions_(std::move(positions)),
          normals_(std::move(normals)),
          uvs_(std::move(uvs)),
          indices_(std::move(indices)) {
        build_internal();
    }

    std::size_t triangle_count() const {
        return indices_.size() / 3;
    }

    const AABB& bounding_box() const {
        return bounds_;
    }

    bool hit(const Ray& r, float t_min, float t_max, HitRecord& rec) const {
        if (nodes_.empty() || prim_indices_.empty()) {
            return false;
        }

        bool hit_anything = false;
        float closest_so_far = t_max;

        std::uint32_t best_tri = 0;
        float best_b1 = 0.0f;
        float best_b2 = 0.0f;

        std::vector<int> stack;
        stack.reserve(64);
        stack.push_back(root_);

        while (!stack.empty()) {
            const int node_index = stack.back();
            stack.pop_back();
            const Node& node = nodes_[node_index];

            if (!node.box.hit(r, t_min, closest_so_far)) {
                continue;
            }

            if (node.is_leaf) {
                for (std::uint32_t i = 0; i < node.count; ++i) {
                    const std::uint32_t tri = prim_indices_[node.start + i];
                    float t = 0.0f;
                    float b1 = 0.0f;
                    float b2 = 0.0f;
                    if (intersect_triangle(tri, r, t_min, closest_so_far, t, b1, b2)) {
                        hit_anything = true;
                        closest_so_far = t;
                        best_tri = tri;
                        best_b1 = b1;
                        best_b2 = b2;
                    }
                }
            } else {
                stack.push_back(node.left);
                stack.push_back(node.right);
            }
        }

        if (!hit_anything) {
            return false;
        }

        const std::size_t i0 = indices_[3 * static_cast<std::size_t>(best_tri) + 0];
        const std::size_t i1 = indices_[3 * static_cast<std::size_t>(best_tri) + 1];
        const std::size_t i2 = indices_[3 * static_cast<std::size_t>(best_tri) + 2];

        const Vec3& p0 = positions_[i0];
        const Vec3& p1 = positions_[i1];
        const Vec3& p2 = positions_[i2];

        const float w0 = 1.0f - best_b1 - best_b2;
        const float w1 = best_b1;
        const float w2 = best_b2;

        rec.t = closest_so_far;
        rec.point = r.at(rec.t);

        Vec3 geom_normal = cross(p1 - p0, p2 - p0);
        if (geom_normal.length_squared() > 1e-20f) {
            geom_normal = normalize(geom_normal);
        } else {
            geom_normal = Vec3(0.0f, 1.0f, 0.0f);
        }

        Vec2 uv = uvs_[i0] * w0 + uvs_[i1] * w1 + uvs_[i2] * w2;
        rec.u = uv.x;
        rec.v = uv.y;

        Vec3 shading_normal = normals_[i0] * w0 + normals_[i1] * w1 + normals_[i2] * w2;
        if (shading_normal.length_squared() > 1e-20f) {
            shading_normal = normalize(shading_normal);
        } else {
            shading_normal = geom_normal;
        }

        if (dot(shading_normal, geom_normal) < 0.0f) {
            shading_normal = -shading_normal;
        }

        Vec3 tangent = tangents_[i0] * w0 + tangents_[i1] * w1 + tangents_[i2] * w2;
        if (tangent.length_squared() > 1e-20f) {
            tangent = normalize(tangent);
        } else {
            tangent = normalize(p1 - p0);
        }

        rec.front_face = dot(r.direction, geom_normal) < 0.0f;
        rec.normal = rec.front_face ? shading_normal : -shading_normal;

        Vec3 ortho_tangent = tangent - dot(tangent, rec.normal) * rec.normal;
        if (ortho_tangent.length_squared() < 1e-16f) {
            ortho_tangent = orthonormal_tangent(rec.normal);
        } else {
            ortho_tangent = normalize(ortho_tangent);
        }
        rec.tangent = rec.front_face ? ortho_tangent : -ortho_tangent;

        const Vec2& uv0 = uvs_[i0];
        const Vec2& uv1 = uvs_[i1];
        const Vec2& uv2 = uvs_[i2];
        const Vec3 e1 = p1 - p0;
        const Vec3 e2 = p2 - p0;
        const Vec2 duv1 = uv1 - uv0;
        const Vec2 duv2 = uv2 - uv0;
        const float det = duv1.x * duv2.y - duv1.y * duv2.x;
        if (std::fabs(det) > 1e-10f) {
            const float inv_det = 1.0f / det;
            rec.dpdu = (e1 * duv2.y - e2 * duv1.y) * inv_det;
            rec.dpdv = (e2 * duv1.x - e1 * duv2.x) * inv_det;
        } else {
            rec.dpdu = rec.tangent;
            rec.dpdv = cross(rec.normal, rec.dpdu);
        }

        rec.material = nullptr;
        rec.object = nullptr;

        return true;
    }

    bool triangle_vertices(std::uint32_t tri,
                           Vec3& out_p0,
                           Vec3& out_p1,
                           Vec3& out_p2,
                           Vec2& out_uv0,
                           Vec2& out_uv1,
                           Vec2& out_uv2) const {
        const std::size_t tri_count = triangle_count();
        if (tri_count == 0 || tri >= tri_count) {
            return false;
        }

        const std::size_t i0 = indices_[3 * static_cast<std::size_t>(tri) + 0];
        const std::size_t i1 = indices_[3 * static_cast<std::size_t>(tri) + 1];
        const std::size_t i2 = indices_[3 * static_cast<std::size_t>(tri) + 2];

        if (i0 >= positions_.size() || i1 >= positions_.size() || i2 >= positions_.size()) {
            return false;
        }
        if (uvs_.empty() || i0 >= uvs_.size() || i1 >= uvs_.size() || i2 >= uvs_.size()) {
            out_uv0 = Vec2(0.0f);
            out_uv1 = Vec2(0.0f);
            out_uv2 = Vec2(0.0f);
        } else {
            out_uv0 = uvs_[i0];
            out_uv1 = uvs_[i1];
            out_uv2 = uvs_[i2];
        }

        out_p0 = positions_[i0];
        out_p1 = positions_[i1];
        out_p2 = positions_[i2];
        return true;
    }

private:
    struct Node {
        AABB box;
        int left = -1;
        int right = -1;
        std::uint32_t start = 0;
        std::uint32_t count = 0;
        bool is_leaf = false;
    };

    static Vec3 orthonormal_tangent(const Vec3& n) {
        Vec3 a = (std::fabs(n.x) < 0.9f) ? Vec3(1.0f, 0.0f, 0.0f) : Vec3(0.0f, 1.0f, 0.0f);
        return normalize(cross(a, n));
    }

    void build_internal() {
        if (indices_.size() % 3 != 0) {
            indices_.resize(indices_.size() - (indices_.size() % 3));
        }

        if (positions_.empty() || indices_.empty()) {
            bounds_ = AABB(Vec3::zero(), Vec3::zero());
            return;
        }

        if (normals_.size() != positions_.size()) {
            compute_vertex_normals();
        }

        if (uvs_.size() != positions_.size()) {
            uvs_.assign(positions_.size(), Vec2::zero());
        }

        compute_vertex_tangents();
        build_triangle_bounds();
        build_bvh();
    }

    void compute_vertex_normals() {
        normals_.assign(positions_.size(), Vec3::zero());

        const std::size_t tri_count = triangle_count();
        for (std::size_t tri = 0; tri < tri_count; ++tri) {
            const std::size_t i0 = indices_[3 * tri + 0];
            const std::size_t i1 = indices_[3 * tri + 1];
            const std::size_t i2 = indices_[3 * tri + 2];

            const Vec3& p0 = positions_[i0];
            const Vec3& p1 = positions_[i1];
            const Vec3& p2 = positions_[i2];

            const Vec3 n = cross(p1 - p0, p2 - p0);
            normals_[i0] += n;
            normals_[i1] += n;
            normals_[i2] += n;
        }

        for (auto& n : normals_) {
            if (n.length_squared() > 1e-20f) {
                n = normalize(n);
            } else {
                n = Vec3(0.0f, 1.0f, 0.0f);
            }
        }
    }

    void compute_vertex_tangents() {
        tangents_.assign(positions_.size(), Vec3::zero());

        const std::size_t tri_count = triangle_count();
        for (std::size_t tri = 0; tri < tri_count; ++tri) {
            const std::size_t i0 = indices_[3 * tri + 0];
            const std::size_t i1 = indices_[3 * tri + 1];
            const std::size_t i2 = indices_[3 * tri + 2];

            const Vec3& p0 = positions_[i0];
            const Vec3& p1 = positions_[i1];
            const Vec3& p2 = positions_[i2];

            const Vec2& uv0 = uvs_[i0];
            const Vec2& uv1 = uvs_[i1];
            const Vec2& uv2 = uvs_[i2];

            const Vec3 dp1 = p1 - p0;
            const Vec3 dp2 = p2 - p0;

            const Vec2 duv1 = uv1 - uv0;
            const Vec2 duv2 = uv2 - uv0;

            const float denom = duv1.x * duv2.y - duv1.y * duv2.x;
            Vec3 tangent;

            if (std::fabs(denom) < 1e-12f) {
                tangent = dp1;
            } else {
                const float r = 1.0f / denom;
                tangent = (dp1 * duv2.y - dp2 * duv1.y) * r;
            }

            tangents_[i0] += tangent;
            tangents_[i1] += tangent;
            tangents_[i2] += tangent;
        }

        for (std::size_t i = 0; i < tangents_.size(); ++i) {
            const Vec3 n = normals_[i];
            Vec3 t = tangents_[i];

            t = t - dot(t, n) * n;
            if (t.length_squared() < 1e-16f) {
                t = orthonormal_tangent(n);
            } else {
                t = normalize(t);
            }

            tangents_[i] = t;
        }
    }

    void build_triangle_bounds() {
        const std::size_t tri_count = triangle_count();
        tri_boxes_.clear();
        tri_centroids_.clear();
        prim_indices_.clear();
        nodes_.clear();

        tri_boxes_.reserve(tri_count);
        tri_centroids_.reserve(tri_count);
        prim_indices_.reserve(tri_count);

        bool first_box = true;

        for (std::size_t tri = 0; tri < tri_count; ++tri) {
            const std::size_t i0 = indices_[3 * tri + 0];
            const std::size_t i1 = indices_[3 * tri + 1];
            const std::size_t i2 = indices_[3 * tri + 2];

            const Vec3& p0 = positions_[i0];
            const Vec3& p1 = positions_[i1];
            const Vec3& p2 = positions_[i2];

            const float min_x = std::min({p0.x, p1.x, p2.x});
            const float min_y = std::min({p0.y, p1.y, p2.y});
            const float min_z = std::min({p0.z, p1.z, p2.z});

            const float max_x = std::max({p0.x, p1.x, p2.x});
            const float max_y = std::max({p0.y, p1.y, p2.y});
            const float max_z = std::max({p0.z, p1.z, p2.z});

            const float epsilon = 1e-5f;
            AABB box(
                Vec3(min_x - epsilon, min_y - epsilon, min_z - epsilon),
                Vec3(max_x + epsilon, max_y + epsilon, max_z + epsilon));

            tri_boxes_.push_back(box);
            tri_centroids_.push_back((p0 + p1 + p2) / 3.0f);
            prim_indices_.push_back(static_cast<std::uint32_t>(tri));

            if (first_box) {
                bounds_ = box;
                first_box = false;
            } else {
                bounds_ = surrounding_box(bounds_, box);
            }
        }
    }

    static int choose_split_axis(const AABB& bounds) {
        const float dx = bounds.max.x - bounds.min.x;
        const float dy = bounds.max.y - bounds.min.y;
        const float dz = bounds.max.z - bounds.min.z;

        if (dx > dy && dx > dz) {
            return 0;
        } else if (dy > dz) {
            return 1;
        }
        return 2;
    }

    static float surface_area(const AABB& b) {
        const float dx = std::max(0.0f, b.max.x - b.min.x);
        const float dy = std::max(0.0f, b.max.y - b.min.y);
        const float dz = std::max(0.0f, b.max.z - b.min.z);
        return 2.0f * (dx * dy + dy * dz + dz * dx);
    }

    struct SAHSplit {
        int axis = -1;
        int bin_split = -1;
    };

    SAHSplit find_sah_split(std::uint32_t start, std::uint32_t end) const {
        constexpr int kBins = 16;
        const std::uint32_t n = end - start;
        if (n < 2) {
            return {};
        }

        AABB parent_bounds;
        bool first_box = true;
        for (std::uint32_t i = start; i < end; ++i) {
            const std::uint32_t tri = prim_indices_[i];
            const AABB& b = tri_boxes_[tri];
            if (first_box) {
                parent_bounds = b;
                first_box = false;
            } else {
                parent_bounds = surrounding_box(parent_bounds, b);
            }
        }

        const float parent_area = surface_area(parent_bounds);
        if (!(parent_area > 0.0f)) {
            return {};
        }

        SAHSplit best;
        float best_cost = FLT_MAX;

        for (int axis = 0; axis < 3; ++axis) {
            float cmin = FLT_MAX;
            float cmax = -FLT_MAX;

            for (std::uint32_t i = start; i < end; ++i) {
                const std::uint32_t tri = prim_indices_[i];
                const Vec3& c = tri_centroids_[tri];
                const float v = (axis == 0) ? c.x : (axis == 1) ? c.y : c.z;
                cmin = std::min(cmin, v);
                cmax = std::max(cmax, v);
            }

            const float extent = cmax - cmin;
            if (!(extent > 1e-6f)) {
                continue;
            }

            struct Bin {
                AABB bounds;
                int count = 0;
                bool init = false;
            };
            std::array<Bin, kBins> bins;

            auto bin_index = [&](float c) {
                const float t = (c - cmin) / extent;
                int idx = static_cast<int>(t * static_cast<float>(kBins));
                if (idx < 0) idx = 0;
                if (idx >= kBins) idx = kBins - 1;
                return idx;
            };

            for (std::uint32_t i = start; i < end; ++i) {
                const std::uint32_t tri = prim_indices_[i];
                const Vec3& c = tri_centroids_[tri];
                const float v = (axis == 0) ? c.x : (axis == 1) ? c.y : c.z;
                const int idx = bin_index(v);
                Bin& bin = bins[static_cast<std::size_t>(idx)];
                const AABB& b = tri_boxes_[tri];
                if (!bin.init) {
                    bin.bounds = b;
                    bin.init = true;
                } else {
                    bin.bounds = surrounding_box(bin.bounds, b);
                }
                bin.count += 1;
            }

            std::array<AABB, kBins> prefix_bounds;
            std::array<int, kBins> prefix_count;
            std::array<AABB, kBins> suffix_bounds;
            std::array<int, kBins> suffix_count;

            bool pref_init = false;
            int cnt = 0;
            AABB bnd;
            for (int i = 0; i < kBins; ++i) {
                const Bin& bin = bins[static_cast<std::size_t>(i)];
                if (bin.count > 0) {
                    if (!pref_init) {
                        bnd = bin.bounds;
                        pref_init = true;
                    } else {
                        bnd = surrounding_box(bnd, bin.bounds);
                    }
                    cnt += bin.count;
                }
                prefix_bounds[static_cast<std::size_t>(i)] = bnd;
                prefix_count[static_cast<std::size_t>(i)] = cnt;
            }

            bool suf_init = false;
            cnt = 0;
            for (int i = kBins - 1; i >= 0; --i) {
                const Bin& bin = bins[static_cast<std::size_t>(i)];
                if (bin.count > 0) {
                    if (!suf_init) {
                        bnd = bin.bounds;
                        suf_init = true;
                    } else {
                        bnd = surrounding_box(bnd, bin.bounds);
                    }
                    cnt += bin.count;
                }
                suffix_bounds[static_cast<std::size_t>(i)] = bnd;
                suffix_count[static_cast<std::size_t>(i)] = cnt;
            }

            for (int split = 0; split < kBins - 1; ++split) {
                const int left_count = prefix_count[static_cast<std::size_t>(split)];
                const int right_count = suffix_count[static_cast<std::size_t>(split + 1)];
                if (left_count == 0 || right_count == 0) {
                    continue;
                }

                const float left_area = surface_area(prefix_bounds[static_cast<std::size_t>(split)]);
                const float right_area = surface_area(suffix_bounds[static_cast<std::size_t>(split + 1)]);
                const float cost =
                    1.0f +
                    (left_area / parent_area) * static_cast<float>(left_count) +
                    (right_area / parent_area) * static_cast<float>(right_count);

                if (cost < best_cost) {
                    best_cost = cost;
                    best.axis = axis;
                    best.bin_split = split;
                }
            }
        }

        return best;
    }

    std::uint32_t partition_by_sah_split(std::uint32_t start, std::uint32_t end, const SAHSplit& split) {
        constexpr int kBins = 16;
        if (split.axis < 0) {
            return start;
        }

        float cmin = FLT_MAX;
        float cmax = -FLT_MAX;
        for (std::uint32_t i = start; i < end; ++i) {
            const std::uint32_t tri = prim_indices_[i];
            const Vec3& c = tri_centroids_[tri];
            const float v = (split.axis == 0) ? c.x : (split.axis == 1) ? c.y : c.z;
            cmin = std::min(cmin, v);
            cmax = std::max(cmax, v);
        }

        const float extent = cmax - cmin;
        if (!(extent > 1e-6f)) {
            return start;
        }

        auto bin_index = [&](std::uint32_t tri) {
            const Vec3& c = tri_centroids_[tri];
            const float v = (split.axis == 0) ? c.x : (split.axis == 1) ? c.y : c.z;
            const float t = (v - cmin) / extent;
            int idx = static_cast<int>(t * static_cast<float>(kBins));
            if (idx < 0) idx = 0;
            if (idx >= kBins) idx = kBins - 1;
            return idx;
        };

        auto mid_it = std::partition(
            prim_indices_.begin() + static_cast<std::ptrdiff_t>(start),
            prim_indices_.begin() + static_cast<std::ptrdiff_t>(end),
            [&](std::uint32_t tri) { return bin_index(tri) <= split.bin_split; });

        return static_cast<std::uint32_t>(mid_it - prim_indices_.begin());
    }

    void build_bvh() {
        if (prim_indices_.empty()) {
            root_ = -1;
            return;
        }

        root_ = build_node(0, static_cast<std::uint32_t>(prim_indices_.size()));
    }

    int build_node(std::uint32_t start, std::uint32_t end) {
        const std::uint32_t span = end - start;

        AABB range_box;
        bool first = true;
        for (std::uint32_t i = start; i < end; ++i) {
            const std::uint32_t tri = prim_indices_[i];
            const AABB& b = tri_boxes_[tri];
            if (first) {
                range_box = b;
                first = false;
            } else {
                range_box = surrounding_box(range_box, b);
            }
        }

        Node node;
        node.box = range_box;

        const int node_index = static_cast<int>(nodes_.size());
        nodes_.push_back(node);

        const std::uint32_t leaf_size = 4;
        if (span <= leaf_size) {
            nodes_[node_index].start = start;
            nodes_[node_index].count = span;
            nodes_[node_index].is_leaf = true;
            return node_index;
        }

        const SAHSplit split = find_sah_split(start, end);
        std::uint32_t mid = start;
        if (split.axis >= 0) {
            mid = partition_by_sah_split(start, end, split);
        }

        if (mid == start || mid == end) {
            const int axis = choose_split_axis(range_box);
            mid = start + span / 2;

            auto centroid_less = [&](std::uint32_t a, std::uint32_t b) {
                const Vec3& ca = tri_centroids_[a];
                const Vec3& cb = tri_centroids_[b];
                if (axis == 0) {
                    return ca.x < cb.x;
                } else if (axis == 1) {
                    return ca.y < cb.y;
                }
                return ca.z < cb.z;
            };

            std::nth_element(
                prim_indices_.begin() + static_cast<std::ptrdiff_t>(start),
                prim_indices_.begin() + static_cast<std::ptrdiff_t>(mid),
                prim_indices_.begin() + static_cast<std::ptrdiff_t>(end),
                centroid_less);
        }

        const int left = build_node(start, mid);
        const int right = build_node(mid, end);

        nodes_[node_index].left = left;
        nodes_[node_index].right = right;
        nodes_[node_index].is_leaf = false;
        nodes_[node_index].count = 0;
        nodes_[node_index].box = surrounding_box(nodes_[left].box, nodes_[right].box);

        return node_index;
    }

    bool intersect_triangle(std::uint32_t tri,
                            const Ray& r,
                            float t_min,
                            float t_max,
                            float& out_t,
                            float& out_b1,
                            float& out_b2) const {
        const std::size_t i0 = indices_[3 * static_cast<std::size_t>(tri) + 0];
        const std::size_t i1 = indices_[3 * static_cast<std::size_t>(tri) + 1];
        const std::size_t i2 = indices_[3 * static_cast<std::size_t>(tri) + 2];

        const Vec3& p0 = positions_[i0];
        const Vec3& p1 = positions_[i1];
        const Vec3& p2 = positions_[i2];

        const Vec3 e1 = p1 - p0;
        const Vec3 e2 = p2 - p0;

        const Vec3 pvec = cross(r.direction, e2);
        const float det = dot(e1, pvec);

        if (std::fabs(det) < 1e-10f) {
            return false;
        }

        const float inv_det = 1.0f / det;
        const Vec3 tvec = r.origin - p0;
        const float u = dot(tvec, pvec) * inv_det;
        if (u < 0.0f || u > 1.0f) {
            return false;
        }

        const Vec3 qvec = cross(tvec, e1);
        const float v = dot(r.direction, qvec) * inv_det;
        if (v < 0.0f || u + v > 1.0f) {
            return false;
        }

        const float t = dot(e2, qvec) * inv_det;
        if (t < t_min || t > t_max) {
            return false;
        }

        out_t = t;
        out_b1 = u;
        out_b2 = v;
        return true;
    }

    std::vector<Vec3> positions_;
    std::vector<Vec3> normals_;
    std::vector<Vec2> uvs_;
    std::vector<std::uint32_t> indices_;
    std::vector<Vec3> tangents_;

    AABB bounds_;
    std::vector<AABB> tri_boxes_;
    std::vector<Vec3> tri_centroids_;
    std::vector<std::uint32_t> prim_indices_;
    std::vector<Node> nodes_;
    int root_ = -1;
};

using MeshDataPtr = std::shared_ptr<const MeshData>;

class Mesh : public Hittable {
public:
    Mesh() = default;

    Mesh(const MeshDataPtr& data, const Transform& transform, const MaterialPtr& material)
        : data_(data),
          transform_(transform),
          material_(material) {
        if (data_) {
            world_bounds_ = transform_aabb(data_->bounding_box(), transform_);
        }
        build_light_distribution();
    }

    Mesh(const MeshDataPtr& data, const MaterialPtr& material)
        : Mesh(data, Transform::identity(), material) {}

    bool hit(const Ray& r, float t_min, float t_max, HitRecord& rec) const override {
        if (!data_) {
            return false;
        }

        Ray local_ray = transform_.apply_inverse(r);
        HitRecord local_rec;
        if (!data_->hit(local_ray, t_min, t_max, local_rec)) {
            return false;
        }

        rec = local_rec;
        rec.point = transform_.apply_point(local_rec.point);
        const float dir_denom = dot(r.direction, r.direction);
        if (dir_denom > 1e-20f) {
            rec.t = dot(rec.point - r.origin, r.direction) / dir_denom;
        }

        Vec3 world_normal = transform_.apply_normal(local_rec.normal);
        if (world_normal.length_squared() > 1e-20f) {
            world_normal = normalize(world_normal);
        } else {
            world_normal = local_rec.normal;
        }
        rec.normal = world_normal;

        Vec3 world_tangent = transform_.apply_vector(local_rec.tangent);
        if (world_tangent.length_squared() > 1e-20f) {
            world_tangent = normalize(world_tangent);
        } else {
            world_tangent = local_rec.tangent;
        }

        world_tangent = world_tangent - dot(world_tangent, rec.normal) * rec.normal;
        if (world_tangent.length_squared() < 1e-16f) {
            world_tangent = normalize(cross(
                (std::fabs(rec.normal.x) < 0.9f) ? Vec3(1.0f, 0.0f, 0.0f) : Vec3(0.0f, 1.0f, 0.0f),
                rec.normal));
        } else {
            world_tangent = normalize(world_tangent);
        }
        rec.tangent = world_tangent;

        rec.dpdu = transform_.apply_vector(local_rec.dpdu);
        rec.dpdv = transform_.apply_vector(local_rec.dpdv);

        const float travel = (rec.point - r.origin).length();
        rec.ray_footprint = r.cone.width + r.cone.spread_angle * travel;

        rec.material = material_.get();
        rec.object = this;

        return true;
    }

    AABB bounding_box() const override {
        return world_bounds_;
    }

    const Material* material() const {
        return material_.get();
    }

    float area() const {
        return total_area_;
    }

    bool sample_surface(RNG& rng,
                        Vec3& out_pos,
                        Vec3& out_normal,
                        float& out_u,
                        float& out_v,
                        float& out_pdf_area) const {
        if (!data_ || total_area_ <= 0.0f || tri_cdf_.empty()) {
            return false;
        }

        const float uu = clamp_float(rng.uniform(), 0.0f, 0.99999994f);
        auto it = std::upper_bound(tri_cdf_.begin() + 1, tri_cdf_.end(), uu);
        std::size_t tri = static_cast<std::size_t>((it - tri_cdf_.begin()) - 1);
        if (tri >= data_->triangle_count()) {
            tri = data_->triangle_count() - 1;
        }

        Vec3 p0, p1, p2;
        Vec2 uv0, uv1, uv2;
        if (!data_->triangle_vertices(static_cast<std::uint32_t>(tri), p0, p1, p2, uv0, uv1, uv2)) {
            return false;
        }

        const float r1 = std::sqrt(rng.uniform());
        const float r2 = rng.uniform();
        const float b0 = 1.0f - r1;
        const float b1 = r1 * (1.0f - r2);
        const float b2 = r1 * r2;

        const Vec3 local_pos = b0 * p0 + b1 * p1 + b2 * p2;
        const Vec2 local_uv = b0 * uv0 + b1 * uv1 + b2 * uv2;

        const Vec3 w0 = transform_.apply_point(p0);
        const Vec3 w1 = transform_.apply_point(p1);
        const Vec3 w2 = transform_.apply_point(p2);
        Vec3 gn = cross(w1 - w0, w2 - w0);
        if (gn.length_squared() > 1e-20f) {
            gn = normalize(gn);
        } else {
            gn = Vec3(0.0f, 1.0f, 0.0f);
        }

        out_pos = transform_.apply_point(local_pos);
        out_normal = gn;
        out_u = local_uv.x;
        out_v = local_uv.y;
        out_pdf_area = 1.0f / total_area_;
        return true;
    }

private:
    void build_light_distribution() {
        total_area_ = 0.0f;
        tri_cdf_.clear();
        if (!data_) {
            return;
        }

        const std::size_t tri_count = data_->triangle_count();
        tri_cdf_.assign(tri_count + 1, 0.0f);
        if (tri_count == 0) {
            return;
        }

        for (std::size_t tri = 0; tri < tri_count; ++tri) {
            Vec3 p0, p1, p2;
            Vec2 uv0, uv1, uv2;
            if (!data_->triangle_vertices(static_cast<std::uint32_t>(tri), p0, p1, p2, uv0, uv1, uv2)) {
                tri_cdf_[tri + 1] = total_area_;
                continue;
            }

            const Vec3 w0 = transform_.apply_point(p0);
            const Vec3 w1 = transform_.apply_point(p1);
            const Vec3 w2 = transform_.apply_point(p2);
            const float a = 0.5f * cross(w1 - w0, w2 - w0).length();
            total_area_ += std::max(0.0f, a);
            tri_cdf_[tri + 1] = total_area_;
        }

        if (total_area_ <= 0.0f) {
            tri_cdf_.clear();
            return;
        }

        const float inv = 1.0f / total_area_;
        for (std::size_t i = 1; i < tri_cdf_.size(); ++i) {
            tri_cdf_[i] *= inv;
        }
    }

    MeshDataPtr data_;
    Transform transform_;
    AABB world_bounds_;
    MaterialPtr material_;
    float total_area_ = 0.0f;
    std::vector<float> tri_cdf_;
};
