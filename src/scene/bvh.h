#pragma once

#include <algorithm>
#include <array>
#include <cfloat>
#include <memory>
#include <vector>

#include "math/aabb.h"
#include "scene/hittable.h"

class BVHNode : public Hittable {
public:
    BVHNode() = default;

    BVHNode(std::vector<HittablePtr>& objects, std::size_t start, std::size_t end) {
        const std::size_t object_span = end - start;
        if (object_span == 0) {
            return;
        }

        constexpr std::size_t kLeafSize = 4;
        if (object_span <= kLeafSize) {
            is_leaf_ = true;
            leaf_objects_.reserve(object_span);

            bool first = true;
            for (std::size_t i = start; i < end; ++i) {
                if (!objects[i]) {
                    continue;
                }
                leaf_objects_.push_back(objects[i]);
                const AABB b = objects[i]->bounding_box();
                if (first) {
                    box_ = b;
                    first = false;
                } else {
                    box_ = surrounding_box(box_, b);
                }
            }
            return;
        }

        const Split split = find_sah_split(objects, start, end);
        std::size_t mid = start;

        if (split.axis >= 0) {
            mid = partition_by_split(objects, start, end, split);
        }

        if (mid == start || mid == end) {
            const int axis = choose_split_axis(objects, start, end);
            const std::size_t fallback_mid = start + object_span / 2;
            auto comparator = [axis](const HittablePtr& a, const HittablePtr& b) {
                return centroid_axis(a->bounding_box(), axis) <
                       centroid_axis(b->bounding_box(), axis);
            };
            std::nth_element(objects.begin() + static_cast<std::ptrdiff_t>(start),
                             objects.begin() + static_cast<std::ptrdiff_t>(fallback_mid),
                             objects.begin() + static_cast<std::ptrdiff_t>(end),
                             comparator);
            mid = fallback_mid;
        }

        left_ = std::make_shared<BVHNode>(objects, start, mid);
        right_ = std::make_shared<BVHNode>(objects, mid, end);

        const AABB box_left = left_ ? left_->bounding_box() : AABB();
        const AABB box_right = right_ ? right_->bounding_box() : AABB();
        box_ = surrounding_box(box_left, box_right);
    }

    bool hit(const Ray& r, float t_min, float t_max, HitRecord& rec) const override {
        if (!box_.hit(r, t_min, t_max)) {
            return false;
        }

        if (is_leaf_) {
            bool hit_anything = false;
            float closest_so_far = t_max;
            HitRecord temp;
            for (const auto& obj : leaf_objects_) {
                if (obj && obj->hit(r, t_min, closest_so_far, temp)) {
                    hit_anything = true;
                    closest_so_far = temp.t;
                    rec = temp;
                }
            }
            return hit_anything;
        }

        HitRecord left_rec;
        HitRecord right_rec;
        const bool hit_left = left_ && left_->hit(r, t_min, t_max, left_rec);
        const bool hit_right = right_ && right_->hit(r, t_min, t_max, right_rec);

        if (hit_left && hit_right) {
            rec = (left_rec.t < right_rec.t) ? left_rec : right_rec;
            return true;
        }
        if (hit_left) {
            rec = left_rec;
            return true;
        }
        if (hit_right) {
            rec = right_rec;
            return true;
        }
        return false;
    }

    AABB bounding_box() const override {
        return box_;
    }

private:
    struct Split {
        int axis = -1;
        int bin_split = -1;
    };

    static float centroid_axis(const AABB& b, int axis) {
        if (axis == 0) {
            return 0.5f * (b.min.x + b.max.x);
        }
        if (axis == 1) {
            return 0.5f * (b.min.y + b.max.y);
        }
        return 0.5f * (b.min.z + b.max.z);
    }

    static float surface_area(const AABB& b) {
        const float dx = std::max(0.0f, b.max.x - b.min.x);
        const float dy = std::max(0.0f, b.max.y - b.min.y);
        const float dz = std::max(0.0f, b.max.z - b.min.z);
        return 2.0f * (dx * dy + dy * dz + dz * dx);
    }

    static int choose_split_axis(const std::vector<HittablePtr>& objects,
                                 std::size_t start,
                                 std::size_t end) {
        AABB bounds;
        bool first_box = true;

        for (std::size_t i = start; i < end; ++i) {
            if (!objects[i]) continue;
            const AABB box = objects[i]->bounding_box();
            if (first_box) {
                bounds = box;
                first_box = false;
            } else {
                bounds = surrounding_box(bounds, box);
            }
        }

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

    static AABB bounds_range(const std::vector<HittablePtr>& objects,
                             std::size_t start,
                             std::size_t end) {
        bool first = true;
        AABB out;
        for (std::size_t i = start; i < end; ++i) {
            if (!objects[i]) continue;
            const AABB b = objects[i]->bounding_box();
            if (first) {
                out = b;
                first = false;
            } else {
                out = surrounding_box(out, b);
            }
        }
        return out;
    }

    static Split find_sah_split(const std::vector<HittablePtr>& objects,
                                std::size_t start,
                                std::size_t end) {
        constexpr int kBins = 16;
        const std::size_t n = end - start;
        if (n < 2) {
            return {};
        }

        const AABB parent_bounds = bounds_range(objects, start, end);
        const float parent_area = surface_area(parent_bounds);
        if (parent_area <= 0.0f) {
            return {};
        }

        Split best;
        float best_cost = FLT_MAX;

        for (int axis = 0; axis < 3; ++axis) {
            float cmin = FLT_MAX;
            float cmax = -FLT_MAX;
            for (std::size_t i = start; i < end; ++i) {
                if (!objects[i]) continue;
                const float c = centroid_axis(objects[i]->bounding_box(), axis);
                cmin = std::min(cmin, c);
                cmax = std::max(cmax, c);
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

            for (std::size_t i = start; i < end; ++i) {
                if (!objects[i]) continue;
                const AABB b = objects[i]->bounding_box();
                const float c = centroid_axis(b, axis);
                const int idx = bin_index(c);
                Bin& bin = bins[static_cast<std::size_t>(idx)];
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
                const float cost = 1.0f + (left_area / parent_area) * static_cast<float>(left_count) +
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

    static std::size_t partition_by_split(std::vector<HittablePtr>& objects,
                                          std::size_t start,
                                          std::size_t end,
                                          const Split& split) {
        constexpr int kBins = 16;
        float cmin = FLT_MAX;
        float cmax = -FLT_MAX;
        for (std::size_t i = start; i < end; ++i) {
            if (!objects[i]) continue;
            const float c = centroid_axis(objects[i]->bounding_box(), split.axis);
            cmin = std::min(cmin, c);
            cmax = std::max(cmax, c);
        }
        const float extent = cmax - cmin;
        if (!(extent > 1e-6f)) {
            return start;
        }

        auto bin_index = [&](const HittablePtr& obj) {
            const float c = centroid_axis(obj->bounding_box(), split.axis);
            const float t = (c - cmin) / extent;
            int idx = static_cast<int>(t * static_cast<float>(kBins));
            if (idx < 0) idx = 0;
            if (idx >= kBins) idx = kBins - 1;
            return idx;
        };

        auto mid_it = std::partition(objects.begin() + static_cast<std::ptrdiff_t>(start),
                                     objects.begin() + static_cast<std::ptrdiff_t>(end),
                                     [&](const HittablePtr& obj) {
                                         if (!obj) return true;
                                         return bin_index(obj) <= split.bin_split;
                                     });
        return static_cast<std::size_t>(mid_it - objects.begin());
    }

    HittablePtr left_;
    HittablePtr right_;
    AABB box_;
    bool is_leaf_ = false;
    std::vector<HittablePtr> leaf_objects_;
};
