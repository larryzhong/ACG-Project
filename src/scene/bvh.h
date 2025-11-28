#pragma once

#include <algorithm>
#include <memory>
#include <vector>

#include "math/aabb.h"
#include "scene/hittable.h"

class BVHNode : public Hittable {
public:
    BVHNode() = default;

    BVHNode(std::vector<HittablePtr>& objects, std::size_t start, std::size_t end) {
        const int axis = choose_split_axis(objects, start, end);

        auto comparator = [axis](const HittablePtr& a, const HittablePtr& b) {
            const AABB box_a = a->bounding_box();
            const AABB box_b = b->bounding_box();

            float center_a = 0.0f;
            float center_b = 0.0f;

            if (axis == 0) {
                center_a = 0.5f * (box_a.min.x + box_a.max.x);
                center_b = 0.5f * (box_b.min.x + box_b.max.x);
            } else if (axis == 1) {
                center_a = 0.5f * (box_a.min.y + box_a.max.y);
                center_b = 0.5f * (box_b.min.y + box_b.max.y);
            } else {
                center_a = 0.5f * (box_a.min.z + box_a.max.z);
                center_b = 0.5f * (box_b.min.z + box_b.max.z);
            }

            return center_a < center_b;
        };

        const std::size_t object_span = end - start;

        if (object_span == 1) {
            left_ = right_ = objects[start];
        } else if (object_span == 2) {
            if (comparator(objects[start], objects[start + 1])) {
                left_ = objects[start];
                right_ = objects[start + 1];
            } else {
                left_ = objects[start + 1];
                right_ = objects[start];
            }
        } else {
            const std::size_t mid = start + object_span / 2;
            std::nth_element(objects.begin() + static_cast<std::ptrdiff_t>(start),
                             objects.begin() + static_cast<std::ptrdiff_t>(mid),
                             objects.begin() + static_cast<std::ptrdiff_t>(end),
                             comparator);

            left_ = std::make_shared<BVHNode>(objects, start, mid);
            right_ = std::make_shared<BVHNode>(objects, mid, end);
        }

        const AABB box_left = left_->bounding_box();
        const AABB box_right = right_->bounding_box();
        box_ = surrounding_box(box_left, box_right);
    }

    bool hit(const Ray& r, float t_min, float t_max, HitRecord& rec) const override {
        if (!box_.hit(r, t_min, t_max)) {
            return false;
        }

        HitRecord left_rec;
        HitRecord right_rec;
        bool hit_left = left_->hit(r, t_min, t_max, left_rec);
        bool hit_right = right_->hit(r, t_min, t_max, right_rec);

        if (hit_left && hit_right) {
            if (left_rec.t < right_rec.t) {
                rec = left_rec;
            } else {
                rec = right_rec;
            }
            return true;
        } else if (hit_left) {
            rec = left_rec;
            return true;
        } else if (hit_right) {
            rec = right_rec;
            return true;
        }

        return false;
    }

    AABB bounding_box() const override {
        return box_;
    }

private:
    static int choose_split_axis(const std::vector<HittablePtr>& objects,
                                 std::size_t start,
                                 std::size_t end) {
        AABB bounds;
        bool first_box = true;

        for (std::size_t i = start; i < end; ++i) {
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

    HittablePtr left_;
    HittablePtr right_;
    AABB box_;
};
