#pragma once

#include <algorithm>

#include "core/ray.h"
#include "core/vec3.h"

class AABB {
public:
    Vec3 min;
    Vec3 max;

    AABB() : min(Vec3::zero()), max(Vec3::zero()) {}

    AABB(const Vec3& a, const Vec3& b) : min(a), max(b) {}

    bool hit(const Ray& r, float t_min, float t_max) const {
        for (int axis = 0; axis < 3; ++axis) {
            float origin = 0.0f;
            float direction = 0.0f;
            float axis_min = 0.0f;
            float axis_max = 0.0f;

            if (axis == 0) {
                origin = r.origin.x;
                direction = r.direction.x;
                axis_min = min.x;
                axis_max = max.x;
            } else if (axis == 1) {
                origin = r.origin.y;
                direction = r.direction.y;
                axis_min = min.y;
                axis_max = max.y;
            } else {
                origin = r.origin.z;
                direction = r.direction.z;
                axis_min = min.z;
                axis_max = max.z;
            }

            const float inv_d = 1.0f / direction;
            float t0 = (axis_min - origin) * inv_d;
            float t1 = (axis_max - origin) * inv_d;

            if (inv_d < 0.0f) {
                std::swap(t0, t1);
            }

            t_min = t0 > t_min ? t0 : t_min;
            t_max = t1 < t_max ? t1 : t_max;

            if (t_max <= t_min) {
                return false;
            }
        }

        return true;
    }
};

inline AABB surrounding_box(const AABB& box0, const AABB& box1) {
    const Vec3 small(
        std::min(box0.min.x, box1.min.x),
        std::min(box0.min.y, box1.min.y),
        std::min(box0.min.z, box1.min.z));

    const Vec3 big(
        std::max(box0.max.x, box1.max.x),
        std::max(box0.max.y, box1.max.y),
        std::max(box0.max.z, box1.max.z));

    return AABB(small, big);
}
