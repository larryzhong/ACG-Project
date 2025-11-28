#pragma once

#include <memory>
#include <vector>

#include "scene/hittable.h"

class Scene {
public:
    std::vector<HittablePtr> objects;

    bool hit(const Ray& r, float t_min, float t_max, HitRecord& rec) const {
        HitRecord temp_rec;
        bool hit_anything = false;
        float closest_so_far = t_max;

        for (const auto& obj : objects) {
            if (obj && obj->hit(r, t_min, closest_so_far, temp_rec)) {
                hit_anything = true;
                closest_so_far = temp_rec.t;
                rec = temp_rec;
            }
        }

        return hit_anything;
    }
};
