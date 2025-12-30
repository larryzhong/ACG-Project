#pragma once

#include <memory>
#include <vector>

#include "scene/bvh.h"
#include "scene/environment.h"
#include "scene/hittable.h"
#include "scene/light.h"

class Medium;

class Scene {
public:
    std::vector<HittablePtr> objects;
    HittablePtr accel;
    LightCollection lights;
    EnvironmentMapPtr environment;
    std::shared_ptr<Medium> global_medium;
    bool hide_environment_background = false;

    void build_bvh() {
        if (objects.empty()) {
            accel.reset();
            return;
        }

        std::vector<HittablePtr> copy = objects;
        accel = std::make_shared<BVHNode>(copy, 0, copy.size());
    }

    bool hit(const Ray& r, float t_min, float t_max, HitRecord& rec) const {
        if (accel) {
            return accel->hit(r, t_min, t_max, rec);
        }

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
