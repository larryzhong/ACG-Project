#pragma once

#include <memory>

#include "scene/hittable.h"
#include "scene/medium.h"

// Wraps a surface boundary and annotates hit records with a medium interface.
// The integrator can use hit.front_face to choose which medium the ray enters.
class MediumBoundary : public Hittable {
public:
    MediumBoundary(HittablePtr boundary,
                  MediumPtr inside,
                  MediumPtr outside = nullptr)
        : boundary_(std::move(boundary)),
          inside_(std::move(inside)),
          outside_(std::move(outside)) {}

    bool hit(const Ray& r, float t_min, float t_max, HitRecord& rec) const override {
        if (!boundary_ || !boundary_->hit(r, t_min, t_max, rec)) {
            return false;
        }

        rec.medium_inside = inside_.get();
        rec.medium_outside = outside_.get();
        return true;
    }

    AABB bounding_box() const override {
        return boundary_ ? boundary_->bounding_box() : AABB();
    }

private:
    HittablePtr boundary_;
    MediumPtr inside_;
    MediumPtr outside_;
};
