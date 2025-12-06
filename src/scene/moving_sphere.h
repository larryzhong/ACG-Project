#pragma once

#include "core/vec3.h"
#include "math/aabb.h"
#include "scene/hittable.h"
#include "scene/material.h"
#include "scene/sphere.h"

class MovingSphere : public Hittable {
public:
    MovingSphere() = default;

    MovingSphere(const Vec3& center0,
                 const Vec3& center1,
                 float time0,
                 float time1,
                 float radius,
                 const MaterialPtr& material)
        : center0_(center0),
          center1_(center1),
          time0_(time0),
          time1_(time1),
          radius_(radius),
          material_(material) {}

    Vec3 center(float time) const {
        if (time1_ == time0_) {
            return center0_;
        }
        const float alpha = (time - time0_) / (time1_ - time0_);
        return center0_ + alpha * (center1_ - center0_);
    }

    bool hit(const Ray& r, float t_min, float t_max, HitRecord& rec) const override {
        const Vec3 c = center(r.time);
        const Vec3 oc = r.origin - c;
        const float a = dot(r.direction, r.direction);
        const float half_b = dot(oc, r.direction);
        const float c_term = dot(oc, oc) - radius_ * radius_;
        const float discriminant = half_b * half_b - a * c_term;
        if (discriminant < 0.0f) {
            return false;
        }
        const float sqrt_d = std::sqrt(discriminant);

        float root = (-half_b - sqrt_d) / a;
        if (root < t_min || root > t_max) {
            root = (-half_b + sqrt_d) / a;
            if (root < t_min || root > t_max) {
                return false;
            }
        }

        rec.t = root;
        rec.point = r.at(rec.t);
        const Vec3 outward_normal = (rec.point - c) / radius_;
        rec.set_face_normal(r, outward_normal);
        get_sphere_uv(outward_normal, rec.u, rec.v);
        rec.material = material_.get();

        rec.tangent = normalize(Vec3(-outward_normal.z, 0.0f, outward_normal.x));
        if (rec.tangent.length_squared() < 1e-6f) {
            rec.tangent = Vec3(1.0f, 0.0f, 0.0f);
        }

        return true;
    }

    AABB bounding_box() const override {
        const Vec3 radius_vec(radius_, radius_, radius_);
        const AABB box0(center0_ - radius_vec, center0_ + radius_vec);
        const AABB box1(center1_ - radius_vec, center1_ + radius_vec);
        return surrounding_box(box0, box1);
    }

private:
    Vec3 center0_;
    Vec3 center1_;
    float time0_;
    float time1_;
    float radius_;
    MaterialPtr material_;
};
