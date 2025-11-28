#pragma once

#include "core/vec3.h"
#include "scene/hittable.h"
#include "scene/material.h"

inline void get_sphere_uv(const Vec3& p, float& u, float& v) {
    const float theta = std::acos(-p.y);
    const float phi = std::atan2(-p.z, p.x) + kPi;
    u = phi / (2.0f * kPi);
    v = theta / kPi;
}

class Sphere : public Hittable {
public:
    Sphere() : center_(0.0f), radius_(0.5f) {}

    Sphere(const Vec3& center, float radius, const MaterialPtr& material)
        : center_(center), radius_(radius), material_(material) {}

    bool hit(const Ray& r, float t_min, float t_max, HitRecord& rec) const override {
        const Vec3 oc = r.origin - center_;
        const float a = dot(r.direction, r.direction);
        const float half_b = dot(oc, r.direction);
        const float c = dot(oc, oc) - radius_ * radius_;
        const float discriminant = half_b * half_b - a * c;
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
        const Vec3 outward_normal = (rec.point - center_) / radius_;
        rec.set_face_normal(r, outward_normal);
        get_sphere_uv(outward_normal, rec.u, rec.v);
        rec.material = material_.get();
        return true;
    }

    AABB bounding_box() const override {
        const Vec3 radius_vec(radius_, radius_, radius_);
        return AABB(center_ - radius_vec, center_ + radius_vec);
    }

private:
    Vec3 center_;
    float radius_;
    MaterialPtr material_;
};
