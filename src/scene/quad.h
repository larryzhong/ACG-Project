#pragma once

#include <memory>

#include "core/vec3.h"
#include "math/aabb.h"
#include "scene/hittable.h"
#include "scene/material.h"

class Quad : public Hittable {
public:
    Quad() = default;

    Quad(const Vec3& p0,
         const Vec3& u,
         const Vec3& v,
         const MaterialPtr& material)
        : p0_(p0),
          u_(u),
          v_(v),
          material_(material) {
        normal_ = normalize(cross(u_, v_));
        d_ = dot(normal_, p0_);
        uu_ = dot(u_, u_);
        uv_ = dot(u_, v_);
        vv_ = dot(v_, v_);
        denom_ = uu_ * vv_ - uv_ * uv_;

        const Vec3 p1 = p0_;
        const Vec3 p2 = p0_ + u_;
        const Vec3 p3 = p0_ + v_;
        const Vec3 p4 = p0_ + u_ + v_;

        const float min_x = std::min(std::min(p1.x, p2.x), std::min(p3.x, p4.x));
        const float min_y = std::min(std::min(p1.y, p2.y), std::min(p3.y, p4.y));
        const float min_z = std::min(std::min(p1.z, p2.z), std::min(p3.z, p4.z));

        const float max_x = std::max(std::max(p1.x, p2.x), std::max(p3.x, p4.x));
        const float max_y = std::max(std::max(p1.y, p2.y), std::max(p3.y, p4.y));
        const float max_z = std::max(std::max(p1.z, p2.z), std::max(p3.z, p4.z));

        // Slightly expand box to avoid zero-thickness issues.
        const float epsilon = 1e-4f;
        box_ = AABB(
            Vec3(min_x - epsilon, min_y - epsilon, min_z - epsilon),
            Vec3(max_x + epsilon, max_y + epsilon, max_z + epsilon));
    }

    const Vec3& p0() const {
        return p0_;
    }

    const Vec3& u() const {
        return u_;
    }

    const Vec3& v() const {
        return v_;
    }

    const Vec3& normal() const {
        return normal_;
    }

    float area() const {
        return cross(u_, v_).length();
    }

    const Material* material() const {
        return material_.get();
    }

    bool hit(const Ray& r, float t_min, float t_max, HitRecord& rec) const override {
        const float denom = dot(normal_, r.direction);
        if (std::fabs(denom) < 1e-8f) {
            return false;
        }

        const float t = (d_ - dot(normal_, r.origin)) / denom;
        if (t < t_min || t > t_max) {
            return false;
        }

        const Vec3 p = r.at(t);
        const Vec3 rel = p - p0_;

        const float wu = dot(rel, u_);
        const float wv = dot(rel, v_);

        if (std::fabs(denom_) < 1e-8f) {
            return false;
        }

        const float a = (wu * vv_ - wv * uv_) / denom_;
        const float b = (wv * uu_ - wu * uv_) / denom_;

        if (a < 0.0f || a > 1.0f || b < 0.0f || b > 1.0f) {
            return false;
        }

        rec.t = t;
        rec.point = p;
        rec.u = a;
        rec.v = b;
        rec.material = material_.get();
        rec.object = this;
        rec.set_face_normal(r, normal_);

        rec.tangent = normalize(u_);

        rec.dpdu = u_;
        rec.dpdv = v_;
        const float travel = (rec.point - r.origin).length();
        rec.ray_footprint = r.cone.width + r.cone.spread_angle * travel;

        return true;
    }

    AABB bounding_box() const override {
        return box_;
    }

private:
    Vec3 p0_;
    Vec3 u_;
    Vec3 v_;
    Vec3 normal_;
    float d_ = 0.0f;
    float uu_ = 0.0f;
    float uv_ = 0.0f;
    float vv_ = 0.0f;
    float denom_ = 0.0f;
    AABB box_;
    MaterialPtr material_;
};
