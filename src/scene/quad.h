#pragma once

#include <memory>

#include "core/vec3.h"
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
        rec.set_face_normal(r, normal_);

        return true;
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
    MaterialPtr material_;
};

