#pragma once

#include <cmath>
#include "core/vec3.h"

class ONB {
public:
    ONB() = default;

    void build_from_w(const Vec3& n) {
        w_ = normalize(n);
        Vec3 a = (std::fabs(w_.x) > 0.9f) ? Vec3(0.0f, 1.0f, 0.0f) : Vec3(1.0f, 0.0f, 0.0f);
        v_ = normalize(cross(w_, a));
        u_ = cross(w_, v_);
    }

    Vec3 local(const Vec3& a) const {
        return a.x * u_ + a.y * v_ + a.z * w_;
    }

    const Vec3& u() const { return u_; }
    const Vec3& v() const { return v_; }
    const Vec3& w() const { return w_; }

private:
    Vec3 u_, v_, w_;
};
