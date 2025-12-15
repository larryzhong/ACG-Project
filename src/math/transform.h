#pragma once

#include <algorithm>
#include <cmath>

#include "core/ray.h"
#include "core/vec3.h"
#include "math/aabb.h"

struct Mat3 {
    Vec3 c0;
    Vec3 c1;
    Vec3 c2;

    Mat3()
        : c0(1.0f, 0.0f, 0.0f),
          c1(0.0f, 1.0f, 0.0f),
          c2(0.0f, 0.0f, 1.0f) {}

    Mat3(const Vec3& c0_, const Vec3& c1_, const Vec3& c2_)
        : c0(c0_), c1(c1_), c2(c2_) {}

    static Mat3 identity() {
        return Mat3();
    }

    static Mat3 diagonal(const Vec3& d) {
        return Mat3(
            Vec3(d.x, 0.0f, 0.0f),
            Vec3(0.0f, d.y, 0.0f),
            Vec3(0.0f, 0.0f, d.z));
    }
};

inline Vec3 operator*(const Mat3& m, const Vec3& v) {
    return m.c0 * v.x + m.c1 * v.y + m.c2 * v.z;
}

inline Mat3 operator*(const Mat3& a, const Mat3& b) {
    return Mat3(a * b.c0, a * b.c1, a * b.c2);
}

inline Mat3 transpose(const Mat3& m) {
    return Mat3(
        Vec3(m.c0.x, m.c1.x, m.c2.x),
        Vec3(m.c0.y, m.c1.y, m.c2.y),
        Vec3(m.c0.z, m.c1.z, m.c2.z));
}

inline float determinant(const Mat3& m) {
    return dot(m.c0, cross(m.c1, m.c2));
}

inline Mat3 inverse(const Mat3& m) {
    const Vec3 r0 = cross(m.c1, m.c2);
    const Vec3 r1 = cross(m.c2, m.c0);
    const Vec3 r2 = cross(m.c0, m.c1);

    const float det = dot(m.c0, r0);
    if (std::fabs(det) < 1e-12f) {
        return Mat3::identity();
    }

    const float inv_det = 1.0f / det;
    const Vec3 rr0 = r0 * inv_det;
    const Vec3 rr1 = r1 * inv_det;
    const Vec3 rr2 = r2 * inv_det;

    return Mat3(
        Vec3(rr0.x, rr1.x, rr2.x),
        Vec3(rr0.y, rr1.y, rr2.y),
        Vec3(rr0.z, rr1.z, rr2.z));
}

struct Transform {
    Mat3 m;
    Vec3 t;
    Mat3 m_inv;
    Vec3 t_inv;

    Transform()
        : m(Mat3::identity()),
          t(Vec3::zero()),
          m_inv(Mat3::identity()),
          t_inv(Vec3::zero()) {}

    Transform(const Mat3& linear, const Vec3& translation)
        : m(linear),
          t(translation),
          m_inv(inverse(linear)),
          t_inv(Vec3::zero()) {
        t_inv = -(m_inv * t);
    }

    Transform(const Mat3& linear, const Vec3& translation, const Mat3& linear_inv)
        : m(linear),
          t(translation),
          m_inv(linear_inv),
          t_inv(Vec3::zero()) {
        t_inv = -(m_inv * t);
    }

    static Transform identity() {
        return Transform();
    }

    static Transform translate(const Vec3& delta) {
        Transform tr;
        tr.t = delta;
        tr.t_inv = -delta;
        return tr;
    }

    static Transform scale(const Vec3& s) {
        Mat3 linear = Mat3::diagonal(s);
        Mat3 linear_inv = Mat3::diagonal(Vec3(1.0f / s.x, 1.0f / s.y, 1.0f / s.z));
        return Transform(linear, Vec3::zero(), linear_inv);
    }

    static Transform uniform_scale(float s) {
        return scale(Vec3(s, s, s));
    }

    static Transform rotate_x(float radians) {
        const float c = std::cos(radians);
        const float s = std::sin(radians);
        Mat3 linear(
            Vec3(1.0f, 0.0f, 0.0f),
            Vec3(0.0f, c, s),
            Vec3(0.0f, -s, c));
        return Transform(linear, Vec3::zero(), transpose(linear));
    }

    static Transform rotate_y(float radians) {
        const float c = std::cos(radians);
        const float s = std::sin(radians);
        Mat3 linear(
            Vec3(c, 0.0f, -s),
            Vec3(0.0f, 1.0f, 0.0f),
            Vec3(s, 0.0f, c));
        return Transform(linear, Vec3::zero(), transpose(linear));
    }

    static Transform rotate_z(float radians) {
        const float c = std::cos(radians);
        const float s = std::sin(radians);
        Mat3 linear(
            Vec3(c, s, 0.0f),
            Vec3(-s, c, 0.0f),
            Vec3(0.0f, 0.0f, 1.0f));
        return Transform(linear, Vec3::zero(), transpose(linear));
    }

    Vec3 apply_point(const Vec3& p) const {
        return m * p + t;
    }

    Vec3 apply_vector(const Vec3& v) const {
        return m * v;
    }

    Vec3 apply_normal(const Vec3& n) const {
        return transpose(m_inv) * n;
    }

    Vec3 apply_point_inverse(const Vec3& p) const {
        return m_inv * p + t_inv;
    }

    Vec3 apply_vector_inverse(const Vec3& v) const {
        return m_inv * v;
    }

    Ray apply_inverse(const Ray& r) const {
        return Ray(apply_point_inverse(r.origin), apply_vector_inverse(r.direction), r.time);
    }
};

inline Transform operator*(const Transform& a, const Transform& b) {
    Transform out;
    out.m = a.m * b.m;
    out.t = a.m * b.t + a.t;

    out.m_inv = b.m_inv * a.m_inv;
    out.t_inv = b.m_inv * a.t_inv + b.t_inv;

    return out;
}

inline AABB transform_aabb(const AABB& box, const Transform& tr) {
    Vec3 p000(box.min.x, box.min.y, box.min.z);
    Vec3 p100(box.max.x, box.min.y, box.min.z);
    Vec3 p010(box.min.x, box.max.y, box.min.z);
    Vec3 p110(box.max.x, box.max.y, box.min.z);
    Vec3 p001(box.min.x, box.min.y, box.max.z);
    Vec3 p101(box.max.x, box.min.y, box.max.z);
    Vec3 p011(box.min.x, box.max.y, box.max.z);
    Vec3 p111(box.max.x, box.max.y, box.max.z);

    Vec3 corners[8] = {
        tr.apply_point(p000),
        tr.apply_point(p100),
        tr.apply_point(p010),
        tr.apply_point(p110),
        tr.apply_point(p001),
        tr.apply_point(p101),
        tr.apply_point(p011),
        tr.apply_point(p111),
    };

    Vec3 minv = corners[0];
    Vec3 maxv = corners[0];

    for (int i = 1; i < 8; ++i) {
        minv.x = std::min(minv.x, corners[i].x);
        minv.y = std::min(minv.y, corners[i].y);
        minv.z = std::min(minv.z, corners[i].z);

        maxv.x = std::max(maxv.x, corners[i].x);
        maxv.y = std::max(maxv.y, corners[i].y);
        maxv.z = std::max(maxv.z, corners[i].z);
    }

    return AABB(minv, maxv);
}

