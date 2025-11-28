#pragma once

#include <cmath>

struct Vec3 {
    float x;
    float y;
    float z;

    Vec3() : x(0.0f), y(0.0f), z(0.0f) {}
    Vec3(float v) : x(v), y(v), z(v) {}
    Vec3(float x_, float y_, float z_) : x(x_), y(y_), z(z_) {}

    Vec3& operator+=(const Vec3& rhs) {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        return *this;
    }

    Vec3& operator-=(const Vec3& rhs) {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        return *this;
    }

    Vec3& operator*=(float s) {
        x *= s;
        y *= s;
        z *= s;
        return *this;
    }

    Vec3& operator/=(float s) {
        const float inv = 1.0f / s;
        x *= inv;
        y *= inv;
        z *= inv;
        return *this;
    }

    float length() const {
        return std::sqrt(length_squared());
    }

    float length_squared() const {
        return x * x + y * y + z * z;
    }

    static Vec3 zero() {
        return Vec3(0.0f);
    }
};

inline Vec3 operator+(Vec3 lhs, const Vec3& rhs) {
    lhs += rhs;
    return lhs;
}

inline Vec3 operator-(const Vec3& v) {
    return Vec3(-v.x, -v.y, -v.z);
}

inline Vec3 operator-(Vec3 lhs, const Vec3& rhs) {
    lhs -= rhs;
    return lhs;
}

inline Vec3 operator*(Vec3 v, float s) {
    v *= s;
    return v;
}

inline Vec3 operator*(float s, Vec3 v) {
    v *= s;
    return v;
}

inline Vec3 operator/(Vec3 v, float s) {
    v /= s;
    return v;
}

inline Vec3 operator*(const Vec3& a, const Vec3& b) {
    return Vec3(a.x * b.x, a.y * b.y, a.z * b.z);
}

inline float dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return Vec3(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    );
}

inline Vec3 normalize(const Vec3& v) {
    const float len = v.length();
    if (len == 0.0f) {
        return Vec3::zero();
    }
    return v / len;
}

inline Vec3 reflect(const Vec3& v, const Vec3& n) {
    return v - 2.0f * dot(v, n) * n;
}
