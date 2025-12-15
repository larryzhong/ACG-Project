#pragma once

#include <cmath>

struct Vec2 {
    float x;
    float y;

    Vec2() : x(0.0f), y(0.0f) {}
    Vec2(float v) : x(v), y(v) {}
    Vec2(float x_, float y_) : x(x_), y(y_) {}

    Vec2& operator+=(const Vec2& rhs) {
        x += rhs.x;
        y += rhs.y;
        return *this;
    }

    Vec2& operator-=(const Vec2& rhs) {
        x -= rhs.x;
        y -= rhs.y;
        return *this;
    }

    Vec2& operator*=(float s) {
        x *= s;
        y *= s;
        return *this;
    }

    Vec2& operator/=(float s) {
        const float inv = 1.0f / s;
        x *= inv;
        y *= inv;
        return *this;
    }

    float length() const {
        return std::sqrt(length_squared());
    }

    float length_squared() const {
        return x * x + y * y;
    }

    static Vec2 zero() {
        return Vec2(0.0f);
    }
};

inline Vec2 operator+(Vec2 lhs, const Vec2& rhs) {
    lhs += rhs;
    return lhs;
}

inline Vec2 operator-(const Vec2& v) {
    return Vec2(-v.x, -v.y);
}

inline Vec2 operator-(Vec2 lhs, const Vec2& rhs) {
    lhs -= rhs;
    return lhs;
}

inline Vec2 operator*(Vec2 v, float s) {
    v *= s;
    return v;
}

inline Vec2 operator*(float s, Vec2 v) {
    v *= s;
    return v;
}

inline Vec2 operator/(Vec2 v, float s) {
    v /= s;
    return v;
}

