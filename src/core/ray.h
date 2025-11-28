#pragma once

#include "core/vec3.h"

struct Ray {
    Vec3 origin;
    Vec3 direction;
    float time;

    Ray() : origin(), direction(0.0f, 0.0f, 1.0f), time(0.0f) {}

    Ray(const Vec3& o, const Vec3& d, float t = 0.0f)
        : origin(o), direction(d), time(t) {}

    Vec3 at(float t) const {
        return origin + t * direction;
    }
};

