#pragma once

#include "core/vec3.h"

class Medium;

struct RayCone {
    float width = 0.0f;        // world-space radius at ray origin
    float spread_angle = 0.0f; // radians per unit distance
};

struct Ray {
    Vec3 origin;
    Vec3 direction;
    float time;
    RayCone cone;
    const Medium* medium = nullptr;

    Ray() : origin(), direction(0.0f, 0.0f, 1.0f), time(0.0f), cone(), medium(nullptr) {}

    Ray(const Vec3& o, const Vec3& d, float t = 0.0f)
        : origin(o), direction(d), time(t), cone(), medium(nullptr) {}

    Ray(const Vec3& o, const Vec3& d, float t, const RayCone& c)
        : origin(o), direction(d), time(t), cone(c), medium(nullptr) {}

    Vec3 at(float t) const {
        return origin + t * direction;
    }
};
