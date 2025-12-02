#pragma once

#include "core/rng.h"
#include "core/vec3.h"

inline Vec3 random_in_unit_sphere(RNG& rng) {
    for (;;) {
        const float x = 2.0f * rng.uniform() - 1.0f;
        const float y = 2.0f * rng.uniform() - 1.0f;
        const float z = 2.0f * rng.uniform() - 1.0f;
        const Vec3 p(x, y, z);
        if (p.length_squared() < 1.0f) {
            return p;
        }
    }
}

inline Vec3 random_unit_vector(RNG& rng) {
    return normalize(random_in_unit_sphere(rng));
}

inline Vec3 random_in_hemisphere(const Vec3& normal, RNG& rng) {
    Vec3 in_sphere = random_in_unit_sphere(rng);
    if (dot(in_sphere, normal) < 0.0f) {
        in_sphere = -in_sphere;
    }
    return in_sphere;
}

inline Vec3 random_cosine_direction(RNG& rng) {
    const float r1 = rng.uniform();
    const float r2 = rng.uniform();

    const float phi = 2.0f * kPi * r1;
    const float sqrt_r2 = std::sqrt(r2);

    const float x = std::cos(phi) * sqrt_r2;
    const float y = std::sin(phi) * sqrt_r2;
    const float z = std::sqrt(1.0f - r2);

    return Vec3(x, y, z);
}
