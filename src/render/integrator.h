#pragma once

#include "core/color.h"
#include "core/ray.h"
#include "core/rng.h"

class Scene;

class Integrator {
public:
    virtual ~Integrator() = default;

    virtual Color Li(const Ray& r,
                     const Scene& scene,
                     RNG& rng,
                     int depth) const = 0;
};

