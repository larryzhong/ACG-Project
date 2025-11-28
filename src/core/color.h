#pragma once

#include <algorithm>
#include "core/vec3.h"

using Color = Vec3;

inline float clamp_float(float x, float min_val, float max_val) {
    return std::max(min_val, std::min(max_val, x));
}

