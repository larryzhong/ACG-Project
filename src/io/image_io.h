#pragma once

#include <cmath>
#include <fstream>
#include <string>

#include "core/color.h"
#include "render/film.h"

inline void write_ppm(const std::string& filename, const Film& film, bool apply_gamma = true) {
    std::ofstream out(filename, std::ios::out | std::ios::binary);
    if (!out) {
        return;
    }

    const int width = film.width();
    const int height = film.height();

    out << "P3\n" << width << ' ' << height << "\n255\n";

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            Color c = film.get_pixel(x, y);

            if (apply_gamma) {
                c.x = std::sqrt(c.x);
                c.y = std::sqrt(c.y);
                c.z = std::sqrt(c.z);
            }

            const float r = clamp_float(c.x, 0.0f, 1.0f);
            const float g = clamp_float(c.y, 0.0f, 1.0f);
            const float b = clamp_float(c.z, 0.0f, 1.0f);

            const int ir = static_cast<int>(255.999f * r);
            const int ig = static_cast<int>(255.999f * g);
            const int ib = static_cast<int>(255.999f * b);

            out << ir << ' ' << ig << ' ' << ib << '\n';
        }
    }
}

