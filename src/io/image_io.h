#pragma once

#include <cmath>
#include <fstream>
#include <string>

#include "core/color.h"
#include "render/film.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "io/stb_image_write.h"

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

inline void write_png(const std::string& filename, const Film& film, bool apply_gamma = true) {
    const int width = film.width();
    const int height = film.height();
    const int channels = 3;

    std::vector<unsigned char> data(static_cast<size_t>(width * height * channels));

    int index = 0;
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

            data[index++] = static_cast<unsigned char>(255.999f * r);
            data[index++] = static_cast<unsigned char>(255.999f * g);
            data[index++] = static_cast<unsigned char>(255.999f * b);
        }
    }

    stbi_write_png(filename.c_str(), width, height, channels, data.data(), width * channels);
}
