#pragma once

#include <cmath>
#include <fstream>
#include <string>

#include "core/color.h"
#include "render/film.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "io/stb_image_write.h"

inline float linear_to_srgb(float x) {
    x = std::max(0.0f, x);
    if (x <= 0.0031308f) {
        return 12.92f * x;
    }
    return 1.055f * std::pow(x, 1.0f / 2.4f) - 0.055f;
}

inline Color tonemap_aces(Color x) {
    // Narkowicz 2015, "ACES Filmic Tone Mapping Curve"
    x.x = std::max(0.0f, x.x);
    x.y = std::max(0.0f, x.y);
    x.z = std::max(0.0f, x.z);

    const float a = 2.51f;
    const float b = 0.03f;
    const float c = 2.43f;
    const float d = 0.59f;
    const float e = 0.14f;

    auto map = [&](float v) {
        return (v * (a * v + b)) / (v * (c * v + d) + e);
    };

    return Color(map(x.x), map(x.y), map(x.z));
}

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
                c = tonemap_aces(c);
                c.x = linear_to_srgb(c.x);
                c.y = linear_to_srgb(c.y);
                c.z = linear_to_srgb(c.z);
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
                c = tonemap_aces(c);
                c.x = linear_to_srgb(c.x);
                c.y = linear_to_srgb(c.y);
                c.z = linear_to_srgb(c.z);
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
