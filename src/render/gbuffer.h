#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

#include "core/color.h"
#include "core/vec3.h"

class GBuffer {
public:
    GBuffer(int width, int height)
        : width_(width),
          height_(height),
          normals_(static_cast<std::size_t>(width) * static_cast<std::size_t>(height), Vec3(0.0f)),
          albedos_(static_cast<std::size_t>(width) * static_cast<std::size_t>(height), Color(0.0f)),
          depths_(static_cast<std::size_t>(width) * static_cast<std::size_t>(height),
                  std::numeric_limits<float>::infinity()),
          valid_(static_cast<std::size_t>(width) * static_cast<std::size_t>(height), 0) {}

    int width() const {
        return width_;
    }

    int height() const {
        return height_;
    }

    void set(int x, int y, const Vec3& n, const Color& a, float d, bool v) {
        if (x < 0 || x >= width_ || y < 0 || y >= height_) {
            return;
        }
        const std::size_t i = index(x, y);
        normals_[i] = n;
        albedos_[i] = a;
        depths_[i] = d;
        valid_[i] = v ? 1u : 0u;
    }

    const Vec3& normal_at(std::size_t idx) const {
        return normals_[idx];
    }

    const Color& albedo_at(std::size_t idx) const {
        return albedos_[idx];
    }

    float depth_at(std::size_t idx) const {
        return depths_[idx];
    }

    std::uint8_t valid_at(std::size_t idx) const {
        return valid_[idx];
    }

private:
    std::size_t index(int x, int y) const {
        return static_cast<std::size_t>(y) * static_cast<std::size_t>(width_) +
               static_cast<std::size_t>(x);
    }

    int width_;
    int height_;
    std::vector<Vec3> normals_;
    std::vector<Color> albedos_;
    std::vector<float> depths_;
    std::vector<std::uint8_t> valid_;
};

