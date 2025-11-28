#pragma once

#include <vector>

#include "core/color.h"

class Film {
public:
    Film(int width, int height)
        : width_(width),
          height_(height),
          pixels_(static_cast<std::size_t>(width) * static_cast<std::size_t>(height), Color(0.0f)) {}

    int width() const {
        return width_;
    }

    int height() const {
        return height_;
    }

    void add_sample(int x, int y, const Color& c) {
        if (x < 0 || x >= width_ || y < 0 || y >= height_) {
            return;
        }
        pixels_[index(x, y)] += c;
    }

    Color get_pixel(int x, int y) const {
        if (x < 0 || x >= width_ || y < 0 || y >= height_) {
            return Color(0.0f);
        }
        return pixels_[index(x, y)];
    }

private:
    std::size_t index(int x, int y) const {
        return static_cast<std::size_t>(y) * static_cast<std::size_t>(width_) +
               static_cast<std::size_t>(x);
    }

    int width_;
    int height_;
    std::vector<Color> pixels_;
};

