#pragma once

#include <memory>
#include <string>
#include <vector>

#include "core/color.h"
#include "core/rng.h"
#include "core/vec3.h"

class EnvironmentMap {
public:
    EnvironmentMap() = default;

    bool load_hdr(const std::string& path, std::string* out_error = nullptr);

    bool valid() const { return width_ > 0 && height_ > 0 && pixels_.size() == static_cast<size_t>(width_ * height_ * 3); }

    void set_intensity(float intensity) { intensity_ = intensity; }
    float intensity() const { return intensity_; }

    // Rotate the environment around +Y (in degrees). Positive values rotate the map in the +phi direction.
    void set_rotation_deg(float rotation_deg) { rotation_y_rad_ = rotation_deg * (kPi / 180.0f); }
    float rotation_deg() const { return rotation_y_rad_ * (180.0f / kPi); }

    Color Le(const Vec3& dir) const;

    Vec3 sample(float& out_pdf, RNG& rng) const;
    float pdf(const Vec3& dir) const;

private:
    void clear();
    void build_distributions();

    Color texel(int x, int y) const;
    Color texel_bilinear(float u, float v) const;

    int width_ = 0;
    int height_ = 0;
    std::vector<float> pixels_;

    float intensity_ = 1.0f;
    float rotation_y_rad_ = 0.0f;

    float total_weight_ = 0.0f;
    std::vector<float> weights_;
    std::vector<float> marginal_cdf_;
    std::vector<float> conditional_cdf_;
};

using EnvironmentMapPtr = std::shared_ptr<EnvironmentMap>;
