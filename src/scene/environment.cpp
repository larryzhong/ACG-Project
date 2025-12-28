#include "scene/environment.h"

#include <algorithm>
#include <cmath>
#include <cstring>

#include "core/color.h"

#include "io/stb_image.h"

namespace {
inline float luminance(const Color& c) {
    return 0.2126f * c.x + 0.7152f * c.y + 0.0722f * c.z;
}

inline float inv_two_pi() {
    return 1.0f / (2.0f * kPi);
}

inline int clamp_int(int v, int lo, int hi) {
    return std::max(lo, std::min(v, hi));
}

inline float wrap01(float x) {
    x -= std::floor(x);
    if (x < 0.0f) x += 1.0f;
    if (x >= 1.0f) x -= 1.0f;
    return x;
}

struct DiscreteSample {
    int index = 0;
    float du = 0.0f;
};

inline DiscreteSample sample_cdf(const float* cdf, int n, float u) {
    const float uu = clamp_float(u, 0.0f, 0.99999994f);
    const float* begin = cdf;
    const float* end = cdf + (n + 1);
    const float* it = std::upper_bound(begin + 1, end, uu);
    int idx = static_cast<int>((it - begin) - 1);
    idx = clamp_int(idx, 0, n - 1);
    const float c0 = cdf[idx];
    const float c1 = cdf[idx + 1];
    float du = 0.0f;
    if (c1 > c0) {
        du = (uu - c0) / (c1 - c0);
    }
    return {idx, clamp_float(du, 0.0f, 0.99999994f)};
}

inline Vec3 uniform_sphere(RNG& rng) {
    const float u1 = rng.uniform();
    const float u2 = rng.uniform();
    const float y = 1.0f - 2.0f * u1;
    const float r = std::sqrt(std::max(0.0f, 1.0f - y * y));
    const float phi = 2.0f * kPi * u2;
    return Vec3(r * std::cos(phi), y, r * std::sin(phi));
}
}  // namespace

void EnvironmentMap::clear() {
    width_ = 0;
    height_ = 0;
    pixels_.clear();
    total_weight_ = 0.0f;
    weights_.clear();
    marginal_cdf_.clear();
    conditional_cdf_.clear();
}

bool EnvironmentMap::load_hdr(const std::string& path, std::string* out_error) {
    clear();
    if (out_error) out_error->clear();

    int w = 0;
    int h = 0;
    int comp = 0;
    float* data = stbi_loadf(path.c_str(), &w, &h, &comp, 3);
    if (!data || w <= 0 || h <= 0) {
        if (out_error) {
            *out_error = "Failed to load HDR environment: " + path;
        }
        if (data) {
            stbi_image_free(data);
        }
        return false;
    }

    width_ = w;
    height_ = h;
    pixels_.assign(data, data + static_cast<size_t>(w * h * 3));
    stbi_image_free(data);

    build_distributions();
    return valid();
}

void EnvironmentMap::build_distributions() {
    total_weight_ = 0.0f;
    weights_.assign(static_cast<size_t>(width_ * height_), 0.0f);
    marginal_cdf_.assign(static_cast<size_t>(height_ + 1), 0.0f);
    conditional_cdf_.assign(static_cast<size_t>(height_ * (width_ + 1)), 0.0f);

    if (!valid()) {
        total_weight_ = 0.0f;
        return;
    }

    std::vector<float> row_sum(static_cast<size_t>(height_), 0.0f);

    for (int y = 0; y < height_; ++y) {
        const float theta = kPi * (static_cast<float>(y) + 0.5f) / static_cast<float>(height_);
        const float sin_theta = std::sin(theta);
        float sum = 0.0f;

        float* cdf = conditional_cdf_.data() + static_cast<size_t>(y) * static_cast<size_t>(width_ + 1);
        cdf[0] = 0.0f;

        for (int x = 0; x < width_; ++x) {
            const Color c = texel(x, y);
            float w = luminance(c);
            if (!std::isfinite(w) || w < 0.0f) {
                w = 0.0f;
            }
            w *= sin_theta;
            weights_[static_cast<size_t>(y) * static_cast<size_t>(width_) + static_cast<size_t>(x)] = w;
            sum += w;
            cdf[x + 1] = sum;
        }

        row_sum[static_cast<size_t>(y)] = sum;
        if (sum > 0.0f) {
            const float inv = 1.0f / sum;
            for (int i = 1; i <= width_; ++i) {
                cdf[i] *= inv;
            }
        } else {
            for (int i = 1; i <= width_; ++i) {
                cdf[i] = static_cast<float>(i) / static_cast<float>(width_);
            }
        }
    }

    marginal_cdf_[0] = 0.0f;
    for (int y = 0; y < height_; ++y) {
        total_weight_ += row_sum[static_cast<size_t>(y)];
        marginal_cdf_[static_cast<size_t>(y + 1)] = total_weight_;
    }

    if (total_weight_ > 0.0f) {
        const float inv = 1.0f / total_weight_;
        for (int y = 1; y <= height_; ++y) {
            marginal_cdf_[static_cast<size_t>(y)] *= inv;
        }
    } else {
        for (int y = 1; y <= height_; ++y) {
            marginal_cdf_[static_cast<size_t>(y)] = static_cast<float>(y) / static_cast<float>(height_);
        }
    }
}

Color EnvironmentMap::texel(int x, int y) const {
    const int ix = clamp_int(x, 0, width_ - 1);
    const int iy = clamp_int(y, 0, height_ - 1);
    const size_t idx = (static_cast<size_t>(iy) * static_cast<size_t>(width_) + static_cast<size_t>(ix)) * 3u;
    return Color(pixels_[idx], pixels_[idx + 1], pixels_[idx + 2]);
}

Color EnvironmentMap::texel_bilinear(float u, float v) const {
    if (!valid()) {
        return Color(0.0f);
    }

    const float uu = wrap01(u);
    const float vv = clamp_float(v, 0.0f, 1.0f);

    const float x = uu * static_cast<float>(width_) - 0.5f;
    const float y = vv * static_cast<float>(height_) - 0.5f;

    const int x0 = static_cast<int>(std::floor(x));
    const int y0 = static_cast<int>(std::floor(y));
    const float tx = x - static_cast<float>(x0);
    const float ty = y - static_cast<float>(y0);

    auto wrap_x = [&](int ix) {
        int r = ix % width_;
        if (r < 0) r += width_;
        return r;
    };

    const int ix0 = wrap_x(x0);
    const int ix1 = wrap_x(x0 + 1);
    const int iy0 = clamp_int(y0, 0, height_ - 1);
    const int iy1 = clamp_int(y0 + 1, 0, height_ - 1);

    const Color c00 = texel(ix0, iy0);
    const Color c10 = texel(ix1, iy0);
    const Color c01 = texel(ix0, iy1);
    const Color c11 = texel(ix1, iy1);

    const Color c0 = (1.0f - tx) * c00 + tx * c10;
    const Color c1 = (1.0f - tx) * c01 + tx * c11;
    return (1.0f - ty) * c0 + ty * c1;
}

Color EnvironmentMap::Le(const Vec3& dir) const {
    if (!valid()) {
        return Color(0.0f);
    }

    const Vec3 d = normalize(dir);
    const float y = clamp_float(d.y, -1.0f, 1.0f);
    const float theta = std::acos(y);
    const float phi = std::atan2(d.z, d.x) - rotation_y_rad_;
    const float u = wrap01(phi * inv_two_pi());
    const float v = theta * (1.0f / kPi);
    return intensity_ * texel_bilinear(u, v);
}

Vec3 EnvironmentMap::sample(float& out_pdf, RNG& rng) const {
    if (!valid() || marginal_cdf_.empty() || conditional_cdf_.empty()) {
        out_pdf = 1.0f / (4.0f * kPi);
        return uniform_sphere(rng);
    }

    const DiscreteSample row = sample_cdf(marginal_cdf_.data(), height_, rng.uniform());
    const float* row_cdf = conditional_cdf_.data() + static_cast<size_t>(row.index) * static_cast<size_t>(width_ + 1);
    const DiscreteSample col = sample_cdf(row_cdf, width_, rng.uniform());

    const float u = (static_cast<float>(col.index) + col.du) / static_cast<float>(width_);
    const float v = (static_cast<float>(row.index) + row.du) / static_cast<float>(height_);

    const float phi = 2.0f * kPi * wrap01(u + rotation_y_rad_ * inv_two_pi());
    const float theta = kPi * v;

    const float sin_theta = std::sin(theta);
    const Vec3 wi(sin_theta * std::cos(phi), std::cos(theta), sin_theta * std::sin(phi));

    out_pdf = pdf(wi);
    return wi;
}

float EnvironmentMap::pdf(const Vec3& dir) const {
    if (!valid() || total_weight_ <= 0.0f) {
        return 1.0f / (4.0f * kPi);
    }

    const Vec3 d = normalize(dir);
    const float y = clamp_float(d.y, -1.0f, 1.0f);
    const float theta = std::acos(y);
    const float phi = std::atan2(d.z, d.x) - rotation_y_rad_;
    const float u = wrap01(phi * inv_two_pi());
    const float v = theta * (1.0f / kPi);

    const int ix = clamp_int(static_cast<int>(u * static_cast<float>(width_)), 0, width_ - 1);
    const int iy = clamp_int(static_cast<int>(v * static_cast<float>(height_)), 0, height_ - 1);

    const float theta_center = kPi * (static_cast<float>(iy) + 0.5f) / static_cast<float>(height_);
    const float sin_theta = std::sin(theta_center);
    if (sin_theta <= 0.0f) {
        return 0.0f;
    }

    const float w = weights_[static_cast<size_t>(iy) * static_cast<size_t>(width_) + static_cast<size_t>(ix)];
    if (w <= 0.0f) {
        return 0.0f;
    }

    const float delta_phi = 2.0f * kPi / static_cast<float>(width_);
    const float delta_theta = kPi / static_cast<float>(height_);
    const float pixel_solid_angle = delta_phi * delta_theta * sin_theta;
    if (pixel_solid_angle <= 0.0f) {
        return 0.0f;
    }

    const float p_pixel = w / total_weight_;
    return p_pixel / pixel_solid_angle;
}
