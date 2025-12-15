#pragma once

#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

#include "core/color.h"
#include "core/vec3.h"
#include "io/stb_image.h"

class Texture {
public:
    virtual ~Texture() = default;

    virtual Color value(float u, float v, const Vec3& p) const = 0;

    virtual float alpha(float u, float v, const Vec3& p) const {
        (void)u;
        (void)v;
        (void)p;
        return 1.0f;
    }
};

using TexturePtr = std::shared_ptr<Texture>;

class SolidColor : public Texture {
public:
    SolidColor() : color_(0.0f) {}
    explicit SolidColor(const Color& c) : color_(c) {}

    Color value(float /*u*/, float /*v*/, const Vec3& /*p*/) const override {
        return color_;
    }

private:
    Color color_;
};

class CheckerTexture : public Texture {
public:
    CheckerTexture() = default;

    CheckerTexture(const TexturePtr& even, const TexturePtr& odd, float scale = 1.0f)
        : even_(even), odd_(odd), scale_(scale) {}

    Color value(float u, float v, const Vec3& p) const override {
        const float sx = std::floor(scale_ * p.x);
        const float sy = std::floor(scale_ * p.y);
        const float sz = std::floor(scale_ * p.z);

        const int sum = static_cast<int>(sx + sy + sz);
        if (sum & 1) {
            return odd_ ? odd_->value(u, v, p) : Color(0.0f);
        }
        return even_ ? even_->value(u, v, p) : Color(1.0f);
    }

private:
    TexturePtr even_;
    TexturePtr odd_;
    float scale_ = 1.0f;
};

class AlphaCheckerTexture : public Texture {
public:
    AlphaCheckerTexture() = default;

    explicit AlphaCheckerTexture(float scale) : scale_(scale) {}

    Color value(float /*u*/, float /*v*/, const Vec3& /*p*/) const override {
        return Color(1.0f);
    }

    float alpha(float /*u*/, float /*v*/, const Vec3& p) const override {
        const float sx = std::floor(scale_ * p.x);
        const float sz = std::floor(scale_ * p.z);
        const int sum = static_cast<int>(sx + sz);
        return (sum & 1) ? 0.0f : 1.0f;
    }

private:
    float scale_ = 1.0f;
};

class ImageTexture : public Texture {
public:
    enum class ColorSpace {
        sRGB,
        Linear
    };

    ImageTexture()
        : data_(nullptr),
          width_(0),
          height_(0),
          channels_(0),
          color_space_(ColorSpace::sRGB),
          channel_(-1),
          owns_memory_(false) {}

    explicit ImageTexture(const std::string& filename,
                          ColorSpace color_space = ColorSpace::sRGB,
                          int channel = -1)
        : data_(nullptr),
          width_(0),
          height_(0),
          channels_(0),
          color_space_(color_space),
          channel_(channel),
          owns_memory_(false) {
        data_ = stbi_load(filename.c_str(), &width_, &height_, &channels_, 0);
        if (!data_) {
            std::cerr << "ERROR: Could not load texture image: " << filename << "\n";
            width_ = height_ = 0;
        }
    }

    ImageTexture(std::vector<unsigned char> bytes,
                 int width,
                 int height,
                 int channels,
                 ColorSpace color_space = ColorSpace::sRGB,
                 int channel = -1)
        : owned_data_(std::move(bytes)),
          data_(nullptr),
          width_(width),
          height_(height),
          channels_(channels),
          color_space_(color_space),
          channel_(channel),
          owns_memory_(true) {
        if (width_ <= 0 || height_ <= 0 || channels_ <= 0) {
            width_ = height_ = channels_ = 0;
            owns_memory_ = false;
            owned_data_.clear();
        }
        if (!owned_data_.empty()) {
            data_ = owned_data_.data();
        } else {
            data_ = nullptr;
        }
    }

    ~ImageTexture() {
        if (data_) {
            if (!owns_memory_) {
                stbi_image_free(data_);
            }
        }
    }

    // Disable copy
    ImageTexture(const ImageTexture&) = delete;
    ImageTexture& operator=(const ImageTexture&) = delete;

    // Enable move
    ImageTexture(ImageTexture&& other) noexcept
        : owned_data_(std::move(other.owned_data_)),
          data_(other.data_),
          width_(other.width_),
          height_(other.height_),
          channels_(other.channels_),
          color_space_(other.color_space_),
          channel_(other.channel_),
          owns_memory_(other.owns_memory_) {
        other.data_ = nullptr;
        other.width_ = other.height_ = other.channels_ = 0;
        other.owns_memory_ = false;
    }

    ImageTexture& operator=(ImageTexture&& other) noexcept {
        if (this != &other) {
            if (data_) {
                if (!owns_memory_) {
                    stbi_image_free(data_);
                }
            }
            owned_data_ = std::move(other.owned_data_);
            data_ = other.data_;
            width_ = other.width_;
            height_ = other.height_;
            channels_ = other.channels_;
            color_space_ = other.color_space_;
            channel_ = other.channel_;
            owns_memory_ = other.owns_memory_;
            other.data_ = nullptr;
            other.width_ = other.height_ = other.channels_ = 0;
            other.owns_memory_ = false;
        }
        return *this;
    }

    Color value(float u, float v, const Vec3& /*p*/) const override {
        if (!data_) {
            return Color(1.0f, 0.0f, 1.0f); // Magenta for missing texture
        }

        // Clamp UV to [0, 1]
        u = clamp_float(u, 0.0f, 1.0f);
        v = 1.0f - clamp_float(v, 0.0f, 1.0f); // Flip V (image origin is top-left)

        // Map to pixel coordinates (for bilinear interpolation)
        float fx = u * (width_ - 1);
        float fy = v * (height_ - 1);

        int x0 = static_cast<int>(fx);
        int y0 = static_cast<int>(fy);
        int x1 = std::min(x0 + 1, width_ - 1);
        int y1 = std::min(y0 + 1, height_ - 1);

        float tx = fx - x0;
        float ty = fy - y0;

        // Bilinear interpolation
        Color c00 = get_pixel(x0, y0);
        Color c10 = get_pixel(x1, y0);
        Color c01 = get_pixel(x0, y1);
        Color c11 = get_pixel(x1, y1);

        Color c0 = c00 * (1.0f - tx) + c10 * tx;
        Color c1 = c01 * (1.0f - tx) + c11 * tx;

        return c0 * (1.0f - ty) + c1 * ty;
    }

    float alpha(float u, float v, const Vec3& /*p*/) const override {
        if (!data_) {
            return 1.0f;
        }

        u = clamp_float(u, 0.0f, 1.0f);
        v = 1.0f - clamp_float(v, 0.0f, 1.0f);

        int x = static_cast<int>(u * (width_ - 1));
        int y = static_cast<int>(v * (height_ - 1));

        const std::size_t idx = (static_cast<std::size_t>(y) * width_ + x) * channels_;
        if (channels_ >= 4) {
            return data_[idx + 3] / 255.0f;
        }
        if (channels_ == 2) {
            return data_[idx + 1] / 255.0f;
        }
        return 1.0f;
    }

    bool valid() const { return data_ != nullptr; }
    int width() const { return width_; }
    int height() const { return height_; }

private:
    Color get_pixel(int x, int y) const {
        const std::size_t idx = (static_cast<std::size_t>(y) * width_ + x) * channels_;

        if (!data_ || width_ <= 0 || height_ <= 0 || channels_ <= 0) {
            return Color(1.0f, 0.0f, 1.0f);
        }

        auto decode = [&](int c, bool apply_srgb) -> float {
            if (c < 0 || c >= channels_) {
                return 0.0f;
            }
            const float v = data_[idx + c] / 255.0f;
            if (apply_srgb) {
                return std::pow(v, 2.2f);
            }
            return v;
        };

        if (channel_ >= 0) {
            const bool apply_srgb = (color_space_ == ColorSpace::sRGB) && (channel_ != 3);
            const float v = decode(channel_, apply_srgb);
            return Color(v, v, v);
        }

        if (channels_ == 1 || channels_ == 2) {
            const bool apply_srgb = (color_space_ == ColorSpace::sRGB);
            const float gray = decode(0, apply_srgb);
            return Color(gray, gray, gray);
        }

        const bool apply_srgb = (color_space_ == ColorSpace::sRGB);
        const float r = decode(0, apply_srgb);
        const float g = decode(1, apply_srgb);
        const float b = decode(2, apply_srgb);
        return Color(r, g, b);
    }

    std::vector<unsigned char> owned_data_;
    unsigned char* data_;
    int width_;
    int height_;
    int channels_;
    ColorSpace color_space_;
    int channel_;
    bool owns_memory_;
};

class NormalMapTexture {
public:
    NormalMapTexture()
        : owned_data_(),
          data_(nullptr),
          width_(0),
          height_(0),
          channels_(0),
          owns_memory_(false) {}

    explicit NormalMapTexture(const std::string& filename) {
        int channels_in_file = 0;
        data_ = stbi_load(filename.c_str(), &width_, &height_, &channels_in_file, 3);
        owns_memory_ = false;
        if (!data_) {
            std::cerr << "ERROR: Could not load normal map: " << filename << "\n";
            width_ = height_ = 0;
            channels_ = 0;
            return;
        }
        channels_ = 3;
    }

    NormalMapTexture(std::vector<unsigned char> bytes,
                     int width,
                     int height,
                     int channels)
        : owned_data_(std::move(bytes)),
          data_(nullptr),
          width_(width),
          height_(height),
          channels_(channels),
          owns_memory_(true) {
        if (width_ <= 0 || height_ <= 0 || channels_ <= 0) {
            width_ = height_ = channels_ = 0;
            owns_memory_ = false;
            owned_data_.clear();
        }
        if (!owned_data_.empty()) {
            data_ = owned_data_.data();
        } else {
            data_ = nullptr;
        }
    }

    ~NormalMapTexture() {
        if (data_) {
            if (!owns_memory_) {
                stbi_image_free(data_);
            }
        }
    }

    NormalMapTexture(const NormalMapTexture&) = delete;
    NormalMapTexture& operator=(const NormalMapTexture&) = delete;

    NormalMapTexture(NormalMapTexture&& other) noexcept
        : owned_data_(std::move(other.owned_data_)),
          data_(other.data_),
          width_(other.width_),
          height_(other.height_),
          channels_(other.channels_),
          owns_memory_(other.owns_memory_) {
        other.data_ = nullptr;
        other.width_ = other.height_ = other.channels_ = 0;
        other.owns_memory_ = false;
    }

    Vec3 get_normal(float u, float v) const {
        if (!data_) {
            return Vec3(0.0f, 0.0f, 1.0f);
        }

        u = clamp_float(u, 0.0f, 1.0f);
        v = 1.0f - clamp_float(v, 0.0f, 1.0f);

        int x = static_cast<int>(u * (width_ - 1));
        int y = static_cast<int>(v * (height_ - 1));

        const int stride = (channels_ >= 3) ? channels_ : 3;
        const std::size_t idx = (static_cast<std::size_t>(y) * width_ + x) * static_cast<std::size_t>(stride);

        float nx = data_[idx + 0] / 255.0f * 2.0f - 1.0f;
        float ny = data_[idx + 1] / 255.0f * 2.0f - 1.0f;
        float nz = data_[idx + 2] / 255.0f * 2.0f - 1.0f;

        return normalize(Vec3(nx, ny, nz));
    }

    static Vec3 apply_normal_map(const Vec3& tangent_normal,
                                  const Vec3& normal,
                                  const Vec3& tangent) {
        Vec3 N = normalize(normal);
        Vec3 T = normalize(tangent - dot(tangent, N) * N);
        Vec3 B = cross(N, T);

        return normalize(
            tangent_normal.x * T +
            tangent_normal.y * B +
            tangent_normal.z * N
        );
    }

    bool valid() const { return data_ != nullptr; }

private:
    std::vector<unsigned char> owned_data_;
    unsigned char* data_;
    int width_;
    int height_;
    int channels_;
    bool owns_memory_;
};

using NormalMapPtr = std::shared_ptr<NormalMapTexture>;

class ColorFactorTexture : public Texture {
public:
    ColorFactorTexture(const TexturePtr& base, const Color& color_factor, float alpha_factor = 1.0f)
        : base_(base), color_factor_(color_factor), alpha_factor_(alpha_factor) {}

    Color value(float u, float v, const Vec3& p) const override {
        Color c = base_ ? base_->value(u, v, p) : Color(1.0f);
        return Color(c.x * color_factor_.x, c.y * color_factor_.y, c.z * color_factor_.z);
    }

    float alpha(float u, float v, const Vec3& p) const override {
        const float a = base_ ? base_->alpha(u, v, p) : 1.0f;
        return clamp_float(a * alpha_factor_, 0.0f, 1.0f);
    }

private:
    TexturePtr base_;
    Color color_factor_;
    float alpha_factor_ = 1.0f;
};

class AlphaCutoffTexture : public Texture {
public:
    AlphaCutoffTexture(const TexturePtr& base, float cutoff)
        : base_(base), cutoff_(cutoff) {}

    Color value(float u, float v, const Vec3& p) const override {
        return base_ ? base_->value(u, v, p) : Color(1.0f);
    }

    float alpha(float u, float v, const Vec3& p) const override {
        const float a = base_ ? base_->alpha(u, v, p) : 1.0f;
        return (a >= cutoff_) ? 1.0f : 0.0f;
    }

private:
    TexturePtr base_;
    float cutoff_ = 0.5f;
};
