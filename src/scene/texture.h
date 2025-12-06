#pragma once

#include <memory>

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
    ImageTexture() : data_(nullptr), width_(0), height_(0), channels_(0) {}

    explicit ImageTexture(const std::string& filename) {
        data_ = stbi_load(filename.c_str(), &width_, &height_, &channels_, 0);
        if (!data_) {
            std::cerr << "ERROR: Could not load texture image: " << filename << "\n";
            width_ = height_ = 0;
        }
    }

    ~ImageTexture() {
        if (data_) {
            stbi_image_free(data_);
        }
    }

    // Disable copy
    ImageTexture(const ImageTexture&) = delete;
    ImageTexture& operator=(const ImageTexture&) = delete;

    // Enable move
    ImageTexture(ImageTexture&& other) noexcept
        : data_(other.data_), width_(other.width_),
          height_(other.height_), channels_(other.channels_) {
        other.data_ = nullptr;
        other.width_ = other.height_ = other.channels_ = 0;
    }

    ImageTexture& operator=(ImageTexture&& other) noexcept {
        if (this != &other) {
            if (data_) stbi_image_free(data_);
            data_ = other.data_;
            width_ = other.width_;
            height_ = other.height_;
            channels_ = other.channels_;
            other.data_ = nullptr;
            other.width_ = other.height_ = other.channels_ = 0;
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
        if (!data_ || channels_ < 4) {
            return 1.0f;
        }

        u = clamp_float(u, 0.0f, 1.0f);
        v = 1.0f - clamp_float(v, 0.0f, 1.0f);

        int x = static_cast<int>(u * (width_ - 1));
        int y = static_cast<int>(v * (height_ - 1));

        const std::size_t idx = (static_cast<std::size_t>(y) * width_ + x) * channels_;
        return data_[idx + 3] / 255.0f;
    }

    bool valid() const { return data_ != nullptr; }
    int width() const { return width_; }
    int height() const { return height_; }

private:
    Color get_pixel(int x, int y) const {
        const std::size_t idx = (static_cast<std::size_t>(y) * width_ + x) * channels_;

        float r, g, b;
        if (channels_ == 1 || channels_ == 2) {
            float gray = std::pow(data_[idx + 0] / 255.0f, 2.2f);
            r = g = b = gray;
        } else {
            r = std::pow(data_[idx + 0] / 255.0f, 2.2f);
            g = std::pow(data_[idx + 1] / 255.0f, 2.2f);
            b = std::pow(data_[idx + 2] / 255.0f, 2.2f);
        }

        return Color(r, g, b);
    }

    unsigned char* data_;
    int width_;
    int height_;
    int channels_;
};

class NormalMapTexture {
public:
    NormalMapTexture() : data_(nullptr), width_(0), height_(0), channels_(0) {}

    explicit NormalMapTexture(const std::string& filename) {
        data_ = stbi_load(filename.c_str(), &width_, &height_, &channels_, 3);
        if (!data_) {
            std::cerr << "ERROR: Could not load normal map: " << filename << "\n";
            width_ = height_ = 0;
        }
    }

    ~NormalMapTexture() {
        if (data_) {
            stbi_image_free(data_);
        }
    }

    NormalMapTexture(const NormalMapTexture&) = delete;
    NormalMapTexture& operator=(const NormalMapTexture&) = delete;

    NormalMapTexture(NormalMapTexture&& other) noexcept
        : data_(other.data_), width_(other.width_),
          height_(other.height_), channels_(other.channels_) {
        other.data_ = nullptr;
        other.width_ = other.height_ = other.channels_ = 0;
    }

    Vec3 get_normal(float u, float v) const {
        if (!data_) {
            return Vec3(0.0f, 0.0f, 1.0f);
        }

        u = clamp_float(u, 0.0f, 1.0f);
        v = 1.0f - clamp_float(v, 0.0f, 1.0f);

        int x = static_cast<int>(u * (width_ - 1));
        int y = static_cast<int>(v * (height_ - 1));

        const std::size_t idx = (static_cast<std::size_t>(y) * width_ + x) * 3;

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
    unsigned char* data_;
    int width_;
    int height_;
    int channels_;
};

using NormalMapPtr = std::shared_ptr<NormalMapTexture>;
