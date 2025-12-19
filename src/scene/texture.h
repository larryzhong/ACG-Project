#pragma once

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <utility>
#include <vector>

#include "core/color.h"
#include "core/vec3.h"
#include "io/stb_image.h"
#include "scene/hittable.h"

class Texture {
public:
    virtual ~Texture() = default;

    virtual Color value(const HitRecord& hit) const = 0;

    virtual float alpha(const HitRecord& /*hit*/) const {
        return 1.0f;
    }
};

using TexturePtr = std::shared_ptr<Texture>;

class SolidColor : public Texture {
public:
    SolidColor() : color_(0.0f) {}
    explicit SolidColor(const Color& c) : color_(c) {}

    Color value(const HitRecord& /*hit*/) const override {
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

    Color value(const HitRecord& hit) const override {
        const Vec3& p = hit.point;
        const float sx = std::floor(scale_ * p.x);
        const float sy = std::floor(scale_ * p.y);
        const float sz = std::floor(scale_ * p.z);

        const int sum = static_cast<int>(sx + sy + sz);
        if (sum & 1) {
            return odd_ ? odd_->value(hit) : Color(0.0f);
        }
        return even_ ? even_->value(hit) : Color(1.0f);
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

    Color value(const HitRecord& /*hit*/) const override {
        return Color(1.0f);
    }

    float alpha(const HitRecord& hit) const override {
        const Vec3& p = hit.point;
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
          flip_v_(true),
          owns_memory_(false) {}

    explicit ImageTexture(const std::string& filename,
                          ColorSpace color_space = ColorSpace::sRGB,
                          int channel = -1,
                          bool flip_v = true)
        : data_(nullptr),
          width_(0),
          height_(0),
          channels_(0),
          color_space_(color_space),
          channel_(channel),
          flip_v_(flip_v),
          owns_memory_(false) {
        data_ = stbi_load(filename.c_str(), &width_, &height_, &channels_, 0);
        if (!data_) {
            std::cerr << "ERROR: Could not load texture image: " << filename << "\n";
            width_ = height_ = 0;
            channels_ = 0;
        } else {
            build_mipmaps();
        }
    }

    ImageTexture(std::vector<unsigned char> bytes,
                 int width,
                 int height,
                 int channels,
                 ColorSpace color_space = ColorSpace::sRGB,
                 int channel = -1,
                 bool flip_v = true)
        : owned_data_(std::move(bytes)),
          data_(nullptr),
          width_(width),
          height_(height),
          channels_(channels),
          color_space_(color_space),
          channel_(channel),
          flip_v_(flip_v),
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
        if (data_) {
            build_mipmaps();
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
          flip_v_(other.flip_v_),
          mip_levels_(std::move(other.mip_levels_)),
          owns_memory_(other.owns_memory_) {
        other.data_ = nullptr;
        other.width_ = other.height_ = other.channels_ = 0;
        other.mip_levels_.clear();
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
            flip_v_ = other.flip_v_;
            mip_levels_ = std::move(other.mip_levels_);
            owns_memory_ = other.owns_memory_;
            other.data_ = nullptr;
            other.width_ = other.height_ = other.channels_ = 0;
            other.mip_levels_.clear();
            other.owns_memory_ = false;
        }
        return *this;
    }

    Color value(const HitRecord& hit) const override {
        if (mip_levels_.empty()) {
            return Color(1.0f, 0.0f, 1.0f); // Magenta for missing texture
        }

        const float u = hit.u;
        const float v = hit.v;
        const float lod = compute_lod(hit);
        return sample_trilinear(u, v, lod);
    }

    float alpha(const HitRecord& hit) const override {
        if (mip_levels_.empty()) {
            return 1.0f;
        }

        const float u = hit.u;
        const float v = hit.v;
        const float lod = compute_lod(hit);
        return sample_alpha_trilinear(u, v, lod);
    }

    bool valid() const { return !mip_levels_.empty(); }
    int width() const { return width_; }
    int height() const { return height_; }

private:
    struct MipLevel {
        int width = 0;
        int height = 0;
        std::vector<unsigned char> rgba; // Always 4 channels, row-major.

        bool valid() const {
            return width > 0 && height > 0 &&
                   rgba.size() == static_cast<std::size_t>(width) * static_cast<std::size_t>(height) * 4u;
        }
    };

    void release_source_data() {
        if (data_ && !owns_memory_) {
            stbi_image_free(data_);
        }
        owned_data_.clear();
        owned_data_.shrink_to_fit();
        data_ = nullptr;
        owns_memory_ = false;
        channels_ = 0;
    }

    static unsigned char float_to_byte(float v) {
        v = clamp_float(v, 0.0f, 1.0f);
        return static_cast<unsigned char>(v * 255.0f + 0.5f);
    }

    static float byte_to_float(unsigned char b) {
        return static_cast<float>(b) / 255.0f;
    }

    static float linear_to_srgb(float v) {
        v = clamp_float(v, 0.0f, 1.0f);
        return std::pow(v, 1.0f / 2.2f);
    }

    static float srgb_to_linear(float v) {
        v = clamp_float(v, 0.0f, 1.0f);
        return std::pow(v, 2.2f);
    }

    unsigned char source_channel(int pixel_index, int c) const {
        if (!data_ || channels_ <= 0 || width_ <= 0 || height_ <= 0) {
            return 0;
        }
        if (c < 0) {
            return 0;
        }
        if (c >= channels_) {
            return 0;
        }
        return data_[static_cast<std::size_t>(pixel_index) * static_cast<std::size_t>(channels_) + static_cast<std::size_t>(c)];
    }

    void build_mipmaps() {
        mip_levels_.clear();
        if (!data_ || width_ <= 0 || height_ <= 0 || channels_ <= 0) {
            return;
        }

        MipLevel base;
        base.width = width_;
        base.height = height_;
        base.rgba.resize(static_cast<std::size_t>(width_) * static_cast<std::size_t>(height_) * 4u);

        for (int y = 0; y < height_; ++y) {
            for (int x = 0; x < width_; ++x) {
                const int pixel_index = y * width_ + x;
                unsigned char r = 0, g = 0, b = 0, a = 255;

                if (channels_ == 1) {
                    r = g = b = source_channel(pixel_index, 0);
                } else if (channels_ == 2) {
                    r = g = b = source_channel(pixel_index, 0);
                    a = source_channel(pixel_index, 1);
                } else if (channels_ == 3) {
                    r = source_channel(pixel_index, 0);
                    g = source_channel(pixel_index, 1);
                    b = source_channel(pixel_index, 2);
                } else {
                    r = source_channel(pixel_index, 0);
                    g = source_channel(pixel_index, 1);
                    b = source_channel(pixel_index, 2);
                    a = source_channel(pixel_index, 3);
                }

                const std::size_t dst = (static_cast<std::size_t>(pixel_index) * 4u);
                base.rgba[dst + 0] = r;
                base.rgba[dst + 1] = g;
                base.rgba[dst + 2] = b;
                base.rgba[dst + 3] = a;
            }
        }

        mip_levels_.push_back(std::move(base));

        while (true) {
            const MipLevel& src = mip_levels_.back();
            if (src.width == 1 && src.height == 1) {
                break;
            }

            MipLevel dst;
            dst.width = std::max(1, src.width / 2);
            dst.height = std::max(1, src.height / 2);
            dst.rgba.resize(static_cast<std::size_t>(dst.width) * static_cast<std::size_t>(dst.height) * 4u);

            for (int y = 0; y < dst.height; ++y) {
                for (int x = 0; x < dst.width; ++x) {
                    float accum[4] = {0.0f, 0.0f, 0.0f, 0.0f};
                    float wsum = 0.0f;

                    for (int oy = 0; oy < 2; ++oy) {
                        for (int ox = 0; ox < 2; ++ox) {
                            const int sx = std::min(src.width - 1, 2 * x + ox);
                            const int sy = std::min(src.height - 1, 2 * y + oy);
                            const std::size_t sidx =
                                (static_cast<std::size_t>(sy) * static_cast<std::size_t>(src.width) +
                                 static_cast<std::size_t>(sx)) *
                                4u;

                            float r = byte_to_float(src.rgba[sidx + 0]);
                            float g = byte_to_float(src.rgba[sidx + 1]);
                            float b = byte_to_float(src.rgba[sidx + 2]);
                            float a = byte_to_float(src.rgba[sidx + 3]);

                            if (color_space_ == ColorSpace::sRGB) {
                                r = srgb_to_linear(r);
                                g = srgb_to_linear(g);
                                b = srgb_to_linear(b);
                            }

                            accum[0] += r;
                            accum[1] += g;
                            accum[2] += b;
                            accum[3] += a;
                            wsum += 1.0f;
                        }
                    }

                    if (wsum <= 0.0f) {
                        wsum = 1.0f;
                    }

                    float r = accum[0] / wsum;
                    float g = accum[1] / wsum;
                    float b = accum[2] / wsum;
                    float a = accum[3] / wsum;

                    if (color_space_ == ColorSpace::sRGB) {
                        r = linear_to_srgb(r);
                        g = linear_to_srgb(g);
                        b = linear_to_srgb(b);
                    }

                    const std::size_t didx =
                        (static_cast<std::size_t>(y) * static_cast<std::size_t>(dst.width) +
                         static_cast<std::size_t>(x)) *
                        4u;
                    dst.rgba[didx + 0] = float_to_byte(r);
                    dst.rgba[didx + 1] = float_to_byte(g);
                    dst.rgba[didx + 2] = float_to_byte(b);
                    dst.rgba[didx + 3] = float_to_byte(a);
                }
            }

            mip_levels_.push_back(std::move(dst));
        }

        release_source_data();
    }

    Color decode_pixel(const MipLevel& level, int x, int y) const {
        if (!level.valid()) {
            return Color(1.0f, 0.0f, 1.0f);
        }

        x = std::min(std::max(x, 0), level.width - 1);
        y = std::min(std::max(y, 0), level.height - 1);

        const std::size_t idx =
            (static_cast<std::size_t>(y) * static_cast<std::size_t>(level.width) +
             static_cast<std::size_t>(x)) *
            4u;

        auto decode_channel = [&](int c, bool apply_srgb) -> float {
            if (c < 0 || c >= 4) {
                return 0.0f;
            }
            float v = byte_to_float(level.rgba[idx + static_cast<std::size_t>(c)]);
            if (apply_srgb) {
                v = srgb_to_linear(v);
            }
            return v;
        };

        if (channel_ >= 0) {
            const bool apply_srgb = (color_space_ == ColorSpace::sRGB) && (channel_ != 3);
            const float v = decode_channel(channel_, apply_srgb);
            return Color(v, v, v);
        }

        const bool apply_srgb = (color_space_ == ColorSpace::sRGB);
        const float r = decode_channel(0, apply_srgb);
        const float g = decode_channel(1, apply_srgb);
        const float b = decode_channel(2, apply_srgb);
        return Color(r, g, b);
    }

    float decode_alpha(const MipLevel& level, int x, int y) const {
        if (!level.valid()) {
            return 1.0f;
        }
        x = std::min(std::max(x, 0), level.width - 1);
        y = std::min(std::max(y, 0), level.height - 1);
        const std::size_t idx =
            (static_cast<std::size_t>(y) * static_cast<std::size_t>(level.width) +
             static_cast<std::size_t>(x)) *
            4u;
        return byte_to_float(level.rgba[idx + 3u]);
    }

    Color sample_bilinear(const MipLevel& level, float u, float v) const {
        if (!level.valid()) {
            return Color(1.0f, 0.0f, 1.0f);
        }

        u = u - std::floor(u);
        if (flip_v_) v = 1.0f - v;
        v = v - std::floor(v);

        const float fx = u * static_cast<float>(level.width - 1);
        const float fy = v * static_cast<float>(level.height - 1);

        const int x0 = static_cast<int>(fx);
        const int y0 = static_cast<int>(fy);
        const int x1 = std::min(x0 + 1, level.width - 1);
        const int y1 = std::min(y0 + 1, level.height - 1);

        const float tx = fx - static_cast<float>(x0);
        const float ty = fy - static_cast<float>(y0);

        const Color c00 = decode_pixel(level, x0, y0);
        const Color c10 = decode_pixel(level, x1, y0);
        const Color c01 = decode_pixel(level, x0, y1);
        const Color c11 = decode_pixel(level, x1, y1);

        const Color c0 = c00 * (1.0f - tx) + c10 * tx;
        const Color c1 = c01 * (1.0f - tx) + c11 * tx;
        return c0 * (1.0f - ty) + c1 * ty;
    }

    float sample_alpha_bilinear(const MipLevel& level, float u, float v) const {
        if (!level.valid()) {
            return 1.0f;
        }

        u = u - std::floor(u);
        if (flip_v_) v = 1.0f - v;
        v = v - std::floor(v);

        const float fx = u * static_cast<float>(level.width - 1);
        const float fy = v * static_cast<float>(level.height - 1);

        const int x0 = static_cast<int>(fx);
        const int y0 = static_cast<int>(fy);
        const int x1 = std::min(x0 + 1, level.width - 1);
        const int y1 = std::min(y0 + 1, level.height - 1);

        const float tx = fx - static_cast<float>(x0);
        const float ty = fy - static_cast<float>(y0);

        const float a00 = decode_alpha(level, x0, y0);
        const float a10 = decode_alpha(level, x1, y0);
        const float a01 = decode_alpha(level, x0, y1);
        const float a11 = decode_alpha(level, x1, y1);

        const float a0 = a00 * (1.0f - tx) + a10 * tx;
        const float a1 = a01 * (1.0f - tx) + a11 * tx;
        return a0 * (1.0f - ty) + a1 * ty;
    }

    Color sample_trilinear(float u, float v, float lod) const {
        if (mip_levels_.empty()) {
            return Color(1.0f, 0.0f, 1.0f);
        }

        const int max_level = static_cast<int>(mip_levels_.size()) - 1;
        if (!(lod >= 0.0f)) {
            lod = 0.0f;
        }
        lod = clamp_float(lod, 0.0f, static_cast<float>(max_level));

        const int l0 = static_cast<int>(std::floor(lod));
        const int l1 = std::min(l0 + 1, max_level);
        const float t = lod - static_cast<float>(l0);

        const Color c0 = sample_bilinear(mip_levels_[static_cast<std::size_t>(l0)], u, v);
        const Color c1 = sample_bilinear(mip_levels_[static_cast<std::size_t>(l1)], u, v);
        return c0 * (1.0f - t) + c1 * t;
    }

    float sample_alpha_trilinear(float u, float v, float lod) const {
        if (mip_levels_.empty()) {
            return 1.0f;
        }

        const int max_level = static_cast<int>(mip_levels_.size()) - 1;
        if (!(lod >= 0.0f)) {
            lod = 0.0f;
        }
        lod = clamp_float(lod, 0.0f, static_cast<float>(max_level));

        const int l0 = static_cast<int>(std::floor(lod));
        const int l1 = std::min(l0 + 1, max_level);
        const float t = lod - static_cast<float>(l0);

        const float a0 = sample_alpha_bilinear(mip_levels_[static_cast<std::size_t>(l0)], u, v);
        const float a1 = sample_alpha_bilinear(mip_levels_[static_cast<std::size_t>(l1)], u, v);
        return a0 * (1.0f - t) + a1 * t;
    }

    float compute_lod(const HitRecord& hit) const {
        if (mip_levels_.size() <= 1) {
            return 0.0f;
        }

        const float footprint = hit.ray_footprint;
        if (!(footprint > 0.0f)) {
            return 0.0f;
        }

        const float dpdu_len = hit.dpdu.length();
        const float dpdv_len = hit.dpdv.length();

        if (!(dpdu_len > 1e-8f) || !(dpdv_len > 1e-8f)) {
            return 0.0f;
        }

        const float u_radius = footprint / dpdu_len;
        const float v_radius = footprint / dpdv_len;
        const float texels =
            std::max(u_radius * static_cast<float>(width_), v_radius * static_cast<float>(height_));

        if (!(texels > 1.0f)) {
            return 0.0f;
        }
        return std::log2(texels);
    }

    std::vector<unsigned char> owned_data_;
    unsigned char* data_;
    int width_;
    int height_;
    int channels_;
    ColorSpace color_space_;
    int channel_;
    bool flip_v_;
    std::vector<MipLevel> mip_levels_;
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
          flip_v_(true),
          owns_memory_(false) {}

    explicit NormalMapTexture(const std::string& filename, bool flip_v = true) {
        int channels_in_file = 0;
        data_ = stbi_load(filename.c_str(), &width_, &height_, &channels_in_file, 3);
        flip_v_ = flip_v;
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
                     int channels,
                     bool flip_v = true)
        : owned_data_(std::move(bytes)),
          data_(nullptr),
          width_(width),
          height_(height),
          channels_(channels),
          flip_v_(flip_v),
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
          flip_v_(other.flip_v_),
          owns_memory_(other.owns_memory_) {
        other.data_ = nullptr;
        other.width_ = other.height_ = other.channels_ = 0;
        other.owns_memory_ = false;
    }

    Vec3 get_normal(float u, float v) const {
        if (!data_) {
            return Vec3(0.0f, 0.0f, 1.0f);
        }

        u = u - std::floor(u);
        if (flip_v_) v = 1.0f - v;
        v = v - std::floor(v);

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
                                 const Vec3& tangent,
                                 float tangent_sign = 1.0f) {
        Vec3 N = normalize(normal);
        Vec3 T = normalize(tangent - dot(tangent, N) * N);
        Vec3 B = tangent_sign * cross(N, T);

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
    bool flip_v_;
    bool owns_memory_;
};

using NormalMapPtr = std::shared_ptr<NormalMapTexture>;

class ColorFactorTexture : public Texture {
public:
    ColorFactorTexture(const TexturePtr& base, const Color& color_factor, float alpha_factor = 1.0f)
        : base_(base), color_factor_(color_factor), alpha_factor_(alpha_factor) {}

    Color value(const HitRecord& hit) const override {
        Color c = base_ ? base_->value(hit) : Color(1.0f);
        return Color(c.x * color_factor_.x, c.y * color_factor_.y, c.z * color_factor_.z);
    }

    float alpha(const HitRecord& hit) const override {
        const float a = base_ ? base_->alpha(hit) : 1.0f;
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

    Color value(const HitRecord& hit) const override {
        return base_ ? base_->value(hit) : Color(1.0f);
    }

    float alpha(const HitRecord& hit) const override {
        const float a = base_ ? base_->alpha(hit) : 1.0f;
        return (a >= cutoff_) ? 1.0f : 0.0f;
    }

private:
    TexturePtr base_;
    float cutoff_ = 0.5f;
};

class ForceAlphaTexture : public Texture {
public:
    ForceAlphaTexture(const TexturePtr& base, float alpha)
        : base_(base), alpha_(clamp_float(alpha, 0.0f, 1.0f)) {}

    Color value(const HitRecord& hit) const override {
        return base_ ? base_->value(hit) : Color(1.0f);
    }

    float alpha(const HitRecord& /*hit*/) const override {
        return alpha_;
    }

private:
    TexturePtr base_;
    float alpha_ = 1.0f;
};
