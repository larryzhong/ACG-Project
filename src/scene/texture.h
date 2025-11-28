#pragma once

#include <memory>

#include "core/color.h"
#include "core/vec3.h"

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
