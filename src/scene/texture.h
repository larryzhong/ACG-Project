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

