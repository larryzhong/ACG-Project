#pragma once

#include <algorithm>
#include <cmath>

#include "core/vec3.h"

class DensityField {
public:
    virtual ~DensityField() = default;

    // Returns non-negative density scale at point p.
    virtual float density(const Vec3& p) const = 0;

    // Returns a global upper bound for density(p).
    virtual float max_density() const = 0;
};

class ConstantDensityField : public DensityField {
public:
    explicit ConstantDensityField(float d) : d_(std::max(0.0f, d)) {}

    float density(const Vec3& /*p*/) const override {
        return d_;
    }

    float max_density() const override {
        return d_;
    }

private:
    float d_ = 0.0f;
};

// Simple analytic inhomogeneous density: product of sines.
// density(p) = clamp(bias + amp * sin(fx*x) * sin(fy*y) * sin(fz*z), 0, bias+amp)
class SineProductDensityField : public DensityField {
public:
    SineProductDensityField(float freq, float bias, float amp)
        : freq_(freq), bias_(bias), amp_(amp) {
        if (freq_ < 0.0f) freq_ = 0.0f;
        if (bias_ < 0.0f) bias_ = 0.0f;
        if (amp_ < 0.0f) amp_ = 0.0f;
    }

    float density(const Vec3& p) const override {
        if (freq_ <= 0.0f) {
            return bias_;
        }
        const float s = std::sin(freq_ * p.x) * std::sin(freq_ * p.y) * std::sin(freq_ * p.z);
        const float d = bias_ + amp_ * s;
        return std::max(0.0f, std::min(d, max_density()));
    }

    float max_density() const override {
        return bias_ + amp_;
    }

private:
    float freq_ = 0.0f;
    float bias_ = 0.0f;
    float amp_ = 0.0f;
};

// Height falloff: density(p) = exp(-k * max(0, p.y - y0))
class HeightFogDensityField : public DensityField {
public:
    HeightFogDensityField(float y0, float k)
        : y0_(y0), k_(std::max(0.0f, k)) {}

    float density(const Vec3& p) const override {
        const float dy = std::max(0.0f, p.y - y0_);
        if (k_ <= 0.0f) {
            return 1.0f;
        }
        return std::exp(-k_ * dy);
    }

    float max_density() const override {
        return 1.0f;
    }

private:
    float y0_ = 0.0f;
    float k_ = 0.0f;
};

// Linear gradient along Y:
// density(p) = lerp(d0, d1, clamp((p.y - y0) / (y1 - y0), 0, 1))
class LinearGradientDensityField : public DensityField {
public:
    LinearGradientDensityField(float y0, float y1, float d0, float d1)
        : y0_(y0), y1_(y1), d0_(std::max(0.0f, d0)), d1_(std::max(0.0f, d1)) {}

    float density(const Vec3& p) const override {
        if (std::fabs(y1_ - y0_) < 1e-8f) {
            return std::max(d0_, d1_);
        }
        float t = (p.y - y0_) / (y1_ - y0_);
        t = std::max(0.0f, std::min(t, 1.0f));
        const float d = d0_ + t * (d1_ - d0_);
        return std::max(0.0f, d);
    }

    float max_density() const override {
        return std::max(d0_, d1_);
    }

private:
    float y0_ = 0.0f;
    float y1_ = 1.0f;
    float d0_ = 0.0f;
    float d1_ = 0.0f;
};
