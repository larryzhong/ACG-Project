#pragma once

#include "core/color.h"
#include "core/rng.h"
#include "core/sampling.h"
#include "core/vec3.h"

class PhaseFunction {
public:
    virtual ~PhaseFunction() = default;

    virtual Color eval(const Vec3& /*wo*/, const Vec3& /*wi*/) const {
        return Color(0.0f);
    }

    virtual float pdf(const Vec3& /*wo*/, const Vec3& /*wi*/) const {
        return 0.0f;
    }

    virtual bool sample(const Vec3& wo,
                        Vec3& wi,
                        float& out_pdf,
                        Color& out_f,
                        RNG& rng) const = 0;
};

class IsotropicPhaseFunction : public PhaseFunction {
public:
    Color eval(const Vec3& /*wo*/, const Vec3& /*wi*/) const override {
        return Color(1.0f) * (1.0f / (4.0f * kPi));
    }

    float pdf(const Vec3& /*wo*/, const Vec3& /*wi*/) const override {
        return 1.0f / (4.0f * kPi);
    }

    bool sample(const Vec3& /*wo*/, Vec3& wi, float& out_pdf, Color& out_f, RNG& rng) const override {
        wi = random_unit_vector(rng);
        out_pdf = pdf(Vec3(0.0f), wi);
        out_f = eval(Vec3(0.0f), wi);
        return out_pdf > 0.0f;
    }
};
