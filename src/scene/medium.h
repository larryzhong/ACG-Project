#pragma once

#include <cmath>
#include <cstdint>
#include <memory>

#include "core/color.h"
#include "core/ray.h"
#include "core/rng.h"
#include "math/aabb.h"
#include "scene/density_field.h"
#include "scene/phase_function.h"

struct MediumSample {
    bool happened = false;   // true if an interaction occurred before t_max
    bool scattered = false;  // true if interaction was scattering; false means absorption
    float t = 0.0f;          // distance along ray to interaction
    Vec3 p = Vec3(0.0f);     // interaction point
    Color weight = Color(1.0f); // throughput multiplier for this segment/event
};

class Medium {
public:
    virtual ~Medium() = default;

    virtual Color Tr(const Ray& r, float t_max) const = 0;

    virtual bool sample(const Ray& r,
                        float t_max,
                        RNG& rng,
                        MediumSample& out) const = 0;

    virtual const PhaseFunction* phase() const = 0;
};

using MediumPtr = std::shared_ptr<Medium>;

// Simple gray homogeneous medium (same extinction for all color channels).
// This is a good first step; later you can generalize to RGB/channel tracking.
class HomogeneousMedium : public Medium {
public:
    HomogeneousMedium(float sigma_s, float sigma_a)
        : sigma_s_(std::max(0.0f, sigma_s)),
          sigma_a_(std::max(0.0f, sigma_a)),
          sigma_t_(std::max(0.0f, sigma_s_) + std::max(0.0f, sigma_a_)) {}

    Color Tr(const Ray& /*r*/, float t_max) const override {
        if (!(t_max > 0.0f) || sigma_t_ <= 0.0f) {
            return Color(1.0f);
        }
        const float tr = std::exp(-sigma_t_ * t_max);
        return Color(tr, tr, tr);
    }

    bool sample(const Ray& r, float t_max, RNG& rng, MediumSample& out) const override {
        out = MediumSample{};

        if (!(t_max > 0.0f) || sigma_t_ <= 0.0f) {
            out.happened = false;
            out.weight = Color(1.0f);
            return false;
        }

        // Sample a distance from the exponential distribution.
        const float u = std::max(1e-7f, 1.0f - rng.uniform());
        const float t = -std::log(u) / sigma_t_;

        if (t >= t_max) {
            out.happened = false;
            out.weight = Tr(r, t_max);
            return false;
        }

        out.happened = true;
        out.t = t;
        out.p = r.at(t);

        const float scatter_prob = (sigma_t_ > 0.0f) ? (sigma_s_ / sigma_t_) : 0.0f;
        if (rng.uniform() < scatter_prob) {
            out.scattered = true;
            // For exponential distance sampling in a homogeneous medium, the event weight
            // simplifies to albedo (sigma_s / sigma_t) when using a phase function for direction.
            out.weight = Color(scatter_prob, scatter_prob, scatter_prob);
        } else {
            out.scattered = false;
            out.weight = Color(0.0f);
        }

        return out.scattered;
    }

    const PhaseFunction* phase() const override {
        return &phase_;
    }

private:
    float sigma_s_ = 0.0f;
    float sigma_a_ = 0.0f;
    float sigma_t_ = 0.0f;
    IsotropicPhaseFunction phase_;
};

// Heterogeneous (inhomogeneous) bounded medium using delta tracking (Woodcock tracking).
//
// Notes:
// - This implementation is scalar/gray (same coefficients for all color channels).
// - It supports a bounded volume region via an AABB.
// - Tr() uses ratio tracking (unbiased but noisy).
class HeterogeneousMedium : public Medium {
public:
    HeterogeneousMedium(const AABB& bounds,
                        std::shared_ptr<DensityField> density,
                        float sigma_s,
                        float sigma_a,
                        int max_steps = 1024)
        : bounds_(bounds),
          density_(std::move(density)),
          sigma_s_(std::max(0.0f, sigma_s)),
          sigma_a_(std::max(0.0f, sigma_a)),
          sigma_t_base_(sigma_s_ + sigma_a_),
          max_steps_(std::max(1, max_steps)) {
        const float dmax = density_ ? std::max(0.0f, density_->max_density()) : 0.0f;
        majorant_ = sigma_t_base_ * dmax;
    }

    Color Tr(const Ray& r, float t_max) const override {
        if (!(t_max > 0.0f) || majorant_ <= 0.0f || !density_) {
            return Color(1.0f);
        }

        float t0 = 0.0f;
        float t1 = 0.0f;
        if (!intersect_bounds(r, t_max, t0, t1)) {
            return Color(1.0f);
        }

        // Ratio tracking estimator for transmittance.
        float T = 1.0f;
        float t = t0;
        for (int iter = 0; iter < max_steps_; ++iter) {
            const float u = std::max(1e-7f, 1.0f - rng_hash(r, iter));
            const float dt = -std::log(u) / majorant_;
            t += dt;
            if (t >= t1) {
                break;
            }

            const Vec3 p = r.at(t);
            const float sigma_t = sigma_t_at(p);
            const float q = clamp01(1.0f - (sigma_t / majorant_));
            T *= q;
            if (T <= 0.0f) {
                T = 0.0f;
                break;
            }
        }

        return Color(T, T, T);
    }

    bool sample(const Ray& r, float t_max, RNG& rng, MediumSample& out) const override {
        out = MediumSample{};

        if (!(t_max > 0.0f) || majorant_ <= 0.0f || !density_) {
            out.happened = false;
            out.weight = Color(1.0f);
            return false;
        }

        float t0 = 0.0f;
        float t1 = 0.0f;
        if (!intersect_bounds(r, t_max, t0, t1)) {
            out.happened = false;
            out.weight = Color(1.0f);
            return false;
        }

        float t = t0;
        for (int iter = 0; iter < max_steps_; ++iter) {
            const float u = std::max(1e-7f, 1.0f - rng.uniform());
            const float dt = -std::log(u) / majorant_;
            t += dt;
            if (t >= t1) {
                out.happened = false;
                return false;
            }

            const Vec3 p = r.at(t);
            const float sigma_t = sigma_t_at(p);
            if (sigma_t <= 0.0f) {
                continue;
            }

            const float accept = std::min(1.0f, sigma_t / majorant_);
            if (rng.uniform() >= accept) {
                continue;
            }

            out.happened = true;
            out.t = t;
            out.p = p;

            const float scatter_prob = (sigma_t > 0.0f) ? (sigma_s_at(p) / sigma_t) : 0.0f;
            if (rng.uniform() < scatter_prob) {
                out.scattered = true;
                out.weight = Color(scatter_prob, scatter_prob, scatter_prob);
            } else {
                out.scattered = false;
                out.weight = Color(0.0f);
            }

            return out.scattered;
        }

        out.happened = false;
        return false;
    }

    const PhaseFunction* phase() const override {
        return &phase_;
    }

private:
    static float clamp01(float v) {
        return std::max(0.0f, std::min(v, 1.0f));
    }

    // Deterministic uniform in (0,1) from ray + iteration, used to make Tr() const.
    static float rng_hash(const Ray& r, int iter) {
        auto hash_u32 = [](std::uint32_t x) {
            x ^= x >> 16;
            x *= 0x7feb352du;
            x ^= x >> 15;
            x *= 0x846ca68bu;
            x ^= x >> 16;
            return x;
        };

        const std::uint32_t hx = static_cast<std::uint32_t>(std::fabs(r.origin.x) * 1000.0f);
        const std::uint32_t hy = static_cast<std::uint32_t>(std::fabs(r.origin.y) * 1000.0f);
        const std::uint32_t hz = static_cast<std::uint32_t>(std::fabs(r.origin.z) * 1000.0f);
        std::uint32_t h = hx;
        h = hash_u32(h ^ (hy + 0x9e3779b9u));
        h = hash_u32(h ^ (hz + 0x7f4a7c15u));
        h = hash_u32(h ^ static_cast<std::uint32_t>(iter + 1));
        // Map to (0,1).
        const float v = (static_cast<float>(h) + 0.5f) / (static_cast<float>(0xffffffffu) + 1.0f);
        return std::max(1e-7f, std::min(v, 1.0f - 1e-7f));
    }

    bool intersect_bounds(const Ray& r, float t_max, float& out_t0, float& out_t1) const {
        float t0 = 0.0f;
        float t1 = t_max;

        const float o[3] = {r.origin.x, r.origin.y, r.origin.z};
        const float d[3] = {r.direction.x, r.direction.y, r.direction.z};
        const float bmin[3] = {bounds_.min.x, bounds_.min.y, bounds_.min.z};
        const float bmax[3] = {bounds_.max.x, bounds_.max.y, bounds_.max.z};

        for (int axis = 0; axis < 3; ++axis) {
            const float dir = d[axis];
            const float orig = o[axis];

            if (std::fabs(dir) < 1e-12f) {
                // Ray parallel to slab: must be within bounds.
                if (orig < bmin[axis] || orig > bmax[axis]) {
                    return false;
                }
                continue;
            }

            const float inv_d = 1.0f / dir;
            float t_near = (bmin[axis] - orig) * inv_d;
            float t_far = (bmax[axis] - orig) * inv_d;
            if (t_near > t_far) {
                std::swap(t_near, t_far);
            }

            t0 = std::max(t0, t_near);
            t1 = std::min(t1, t_far);
            if (t1 <= t0) {
                return false;
            }
        }

        out_t0 = std::max(0.0f, t0);
        out_t1 = t1;
        return out_t1 > out_t0;
    }

    float sigma_t_at(const Vec3& p) const {
        const float d = std::max(0.0f, density_->density(p));
        return sigma_t_base_ * d;
    }

    float sigma_s_at(const Vec3& p) const {
        const float d = std::max(0.0f, density_->density(p));
        return sigma_s_ * d;
    }

    AABB bounds_;
    std::shared_ptr<DensityField> density_;
    float sigma_s_ = 0.0f;
    float sigma_a_ = 0.0f;
    float sigma_t_base_ = 0.0f;
    float majorant_ = 0.0f;
    int max_steps_ = 1024;
    IsotropicPhaseFunction phase_;
};
