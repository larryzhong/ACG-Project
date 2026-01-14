#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>

#include "core/color.h"
#include "core/rng.h"
#include "core/sampling.h"
#include "core/vec3.h"
#include "scene/hittable.h"
#include "scene/material.h"

struct SubsurfaceSample {
    bool ok = false;
    HitRecord exit_hit;
    Vec3 dir_in_medium = Vec3(0.0f); // direction of the last segment inside the medium (towards the exit)
    Color weight = Color(1.0f);      // throughput from entry to exit
};

// Random-walk subsurface scattering: treat the interior of the object as a homogeneous medium.
// This is integrator-driven (B1): the integrator detects this material and runs a subsurface branch.
class SubsurfaceRandomWalkMaterial : public Material {
public:
    SubsurfaceRandomWalkMaterial(const TexturePtr& albedo,
                                 float sigma_s,
                                 float sigma_a,
                                 int max_steps = 64,
                                 float max_distance = 1e30f)
        : albedo_(albedo),
                    sigma_s_rgb_(Color(std::max(0.0f, sigma_s))),
                    sigma_a_rgb_(Color(std::max(0.0f, sigma_a))),
                    sigma_t_scalar_(compute_sigma_t_scalar(sigma_s_rgb_, sigma_a_rgb_)),
          max_steps_(std::max(1, max_steps)),
          max_distance_(max_distance) {}

        SubsurfaceRandomWalkMaterial(const TexturePtr& albedo,
                                                                 const Color& sigma_s_rgb,
                                                                 const Color& sigma_a_rgb,
                                                                 int max_steps = 64,
                                                                 float max_distance = 1e30f)
                : albedo_(albedo),
                    sigma_s_rgb_(clamp_nonneg(sigma_s_rgb)),
                    sigma_a_rgb_(clamp_nonneg(sigma_a_rgb)),
                    sigma_t_scalar_(compute_sigma_t_scalar(sigma_s_rgb_, sigma_a_rgb_)),
                    max_steps_(std::max(1, max_steps)),
                    max_distance_(max_distance) {}

    // Surface BSDF: keep it simple (diffuse). The integrator will *not* use these for the entry event,
    // but they can be used if you ever hit the material from the outside and choose to shade as a surface.
    Color eval(const Vec3& /*wo*/, const Vec3& wi, const HitRecord& hit) const override {
        const Vec3 n = get_shading_normal(hit);
        if (dot(n, wi) <= 0.0f) {
            return Color(0.0f);
        }
        const Color a = albedo_ ? albedo_->value(hit) : Color(1.0f);
        return a * (1.0f / kPi);
    }

    float pdf(const Vec3& /*wo*/, const Vec3& wi, const HitRecord& hit) const override {
        const Vec3 n = get_shading_normal(hit);
        const float cosine = dot(n, wi);
        return (cosine <= 0.0f) ? 0.0f : cosine * (1.0f / kPi);
    }

    bool sample(const Vec3& /*wo*/,
                const HitRecord& hit,
                Vec3& wi,
                float& out_pdf,
                Color& out_f,
                bool& out_is_delta,
                RNG& rng) const override {
        // Default to diffuse; the integrator overrides entry behavior for SSS.
        ONB uvw;
        uvw.build_from_w(get_shading_normal(hit));
        wi = uvw.local(random_cosine_direction(rng));
        out_pdf = pdf(Vec3(0.0f), wi, hit);
        out_f = eval(Vec3(0.0f), wi, hit);
        out_is_delta = false;
        return out_pdf > 0.0f;
    }

    float cone_roughness(const HitRecord& /*hit*/) const override {
        return 1.0f;
    }

    // Scalar coefficients used for distance sampling (channel-independent transport).
    float sigma_s() const { return luminance(sigma_s_rgb_); }
    float sigma_a() const { return luminance(sigma_a_rgb_); }
    float sigma_t() const { return sigma_t_scalar_; }

    // RGB coefficients used for per-channel throughput.
    Color sigma_s_rgb() const { return sigma_s_rgb_; }
    Color sigma_a_rgb() const { return sigma_a_rgb_; }
    Color sigma_t_rgb() const { return sigma_s_rgb_ + sigma_a_rgb_; }

    Color sss_albedo_rgb() const {
        const Color st = sigma_t_rgb();
        return safe_div(sigma_s_rgb_, st);
    }

    Color albedo(const HitRecord& hit) const {
        return albedo_ ? albedo_->value(hit) : Color(1.0f);
    }

    // Entry->exit random walk inside the same hittable.
    bool sample_subsurface(const HitRecord& entry_hit,
                           const Hittable& boundary,
                           RNG& rng,
                           SubsurfaceSample& out) const {
        out = SubsurfaceSample{};

        const float st = sigma_t();
        if (!(st > 0.0f) || max_component(sigma_s_rgb_) <= 0.0f) {
            return false;
        }

        // Start just inside the surface.
        const Vec3 outward_n = entry_hit.normal;
        Vec3 p = entry_hit.point - 0.001f * outward_n;

        // Initial direction in the medium: isotropic.
        Vec3 dir = random_unit_vector(rng);

        Color throughput(1.0f);
        float traveled = 0.0f;

        const Color sigma_t_rgb_local = sigma_t_rgb();
        const Color albedo_rgb_local = sss_albedo_rgb();
        const float scatter_q = clamp_float(max_component(albedo_rgb_local), 0.0f, 0.999f);

        for (int step = 0; step < max_steps_; ++step) {
            // Sample free-flight distance.
            const float u = std::max(1e-7f, 1.0f - rng.uniform());
            const float s = -std::log(u) / st;

            // Find distance to the boundary in this direction.
            Ray ray(p, dir, 0.0f);
            // We need the next intersection with the same object boundary.
            HitRecord hit;
            if (!boundary.hit(ray, 0.001f, std::min(max_distance_, 1e30f), hit)) {
                return false;
            }

            const float t_boundary = hit.t;

            // Exit if boundary is closer than the sampled collision distance.
            if (s >= t_boundary) {
                // Apply RGB transmittance to the boundary.
                throughput = throughput * exp_rgb(-sigma_t_rgb_local * t_boundary);
                traveled += t_boundary;

                // If we accumulated too much distance, stop.
                if (traveled > max_distance_) {
                    return false;
                }

                // We exited the object.
                out.ok = true;
                out.exit_hit = hit;
                out.dir_in_medium = dir;

                // Apply a final albedo tint (channel-independent transport, channel-dependent color).
                if (scatter_q > 0.0f) {
                    throughput = throughput * (albedo_rgb_local / scatter_q);
                }

                // Tint by surface albedo at entry (simple artistic control).
                if (albedo_) {
                    throughput = throughput * albedo_->value(entry_hit);
                }

                out.weight = throughput;
                return true;
            }

            // We scatter inside before reaching boundary.
            throughput = throughput * exp_rgb(-sigma_t_rgb_local * s);
            traveled += s;
            if (traveled > max_distance_) {
                return false;
            }

            // Absorption vs scattering: roulette based on max-channel albedo.
            if (scatter_q <= 0.0f || rng.uniform() > scatter_q) {
                return false;
            }

            throughput = throughput * (albedo_rgb_local / scatter_q);

            // Scatter: move to new point and pick a new direction.
            p = ray.at(s);
            dir = random_unit_vector(rng);

            // Russian roulette after a few steps.
            if (step > 8) {
                const float q = std::max(0.05f, std::min(0.95f, std::max({throughput.x, throughput.y, throughput.z})));
                if (rng.uniform() > q) {
                    return false;
                }
                throughput /= q;
            }
        }

        return false;
    }

private:
    static Color clamp_nonneg(const Color& c) {
        return Color(std::max(0.0f, c.x), std::max(0.0f, c.y), std::max(0.0f, c.z));
    }

    static float luminance(const Color& c) {
        return 0.2126f * c.x + 0.7152f * c.y + 0.0722f * c.z;
    }

    static float max_component(const Color& c) {
        return std::max({c.x, c.y, c.z});
    }

    static Color safe_div(const Color& a, const Color& b) {
        const float ex = (b.x > 1e-8f) ? (a.x / b.x) : 0.0f;
        const float ey = (b.y > 1e-8f) ? (a.y / b.y) : 0.0f;
        const float ez = (b.z > 1e-8f) ? (a.z / b.z) : 0.0f;
        return Color(ex, ey, ez);
    }

    static Color exp_rgb(const Color& v) {
        return Color(std::exp(v.x), std::exp(v.y), std::exp(v.z));
    }

    static float compute_sigma_t_scalar(const Color& sigma_s_rgb, const Color& sigma_a_rgb) {
        return std::max(0.0f, luminance(sigma_s_rgb + sigma_a_rgb));
    }

    TexturePtr albedo_;
    Color sigma_s_rgb_ = Color(0.0f);
    Color sigma_a_rgb_ = Color(0.0f);
    float sigma_t_scalar_ = 0.0f;
    int max_steps_ = 64;
    float max_distance_ = 1e30f;
};
