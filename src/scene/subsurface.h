#pragma once

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
          sigma_s_(std::max(0.0f, sigma_s)),
          sigma_a_(std::max(0.0f, sigma_a)),
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

    float sigma_s() const { return sigma_s_; }
    float sigma_a() const { return sigma_a_; }
    float sigma_t() const { return sigma_s_ + sigma_a_; }

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
        if (!(st > 0.0f) || sigma_s_ <= 0.0f) {
            return false;
        }

        // Start just inside the surface.
        const Vec3 outward_n = entry_hit.normal;
        Vec3 p = entry_hit.point - 0.001f * outward_n;

        // Initial direction in the medium: isotropic.
        Vec3 dir = random_unit_vector(rng);

        Color throughput(1.0f);
        float traveled = 0.0f;

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
                // Apply transmittance to the boundary.
                const float seg = std::exp(-st * t_boundary);
                throughput *= seg;
                traveled += t_boundary;

                // If we accumulated too much distance, stop.
                if (traveled > max_distance_) {
                    return false;
                }

                // We exited the object.
                out.ok = true;
                out.exit_hit = hit;
                out.dir_in_medium = dir;

                // Apply scattering albedo to approximate scattering survival.
                const float albedo = sigma_s_ / st;
                throughput = throughput * Color(albedo, albedo, albedo);

                // Tint by surface albedo at entry (simple artistic control).
                if (albedo_) {
                    throughput = throughput * albedo_->value(entry_hit);
                }

                out.weight = throughput;
                return true;
            }

            // We scatter inside before reaching boundary.
            const float seg = std::exp(-st * s);
            throughput *= seg;
            traveled += s;
            if (traveled > max_distance_) {
                return false;
            }

            // Convert absorption vs scattering: roulette with single-scattering albedo.
            const float scatter_prob = sigma_s_ / st;
            if (rng.uniform() > scatter_prob) {
                return false;
            }

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
    TexturePtr albedo_;
    float sigma_s_ = 0.0f;
    float sigma_a_ = 0.0f;
    int max_steps_ = 64;
    float max_distance_ = 1e30f;
};
