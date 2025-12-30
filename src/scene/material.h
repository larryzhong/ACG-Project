#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>

#include "core/color.h"
#include "core/ray.h"
#include "core/rng.h"
#include "core/sampling.h"
#include "scene/hittable.h"
#include "scene/texture.h"
#include "math/onb.h"

namespace {
inline float luminance(const Color& c) {
    return 0.2126f * c.x + 0.7152f * c.y + 0.0722f * c.z;
}

inline float abs_dot(const Vec3& a, const Vec3& b) {
    return std::fabs(dot(a, b));
}

inline Color schlick_fresnel(const Color& f0, float cos_theta) {
    const float m = clamp_float(1.0f - cos_theta, 0.0f, 1.0f);
    const float m2 = m * m;
    const float m5 = m2 * m2 * m;
    return f0 + (Color(1.0f) - f0) * m5;
}

inline float ggx_D(float alpha, float n_dot_h) {
    const float a2 = alpha * alpha;
    const float nh2 = n_dot_h * n_dot_h;
    const float denom = nh2 * (a2 - 1.0f) + 1.0f;
    return a2 / (kPi * denom * denom);
}

inline float ggx_G1(float alpha, float n_dot_v) {
    const float a2 = alpha * alpha;
    const float nv = clamp_float(n_dot_v, 0.0f, 1.0f);
    const float denom = nv + std::sqrt(a2 + (1.0f - a2) * nv * nv);
    if (denom <= 0.0f) {
        return 0.0f;
    }
    return 2.0f * nv / denom;
}

inline float ggx_G(float alpha, float n_dot_v, float n_dot_l) {
    return ggx_G1(alpha, n_dot_v) * ggx_G1(alpha, n_dot_l);
}

inline Vec3 sample_ggx_half_vector(float alpha, RNG& rng) {
    const float u1 = rng.uniform();
    const float u2 = rng.uniform();

    const float phi = 2.0f * kPi * u1;
    const float a2 = alpha * alpha;
    const float cos_theta = std::sqrt((1.0f - u2) / (1.0f + (a2 - 1.0f) * u2));
    const float sin_theta = std::sqrt(std::max(0.0f, 1.0f - cos_theta * cos_theta));

    return Vec3(std::cos(phi) * sin_theta, std::sin(phi) * sin_theta, cos_theta);
}
}  // namespace

class Material {
public:
    virtual ~Material() = default;

    virtual Color eval(const Vec3& /*wo*/,
                       const Vec3& /*wi*/,
                       const HitRecord& /*hit*/) const {
        return Color(0.0f);
    }

    virtual float pdf(const Vec3& /*wo*/,
                      const Vec3& /*wi*/,
                      const HitRecord& /*hit*/) const {
        return 0.0f;
    }

    virtual bool sample(const Vec3& wo,
                        const HitRecord& hit,
                        Vec3& wi,
                        float& out_pdf,
                        Color& out_f,
                        bool& out_is_delta,
                        RNG& rng) const = 0;

    virtual Color emitted(const HitRecord& /*hit*/) const {
        return Color(0.0f);
    }

    virtual float opacity(const HitRecord& /*hit*/) const {
        return 1.0f;
    }

    virtual Vec3 get_shading_normal(const HitRecord& hit) const {
        return hit.normal;
    }

    virtual float cone_roughness(const HitRecord& /*hit*/) const {
        return 1.0f;
    }
};

using MaterialPtr = std::shared_ptr<Material>;

class Lambertian : public Material {
public:
    explicit Lambertian(const TexturePtr& albedo) : albedo_(albedo) {}

    Color eval(const Vec3& /*wo*/,
               const Vec3& wi,
               const HitRecord& hit) const override {
        const Vec3 n = get_shading_normal(hit);
        if (dot(n, wi) <= 0.0f) {
            return Color(0.0f);
        }
        const Color albedo =
            albedo_ ? albedo_->value(hit) : Color(1.0f);
        return albedo * (1.0f / kPi);
    }

    float pdf(const Vec3& /*wo*/,
              const Vec3& wi,
              const HitRecord& hit) const override {
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
        const Vec3 n = get_shading_normal(hit);
        ONB uvw;
        uvw.build_from_w(n);
        wi = uvw.local(random_cosine_direction(rng));

        out_pdf = pdf(Vec3(0.0f), wi, hit);
        out_f = eval(Vec3(0.0f), wi, hit);
        out_is_delta = false;
        return out_pdf > 0.0f;
    }

    float opacity(const HitRecord& hit) const override {
        if (!albedo_) {
            return 1.0f;
        }
        return albedo_->alpha(hit);
    }

private:
    TexturePtr albedo_;
};

class Metal : public Material {
public:
    Metal(const TexturePtr& albedo, float fuzz)
        : albedo_(albedo), fuzz_(fuzz < 1.0f ? fuzz : 1.0f) {}

    Color eval(const Vec3& wo,
               const Vec3& wi,
               const HitRecord& hit) const override {
        if (fuzz_ <= 1e-4f) {
            return Color(0.0f);
        }

        const Vec3 n = get_shading_normal(hit);
        const float n_dot_v = dot(n, wo);
        const float n_dot_l = dot(n, wi);
        if (n_dot_v <= 0.0f || n_dot_l <= 0.0f) {
            return Color(0.0f);
        }

        const Vec3 h = normalize(wo + wi);
        const float n_dot_h = dot(n, h);
        const float v_dot_h = dot(wo, h);
        if (n_dot_h <= 0.0f || v_dot_h <= 0.0f) {
            return Color(0.0f);
        }

        const float alpha = std::max(0.001f, fuzz_ * fuzz_);
        const float D = ggx_D(alpha, n_dot_h);
        const float G = ggx_G(alpha, n_dot_v, n_dot_l);

        const Color f0 =
            albedo_ ? albedo_->value(hit) : Color(1.0f);
        const Color F = schlick_fresnel(f0, v_dot_h);

        return (D * G) * F / (4.0f * n_dot_v * n_dot_l);
    }

    float pdf(const Vec3& wo,
              const Vec3& wi,
              const HitRecord& hit) const override {
        if (fuzz_ <= 1e-4f) {
            return 0.0f;
        }

        const Vec3 n = get_shading_normal(hit);
        const float n_dot_v = dot(n, wo);
        const float n_dot_l = dot(n, wi);
        if (n_dot_v <= 0.0f || n_dot_l <= 0.0f) {
            return 0.0f;
        }

        const Vec3 h = normalize(wo + wi);
        const float n_dot_h = dot(n, h);
        const float v_dot_h = dot(wo, h);
        if (n_dot_h <= 0.0f || v_dot_h <= 0.0f) {
            return 0.0f;
        }

        const float alpha = std::max(0.001f, fuzz_ * fuzz_);
        const float D = ggx_D(alpha, n_dot_h);
        const float pdf_h = D * n_dot_h;
        return pdf_h / (4.0f * v_dot_h);
    }

    bool sample(const Vec3& wo,
                const HitRecord& hit,
                Vec3& wi,
                float& out_pdf,
                Color& out_f,
                bool& out_is_delta,
                RNG& rng) const override {
        const Vec3 n = get_shading_normal(hit);
        const Color f0 =
            albedo_ ? albedo_->value(hit) : Color(1.0f);

        if (fuzz_ <= 1e-4f) {
            wi = reflect(-wo, n);
            if (dot(wi, n) <= 0.0f) {
                return false;
            }
            out_pdf = 1.0f;
            out_f = f0 * (1.0f / std::max(1e-6f, abs_dot(n, wi)));
            out_is_delta = true;
            return true;
        }

        const float alpha = std::max(0.001f, fuzz_ * fuzz_);
        ONB uvw;
        uvw.build_from_w(n);
        const Vec3 h = normalize(uvw.local(sample_ggx_half_vector(alpha, rng)));

        wi = reflect(-wo, h);
        if (dot(wi, n) <= 0.0f) {
            return false;
        }

        out_pdf = pdf(wo, wi, hit);
        out_f = eval(wo, wi, hit);
        out_is_delta = false;
        return out_pdf > 0.0f;
    }

    float opacity(const HitRecord& hit) const override {
        if (!albedo_) {
            return 1.0f;
        }
        return albedo_->alpha(hit);
    }

    float cone_roughness(const HitRecord& /*hit*/) const override {
        return clamp_float(fuzz_, 0.0f, 1.0f);
    }

private:
    TexturePtr albedo_;
    float fuzz_;
};

class DiffuseLight : public Material {
public:
    explicit DiffuseLight(const TexturePtr& emit) : emit_(emit) {}

    bool sample(const Vec3& /*wo*/,
                const HitRecord& /*hit*/,
                Vec3& /*wi*/,
                float& /*out_pdf*/,
                Color& /*out_f*/,
                bool& /*out_is_delta*/,
                RNG& /*rng*/) const override {
        return false;
    }

    Color emitted(const HitRecord& hit) const override {
        if (!emit_ || !hit.front_face) {
            return Color(0.0f);
        }
        return emit_->value(hit);
    }

private:
    TexturePtr emit_;
};

class EmissiveMaterial : public Material {
public:
    EmissiveMaterial(const MaterialPtr& base, const TexturePtr& emission, bool double_sided = false)
        : base_(base), emission_(emission), double_sided_(double_sided) {}

    Color eval(const Vec3& wo,
               const Vec3& wi,
               const HitRecord& hit) const override {
        if (!base_) {
            return Color(0.0f);
        }
        return base_->eval(wo, wi, hit);
    }

    float pdf(const Vec3& wo,
              const Vec3& wi,
              const HitRecord& hit) const override {
        if (!base_) {
            return 0.0f;
        }
        return base_->pdf(wo, wi, hit);
    }

    bool sample(const Vec3& wo,
                const HitRecord& hit,
                Vec3& wi,
                float& out_pdf,
                Color& out_f,
                bool& out_is_delta,
                RNG& rng) const override {
        if (!base_) {
            return false;
        }
        return base_->sample(wo, hit, wi, out_pdf, out_f, out_is_delta, rng);
    }

    Color emitted(const HitRecord& hit) const override {
        if (!emission_) {
            return Color(0.0f);
        }
        if (!double_sided_ && !hit.front_face) {
            return Color(0.0f);
        }
        return emission_->value(hit);
    }

    float opacity(const HitRecord& hit) const override {
        if (!base_) {
            return 1.0f;
        }
        return base_->opacity(hit);
    }

    Vec3 get_shading_normal(const HitRecord& hit) const override {
        if (!base_) {
            return hit.normal;
        }
        return base_->get_shading_normal(hit);
    }

private:
    MaterialPtr base_;
    TexturePtr emission_;
    bool double_sided_ = false;
};

class Dielectric : public Material {
public:
    explicit Dielectric(float index_of_refraction)
        : ir_(index_of_refraction) {}

    bool sample(const Vec3& wo,
                const HitRecord& hit,
                Vec3& wi,
                float& out_pdf,
                Color& out_f,
                bool& out_is_delta,
                RNG& rng) const override {
        const Vec3 n = get_shading_normal(hit);
        const Vec3 incident = normalize(-wo);

        const float refraction_ratio = hit.front_face ? (1.0f / ir_) : ir_;

        const float cos_theta = std::fmin(-dot(incident, n), 1.0f);
        const float sin_theta = std::sqrt(std::max(0.0f, 1.0f - cos_theta * cos_theta));

        const bool cannot_refract = refraction_ratio * sin_theta > 1.0f;
        Vec3 direction;

        const float reflect_prob = schlick(cos_theta, ir_);

        if (cannot_refract || rng.uniform() < reflect_prob) {
            direction = reflect(incident, n);
        } else {
            direction = refract(incident, n, refraction_ratio);
            if (direction.length_squared() < 1e-10f) {
                direction = reflect(incident, n);
            }
        }

        wi = normalize(direction);
        out_pdf = 1.0f;
        out_f = Color(1.0f) * (1.0f / std::max(1e-6f, abs_dot(n, wi)));
        out_is_delta = true;
        return true;
    }

    float cone_roughness(const HitRecord& /*hit*/) const override {
        return 0.0f;
    }

private:
    static float schlick(float cosine, float ref_idx) {
        float r0 = (1.0f - ref_idx) / (1.0f + ref_idx);
        r0 = r0 * r0;
        return r0 + (1.0f - r0) * std::pow(1.0f - cosine, 5.0f);
    }

    float ir_;
};

class DielectricTransmissionBSDF : public Material {
public:
    DielectricTransmissionBSDF(float index_of_refraction,
                               const TexturePtr& transmittance,
                               const NormalMapPtr& normal_map = nullptr,
                               float normal_strength = 1.0f,
                               float transmission_factor = 1.0f,
                               const TexturePtr& transmission_tex = nullptr,
                               bool use_base_alpha_as_transmission = false,
                               float thickness_factor = 0.0f,
                               const TexturePtr& thickness_tex = nullptr,
                               const Color& attenuation_color = Color(1.0f),
                               float attenuation_distance = std::numeric_limits<float>::infinity())
        : ir_(index_of_refraction),
          transmittance_(transmittance),
          normal_map_(normal_map),
          normal_strength_(normal_strength),
          transmission_factor_(transmission_factor),
          transmission_tex_(transmission_tex),
          use_base_alpha_as_transmission_(use_base_alpha_as_transmission),
          thickness_factor_(thickness_factor),
          thickness_tex_(thickness_tex),
          attenuation_color_(attenuation_color),
          attenuation_distance_(attenuation_distance) {}

    Vec3 get_shading_normal(const HitRecord& hit) const override {
        if (normal_map_ && normal_map_->valid()) {
            Vec3 tangent_normal = normal_map_->get_normal(hit.u, hit.v);
            tangent_normal.x *= normal_strength_;
            tangent_normal.y *= normal_strength_;
            tangent_normal = normalize(tangent_normal);
            return NormalMapTexture::apply_normal_map(
                tangent_normal, hit.normal, hit.tangent, hit.tangent_sign);
        }
        return hit.normal;
    }

    bool sample(const Vec3& wo,
                const HitRecord& hit,
                Vec3& wi,
                float& out_pdf,
                Color& out_f,
                bool& out_is_delta,
                RNG& rng) const override {
        const Vec3 n = get_shading_normal(hit);
        const Vec3 incident = normalize(-wo);

        const float refraction_ratio = hit.front_face ? (1.0f / ir_) : ir_;

        const float cos_theta = std::fmin(-dot(incident, n), 1.0f);
        const float sin_theta = std::sqrt(std::max(0.0f, 1.0f - cos_theta * cos_theta));

        const float t = transmission(hit);
        const bool cannot_refract = (refraction_ratio * sin_theta > 1.0f);
        Vec3 direction;

        const float reflect_prob = schlick(cos_theta, ir_);

        if (cannot_refract || rng.uniform() < reflect_prob) {
            direction = reflect(incident, n);
            wi = normalize(direction);
            out_pdf = 1.0f;
            out_f = Color(1.0f) * (1.0f / std::max(1e-6f, abs_dot(n, wi)));
            out_is_delta = true;
            return true;
        }

        if (t <= 1e-6f || rng.uniform() > t) {
            return false;
        }

        direction = refract(incident, n, refraction_ratio);
        if (direction.length_squared() < 1e-10f) {
            direction = reflect(incident, n);
            wi = normalize(direction);
            out_pdf = 1.0f;
            out_f = Color(1.0f) * (1.0f / std::max(1e-6f, abs_dot(n, wi)));
            out_is_delta = true;
            return true;
        }

        wi = normalize(direction);
        const float abs_cos = std::max(1e-6f, abs_dot(n, wi));
        const Color Tr = transmission_color(hit, abs_cos);
        out_pdf = 1.0f;
        out_f = Tr * (1.0f / abs_cos);
        out_is_delta = true;
        return true;
    }

    float cone_roughness(const HitRecord& /*hit*/) const override {
        return 0.0f;
    }

    float opacity(const HitRecord& /*hit*/) const override {
        return 1.0f;
    }

private:
    static float schlick(float cosine, float ref_idx) {
        float r0 = (1.0f - ref_idx) / (1.0f + ref_idx);
        r0 = r0 * r0;
        return r0 + (1.0f - r0) * std::pow(1.0f - cosine, 5.0f);
    }

    float transmission(const HitRecord& hit) const {
        float t = clamp_float(transmission_factor_, 0.0f, 1.0f);
        if (transmission_tex_) {
            t *= clamp_float(transmission_tex_->value(hit).x, 0.0f, 1.0f);
        }
        if (use_base_alpha_as_transmission_ && transmittance_) {
            t *= clamp_float(1.0f - transmittance_->alpha(hit), 0.0f, 1.0f);
        }
        return clamp_float(t, 0.0f, 1.0f);
    }

    static float clamp01(float v) {
        return clamp_float(v, 0.0f, 1.0f);
    }

    Color transmission_color(const HitRecord& hit, float abs_cos_theta) const {
        Color color(1.0f);
        if (transmittance_) {
            const Color c = transmittance_->value(hit);
            color = Color(clamp01(c.x), clamp01(c.y), clamp01(c.z));
        }

        const float thickness = thickness_at(hit);
        if (!(thickness > 0.0f)) {
            return color;
        }

        if (!(attenuation_distance_ > 0.0f) || !std::isfinite(attenuation_distance_)) {
            return color;
        }

        const float distance = thickness / std::max(1e-6f, abs_cos_theta);
        const float exponent = distance / attenuation_distance_;
        if (!(exponent > 0.0f) || !std::isfinite(exponent)) {
            return color;
        }

        const Color a = Color(clamp01(attenuation_color_.x),
                              clamp01(attenuation_color_.y),
                              clamp01(attenuation_color_.z));

        auto pow_c = [&](float v) {
            if (!(v > 0.0f)) {
                return 0.0f;
            }
            return std::pow(v, exponent);
        };

        const Color atten(pow_c(a.x), pow_c(a.y), pow_c(a.z));
        return color * atten;
    }

    float thickness_at(const HitRecord& hit) const {
        float thickness = std::max(0.0f, thickness_factor_);
        if (thickness_tex_) {
            thickness *= std::max(0.0f, thickness_tex_->value(hit).x);
        }
        return thickness;
    }

    float ir_ = 1.5f;
    TexturePtr transmittance_;
    NormalMapPtr normal_map_;
    float normal_strength_ = 1.0f;
    float transmission_factor_ = 1.0f;
    TexturePtr transmission_tex_;
    bool use_base_alpha_as_transmission_ = false;
    float thickness_factor_ = 0.0f;
    TexturePtr thickness_tex_;
    Color attenuation_color_ = Color(1.0f);
    float attenuation_distance_ = std::numeric_limits<float>::infinity();
};

// Mix of diffuse reflection and dielectric reflection/refraction.
// This is a simple reference material (not physically exact), intended for comparisons.
class DiffuseDielectricMix : public Material {
public:
    DiffuseDielectricMix(float index_of_refraction,
                         const TexturePtr& diffuse_albedo,
                         const TexturePtr& transmittance,
                         float transmission_weight)
        : ir_(index_of_refraction),
          diffuse_albedo_(diffuse_albedo),
          transmittance_(transmittance),
          transmission_weight_(clamp_float(transmission_weight, 0.0f, 1.0f)) {}

    Color eval(const Vec3& /*wo*/, const Vec3& wi, const HitRecord& hit) const override {
        const Vec3 n = get_shading_normal(hit);
        if (dot(n, wi) <= 0.0f) {
            return Color(0.0f);
        }
        const float w_diff = 1.0f - transmission_weight_;
        const Color a = diffuse_albedo_ ? diffuse_albedo_->value(hit) : Color(1.0f);
        return w_diff * a * (1.0f / kPi);
    }

    float pdf(const Vec3& /*wo*/, const Vec3& wi, const HitRecord& hit) const override {
        const Vec3 n = get_shading_normal(hit);
        const float cosine = dot(n, wi);
        if (cosine <= 0.0f) {
            return 0.0f;
        }
        const float w_diff = 1.0f - transmission_weight_;
        return w_diff * cosine * (1.0f / kPi);
    }

    bool sample(const Vec3& wo,
                const HitRecord& hit,
                Vec3& wi,
                float& out_pdf,
                Color& out_f,
                bool& out_is_delta,
                RNG& rng) const override {
        const Vec3 n = get_shading_normal(hit);

        // Sample dielectric lobe (delta reflection/refraction).
        if (transmission_weight_ > 0.0f && rng.uniform() < transmission_weight_) {
            const Vec3 incident = normalize(-wo);
            const float refraction_ratio = hit.front_face ? (1.0f / ir_) : ir_;

            const float cos_theta = std::fmin(-dot(incident, n), 1.0f);
            const float sin_theta = std::sqrt(std::max(0.0f, 1.0f - cos_theta * cos_theta));

            const bool cannot_refract = (refraction_ratio * sin_theta > 1.0f);
            const float reflect_prob = schlick(cos_theta, ir_);

            Vec3 direction;
            bool did_refract = false;
            if (cannot_refract || rng.uniform() < reflect_prob) {
                direction = reflect(incident, n);
            } else {
                direction = refract(incident, n, refraction_ratio);
                if (direction.length_squared() < 1e-10f) {
                    direction = reflect(incident, n);
                } else {
                    did_refract = true;
                }
            }

            wi = normalize(direction);
            out_is_delta = true;
            out_pdf = transmission_weight_;

            const float abs_cos = std::max(1e-6f, abs_dot(n, wi));
            Color Tr(1.0f);
            if (did_refract && transmittance_) {
                const Color t = transmittance_->value(hit);
                Tr = Color(clamp_float(t.x, 0.0f, 1.0f),
                           clamp_float(t.y, 0.0f, 1.0f),
                           clamp_float(t.z, 0.0f, 1.0f));
            }

            out_f = transmission_weight_ * Tr * (1.0f / abs_cos);
            return true;
        }

        // Sample diffuse lobe.
        const float w_diff = 1.0f - transmission_weight_;
        if (w_diff <= 0.0f) {
            return false;
        }

        ONB uvw;
        uvw.build_from_w(n);
        wi = uvw.local(random_cosine_direction(rng));

        const float cosine = dot(n, wi);
        if (cosine <= 0.0f) {
            return false;
        }

        out_is_delta = false;
        out_pdf = w_diff * cosine * (1.0f / kPi);
        const Color a = diffuse_albedo_ ? diffuse_albedo_->value(hit) : Color(1.0f);
        out_f = w_diff * a * (1.0f / kPi);
        return out_pdf > 0.0f;
    }

    float cone_roughness(const HitRecord& /*hit*/) const override {
        return 1.0f;
    }

private:
    static float schlick(float cosine, float ref_idx) {
        float r0 = (1.0f - ref_idx) / (1.0f + ref_idx);
        r0 = r0 * r0;
        return r0 + (1.0f - r0) * std::pow(1.0f - cosine, 5.0f);
    }

    float ir_ = 1.5f;
    TexturePtr diffuse_albedo_;
    TexturePtr transmittance_;
    float transmission_weight_ = 0.0f;
};

// Mix of diffuse reflection and refraction-only transmission (NO specular reflection).
// Intentionally non-physical: real dielectrics always have Fresnel reflection.
// Useful as a reference material when you want transmission without mirror-like highlights.
class DiffuseRefractionOnlyMix : public Material {
public:
    DiffuseRefractionOnlyMix(float index_of_refraction,
                             const TexturePtr& diffuse_albedo,
                             const TexturePtr& transmittance,
                             float transmission_weight)
        : ir_(index_of_refraction),
          diffuse_albedo_(diffuse_albedo),
          transmittance_(transmittance),
          transmission_weight_(clamp_float(transmission_weight, 0.0f, 1.0f)) {}

    Color eval(const Vec3& /*wo*/, const Vec3& wi, const HitRecord& hit) const override {
        const Vec3 n = get_shading_normal(hit);
        if (dot(n, wi) <= 0.0f) {
            return Color(0.0f);
        }
        const float w_diff = 1.0f - transmission_weight_;
        const Color a = diffuse_albedo_ ? diffuse_albedo_->value(hit) : Color(1.0f);
        return w_diff * a * (1.0f / kPi);
    }

    float pdf(const Vec3& /*wo*/, const Vec3& wi, const HitRecord& hit) const override {
        const Vec3 n = get_shading_normal(hit);
        const float cosine = dot(n, wi);
        if (cosine <= 0.0f) {
            return 0.0f;
        }
        const float w_diff = 1.0f - transmission_weight_;
        return w_diff * cosine * (1.0f / kPi);
    }

    bool sample(const Vec3& wo,
                const HitRecord& hit,
                Vec3& wi,
                float& out_pdf,
                Color& out_f,
                bool& out_is_delta,
                RNG& rng) const override {
        const Vec3 n = get_shading_normal(hit);

        // Transmission lobe (delta): refraction only.
        if (transmission_weight_ > 0.0f && rng.uniform() < transmission_weight_) {
            const Vec3 incident = normalize(-wo);
            const float refraction_ratio = hit.front_face ? (1.0f / ir_) : ir_;

            const float cos_theta = std::fmin(-dot(incident, n), 1.0f);
            const float sin_theta = std::sqrt(std::max(0.0f, 1.0f - cos_theta * cos_theta));
            const bool cannot_refract = (refraction_ratio * sin_theta > 1.0f);

            if (!cannot_refract) {
                Vec3 direction = refract(incident, n, refraction_ratio);
                if (direction.length_squared() > 1e-10f) {
                    wi = normalize(direction);
                    out_is_delta = true;
                    out_pdf = transmission_weight_;

                    const float abs_cos = std::max(1e-6f, abs_dot(n, wi));
                    Color Tr(1.0f);
                    if (transmittance_) {
                        const Color t = transmittance_->value(hit);
                        Tr = Color(clamp_float(t.x, 0.0f, 1.0f),
                                   clamp_float(t.y, 0.0f, 1.0f),
                                   clamp_float(t.z, 0.0f, 1.0f));
                    }

                    out_f = transmission_weight_ * Tr * (1.0f / abs_cos);
                    return true;
                }
            }

            // If refraction is impossible (e.g. TIR), fall back to diffuse (no reflection).
        }

        // Diffuse lobe.
        const float w_diff = 1.0f - transmission_weight_;
        if (w_diff <= 0.0f) {
            return false;
        }

        ONB uvw;
        uvw.build_from_w(n);
        wi = uvw.local(random_cosine_direction(rng));

        const float cosine = dot(n, wi);
        if (cosine <= 0.0f) {
            return false;
        }

        out_is_delta = false;
        out_pdf = w_diff * cosine * (1.0f / kPi);
        const Color a = diffuse_albedo_ ? diffuse_albedo_->value(hit) : Color(1.0f);
        out_f = w_diff * a * (1.0f / kPi);
        return out_pdf > 0.0f;
    }

    float cone_roughness(const HitRecord& /*hit*/) const override {
        return 1.0f;
    }

private:
    float ir_ = 1.5f;
    TexturePtr diffuse_albedo_;
    TexturePtr transmittance_;
    float transmission_weight_ = 0.0f;
};

class NormalMappedLambertian : public Material {
public:
    NormalMappedLambertian(const TexturePtr& albedo, const NormalMapPtr& normal_map, float strength = 1.0f)
        : albedo_(albedo), normal_map_(normal_map), strength_(strength) {}

    Vec3 get_shading_normal(const HitRecord& hit) const override {
        if (normal_map_ && normal_map_->valid()) {
            Vec3 tangent_normal = normal_map_->get_normal(hit.u, hit.v);
            tangent_normal.x *= strength_;
            tangent_normal.y *= strength_;
            tangent_normal = normalize(tangent_normal);
            return NormalMapTexture::apply_normal_map(
                tangent_normal, hit.normal, hit.tangent, hit.tangent_sign);
        }
        return hit.normal;
    }

    Color eval(const Vec3& /*wo*/,
               const Vec3& wi,
               const HitRecord& hit) const override {
        const Vec3 n = get_shading_normal(hit);
        if (dot(n, wi) <= 0.0f) {
            return Color(0.0f);
        }
        const Color albedo =
            albedo_ ? albedo_->value(hit) : Color(1.0f);
        return albedo * (1.0f / kPi);
    }

    float pdf(const Vec3& /*wo*/,
              const Vec3& wi,
              const HitRecord& hit) const override {
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
        const Vec3 n = get_shading_normal(hit);
        ONB uvw;
        uvw.build_from_w(n);
        wi = uvw.local(random_cosine_direction(rng));

        out_pdf = pdf(Vec3(0.0f), wi, hit);
        out_f = eval(Vec3(0.0f), wi, hit);
        out_is_delta = false;
        return out_pdf > 0.0f;
    }

    float opacity(const HitRecord& hit) const override {
        if (!albedo_) {
            return 1.0f;
        }
        return albedo_->alpha(hit);
    }

private:
    TexturePtr albedo_;
    NormalMapPtr normal_map_;
    float strength_;
};

class PrincipledBSDF : public Material {
public:
    PrincipledBSDF(const TexturePtr& base_color,
                   float metallic_factor,
                   const TexturePtr& metallic_tex,
                   float roughness_factor,
                   const TexturePtr& roughness_tex,
                   const NormalMapPtr& normal_map,
                   float normal_strength,
                   const TexturePtr& occlusion_tex = nullptr,
                   float occlusion_strength = 1.0f)
        : base_color_(base_color),
          metallic_factor_(metallic_factor),
          metallic_tex_(metallic_tex),
          roughness_factor_(roughness_factor),
          roughness_tex_(roughness_tex),
          normal_map_(normal_map),
          normal_strength_(normal_strength),
          occlusion_tex_(occlusion_tex),
          occlusion_strength_(occlusion_strength) {}

    Vec3 get_shading_normal(const HitRecord& hit) const override {
        if (normal_map_ && normal_map_->valid()) {
            Vec3 tangent_normal = normal_map_->get_normal(hit.u, hit.v);
            tangent_normal.x *= normal_strength_;
            tangent_normal.y *= normal_strength_;
            tangent_normal = normalize(tangent_normal);
            return NormalMapTexture::apply_normal_map(
                tangent_normal, hit.normal, hit.tangent, hit.tangent_sign);
        }
        return hit.normal;
    }

    float opacity(const HitRecord& hit) const override {
        if (!base_color_) {
            return 1.0f;
        }
        return base_color_->alpha(hit);
    }

    float cone_roughness(const HitRecord& hit) const override {
        return clamp_float(roughness(hit), 0.0f, 1.0f);
    }

    Color eval(const Vec3& wo,
               const Vec3& wi,
               const HitRecord& hit) const override {
        const Vec3 n = get_shading_normal(hit);
        const float n_dot_v = dot(n, wo);
        const float n_dot_l = dot(n, wi);
        if (n_dot_v <= 0.0f || n_dot_l <= 0.0f) {
            return Color(0.0f);
        }

        const Color base_color = base_color_ ? base_color_->value(hit) : Color(1.0f);
        const float metallic_value = metallic(hit);
        const float alpha = roughness_to_alpha(roughness(hit));

        const Color dielectric_f0(0.04f);
        const Color f0 = (1.0f - metallic_value) * dielectric_f0 + metallic_value * base_color;
        const Color diffuse_color = (1.0f - metallic_value) * base_color;

        const float ao = occlusion(hit);
        Color result = diffuse_color * (ao * (1.0f / kPi));

        const Vec3 h_un = wo + wi;
        if (h_un.length_squared() <= 1e-20f) {
            return result;
        }

        const Vec3 h = normalize(h_un);
        const float n_dot_h = dot(n, h);
        const float v_dot_h = dot(wo, h);
        if (n_dot_h <= 0.0f || v_dot_h <= 0.0f) {
            return result;
        }

        const float D = ggx_D(alpha, n_dot_h);
        const float G = ggx_G(alpha, n_dot_v, n_dot_l);
        const Color F = schlick_fresnel(f0, v_dot_h);

        result += (D * G) * F / (4.0f * n_dot_v * n_dot_l);

        return result;
    }

    float pdf(const Vec3& wo,
              const Vec3& wi,
              const HitRecord& hit) const override {
        const Vec3 n = get_shading_normal(hit);
        const float n_dot_v = dot(n, wo);
        const float n_dot_l = dot(n, wi);
        if (n_dot_v <= 0.0f || n_dot_l <= 0.0f) {
            return 0.0f;
        }

        const float spec_prob = specular_probability(hit);
        const float pdf_diff = n_dot_l * (1.0f / kPi);

        const Vec3 h_un = wo + wi;
        float pdf_spec = 0.0f;
        if (h_un.length_squared() > 1e-20f) {
            const Vec3 h = normalize(h_un);
            const float n_dot_h = dot(n, h);
            const float v_dot_h = dot(wo, h);
            if (n_dot_h > 0.0f && v_dot_h > 0.0f) {
                const float alpha = roughness_to_alpha(roughness(hit));
                const float D = ggx_D(alpha, n_dot_h);
                const float pdf_h = D * n_dot_h;
                pdf_spec = pdf_h / (4.0f * v_dot_h);
            }
        }

        return spec_prob * pdf_spec + (1.0f - spec_prob) * pdf_diff;
    }

    bool sample(const Vec3& wo,
                const HitRecord& hit,
                Vec3& wi,
                float& out_pdf,
                Color& out_f,
                bool& out_is_delta,
                RNG& rng) const override {
        const Vec3 n = get_shading_normal(hit);
        const float n_dot_v = dot(n, wo);
        if (n_dot_v <= 0.0f) {
            return false;
        }

        const float spec_prob = specular_probability(hit);
        const float u = rng.uniform();

        if (u < spec_prob) {
            const float alpha = roughness_to_alpha(roughness(hit));
            ONB uvw;
            uvw.build_from_w(n);
            const Vec3 h = normalize(uvw.local(sample_ggx_half_vector(alpha, rng)));
            wi = reflect(-wo, h);
            if (dot(wi, n) <= 0.0f) {
                return false;
            }
        } else {
            ONB uvw;
            uvw.build_from_w(n);
            wi = uvw.local(random_cosine_direction(rng));
            if (dot(wi, n) <= 0.0f) {
                return false;
            }
        }

        out_pdf = pdf(wo, wi, hit);
        out_f = eval(wo, wi, hit);
        out_is_delta = false;
        return out_pdf > 0.0f;
    }

private:
    float metallic(const HitRecord& hit) const {
        const float tex = metallic_tex_ ? metallic_tex_->value(hit).x : 1.0f;
        return clamp_float(metallic_factor_ * tex, 0.0f, 1.0f);
    }

    float roughness(const HitRecord& hit) const {
        const float tex = roughness_tex_ ? roughness_tex_->value(hit).x : 1.0f;
        return clamp_float(roughness_factor_ * tex, 0.0f, 1.0f);
    }

    static float roughness_to_alpha(float r) {
        const float rr = std::max(0.02f, r);
        return std::max(0.001f, rr * rr);
    }

    float specular_probability(const HitRecord& hit) const {
        const Color base_color = base_color_ ? base_color_->value(hit) : Color(1.0f);
        const float m = metallic(hit);
        const Color dielectric_f0(0.04f);
        const Color f0 = (1.0f - m) * dielectric_f0 + m * base_color;
        const float spec_w = std::max(0.0f, luminance(f0));
        const float diff_w = std::max(0.0f, luminance((1.0f - m) * base_color));
        const float sum = spec_w + diff_w;
        if (sum <= 1e-6f) {
            return 0.5f;
        }
        return clamp_float(spec_w / sum, 0.05f, 0.95f);
    }

    TexturePtr base_color_;
    float metallic_factor_ = 0.0f;
    TexturePtr metallic_tex_;
    float roughness_factor_ = 1.0f;
    TexturePtr roughness_tex_;
    NormalMapPtr normal_map_;
    float normal_strength_ = 1.0f;
    TexturePtr occlusion_tex_;
    float occlusion_strength_ = 1.0f;

    float occlusion(const HitRecord& hit) const {
        if (!occlusion_tex_) {
            return 1.0f;
        }
        const float tex = occlusion_tex_->value(hit).x;
        // glTF: occlusion = mix(1.0, tex, strength)
        const float ao = 1.0f + occlusion_strength_ * (tex - 1.0f);
        return clamp_float(ao, 0.0f, 1.0f);
    }
};
