#pragma once

#include <cmath>
#include <memory>

#include "core/color.h"
#include "core/ray.h"
#include "core/rng.h"
#include "core/sampling.h"
#include "scene/hittable.h"
#include "scene/texture.h"
#include "math/onb.h"

struct ScatterRecord {
    Ray scattered;
    Color attenuation;
    bool is_specular;
};

class Material {
public:
    virtual ~Material() = default;

    virtual bool scatter(const Ray& r_in,
                         const HitRecord& hit,
                         ScatterRecord& srec,
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
};

using MaterialPtr = std::shared_ptr<Material>;

class Lambertian : public Material {
public:
    explicit Lambertian(const TexturePtr& albedo) : albedo_(albedo) {}

    bool scatter(const Ray& r_in,
                const HitRecord& hit,
                ScatterRecord& srec,
                RNG& rng) const override {
        ONB uvw;
        uvw.build_from_w(hit.normal);
        Vec3 scatter_direction = uvw.local(random_cosine_direction(rng));

        srec.scattered = Ray(hit.point, scatter_direction, r_in.time);
        srec.attenuation =
            albedo_ ? albedo_->value(hit.u, hit.v, hit.point) : Color(1.0f);
        srec.is_specular = false;
        return true;
    }

    float opacity(const HitRecord& hit) const override {
        if (!albedo_) {
            return 1.0f;
        }
        return albedo_->alpha(hit.u, hit.v, hit.point);
    }

private:
    TexturePtr albedo_;
};

class Metal : public Material {
public:
    Metal(const TexturePtr& albedo, float fuzz)
        : albedo_(albedo), fuzz_(fuzz < 1.0f ? fuzz : 1.0f) {}

    bool scatter(const Ray& r_in,
                 const HitRecord& hit,
                 ScatterRecord& srec,
                 RNG& rng) const override {
        const Vec3 reflected =
            reflect(normalize(r_in.direction), hit.normal);
        const Vec3 perturbed =
            reflected + fuzz_ * random_in_unit_sphere(rng);

        srec.scattered = Ray(hit.point, perturbed, r_in.time);
        srec.attenuation =
            albedo_ ? albedo_->value(hit.u, hit.v, hit.point) : Color(1.0f);
        srec.is_specular = true;

        return dot(srec.scattered.direction, hit.normal) > 0.0f;
    }

    float opacity(const HitRecord& hit) const override {
        if (!albedo_) {
            return 1.0f;
        }
        return albedo_->alpha(hit.u, hit.v, hit.point);
    }

private:
    TexturePtr albedo_;
    float fuzz_;
};

class DiffuseLight : public Material {
public:
    explicit DiffuseLight(const TexturePtr& emit) : emit_(emit) {}

    bool scatter(const Ray& /*r_in*/,
                 const HitRecord& /*hit*/,
                 ScatterRecord& /*srec*/,
                 RNG& /*rng*/) const override {
        return false;
    }

    Color emitted(const HitRecord& hit) const override {
        if (!emit_ || !hit.front_face) {
            return Color(0.0f);
        }
        return emit_->value(hit.u, hit.v, hit.point);
    }

private:
    TexturePtr emit_;
};

class EmissiveMaterial : public Material {
public:
    EmissiveMaterial(const MaterialPtr& base, const TexturePtr& emission, bool double_sided = false)
        : base_(base), emission_(emission), double_sided_(double_sided) {}

    bool scatter(const Ray& r_in,
                 const HitRecord& hit,
                 ScatterRecord& srec,
                 RNG& rng) const override {
        if (!base_) {
            return false;
        }
        return base_->scatter(r_in, hit, srec, rng);
    }

    Color emitted(const HitRecord& hit) const override {
        if (!emission_) {
            return Color(0.0f);
        }
        if (!double_sided_ && !hit.front_face) {
            return Color(0.0f);
        }
        return emission_->value(hit.u, hit.v, hit.point);
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

    bool scatter(const Ray& r_in,
                 const HitRecord& hit,
                 ScatterRecord& srec,
                 RNG& rng) const override {
        srec.attenuation = Color(1.0f, 1.0f, 1.0f);
        srec.is_specular = true;

        const float refraction_ratio = hit.front_face ? (1.0f / ir_) : ir_;

        const Vec3 unit_direction = normalize(r_in.direction);
        const float cos_theta =
            std::fmin(-dot(unit_direction, hit.normal), 1.0f);
        const float sin_theta = std::sqrt(1.0f - cos_theta * cos_theta);

        const bool cannot_refract = refraction_ratio * sin_theta > 1.0f;
        Vec3 direction;

        const float reflect_prob = schlick(cos_theta, ir_);

        if (cannot_refract || rng.uniform() < reflect_prob) {
            direction = reflect(unit_direction, hit.normal);
        } else {
            direction = refract(unit_direction, hit.normal, refraction_ratio);
            if (direction.length_squared() < 1e-10f) {
                direction = reflect(unit_direction, hit.normal);
            }
        }

        srec.scattered = Ray(hit.point, direction, r_in.time);
        return true;
    }

private:
    static float schlick(float cosine, float ref_idx) {
        float r0 = (1.0f - ref_idx) / (1.0f + ref_idx);
        r0 = r0 * r0;
        return r0 + (1.0f - r0) * std::pow(1.0f - cosine, 5.0f);
    }

    float ir_;
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
                tangent_normal, hit.normal, hit.tangent);
        }
        return hit.normal;
    }

    bool scatter(const Ray& r_in,
                const HitRecord& hit,
                ScatterRecord& srec,
                RNG& rng) const override {
        Vec3 shading_normal = get_shading_normal(hit);

        ONB uvw;
        uvw.build_from_w(shading_normal);
        Vec3 scatter_direction = uvw.local(random_cosine_direction(rng));

        srec.scattered = Ray(hit.point, scatter_direction, r_in.time);
        srec.attenuation =
            albedo_ ? albedo_->value(hit.u, hit.v, hit.point) : Color(1.0f);
        srec.is_specular = false;
        return true;
    }

    float opacity(const HitRecord& hit) const override {
        if (!albedo_) {
            return 1.0f;
        }
        return albedo_->alpha(hit.u, hit.v, hit.point);
    }

private:
    TexturePtr albedo_;
    NormalMapPtr normal_map_;
    float strength_;
};
