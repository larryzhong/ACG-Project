#pragma once

#include <algorithm>
#include <cmath>

#include "camera/camera.h"
#include "core/color.h"
#include "core/ray.h"
#include "core/rng.h"
#include "core/sampling.h"
#include "render/film.h"
#include "render/integrator.h"
#include "scene/material.h"
#include "scene/medium.h"
#include "scene/mesh.h"
#include "scene/scene.h"
#include "scene/quad.h"
#include "scene/sphere.h"
#include "scene/subsurface.h"

class PathTracer : public Integrator {
public:
    explicit PathTracer(int max_depth) : max_depth_(max_depth) {}

    static RayCone advance_cone(const RayCone& cone, float distance) {
        RayCone out = cone;
        out.width = cone.width + cone.spread_angle * std::max(0.0f, distance);
        return out;
    }

    static RayCone scatter_cone(const RayCone& cone,
                                float distance,
                                float roughness,
                                bool is_delta) {
        RayCone out = advance_cone(cone, distance);
        roughness = clamp_float(roughness, 0.0f, 1.0f);

        const float min_spread =
            is_delta ? (0.02f + 0.15f * roughness)
                     : (0.15f + 0.35f * roughness);

        out.spread_angle = std::max(out.spread_angle, min_spread);
        return out;
    }

    static float mis_weight(float pdf_a, float pdf_b) {
        const float denom = pdf_a + pdf_b;
        if (denom <= 0.0f) {
            return 0.0f;
        }
        return pdf_a / denom;
    }

    Color estimate_direct_lighting(const HitRecord& rec,
                                   const Scene& scene,
                                   const Vec3& wo,
                                   const Ray& in_ray,
                                   RNG& rng) const {
        Color result(0.0f);

        const auto& lights = scene.lights.lights();
        const bool has_env = scene.environment && scene.environment->valid();
        const int light_count = static_cast<int>(lights.size()) + (has_env ? 1 : 0);
        if (light_count <= 0) {
            return result;
        }

        const Vec3 p = rec.point;
        const Vec3 n = rec.material ? rec.material->get_shading_normal(rec) : rec.normal;
        const Material* bsdf = rec.material;
        if (!bsdf) {
            return result;
        }

        const float choose = rng.uniform() * static_cast<float>(light_count);
        const int picked = std::min(static_cast<int>(choose), light_count - 1);

        if (has_env && picked == static_cast<int>(lights.size())) {
            float pdf_env = 0.0f;
            const Vec3 wi = scene.environment->sample(pdf_env, rng);
            if (pdf_env <= 0.0f) {
                return result;
            }

            const float cos_surface = dot(n, wi);
            if (cos_surface <= 0.0f) {
                return result;
            }

            const float shadow_epsilon = 0.001f;
            Ray shadow_ray(p + shadow_epsilon * wi, wi, in_ray.time);
            shadow_ray.medium = in_ray.medium;
            HitRecord shadow_rec;
            if (trace_first_hit(scene, shadow_ray, shadow_epsilon, 1e30f, rng, shadow_rec)) {
                return result;
            }

            const Color Tr = medium_tr(scene, shadow_ray, 1e30f);
            if (is_black(Tr)) {
                return result;
            }

            const Color Le = scene.environment->Le(wi);
            if (Le.x <= 0.0f && Le.y <= 0.0f && Le.z <= 0.0f) {
                return result;
            }

            const Color f = bsdf->eval(wo, wi, rec);
            if (f.x <= 0.0f && f.y <= 0.0f && f.z <= 0.0f) {
                return result;
            }

            const float pdf_light = (pdf_env / static_cast<float>(light_count));
            const float pdf_bsdf = bsdf->pdf(wo, wi, rec);
            const float w = mis_weight(pdf_light, pdf_bsdf);
            result += Tr * Le * f * (cos_surface / pdf_light) * w;
            return result;
        }

        if (picked < 0 || picked >= static_cast<int>(lights.size())) {
            return result;
        }

        const Light& light = lights[static_cast<std::size_t>(picked)];
        if (!light.shape) {
            return result;
        }

        const Hittable* shape = light.shape.get();
        const bool is_point_light = (light.type == LightType::Point);

        Vec3 light_pos;
        Vec3 light_normal;
        float light_u = 0.0f;
        float light_v = 0.0f;
        float pdf_area = 0.0f;
        const Material* light_material = nullptr;

        if (const auto* sphere = dynamic_cast<const Sphere*>(shape)) {
            const float radius = sphere->radius();
            light_material = sphere->material();

            if (is_point_light) {
                light_pos = sphere->center();
                light_normal = Vec3(0.0f, 1.0f, 0.0f);
            } else {
                if (radius <= 0.0f) {
                    return result;
                }

                const Vec3 d = random_unit_vector(rng);
                light_pos = sphere->center() + radius * d;
                light_normal = d;

                const float area = 4.0f * kPi * radius * radius;
                pdf_area = area > 0.0f ? 1.0f / area : 0.0f;
            }
        } else if (const auto* quad = dynamic_cast<const Quad*>(shape)) {
            if (is_point_light) {
                return result;
            }

            const float su = rng.uniform();
            const float sv = rng.uniform();

            light_pos = quad->p0() + su * quad->u() + sv * quad->v();
            light_normal = quad->normal();
            light_u = su;
            light_v = sv;

            const float area = quad->area();
            pdf_area = area > 0.0f ? 1.0f / area : 0.0f;
            light_material = quad->material();
        } else if (const auto* mesh = dynamic_cast<const Mesh*>(shape)) {
            if (is_point_light) {
                return result;
            }
            light_material = mesh->material();
            if (!mesh->sample_surface(rng, light_pos, light_normal, light_u, light_v, pdf_area)) {
                return result;
            }
        } else {
            return result;
        }

        if (!light_material) {
            return result;
        }

        Vec3 to_light = light_pos - p;
        const float dist_squared = to_light.length_squared();
        if (dist_squared <= 0.000001f) {
            return result;
        }

        const float dist = std::sqrt(dist_squared);
        const Vec3 wi = to_light / dist;

        const float cos_surface = dot(n, wi);
        if (cos_surface <= 0.0f) {
            return result;
        }

        float pdf_dir_given_light = 0.0f;
        if (!is_point_light) {
            if (pdf_area <= 0.0f) {
                return result;
            }
            const float cos_light = std::fabs(dot(light_normal, -wi));
            if (cos_light <= 0.0f) {
                return result;
            }
            pdf_dir_given_light = pdf_area * dist_squared / cos_light;
        }

        const float shadow_epsilon = 0.001f;
        Ray shadow_ray(p + shadow_epsilon * wi, wi, in_ray.time);
        shadow_ray.medium = in_ray.medium;
        HitRecord shadow_rec;
        if (trace_first_hit(scene, shadow_ray, shadow_epsilon, dist - 2.0f * shadow_epsilon, rng, shadow_rec)) {
            return result;
        }

        const Color Tr = medium_tr(scene, shadow_ray, dist - 2.0f * shadow_epsilon);
        if (is_black(Tr)) {
            return result;
        }

        HitRecord light_rec;
        light_rec.point = light_pos;
        light_rec.t = dist;
        light_rec.u = light_u;
        light_rec.v = light_v;
        light_rec.material = light_material;
        light_rec.object = shape;
        light_rec.front_face = (dot(wi, light_normal) < 0.0f);
        light_rec.normal = light_rec.front_face ? light_normal : -light_normal;
        light_rec.tangent = Vec3(1.0f, 0.0f, 0.0f);

        const Color Le = light_material->emitted(light_rec);
        if (Le.x <= 0.0f && Le.y <= 0.0f && Le.z <= 0.0f) {
            return result;
        }

        const Color f = bsdf->eval(wo, wi, rec);
        if (f.x <= 0.0f && f.y <= 0.0f && f.z <= 0.0f) {
            return result;
        }

        if (is_point_light) {
            const float scale = static_cast<float>(light_count);
            result += Tr * scale * Le * f * (cos_surface / dist_squared);
            return result;
        }

        const float pdf_light = (pdf_dir_given_light / static_cast<float>(light_count));
        const float pdf_bsdf = bsdf->pdf(wo, wi, rec);
        const float w = mis_weight(pdf_light, pdf_bsdf);
        result += Tr * Le * f * (cos_surface / pdf_light) * w;

        return result;
    }

    Color Li(const Ray& r,
             const Scene& scene,
             RNG& rng,
             int depth) const override {
        return Li_internal(r, scene, rng, depth, true);
    }

private:
    static const Medium* current_medium(const Scene& scene, const Ray& r) {
        if (r.medium) {
            return r.medium;
        }
        return scene.global_medium.get();
    }

    static bool is_black(const Color& c) {
        return c.x <= 0.0f && c.y <= 0.0f && c.z <= 0.0f;
    }

    static Color medium_tr(const Scene& scene, const Ray& r, float t_max) {
        const Medium* m = current_medium(scene, r);
        if (!m) {
            return Color(1.0f);
        }
        return m->Tr(r, t_max);
    }

    static Color occlusion_transmittance(const Scene& scene,
                                         Ray ray,
                                         float t_min,
                                         float t_max,
                                         RNG& rng,
                                         const Hittable* ignore_object = nullptr) {
        if (t_max <= t_min) {
            return Color(1.0f);
        }

        // If a (sufficiently opaque) surface blocks the segment, return 0.
        // For cutout alpha this is stochastic (matches existing behavior).
        HitRecord hit;
        if (trace_first_hit(scene, ray, t_min, t_max, rng, hit, ignore_object)) {
            return Color(0.0f);
        }

        return medium_tr(scene, ray, t_max);
    }

    Color estimate_direct_lighting_medium(const Vec3& p,
                                          const Vec3& wo,
                                          const Ray& in_ray,
                                          const PhaseFunction& phase,
                                          const Scene& scene,
                                          RNG& rng,
                                          const Hittable* ignore_object = nullptr) const {
        Color result(0.0f);

        const auto& lights = scene.lights.lights();
        const bool has_env = scene.environment && scene.environment->valid();
        const int light_count = static_cast<int>(lights.size()) + (has_env ? 1 : 0);
        if (light_count <= 0) {
            return result;
        }

        const float choose = rng.uniform() * static_cast<float>(light_count);
        const int picked = std::min(static_cast<int>(choose), light_count - 1);

        if (has_env && picked == static_cast<int>(lights.size())) {
            float pdf_env = 0.0f;
            const Vec3 wi = scene.environment->sample(pdf_env, rng);
            if (pdf_env <= 0.0f) {
                return result;
            }

            Ray shadow_ray(p + 0.001f * wi, wi, in_ray.time);
            shadow_ray.medium = in_ray.medium;
            const Color Tr = occlusion_transmittance(scene, shadow_ray, 0.001f, 1e30f, rng, ignore_object);
            if (is_black(Tr)) {
                return result;
            }

            const Color Le = scene.environment->Le(wi);
            if (is_black(Le)) {
                return result;
            }

            const Color ph = phase.eval(wo, wi);
            if (is_black(ph)) {
                return result;
            }

            const float pdf_light = (pdf_env / static_cast<float>(light_count));
            if (pdf_light <= 0.0f) {
                return result;
            }

            result += Tr * Le * ph * (1.0f / pdf_light);
            return result;
        }

        if (picked < 0 || picked >= static_cast<int>(lights.size())) {
            return result;
        }

        const Light& light = lights[static_cast<std::size_t>(picked)];
        if (!light.shape) {
            return result;
        }

        const Hittable* shape = light.shape.get();
        const bool is_point_light = (light.type == LightType::Point);

        Vec3 light_pos;
        Vec3 light_normal;
        float light_u = 0.0f;
        float light_v = 0.0f;
        float pdf_area = 0.0f;
        const Material* light_material = nullptr;

        if (const auto* sphere = dynamic_cast<const Sphere*>(shape)) {
            const float radius = sphere->radius();
            light_material = sphere->material();

            if (is_point_light) {
                light_pos = sphere->center();
                light_normal = Vec3(0.0f, 1.0f, 0.0f);
            } else {
                if (radius <= 0.0f) {
                    return result;
                }

                const Vec3 d = random_unit_vector(rng);
                light_pos = sphere->center() + radius * d;
                light_normal = d;

                const float area = 4.0f * kPi * radius * radius;
                pdf_area = area > 0.0f ? 1.0f / area : 0.0f;
            }
        } else if (const auto* quad = dynamic_cast<const Quad*>(shape)) {
            if (is_point_light) {
                return result;
            }

            const float su = rng.uniform();
            const float sv = rng.uniform();

            light_pos = quad->p0() + su * quad->u() + sv * quad->v();
            light_normal = quad->normal();
            light_u = su;
            light_v = sv;

            const float area = quad->area();
            pdf_area = area > 0.0f ? 1.0f / area : 0.0f;
            light_material = quad->material();
        } else if (const auto* mesh = dynamic_cast<const Mesh*>(shape)) {
            if (is_point_light) {
                return result;
            }
            light_material = mesh->material();
            if (!mesh->sample_surface(rng, light_pos, light_normal, light_u, light_v, pdf_area)) {
                return result;
            }
        } else {
            return result;
        }

        if (!light_material) {
            return result;
        }

        Vec3 to_light = light_pos - p;
        const float dist_squared = to_light.length_squared();
        if (dist_squared <= 0.000001f) {
            return result;
        }

        const float dist = std::sqrt(dist_squared);
        const Vec3 wi = to_light / dist;

        float pdf_dir_given_light = 0.0f;
        if (!is_point_light) {
            if (pdf_area <= 0.0f) {
                return result;
            }
            const float cos_light = std::fabs(dot(light_normal, -wi));
            if (cos_light <= 0.0f) {
                return result;
            }
            pdf_dir_given_light = pdf_area * dist_squared / cos_light;
        }

        const float shadow_epsilon = 0.001f;
        Ray shadow_ray(p + shadow_epsilon * wi, wi, in_ray.time);
        shadow_ray.medium = in_ray.medium;
        const Color Tr = occlusion_transmittance(scene, shadow_ray, shadow_epsilon, dist - 2.0f * shadow_epsilon, rng, ignore_object);
        if (is_black(Tr)) {
            return result;
        }

        HitRecord light_rec;
        light_rec.point = light_pos;
        light_rec.t = dist;
        light_rec.u = light_u;
        light_rec.v = light_v;
        light_rec.material = light_material;
        light_rec.object = shape;
        light_rec.front_face = (dot(wi, light_normal) < 0.0f);
        light_rec.normal = light_rec.front_face ? light_normal : -light_normal;
        light_rec.tangent = Vec3(1.0f, 0.0f, 0.0f);

        const Color Le = light_material->emitted(light_rec);
        if (is_black(Le)) {
            return result;
        }

        const Color ph = phase.eval(wo, wi);
        if (is_black(ph)) {
            return result;
        }

        if (is_point_light) {
            const float scale = static_cast<float>(light_count);
            result += Tr * scale * Le * ph * (1.0f / dist_squared);
            return result;
        }

        const float pdf_light = (pdf_dir_given_light / static_cast<float>(light_count));
        if (pdf_light <= 0.0f) {
            return result;
        }

        result += Tr * Le * ph * (1.0f / pdf_light);
        return result;
    }

    static bool is_light_surface(const Scene& scene, const Hittable* obj) {
        if (!obj) {
            return false;
        }
        for (const auto& light : scene.lights.lights()) {
            if (light.shape && light.shape.get() == obj) {
                return true;
            }
        }
        return false;
    }

    static bool trace_first_hit(const Scene& scene,
                                Ray ray,
                                float t_min,
                                float t_max,
                                RNG& rng,
                                HitRecord& out_rec,
                                const Hittable* ignore_object = nullptr) {
        float current_t_max = t_max;
        const float eps = 0.001f;
        for (int iter = 0; iter < 256; ++iter) {
            HitRecord rec;
            if (!scene.hit(ray, t_min, current_t_max, rec)) {
                return false;
            }

            // Treat the ignored object as fully transparent (used for subsurface shadow rays).
            if (ignore_object && rec.object == ignore_object) {
                Ray next_ray(rec.point + eps * ray.direction, ray.direction, ray.time);
                next_ray.medium = ray.medium;
                ray = next_ray;
                current_t_max -= (rec.t + eps);
                if (current_t_max <= 0.0f) {
                    return false;
                }
                continue;
            }

            float op = 1.0f;
            if (rec.material) {
                op = rec.material->opacity(rec);
            }

            if (rng.uniform() < op) {
                out_rec = rec;
                return true;
            }

            Ray next_ray(rec.point + eps * ray.direction, ray.direction, ray.time);
            next_ray.medium = ray.medium;
            ray = next_ray;
            current_t_max -= (rec.t + eps);
            if (current_t_max <= 0.0f) {
                return false;
            }
        }
        return false;
    }

    static float light_pdf_direction(const Scene& scene,
                                     const Vec3& p,
                                     const Vec3& wi,
                                     const HitRecord& light_hit) {
        if (!is_light_surface(scene, light_hit.object)) {
            return 0.0f;
        }

        const Vec3 diff = light_hit.point - p;
        const float dist_squared = diff.length_squared();
        if (dist_squared <= 0.0f) {
            return 0.0f;
        }

        if (const auto* quad = dynamic_cast<const Quad*>(light_hit.object)) {
            const float cos_light = std::fabs(dot(light_hit.normal, -wi));
            if (cos_light <= 0.0f) {
                return 0.0f;
            }

            const float area = quad->area();
            if (area <= 0.0f) {
                return 0.0f;
            }

            return dist_squared / (cos_light * area);
        } else if (const auto* sphere = dynamic_cast<const Sphere*>(light_hit.object)) {
            const float cos_light = std::fabs(dot(light_hit.normal, -wi));
            if (cos_light <= 0.0f) {
                return 0.0f;
            }

            const float r = sphere->radius();
            const float area = 4.0f * kPi * r * r;
            if (area <= 0.0f) {
                return 0.0f;
            }

            return dist_squared / (cos_light * area);
        } else if (const auto* mesh = dynamic_cast<const Mesh*>(light_hit.object)) {
            const float pdf_area = mesh->pdf_area(light_hit);
            if (pdf_area <= 0.0f) {
                return 0.0f;
            }

            const Vec3 gn = mesh->geometric_normal(light_hit.primitive_id);
            const float cos_light = std::fabs(dot(gn, -wi));
            if (cos_light <= 0.0f) {
                return 0.0f;
            }

            return pdf_area * dist_squared / cos_light;
        } else {
            return 0.0f;
        }
    }

    Color bsdf_sample_emitter(const Scene& scene,
                              const HitRecord& rec,
                              const Ray& in_ray,
                              const Vec3& wi,
                              const Color& f,
                              float pdf_bsdf,
                              RNG& rng) const {
        const auto& lights = scene.lights.lights();
        const bool has_env = scene.environment && scene.environment->valid();
        const int light_count = static_cast<int>(lights.size()) + (has_env ? 1 : 0);
        if (light_count <= 0 || pdf_bsdf <= 0.0f) {
            return Color(0.0f);
        }

        const Vec3 n = rec.material ? rec.material->get_shading_normal(rec) : rec.normal;
        const float cos_surface = dot(n, wi);
        if (cos_surface <= 0.0f) {
            return Color(0.0f);
        }

        const float eps = 0.001f;
        Ray ray(rec.point + eps * wi, wi, in_ray.time);
        ray.medium = in_ray.medium;

        HitRecord hit;
        if (!trace_first_hit(scene, ray, eps, 1e30f, rng, hit)) {
            if (!has_env) {
                return Color(0.0f);
            }

            const Color Tr = medium_tr(scene, ray, 1e30f);
            if (is_black(Tr)) {
                return Color(0.0f);
            }

            const Color Le = scene.environment->Le(wi);
            if (Le.x <= 0.0f && Le.y <= 0.0f && Le.z <= 0.0f) {
                return Color(0.0f);
            }

            const float pdf_light = scene.environment->pdf(wi) / static_cast<float>(light_count);
            const float w = mis_weight(pdf_bsdf, pdf_light);
            return Tr * Le * f * (cos_surface / pdf_bsdf) * w;
        }

        if (!hit.material || !is_light_surface(scene, hit.object)) {
            return Color(0.0f);
        }

        const Color Le = hit.material->emitted(hit);
        if (Le.x <= 0.0f && Le.y <= 0.0f && Le.z <= 0.0f) {
            return Color(0.0f);
        }

        const Color Tr = medium_tr(scene, ray, hit.t);
        if (is_black(Tr)) {
            return Color(0.0f);
        }

        const float pdf_dir = light_pdf_direction(scene, rec.point, wi, hit);
        if (pdf_dir <= 0.0f) {
            return Color(0.0f);
        }

        const float pdf_light = pdf_dir / static_cast<float>(light_count);
        const float w = mis_weight(pdf_bsdf, pdf_light);
        return Tr * Le * f * (cos_surface / pdf_bsdf) * w;
    }

    Color Li_internal(const Ray& r,
                    const Scene& scene,
                    RNG& rng,
                    int depth,
                    bool count_emitted) const {
        if (depth >= max_depth_) {
            return Color(0.0f);
        }

        HitRecord rec;
        const bool hit_surface = scene.hit(r, 0.001f, 1e30f, rec);
        const float t_surface = hit_surface ? rec.t : 1e30f;

        // Sample a medium interaction before the first surface hit.
        const Medium* m = current_medium(scene, r);
        if (m) {
            MediumSample ms;
            m->sample(r, t_surface, rng, ms);
            if (ms.happened) {
                if (!ms.scattered) {
                    return Color(0.0f);
                }

                const PhaseFunction* phase = m->phase();
                if (!phase) {
                    return Color(0.0f);
                }

                const Vec3 wo = normalize(-r.direction);
                Vec3 wi;
                float pdf_phase = 0.0f;
                Color f_phase(0.0f);
                if (!phase->sample(wo, wi, pdf_phase, f_phase, rng) || pdf_phase <= 0.0f) {
                    return Color(0.0f);
                }

                Color direct = estimate_direct_lighting_medium(ms.p, wo, r, *phase, scene, rng);

                float rr_prob = 1.0f;
                if (depth > 3) {
                    const float max_w = std::max({ms.weight.x, ms.weight.y, ms.weight.z});
                    rr_prob = std::max(0.05f, std::min(max_w, 0.95f));
                    if (rng.uniform() >= rr_prob) {
                        return direct;
                    }
                }

                Ray scattered(ms.p, wi, r.time, scatter_cone(r.cone, ms.t, 1.0f, false));
                scattered.medium = r.medium;
                const Color phase_weight = f_phase * (1.0f / pdf_phase);
                const Color indirect = ms.weight * phase_weight *
                    Li_internal(scattered, scene, rng, depth + 1, true);
                return direct + indirect / rr_prob;
            }
        }

        if (!hit_surface) {
            if (scene.environment && scene.environment->valid()) {
                if (!count_emitted) {
                    return Color(0.0f);
                }
                if (depth == 0 && scene.hide_environment_background) {
                    return Color(0.0f);
                }
                const Color Tr = medium_tr(scene, r, 1e30f);
                return Tr * scene.environment->Le(r.direction);
            }
            return Color(0.0f);
        }

        const Color segment_Tr = medium_tr(scene, r, t_surface);

        const Material* material = rec.material;

        if (material) {
            const float op = material->opacity(rec);
            if (op < 1.0f) {
                const float sample = rng.uniform();
                if (sample > op) {
                    const float travel = (rec.point - r.origin).length();
                    Ray continued_ray(rec.point + 0.001f * r.direction,
                                      r.direction,
                                      r.time,
                                      advance_cone(r.cone, travel));
                    continued_ray.medium = r.medium;
                    return Li_internal(continued_ray, scene, rng, depth, count_emitted);
                }
            }
        }

        if (!material) {
            return Color(0.0f);
        }

        // -------------------------
        // Subsurface scattering (integrator-driven random walk)
        // -------------------------
        if (rec.front_face) {
            if (const auto* sss = dynamic_cast<const SubsurfaceRandomWalkMaterial*>(material)) {
                if (rec.object) {
                    // Integrator-driven subsurface as a bounded homogeneous medium.
                    // Key difference vs the previous implementation:
                    // - We do NEE from scattering events *inside* the object, so a pure backlight can still glow
                    //   even if there is no ground/walls.
                    // - Shadow rays from inside ignore hits against the same boundary object.

                    const float st = sss->sigma_t();
                    if (!(st > 0.0f) || sss->sigma_s() <= 0.0f) {
                        return Color(0.0f);
                    }

                    const float eps = 0.001f;
                    Vec3 p = rec.point - eps * rec.normal; // just inside
                    Vec3 dir = random_unit_vector(rng);

                    Color throughput = sss->albedo(rec);
                    Color accum(0.0f);
                    float traveled = 0.0f;

                    IsotropicPhaseFunction phase;

                    int local_depth = depth;
                    for (int step = 0; step < 512; ++step) {
                        if (local_depth >= max_depth_) {
                            break;
                        }

                        // Sample free-flight distance.
                        const float u = std::max(1e-7f, 1.0f - rng.uniform());
                        const float s = -std::log(u) / st;

                        // Find the distance to the boundary along this direction.
                        Ray boundary_ray(p, dir, r.time);
                        HitRecord boundary_hit;
                        if (!rec.object->hit(boundary_ray, eps, 1e30f, boundary_hit)) {
                            break;
                        }

                        const float t_boundary = boundary_hit.t;

                        // If we reach the boundary before the next collision, exit the medium.
                        if (s >= t_boundary) {
                            const float seg = std::exp(-st * t_boundary);
                            throughput *= seg;
                            traveled += t_boundary;

                            const Vec3 exit_p = boundary_hit.point + eps * dir;
                            const float travel_total = (exit_p - r.origin).length();
                            Ray exit_ray(exit_p, dir, r.time,
                                         scatter_cone(r.cone, travel_total, 1.0f, false));
                            exit_ray.medium = r.medium;

                            accum += throughput * Li_internal(exit_ray, scene, rng, local_depth + 1, true);
                            break;
                        }

                        // We collide inside the medium.
                        const float seg = std::exp(-st * s);
                        throughput *= seg;
                        traveled += s;

                        // Absorption vs scattering: terminate if absorbed.
                        const float scatter_prob = sss->sigma_s() / st;
                        if (rng.uniform() > scatter_prob) {
                            break;
                        }

                        // Move to the scattering point.
                        p = boundary_ray.at(s);

                        // Direct lighting from this scattering event.
                        {
                            const Vec3 wo = normalize(-dir);
                            Ray inside_context(p, dir, r.time);
                            inside_context.medium = r.medium;
                            const Color direct = estimate_direct_lighting_medium(
                                p, wo, inside_context, phase, scene, rng, rec.object);
                            accum += throughput * direct;
                        }

                        // Sample new direction (phase function) and continue.
                        {
                            const Vec3 wo = normalize(-dir);
                            Vec3 wi;
                            float pdf_phase = 0.0f;
                            Color f_phase(0.0f);
                            if (!phase.sample(wo, wi, pdf_phase, f_phase, rng) || pdf_phase <= 0.0f) {
                                break;
                            }

                            throughput = throughput * (f_phase * (1.0f / pdf_phase));
                            dir = wi;
                        }

                        // Russian roulette to avoid long random walks.
                        if (step > 8) {
                            const float max_w = std::max({throughput.x, throughput.y, throughput.z});
                            const float q = std::max(0.05f, std::min(0.95f, max_w));
                            if (rng.uniform() > q) {
                                break;
                            }
                            throughput /= q;
                        }

                        ++local_depth;
                    }

                    return segment_Tr * accum;
                }
            }
        }

        const Vec3 wo = normalize(-r.direction);

        Vec3 wi;
        float pdf = 0.0f;
        Color f(0.0f);
        bool is_delta = false;
        const bool did_sample = material->sample(wo, rec, wi, pdf, f, is_delta, rng);

        Color emitted = material->emitted(rec);
        if (!count_emitted && is_light_surface(scene, rec.object)) {
            emitted = Color(0.0f);
        }

        if (!did_sample) {
            return segment_Tr * emitted;
        }

        Color direct(0.0f);
        if (!is_delta) {
            direct = estimate_direct_lighting(rec, scene, wo, r, rng);
            direct += bsdf_sample_emitter(scene, rec, r, wi, f, pdf, rng);
        }

        const Vec3 n = material->get_shading_normal(rec);
        const float cos_term = std::fabs(dot(n, wi));
        if (pdf <= 0.0f || cos_term <= 0.0f) {
            return segment_Tr * (emitted + direct);
        }

        Color weight = f * (cos_term / pdf);

        float rr_prob = 1.0f;
        if (depth > 3) {
            float max_w = std::max({weight.x, weight.y, weight.z});
            rr_prob = std::max(0.05f, std::min(max_w, 0.95f));
            if (rng.uniform() >= rr_prob) {
                return segment_Tr * (emitted + direct);
            }
        }

        const bool next_count_emitted = is_delta;
        const float travel = (rec.point - r.origin).length();
        const float roughness = material->cone_roughness(rec);
        Ray scattered(rec.point, wi, r.time, scatter_cone(r.cone, travel, roughness, is_delta));
        scattered.medium = r.medium;
        const Color indirect = weight * Li_internal(scattered, scene, rng, depth + 1, next_count_emitted);
        return segment_Tr * (emitted + direct + indirect / rr_prob);
    }

    int max_depth_;
};

#include <atomic>
#include <thread>
inline void render_image(const Scene& scene,
                         const Camera& camera,
                         const Integrator& integrator,
                         Film& film,
                         int samples_per_pixel) {
    const int width = film.width();
    const int height = film.height();
    
    std::atomic<int> next_row{0};
    std::atomic<int> completed_rows{0};
    
    auto worker = [&]() {
        while (true) {
            int y = next_row.fetch_add(1);
            if (y >= height) break;
            
            RNG rng(42u + static_cast<std::uint64_t>(y) * 1000000u);
            
            for (int x = 0; x < width; ++x) {
                Color pixel_color(0.0f);
                for (int s = 0; s < samples_per_pixel; ++s) {
                    const float u = (static_cast<float>(x) + rng.uniform()) /
                                    static_cast<float>(width - 1);
                    const float v = 1.0f - ((static_cast<float>(y) + rng.uniform()) /
                                            static_cast<float>(height - 1));

                    Ray r = camera.generate_ray(u, v, rng);
                    r.medium = scene.global_medium.get();
                    pixel_color += integrator.Li(r, scene, rng, 0);
                }
                pixel_color /= static_cast<float>(samples_per_pixel);
                film.add_sample(x, y, pixel_color);
            }
            
            int done = completed_rows.fetch_add(1) + 1;
            if (done % 10 == 0) {
                std::cerr << "\rProgress: " << (100 * done / height) << "%" << std::flush;
            }
        }
    };
    
    unsigned int num_threads = std::thread::hardware_concurrency();
    if (num_threads == 0) num_threads = 4;
    
    std::vector<std::thread> threads;
    for (unsigned int i = 0; i < num_threads; ++i) {
        threads.emplace_back(worker);
    }
    for (auto& t : threads) {
        t.join();
    }
    
    std::cerr << "\rProgress: 100%\n";
}
