#pragma once

#include "camera/camera.h"
#include "core/color.h"
#include "core/ray.h"
#include "core/rng.h"
#include "core/sampling.h"
#include "render/film.h"
#include "render/integrator.h"
#include "scene/material.h"
#include "scene/scene.h"
#include "scene/quad.h"
#include "scene/sphere.h"

class PathTracer : public Integrator {
public:
    explicit PathTracer(int max_depth) : max_depth_(max_depth) {}

    Color estimate_direct_lighting(const HitRecord& rec,
                                   const Scene& scene,
                                   const Ray& in_ray,
                                   RNG& rng) const {
        Color result(0.0f);

        const auto& lights = scene.lights.lights();
        if (lights.empty()) {
            return result;
        }

        const Vec3 p = rec.point;
        const Vec3 n = rec.normal;

        for (const auto& light : lights) {
            if (!light.shape) {
                continue;
            }

            const Hittable* shape = light.shape.get();

            Vec3 light_pos;
            Vec3 light_normal;
            float pdf_area = 0.0f;
            const Material* light_material = nullptr;

            if (const auto* sphere = dynamic_cast<const Sphere*>(shape)) {
                const float radius = sphere->radius();
                if (radius <= 0.0f) {
                    continue;
                }

                const Vec3 d = random_unit_vector(rng);
                light_pos = sphere->center() + radius * d;
                light_normal = d;

                const float area = 4.0f * kPi * radius * radius;
                pdf_area = area > 0.0f ? 1.0f / area : 0.0f;
                light_material = sphere->material();
            } else if (const auto* quad = dynamic_cast<const Quad*>(shape)) {
                const float su = rng.uniform();
                const float sv = rng.uniform();

                light_pos = quad->p0() + su * quad->u() + sv * quad->v();
                light_normal = quad->normal();

                const float area = quad->area();
                pdf_area = area > 0.0f ? 1.0f / area : 0.0f;
                light_material = quad->material();
            } else {
                continue;
            }

            if (!light_material || pdf_area <= 0.0f) {
                continue;
            }

            Vec3 to_light = light_pos - p;
            const float dist_squared = to_light.length_squared();
            if (dist_squared <= 0.0f) {
                continue;
            }

            const float dist = std::sqrt(dist_squared);
            const Vec3 wi = to_light / dist;

            const float cos_surface = dot(n, wi);
            const float cos_light = dot(light_normal, -wi);
            if (cos_surface <= 0.0f || cos_light <= 0.0f) {
                continue;
            }

            const float shadow_epsilon = 0.001f;
            Ray shadow_ray(p + shadow_epsilon * wi, wi, in_ray.time);
            HitRecord shadow_rec;
            if (scene.hit(shadow_ray, 0.001f, dist - shadow_epsilon, shadow_rec)) {
                continue;
            }

            HitRecord light_rec;
            light_rec.point = light_pos;
            light_rec.normal = light_normal;
            light_rec.front_face = true;
            light_rec.u = 0.0f;
            light_rec.v = 0.0f;

            const Color emitted = light_material->emitted(light_rec);
            if (emitted.x <= 0.0f && emitted.y <= 0.0f && emitted.z <= 0.0f) {
                continue;
            }

            const float geometry = (cos_surface * cos_light) / dist_squared;
            const float weight = geometry / pdf_area;

            result += emitted * weight;
        }

        // Lambertian BRDF factor (albedo / pi); albedo is applied outside via
        // srec.attenuation, so we include only the 1 / pi here.
        return (1.0f / kPi) * result;
    }

    Color Li(const Ray& r,
             const Scene& scene,
             RNG& rng,
             int depth) const override {
        if (depth >= max_depth_) {
            return Color(0.0f);
        }

        HitRecord rec;
        if (!scene.hit(r, 0.001f, 1e30f, rec)) {
            return Color(0.0f); 
        }

        const Material* material = rec.material;

        if (material) {
            const float op = material->opacity(rec);
            if (op < 1.0f) {
                const float sample = rng.uniform();
                if (sample > op) {
                    Ray continued_ray(rec.point + 0.001f * r.direction,
                                      r.direction,
                                      r.time);
                    return Li(continued_ray, scene, rng, depth + 1);
                }
            }
        }

        Color emitted(0.0f);
        if (material) {
            emitted = material->emitted(rec);
        }

        if (!material) {
            return emitted;
        }

        ScatterRecord srec;
        if (!material->scatter(r, rec, srec, rng)) {
            return emitted;
        }

        Color direct(0.0f);
        if (!srec.is_specular) {
            direct = srec.attenuation *
                     estimate_direct_lighting(rec, scene, r, rng);
        }

        const Color indirect =
            srec.attenuation *
            Li(srec.scattered, scene, rng, depth + 1);

        return emitted + direct + indirect;
    }

private:
    int max_depth_;
};

inline void render_image(const Scene& scene,
                         const Camera& camera,
                         const Integrator& integrator,
                         Film& film,
                         int samples_per_pixel) {
    RNG rng(42u);

    const int width = film.width();
    const int height = film.height();

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            Color pixel_color(0.0f);
            for (int s = 0; s < samples_per_pixel; ++s) {
                const float u = (static_cast<float>(x) + rng.uniform()) /
                                static_cast<float>(width - 1);
                const float v = 1.0f - ((static_cast<float>(y) + rng.uniform()) /
                                        static_cast<float>(height - 1));

                Ray r = camera.generate_ray(u, v, rng);
                pixel_color += integrator.Li(r, scene, rng, 0);
            }

            pixel_color /= static_cast<float>(samples_per_pixel);
            film.add_sample(x, y, pixel_color);
        }
    }
}
