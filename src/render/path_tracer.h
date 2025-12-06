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
        const Vec3 n = rec.material ? rec.material->get_shading_normal(rec) : rec.normal;

        for (const auto& light : lights) {
            if (!light.shape) {
                continue;
            }

            const Hittable* shape = light.shape.get();

            Vec3 light_pos;
            Vec3 light_normal;
            float pdf_area = 0.0f;
            const Material* light_material = nullptr;
            bool is_point_light = (light.type == LightType::Point);

            if (const auto* sphere = dynamic_cast<const Sphere*>(shape)) {
                const float radius = sphere->radius();
                light_material = sphere->material();

                if (is_point_light) {
                    light_pos = sphere->center();
                    light_normal = Vec3(0.0f, 1.0f, 0.0f);
                } else {
                    if (radius <= 0.0f) continue;

                    const Vec3 d = random_unit_vector(rng);
                    light_pos = sphere->center() + radius * d;
                    light_normal = d;

                    const float area = 4.0f * kPi * radius * radius;
                    pdf_area = area > 0.0f ? 1.0f / area : 0.0f;
                }
            } else if (const auto* quad = dynamic_cast<const Quad*>(shape)) {
                if (is_point_light) continue;

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

            if (!light_material || !is_point_light && pdf_area <= 0.0f) {
                continue;
            }

            Vec3 to_light = light_pos - p;
            const float dist_squared = to_light.length_squared();
            if (dist_squared <= 0.000001f) continue;

            const float dist = std::sqrt(dist_squared);
            const Vec3 wi = to_light / dist;

            const float cos_surface = dot(n, wi);

            if (cos_surface <= 0.0f) {
                continue;
            }

            if (!is_point_light) {
                const float cos_light = dot(light_normal, -wi);
                if (cos_light <= 0.0f) {
                    continue;
                }
            }

            bool occluded = false;
            const float shadow_epsilon = 0.001f;
            
            Ray shadow_ray(p + shadow_epsilon * wi, wi, in_ray.time);
            
            float current_t_max = dist - 2.0f * shadow_epsilon;

            while (true) {
                HitRecord shadow_rec;
                if (!scene.hit(shadow_ray, shadow_epsilon, current_t_max, shadow_rec)) {
                    break;
                }

                float op = 1.0f;
                if (shadow_rec.material) {
                    op = shadow_rec.material->opacity(shadow_rec);
                }

                if (rng.uniform() < op) {
                    occluded = true;
                    break;
                }

                shadow_ray = Ray(shadow_rec.point + shadow_epsilon * wi, wi, in_ray.time);
                
                current_t_max -= (shadow_rec.t + shadow_epsilon);

                if (current_t_max <= 0.0f) {
                    break; 
                }
            }

            if (occluded) {
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

            float weight;
            if (is_point_light) {
                weight = cos_surface / dist_squared;
            } else {
                const float cos_light = dot(light_normal, -wi);
                const float geometry = (cos_surface * cos_light) / dist_squared;
                weight = geometry / pdf_area;
            }

            result += emitted * weight;
        }

        return result * (1.0f / kPi);
    }

    Color Li(const Ray& r,
             const Scene& scene,
             RNG& rng,
             int depth) const override {
        return Li_internal(r, scene, rng, depth, true);
    }

private:
    Color Li_internal(const Ray& r,
                    const Scene& scene,
                    RNG& rng,
                    int depth,
                    bool count_emitted) const {
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
                    return Li_internal(continued_ray, scene, rng, depth, count_emitted);
                }
            }
        }

        Color emitted(0.0f);
        if (material && count_emitted) {
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

        float rr_prob = 1.0f;
        if (depth > 3) {
            float max_attenuation = std::max({srec.attenuation.x, srec.attenuation.y, srec.attenuation.z});
            rr_prob = std::max(0.05f, std::min(max_attenuation, 0.95f));
            if (rng.uniform() >= rr_prob) {
                return emitted + direct;
            }
        }

        const bool next_count_emitted = srec.is_specular;

        const Color indirect =
            srec.attenuation *
            Li_internal(srec.scattered, scene, rng, depth + 1, next_count_emitted);

        return emitted + direct + indirect / rr_prob;
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
