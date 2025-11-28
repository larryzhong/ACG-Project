#pragma once

#include "camera/camera.h"
#include "core/color.h"
#include "core/ray.h"
#include "core/rng.h"
#include "render/film.h"
#include "render/integrator.h"
#include "scene/material.h"
#include "scene/scene.h"

class PathTracer : public Integrator {
public:
    explicit PathTracer(int max_depth) : max_depth_(max_depth) {}

    Color Li(const Ray& r,
             const Scene& scene,
             RNG& rng,
             int depth) const override {
        if (depth >= max_depth_) {
            return Color(0.0f);
        }

        HitRecord rec;
        if (!scene.hit(r, 0.001f, 1e30f, rec)) {
            const Vec3 unit_direction = normalize(r.direction);
            const float t = 0.5f * (unit_direction.y + 1.0f);
            const Color white(1.0f, 1.0f, 1.0f);
            const Color blue(0.5f, 0.7f, 1.0f);
            return (1.0f - t) * white + t * blue;
        }

        Color emitted(0.0f);
        if (rec.material) {
            emitted = rec.material->emitted(rec);
        }

        if (!rec.material) {
            return emitted;
        }

        ScatterRecord srec;
        if (!rec.material->scatter(r, rec, srec, rng)) {
            return emitted;
        }

        return emitted + srec.attenuation *
                              Li(srec.scattered, scene, rng, depth + 1);
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
                const float u =
                    static_cast<float>(x) / static_cast<float>(width - 1);
                const float v =
                    static_cast<float>(y) / static_cast<float>(height - 1);

                Ray r = camera.generate_ray(u, v, rng);
                pixel_color += integrator.Li(r, scene, rng, 0);
            }

            pixel_color /= static_cast<float>(samples_per_pixel);
            film.add_sample(x, y, pixel_color);
        }
    }
}

