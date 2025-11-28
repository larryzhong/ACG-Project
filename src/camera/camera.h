#pragma once

#include "core/ray.h"
#include "core/rng.h"
#include "core/vec3.h"

struct CameraSettings {
    Vec3 look_from;
    Vec3 look_at;
    Vec3 up;
    float vertical_fov_deg;
    float aspect_ratio;
    float aperture;
    float focus_dist;
    float t0;
    float t1;

    CameraSettings()
        : look_from(0.0f, 0.0f, 1.0f),
          look_at(0.0f, 0.0f, 0.0f),
          up(0.0f, 1.0f, 0.0f),
          vertical_fov_deg(40.0f),
          aspect_ratio(16.0f / 9.0f),
          aperture(0.0f),
          focus_dist(1.0f),
          t0(0.0f),
          t1(1.0f) {}
};

class Camera {
public:
    explicit Camera(const CameraSettings& settings);

    Ray generate_ray(float s, float t, RNG& rng) const;

private:
    Vec3 origin_;
    Vec3 lower_left_corner_;
    Vec3 horizontal_;
    Vec3 vertical_;
};

inline Camera::Camera(const CameraSettings& settings) {
    const float theta = settings.vertical_fov_deg * static_cast<float>(3.14159265358979323846) / 180.0f;
    const float h = std::tan(theta / 2.0f);
    const float viewport_height = 2.0f * h;
    const float viewport_width = settings.aspect_ratio * viewport_height;

    const Vec3 w = normalize(settings.look_from - settings.look_at);
    const Vec3 u = normalize(cross(settings.up, w));
    const Vec3 v = cross(w, u);

    origin_ = settings.look_from;
    const float focus_dist = settings.focus_dist;

    horizontal_ = focus_dist * viewport_width * u;
    vertical_ = focus_dist * viewport_height * v;
    lower_left_corner_ =
        origin_ - horizontal_ / 2.0f - vertical_ / 2.0f - focus_dist * w;
}

inline Ray Camera::generate_ray(float s, float t, RNG& /*rng*/) const {
    const Vec3 direction =
        lower_left_corner_ + s * horizontal_ + t * vertical_ - origin_;
    return Ray(origin_, direction);
}

