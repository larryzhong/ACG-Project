#pragma once

#include <memory>

#include "core/color.h"
#include "scene/material.h"
#include "scene/scene.h"
#include "scene/sphere.h"
#include "scene/texture.h"

inline Scene build_simple_scene_basic() {
    Scene scene;

    // Ground
    auto ground_tex = std::make_shared<SolidColor>(Color(0.8f, 0.8f, 0.0f));
    auto ground_mat = std::make_shared<Lambertian>(ground_tex);
    scene.objects.push_back(
        std::make_shared<Sphere>(Vec3(0.0f, -1000.5f, -1.0f), 1000.0f, ground_mat));

    // Diffuse sphere
    auto diffuse_tex = std::make_shared<SolidColor>(Color(0.7f, 0.3f, 0.3f));
    auto diffuse_mat = std::make_shared<Lambertian>(diffuse_tex);
    scene.objects.push_back(
        std::make_shared<Sphere>(Vec3(-0.6f, 0.0f, -1.0f), 0.5f, diffuse_mat));

    // Metal sphere
    auto metal_tex = std::make_shared<SolidColor>(Color(0.8f, 0.8f, 0.8f));
    auto metal_mat = std::make_shared<Metal>(metal_tex, 0.1f);
    scene.objects.push_back(
        std::make_shared<Sphere>(Vec3(0.6f, 0.0f, -1.0f), 0.5f, metal_mat));

    // Emissive light sphere above
    auto light_tex = std::make_shared<SolidColor>(Color(4.0f, 4.0f, 4.0f));
    auto light_mat = std::make_shared<DiffuseLight>(light_tex);
    scene.objects.push_back(
        std::make_shared<Sphere>(Vec3(0.0f, 2.0f, -1.0f), 0.3f, light_mat));

    return scene;
}
