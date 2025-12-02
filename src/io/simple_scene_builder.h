#pragma once

#include <memory>

#include "core/color.h"
#include "scene/material.h"
#include "scene/moving_sphere.h"
#include "scene/quad.h"
#include "scene/scene.h"
#include "scene/sphere.h"
#include "scene/texture.h"

inline Scene build_simple_scene_basic() {
    Scene scene;

    // --- Materials ---
    auto white_tex = std::make_shared<SolidColor>(Color(0.73f, 0.73f, 0.73f));
    auto red_tex = std::make_shared<SolidColor>(Color(0.65f, 0.05f, 0.05f));
    auto green_tex = std::make_shared<SolidColor>(Color(0.12f, 0.45f, 0.15f));

    auto white = std::make_shared<Lambertian>(white_tex);
    auto red = std::make_shared<Lambertian>(red_tex);
    auto green = std::make_shared<Lambertian>(green_tex);

    // Main light (bright)
    auto light_tex = std::make_shared<SolidColor>(Color(15.0f, 15.0f, 15.0f));
    auto light_mat = std::make_shared<DiffuseLight>(light_tex);

    // Fill light (dimmer) for the area behind the camera
    auto fill_light_tex = std::make_shared<SolidColor>(Color(5.0f, 5.0f, 5.0f));
    auto fill_light_mat = std::make_shared<DiffuseLight>(fill_light_tex);

    // --- Geometry ---
    // Scene extends from z = -2 (back of box) to z = 4 (behind camera)
    
    const float x0 = -1.0f, x1 = 1.0f;
    const float y0 = 0.0f,  y1 = 2.0f;
    const float z_box_back = -2.0f;
    const float z_studio_back = 4.0f; // Wall behind camera

    // 1. FLOOR (Continuous checkerboard)
    auto floor_even = std::make_shared<SolidColor>(Color(0.8f, 0.8f, 0.8f));
    auto floor_odd = std::make_shared<SolidColor>(Color(0.3f, 0.3f, 0.3f));
    auto floor_checker = std::make_shared<CheckerTexture>(floor_even, floor_odd, 2.0f);
    auto floor_mat = std::make_shared<Lambertian>(floor_checker);

    scene.objects.push_back(std::make_shared<Quad>(
        Vec3(x0, y0, z_box_back),
        Vec3(x1 - x0, 0.0f, 0.0f),
        Vec3(0.0f, 0.0f, z_studio_back - z_box_back), // Covers entire depth
        floor_mat));

    // 2. CEILING (Continuous white)
    scene.objects.push_back(std::make_shared<Quad>(
        Vec3(x0, y1, z_box_back),
        Vec3(x1 - x0, 0.0f, 0.0f),
        Vec3(0.0f, 0.0f, z_studio_back - z_box_back),
        white));

    // 3. WALLS
    
    // Back wall (White)
    scene.objects.push_back(std::make_shared<Quad>(
        Vec3(x0, y0, z_box_back),
        Vec3(x1 - x0, 0.0f, 0.0f),
        Vec3(0.0f, y1 - y0, 0.0f),
        white));

    // [FIX] Left wall (Red) - EXTENDED
    // Now stretches from the back of the box (-2) all the way to behind camera (4)
    scene.objects.push_back(std::make_shared<Quad>(
        Vec3(x0, y0, z_studio_back),
        Vec3(0.0f, 0.0f, z_box_back - z_studio_back), // Vector pointing to -Z
        Vec3(0.0f, y1 - y0, 0.0f),
        red));

    // [FIX] Right wall (Green) - EXTENDED
    // Now stretches from the back of the box (-2) all the way to behind camera (4)
    scene.objects.push_back(std::make_shared<Quad>(
        Vec3(x1, y0, z_box_back),
        Vec3(0.0f, 0.0f, z_studio_back - z_box_back), // Vector pointing to +Z
        Vec3(0.0f, y1 - y0, 0.0f),
        green));

    // Studio Back Wall (Behind Camera - White)
    scene.objects.push_back(std::make_shared<Quad>(
        Vec3(x0, y0, z_studio_back),
        Vec3(x1 - x0, 0.0f, 0.0f),
        Vec3(0.0f, y1 - y0, 0.0f), // Normal points to -Z (into scene)
        white));

    // 4. LIGHTS
    // Main Ceiling Light (Inside Box)
    const float lx0 = -0.3f, lx1 = 0.3f;
    const float lz0 = -1.3f, lz1 = -0.7f;
    const float ly = y1 - 0.001f;

    auto ceiling_light = std::make_shared<Quad>(
        Vec3(lx0, ly, lz0),
        Vec3(lx1 - lx0, 0.0f, 0.0f),
        Vec3(0.0f, 0.0f, lz1 - lz0),
        light_mat);
    scene.objects.push_back(ceiling_light);
    scene.lights.add_area_light(ceiling_light);

    // Fill Light (Behind Camera)
    auto studio_light = std::make_shared<Quad>(
        Vec3(lx0, ly, 1.8f),
        Vec3(lx1 - lx0, 0.0f, 0.0f),
        Vec3(0.0f, 0.0f, 0.6f),
        fill_light_mat);
    scene.objects.push_back(studio_light);
    scene.lights.add_area_light(studio_light);

    // --- Objects ---
    // Metal sphere (left side)
    auto metal_tex = std::make_shared<SolidColor>(Color(0.9f, 0.9f, 0.9f));
    auto metal = std::make_shared<Metal>(metal_tex, 0.05f);
    scene.objects.push_back(
        std::make_shared<Sphere>(Vec3(-0.5f, 0.35f, -0.7f), 0.35f, metal));

    // Glass sphere (right side)
    auto glass = std::make_shared<Dielectric>(1.5f);
    scene.objects.push_back(
        std::make_shared<Sphere>(Vec3(0.5f, 0.35f, -1.0f), 0.35f, glass));

    // Small diffuse sphere (center back)
    auto blue_tex = std::make_shared<SolidColor>(Color(0.2f, 0.3f, 0.7f));
    auto blue_mat = std::make_shared<Lambertian>(blue_tex);
    scene.objects.push_back(
        std::make_shared<Sphere>(Vec3(0.0f, 0.25f, -1.5f), 0.25f, blue_mat));

    return scene;
}
