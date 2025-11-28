#pragma once

#include <memory>

#include "core/color.h"
#include "scene/material.h"
#include "scene/quad.h"
#include "scene/scene.h"
#include "scene/sphere.h"
#include "scene/texture.h"

inline Scene build_simple_scene_basic() {
    Scene scene;

    // Materials
    auto white_tex = std::make_shared<SolidColor>(Color(0.73f, 0.73f, 0.73f));
    auto red_tex = std::make_shared<SolidColor>(Color(0.65f, 0.05f, 0.05f));
    auto green_tex = std::make_shared<SolidColor>(Color(0.12f, 0.45f, 0.15f));

    auto white = std::make_shared<Lambertian>(white_tex);
    auto red = std::make_shared<Lambertian>(red_tex);
    auto green = std::make_shared<Lambertian>(green_tex);

    auto light_tex = std::make_shared<SolidColor>(Color(12.0f, 12.0f, 12.0f));
    auto light = std::make_shared<DiffuseLight>(light_tex);

    // Room dimensions (Cornell-style box)
    const float x0 = -1.0f;
    const float x1 = 1.0f;
    const float y0 = 0.0f;
    const float y1 = 2.0f;
    const float z0 = -3.0f;
    const float z1 = -1.0f;

    // Checker floor (y = y0)
    auto floor_even = std::make_shared<SolidColor>(Color(0.8f, 0.8f, 0.8f));
    auto floor_odd = std::make_shared<SolidColor>(Color(0.2f, 0.2f, 0.2f));
    auto floor_checker = std::make_shared<CheckerTexture>(floor_even, floor_odd, 4.0f);
    auto floor_mat = std::make_shared<Lambertian>(floor_checker);

    scene.objects.push_back(std::make_shared<Quad>(
        Vec3(x0, y0, z0),
        Vec3(x1 - x0, 0.0f, 0.0f),
        Vec3(0.0f, 0.0f, z1 - z0),
        floor_mat));

    // Ceiling (y = y1)
    scene.objects.push_back(std::make_shared<Quad>(
        Vec3(x0, y1, z0),
        Vec3(x1 - x0, 0.0f, 0.0f),
        Vec3(0.0f, 0.0f, z1 - z0),
        white));

    // Back wall (z = z0)
    scene.objects.push_back(std::make_shared<Quad>(
        Vec3(x0, y0, z0),
        Vec3(x1 - x0, 0.0f, 0.0f),
        Vec3(0.0f, y1 - y0, 0.0f),
        white));

    // Left wall (x = x0) - red
    scene.objects.push_back(std::make_shared<Quad>(
        Vec3(x0, y0, z1),
        Vec3(0.0f, 0.0f, z0 - z1),
        Vec3(0.0f, y1 - y0, 0.0f),
        red));

    // Right wall (x = x1) - green
    scene.objects.push_back(std::make_shared<Quad>(
        Vec3(x1, y0, z0),
        Vec3(0.0f, 0.0f, z1 - z0),
        Vec3(0.0f, y1 - y0, 0.0f),
        green));

    // Ceiling light quad (small emissive patch)
    const float lx0 = -0.5f;
    const float lx1 = 0.5f;
    const float lz0 = -2.5f;
    const float lz1 = -1.5f;
    const float ly = y1 - 0.01f;

    scene.objects.push_back(std::make_shared<Quad>(
        Vec3(lx0, ly, lz0),
        Vec3(lx1 - lx0, 0.0f, 0.0f),
        Vec3(0.0f, 0.0f, lz1 - lz0),
        light));

    // Objects in the box: diffuse, metal, checker, and glass spheres.
    auto diffuse_tex = std::make_shared<SolidColor>(Color(0.7f, 0.3f, 0.3f));
    auto diffuse = std::make_shared<Lambertian>(diffuse_tex);
    scene.objects.push_back(
        std::make_shared<Sphere>(Vec3(-0.5f, 0.5f, -2.2f), 0.5f, diffuse));

    auto metal_tex = std::make_shared<SolidColor>(Color(0.8f, 0.8f, 0.8f));
    auto metal = std::make_shared<Metal>(metal_tex, 0.1f);
    scene.objects.push_back(
        std::make_shared<Sphere>(Vec3(0.5f, 0.5f, -2.0f), 0.5f, metal));

    // Textured sphere using the same checker texture for visual verification.
    auto checker_sphere_mat = std::make_shared<Lambertian>(floor_checker);
    scene.objects.push_back(
        std::make_shared<Sphere>(Vec3(0.0f, 1.0f, -2.3f), 0.3f, checker_sphere_mat));

    // Glass sphere (Dielectric) to verify refraction.
    auto glass = std::make_shared<Dielectric>(1.5f);
    scene.objects.push_back(
        std::make_shared<Sphere>(Vec3(0.0f, 0.5f, -1.4f), 0.5f, glass));

    return scene;
}
