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

// ==========================================
// Scene 2: Depth of Field Test (Improved)
// ==========================================
inline Scene build_dof_scene() {
    Scene scene;

    auto checker = std::make_shared<CheckerTexture>(
        std::make_shared<SolidColor>(Color(0.2f, 0.3f, 0.1f)),
        std::make_shared<SolidColor>(Color(0.9f, 0.9f, 0.9f)),
        2.0f
    );
    auto ground_mat = std::make_shared<Lambertian>(checker);
    
    scene.objects.push_back(std::make_shared<Quad>(
        Vec3(-50, 0, -50), Vec3(100, 0, 0), Vec3(0, 0, 100), ground_mat
    ));

    auto material_left   = std::make_shared<Dielectric>(1.5f);
    auto material_center = std::make_shared<Lambertian>(std::make_shared<SolidColor>(Color(0.1f, 0.2f, 0.5f)));
    auto material_right  = std::make_shared<Metal>(std::make_shared<SolidColor>(Color(0.8f, 0.6f, 0.2f)), 0.0f);

    scene.objects.push_back(std::make_shared<Sphere>(Vec3(-1.0f, 0.5f, -1.0f), 0.5f, material_left));
    
    scene.objects.push_back(std::make_shared<Sphere>(Vec3(0.0f, 0.5f, -2.0f), 0.5f, material_center));
    
    scene.objects.push_back(std::make_shared<Sphere>(Vec3(1.0f, 0.5f, -3.0f), 0.5f, material_right));

    auto light_mat = std::make_shared<DiffuseLight>(std::make_shared<SolidColor>(Color(10, 10, 10)));
    auto light = std::make_shared<Quad>(Vec3(-2, 4, -1), Vec3(4, 0, 0), Vec3(0, 0, 4), light_mat);
    scene.objects.push_back(light);
    scene.lights.add_area_light(light);

    return scene;
}

// ==========================================
// Scene 3: Motion Blur Showcase
// ==========================================
inline Scene build_motion_blur_scene() {
    Scene scene;

    auto orange = std::make_shared<Lambertian>(std::make_shared<SolidColor>(Color(0.8f, 0.3f, 0.1f)));
    auto white  = std::make_shared<Lambertian>(std::make_shared<SolidColor>(Color(0.8f, 0.8f, 0.8f)));

    scene.objects.push_back(std::make_shared<Sphere>(Vec3(0.0f, -1000.0f, 0.0f), 1000.0f, white));

    scene.objects.push_back(std::make_shared<Sphere>(Vec3(-2.0f, 1.0f, 0.0f), 1.0f, orange));

    auto moving_mat = std::make_shared<Lambertian>(std::make_shared<SolidColor>(Color(0.2f, 0.4f, 0.8f)));
    Vec3 center0(2.0f, 1.0f, 0.0f);
    Vec3 center1(2.0f, 2.5f, 0.0f);
    scene.objects.push_back(std::make_shared<MovingSphere>(center0, center1, 0.0f, 1.0f, 1.0f, moving_mat));

    auto light_mat = std::make_shared<DiffuseLight>(std::make_shared<SolidColor>(Color(7, 7, 7)));
    auto light = std::make_shared<Quad>(Vec3(-3, 5, -3), Vec3(6, 0, 0), Vec3(0, 0, 6), light_mat);
    scene.objects.push_back(light);
    scene.lights.add_area_light(light);

    return scene;
}

// ==========================================
// Scene 4: Texture & Shadows
// ==========================================
inline Scene build_texture_scene() {
    Scene scene;

    auto checker = std::make_shared<CheckerTexture>(
        std::make_shared<SolidColor>(Color(0.2f, 0.3f, 0.1f)),
        std::make_shared<SolidColor>(Color(0.9f, 0.9f, 0.9f)),
        10.0f // Scale
    );
    auto ground_mat = std::make_shared<Lambertian>(checker);

    scene.objects.push_back(std::make_shared<Sphere>(Vec3(0.0f, -1000.0f, 0.0f), 1000.0f, ground_mat));

    auto metal = std::make_shared<Metal>(std::make_shared<SolidColor>(Color(0.8f, 0.8f, 0.8f)), 0.1f);
    scene.objects.push_back(std::make_shared<Sphere>(Vec3(0.0f, 2.0f, 0.0f), 2.0f, metal));

    auto light_mat = std::make_shared<DiffuseLight>(std::make_shared<SolidColor>(Color(15, 15, 15)));
    
    auto l1 = std::make_shared<Sphere>(Vec3(-3.0f, 5.0f, 3.0f), 1.0f, light_mat);
    scene.objects.push_back(l1);
    scene.lights.add_area_light(l1);
    
    auto blue_light = std::make_shared<DiffuseLight>(std::make_shared<SolidColor>(Color(5, 5, 20)));
    auto l2 = std::make_shared<Sphere>(Vec3(3.0f, 4.0f, 1.0f), 0.8f, blue_light);
    scene.objects.push_back(l2);
    scene.lights.add_area_light(l2);

    return scene;
}

// ==========================================
// Scene 5: Random Balls (BVH Test)
// ==========================================
inline Scene build_random_scene() {
    Scene scene;
    RNG rng(12345);

    // Ground
    auto checker = std::make_shared<CheckerTexture>(
        std::make_shared<SolidColor>(Color(0.2f, 0.3f, 0.1f)),
        std::make_shared<SolidColor>(Color(0.9f, 0.9f, 0.9f)), 
        2.0f
    );
    scene.objects.push_back(std::make_shared<Sphere>(Vec3(0.0f, -1000.0f, 0.0f), 1000.0f, std::make_shared<Lambertian>(checker)));

    // Random small spheres
    for (int a = -11; a < 11; a++) {
        for (int b = -11; b < 11; b++) {
            float choose_mat = rng.uniform();
            Vec3 center(a + 0.9f * rng.uniform(), 0.2f, b + 0.9f * rng.uniform());

            if ((center - Vec3(4, 0.2f, 0)).length() > 0.9f) {
                std::shared_ptr<Material> sphere_material;

                if (choose_mat < 0.7f) {
                    // diffuse
                    Color albedo = random_in_unit_sphere(rng) * random_in_unit_sphere(rng);
                    sphere_material = std::make_shared<Lambertian>(std::make_shared<SolidColor>(albedo));
                    scene.objects.push_back(std::make_shared<Sphere>(center, 0.2f, sphere_material));
                } else if (choose_mat < 0.90f) {
                    // metal
                    Color albedo = Color(0.5f, 0.5f, 0.5f) + 0.5f * random_in_unit_sphere(rng);
                    float fuzz = rng.uniform() * 0.5f;
                    sphere_material = std::make_shared<Metal>(std::make_shared<SolidColor>(albedo), fuzz);
                    scene.objects.push_back(std::make_shared<Sphere>(center, 0.2f, sphere_material));
                } else {
                    // glass
                    sphere_material = std::make_shared<Dielectric>(1.5f);
                    scene.objects.push_back(std::make_shared<Sphere>(center, 0.2f, sphere_material));
                }
            }
        }
    }

    auto mat1 = std::make_shared<Dielectric>(1.5f);
    scene.objects.push_back(std::make_shared<Sphere>(Vec3(0.0f, 1.0f, 0.0f), 1.0f, mat1));

    auto mat2 = std::make_shared<Lambertian>(std::make_shared<SolidColor>(Color(0.4f, 0.2f, 0.1f)));
    scene.objects.push_back(std::make_shared<Sphere>(Vec3(-4.0f, 1.0f, 0.0f), 1.0f, mat2));

    auto mat3 = std::make_shared<Metal>(std::make_shared<SolidColor>(Color(0.7f, 0.6f, 0.5f)), 0.0f);
    scene.objects.push_back(std::make_shared<Sphere>(Vec3(4.0f, 1.0f, 0.0f), 1.0f, mat3));

    // Sun light (Area light to simulate sun)
    auto light_mat = std::make_shared<DiffuseLight>(std::make_shared<SolidColor>(Color(10, 10, 10)));
    auto light = std::make_shared<Sphere>(Vec3(0.0f, 8.0f, 0.0f), 1.0f, light_mat);
    scene.objects.push_back(light);
    scene.lights.add_area_light(light);

    return scene;
}

inline Scene build_solar_system_scene() {
    Scene scene;

    auto key_light_mat = std::make_shared<DiffuseLight>(std::make_shared<SolidColor>(Color(200.0f, 195.0f, 190.0f)));
    auto key_light = std::make_shared<Sphere>(Vec3(-30.0f, 10.0f, 30.0f), 4.0f, key_light_mat);
    scene.objects.push_back(key_light);
    scene.lights.add_area_light(key_light);

    // auto sun_tex = std::make_shared<ImageTexture>("../assets/textures/2k_sun.jpg");
    // auto sun_visual_mat = std::make_shared<DiffuseLight>(sun_tex);
    // scene.objects.push_back(std::make_shared<Sphere>(Vec3(-8.0f, 5.0f, -15.0f), 2.0f, sun_visual_mat));

    auto fill_light_mat = std::make_shared<DiffuseLight>(std::make_shared<SolidColor>(Color(0.05f, 0.05f, 0.1f)));
    scene.objects.push_back(std::make_shared<Quad>(Vec3(-50, -50, 50), Vec3(100, 0, 0), Vec3(0, 100, 0), fill_light_mat));

    auto earth_tex = std::make_shared<ImageTexture>("../assets/textures/2k_earth_daymap.jpg");
    auto earth_normal = std::make_shared<NormalMapTexture>("../assets/textures/2k_earth_normal.png");
    auto earth_mat = std::make_shared<NormalMappedLambertian>(earth_tex, earth_normal, 10.0f);
    scene.objects.push_back(std::make_shared<Sphere>(Vec3(0.0f, 0.0f, 0.0f), 1.0f, earth_mat));

    auto jupiter_tex = std::make_shared<ImageTexture>("../assets/textures/2k_jupiter.jpg");
    auto jupiter_mat = std::make_shared<Lambertian>(jupiter_tex);
    scene.objects.push_back(std::make_shared<Sphere>(Vec3(15.0f, 2.0f, -25.0f), 9.0f, jupiter_mat));

    auto moon_tex = std::make_shared<ImageTexture>("../assets/textures/2k_moon.jpg");
    auto moon_mat = std::make_shared<Lambertian>(moon_tex);
    scene.objects.push_back(std::make_shared<Sphere>(Vec3(-2.2f, -0.8f, 3.5f), 0.28f, moon_mat));

    RNG rng(999);
    auto star_mat = std::make_shared<DiffuseLight>(std::make_shared<SolidColor>(Color(2.0f, 2.0f, 2.0f))); 
    for (int i = 0; i < 300; ++i) {
        Vec3 pos = random_in_unit_sphere(rng);
        if (pos.length() < 0.2f) continue;
        pos = normalize(pos) * 80.0f;
        float radius = 0.05f + rng.uniform() * 0.05f; 
        scene.objects.push_back(std::make_shared<Sphere>(pos, radius, star_mat));
    }

    return scene;
}
