#pragma once

#include <memory>
#include <string>

#include "core/color.h"
#include "scene/material.h"
#include "scene/mesh.h"
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

inline Scene build_alpha_shadow_scene() {
    Scene scene;

    // Ground
    auto ground_mat = std::make_shared<Lambertian>(
        std::make_shared<SolidColor>(Color(0.8f, 0.8f, 0.8f)));
    scene.objects.push_back(std::make_shared<Quad>(
        Vec3(-5, 0, -5), Vec3(10, 0, 0), Vec3(0, 0, 10), ground_mat));

    // Semi-transparent checker sphere using AlphaCheckerTexture
    auto alpha_checker = std::make_shared<AlphaCheckerTexture>(8.0f);
    auto alpha_mat = std::make_shared<Lambertian>(alpha_checker);
    scene.objects.push_back(
        std::make_shared<Sphere>(Vec3(0.0f, 1.0f, 0.0f), 1.0f, alpha_mat));

    // Solid sphere behind to show shadow
    auto red_mat = std::make_shared<Lambertian>(
        std::make_shared<SolidColor>(Color(0.8f, 0.2f, 0.2f)));
    scene.objects.push_back(
        std::make_shared<Sphere>(Vec3(2.0f, 0.5f, 1.0f), 0.5f, red_mat));

    // Area light from above-left
    auto light_mat = std::make_shared<DiffuseLight>(
        std::make_shared<SolidColor>(Color(15, 15, 15)));
    auto light = std::make_shared<Quad>(
        Vec3(-2, 5, -2), Vec3(2, 0, 0), Vec3(0, 0, 2), light_mat);
    scene.objects.push_back(light);
    scene.lights.add_area_light(light);

    return scene;
}

inline Scene build_mesh_scene() {
    Scene scene;

    auto ground_mat = std::make_shared<Lambertian>(
        std::make_shared<SolidColor>(Color(0.75f, 0.75f, 0.75f)));
    scene.objects.push_back(std::make_shared<Quad>(
        Vec3(-6.0f, 0.0f, -6.0f), Vec3(12.0f, 0.0f, 0.0f), Vec3(0.0f, 0.0f, 12.0f), ground_mat));

    auto light_mat = std::make_shared<DiffuseLight>(
        std::make_shared<SolidColor>(Color(20.0f, 20.0f, 20.0f)));
    auto light = std::make_shared<Quad>(
        Vec3(-1.2f, 3.5f, -2.2f), Vec3(2.4f, 0.0f, 0.0f), Vec3(0.0f, 0.0f, 2.4f), light_mat);
    scene.objects.push_back(light);
    scene.lights.add_area_light(light);

    const float x = 0.5f;
    const float y = 0.5f;
    const float z = 0.5f;

    std::vector<Vec3> positions = {
        // +Z (front)
        Vec3(-x, -y, z), Vec3(x, -y, z), Vec3(x, y, z), Vec3(-x, y, z),
        // -Z (back)
        Vec3(-x, -y, -z), Vec3(-x, y, -z), Vec3(x, y, -z), Vec3(x, -y, -z),
        // +X (right)
        Vec3(x, -y, -z), Vec3(x, y, -z), Vec3(x, y, z), Vec3(x, -y, z),
        // -X (left)
        Vec3(-x, -y, z), Vec3(-x, y, z), Vec3(-x, y, -z), Vec3(-x, -y, -z),
        // +Y (top)
        Vec3(-x, y, z), Vec3(x, y, z), Vec3(x, y, -z), Vec3(-x, y, -z),
        // -Y (bottom)
        Vec3(-x, -y, -z), Vec3(x, -y, -z), Vec3(x, -y, z), Vec3(-x, -y, z),
    };

    std::vector<Vec3> normals = {
        // +Z
        Vec3(0.0f, 0.0f, 1.0f), Vec3(0.0f, 0.0f, 1.0f), Vec3(0.0f, 0.0f, 1.0f), Vec3(0.0f, 0.0f, 1.0f),
        // -Z
        Vec3(0.0f, 0.0f, -1.0f), Vec3(0.0f, 0.0f, -1.0f), Vec3(0.0f, 0.0f, -1.0f), Vec3(0.0f, 0.0f, -1.0f),
        // +X
        Vec3(1.0f, 0.0f, 0.0f), Vec3(1.0f, 0.0f, 0.0f), Vec3(1.0f, 0.0f, 0.0f), Vec3(1.0f, 0.0f, 0.0f),
        // -X
        Vec3(-1.0f, 0.0f, 0.0f), Vec3(-1.0f, 0.0f, 0.0f), Vec3(-1.0f, 0.0f, 0.0f), Vec3(-1.0f, 0.0f, 0.0f),
        // +Y
        Vec3(0.0f, 1.0f, 0.0f), Vec3(0.0f, 1.0f, 0.0f), Vec3(0.0f, 1.0f, 0.0f), Vec3(0.0f, 1.0f, 0.0f),
        // -Y
        Vec3(0.0f, -1.0f, 0.0f), Vec3(0.0f, -1.0f, 0.0f), Vec3(0.0f, -1.0f, 0.0f), Vec3(0.0f, -1.0f, 0.0f),
    };

    std::vector<Vec2> uvs = {
        // +Z
        Vec2(0.0f, 0.0f), Vec2(1.0f, 0.0f), Vec2(1.0f, 1.0f), Vec2(0.0f, 1.0f),
        // -Z
        Vec2(0.0f, 0.0f), Vec2(0.0f, 1.0f), Vec2(1.0f, 1.0f), Vec2(1.0f, 0.0f),
        // +X
        Vec2(0.0f, 0.0f), Vec2(0.0f, 1.0f), Vec2(1.0f, 1.0f), Vec2(1.0f, 0.0f),
        // -X
        Vec2(0.0f, 0.0f), Vec2(0.0f, 1.0f), Vec2(1.0f, 1.0f), Vec2(1.0f, 0.0f),
        // +Y
        Vec2(0.0f, 0.0f), Vec2(1.0f, 0.0f), Vec2(1.0f, 1.0f), Vec2(0.0f, 1.0f),
        // -Y
        Vec2(0.0f, 0.0f), Vec2(1.0f, 0.0f), Vec2(1.0f, 1.0f), Vec2(0.0f, 1.0f),
    };

    std::vector<std::uint32_t> indices;
    indices.reserve(6 * 6);
    for (std::uint32_t face = 0; face < 6; ++face) {
        const std::uint32_t base = face * 4;
        indices.push_back(base + 0);
        indices.push_back(base + 1);
        indices.push_back(base + 2);
        indices.push_back(base + 0);
        indices.push_back(base + 2);
        indices.push_back(base + 3);
    }

    MeshDataPtr cube_data =
        std::make_shared<MeshData>(std::move(positions), std::move(normals), std::move(uvs), std::move(indices));

    auto cube_tex = std::make_shared<ImageTexture>("../assets/textures/2k_earth_daymap.jpg");
    auto cube_normal = std::make_shared<NormalMapTexture>("../assets/textures/2k_earth_normal.png");
    auto cube_mat = std::make_shared<NormalMappedLambertian>(cube_tex, cube_normal, 2.0f);

    Transform t1 =
        Transform::translate(Vec3(-0.9f, 0.7f, -1.2f)) *
        Transform::rotate_y(0.6f) *
        Transform::uniform_scale(0.8f);
    Transform t2 =
        Transform::translate(Vec3(0.9f, 0.6f, -0.6f)) *
        Transform::rotate_y(-0.4f) *
        Transform::uniform_scale(0.65f);

    scene.objects.push_back(std::make_shared<Mesh>(cube_data, t1, cube_mat));
    scene.objects.push_back(std::make_shared<Mesh>(cube_data, t2, cube_mat));

    auto glass = std::make_shared<Dielectric>(1.5f);
    scene.objects.push_back(std::make_shared<Sphere>(Vec3(0.0f, 0.5f, 0.8f), 0.5f, glass));

    return scene;
}

inline Scene build_hotel_room_scene(const std::string& mural_texture_path = "../assets/textures/Starry_Night.jpg") {
    Scene scene;

    auto add_box = [&](const Vec3& bmin, const Vec3& bmax, const MaterialPtr& mat) {
        const float dx = bmax.x - bmin.x;
        const float dy = bmax.y - bmin.y;
        const float dz = bmax.z - bmin.z;

        scene.objects.push_back(std::make_shared<Quad>(
            Vec3(bmax.x, bmin.y, bmin.z), Vec3(0.0f, dy, 0.0f), Vec3(0.0f, 0.0f, dz), mat));  // +X
        scene.objects.push_back(std::make_shared<Quad>(
            Vec3(bmin.x, bmin.y, bmin.z), Vec3(0.0f, 0.0f, dz), Vec3(0.0f, dy, 0.0f), mat));  // -X
        scene.objects.push_back(std::make_shared<Quad>(
            Vec3(bmin.x, bmax.y, bmin.z), Vec3(0.0f, 0.0f, dz), Vec3(dx, 0.0f, 0.0f), mat));  // +Y
        scene.objects.push_back(std::make_shared<Quad>(
            Vec3(bmin.x, bmin.y, bmin.z), Vec3(dx, 0.0f, 0.0f), Vec3(0.0f, 0.0f, dz), mat));  // -Y
        scene.objects.push_back(std::make_shared<Quad>(
            Vec3(bmin.x, bmin.y, bmax.z), Vec3(dx, 0.0f, 0.0f), Vec3(0.0f, dy, 0.0f), mat));  // +Z
        scene.objects.push_back(std::make_shared<Quad>(
            Vec3(bmin.x, bmin.y, bmin.z), Vec3(0.0f, dy, 0.0f), Vec3(dx, 0.0f, 0.0f), mat));  // -Z
    };

    // -------------------------
    // Materials
    // -------------------------
    auto wall_mat = std::make_shared<Lambertian>(
        std::make_shared<SolidColor>(Color(0.46f, 0.55f, 0.66f)));

    auto ceiling_mat = std::make_shared<Lambertian>(
        std::make_shared<SolidColor>(Color(0.86f, 0.87f, 0.88f)));

    auto floor_tile_mat = std::make_shared<PrincipledBSDF>(
        std::make_shared<SolidColor>(Color(0.26f, 0.27f, 0.29f)),
        /*metallic_factor=*/0.0f,
        /*metallic_tex=*/nullptr,
        /*roughness_factor=*/0.18f,
        /*roughness_tex=*/nullptr,
        /*normal_map=*/nullptr,
        /*normal_strength=*/1.0f);

    auto wood_mat = std::make_shared<PrincipledBSDF>(
        std::make_shared<SolidColor>(Color(0.52f, 0.39f, 0.24f)),
        /*metallic_factor=*/0.0f,
        /*metallic_tex=*/nullptr,
        /*roughness_factor=*/0.32f,
        /*roughness_tex=*/nullptr,
        /*normal_map=*/nullptr,
        /*normal_strength=*/1.0f);

    auto bed_frame_mat = std::make_shared<PrincipledBSDF>(
        std::make_shared<SolidColor>(Color(0.14f, 0.14f, 0.15f)),
        /*metallic_factor=*/0.0f,
        /*metallic_tex=*/nullptr,
        /*roughness_factor=*/0.60f,
        /*roughness_tex=*/nullptr,
        /*normal_map=*/nullptr,
        /*normal_strength=*/1.0f);

    auto headboard_mat = std::make_shared<PrincipledBSDF>(
        std::make_shared<SolidColor>(Color(0.18f, 0.19f, 0.20f)),
        /*metallic_factor=*/0.0f,
        /*metallic_tex=*/nullptr,
        /*roughness_factor=*/0.75f,
        /*roughness_tex=*/nullptr,
        /*normal_map=*/nullptr,
        /*normal_strength=*/1.0f);

    auto bedding_mat = std::make_shared<PrincipledBSDF>(
        std::make_shared<SolidColor>(Color(0.90f, 0.90f, 0.90f)),
        /*metallic_factor=*/0.0f,
        /*metallic_tex=*/nullptr,
        /*roughness_factor=*/0.90f,
        /*roughness_tex=*/nullptr,
        /*normal_map=*/nullptr,
        /*normal_strength=*/1.0f);

    auto tv_body_mat = std::make_shared<PrincipledBSDF>(
        std::make_shared<SolidColor>(Color(0.07f, 0.07f, 0.07f)),
        /*metallic_factor=*/0.0f,
        /*metallic_tex=*/nullptr,
        /*roughness_factor=*/0.35f,
        /*roughness_tex=*/nullptr,
        /*normal_map=*/nullptr,
        /*normal_strength=*/1.0f);

    auto tv_screen_mat = std::make_shared<Lambertian>(
        std::make_shared<SolidColor>(Color(0.02f, 0.02f, 0.02f)));

    auto cabinet_blue_mat = std::make_shared<PrincipledBSDF>(
        std::make_shared<SolidColor>(Color(0.24f, 0.34f, 0.48f)),
        /*metallic_factor=*/0.0f,
        /*metallic_tex=*/nullptr,
        /*roughness_factor=*/0.55f,
        /*roughness_tex=*/nullptr,
        /*normal_map=*/nullptr,
        /*normal_strength=*/1.0f);

    auto metal_dark_mat = std::make_shared<Metal>(
        std::make_shared<SolidColor>(Color(0.06f, 0.06f, 0.065f)), 0.45f);

    auto glass_mat = std::make_shared<Dielectric>(1.5f);

    auto mural_tex = std::make_shared<ImageTexture>(
        mural_texture_path,
        ImageTexture::ColorSpace::sRGB,
        /*channel=*/-1,
        /*flip_v=*/true,
        ImageTexture::WrapMode::ClampToEdge,
        ImageTexture::WrapMode::ClampToEdge);
    auto mural_mat = std::make_shared<Lambertian>(mural_tex);

    // -------------------------
    // Room dimensions
    // -------------------------
    const float x0 = -4.0f;
    const float x1 =  4.0f;
    const float y0 =  0.0f;
    const float y1 =  3.0f;
    const float z0 = -9.0f;   // window wall
    const float z1 =  1.0f;   // camera side wall

    // Shell: floor / ceiling / left / right / entrance.
    scene.objects.push_back(std::make_shared<Quad>(
        Vec3(x0, y0, z0), Vec3(x1 - x0, 0.0f, 0.0f), Vec3(0.0f, 0.0f, z1 - z0), floor_tile_mat));
    scene.objects.push_back(std::make_shared<Quad>(
        Vec3(x0, y1, z0), Vec3(x1 - x0, 0.0f, 0.0f), Vec3(0.0f, 0.0f, z1 - z0), ceiling_mat));
    scene.objects.push_back(std::make_shared<Quad>(
        Vec3(x0, y0, z1), Vec3(0.0f, 0.0f, z0 - z1), Vec3(0.0f, y1 - y0, 0.0f), wall_mat));  // left
    scene.objects.push_back(std::make_shared<Quad>(
        Vec3(x1, y0, z0), Vec3(0.0f, 0.0f, z1 - z0), Vec3(0.0f, y1 - y0, 0.0f), wall_mat));  // right
    scene.objects.push_back(std::make_shared<Quad>(
        Vec3(x0, y0, z1), Vec3(x1 - x0, 0.0f, 0.0f), Vec3(0.0f, y1 - y0, 0.0f), wall_mat));  // entrance

    // -------------------------
    // Window wall (large opening + grid frame + glass)
    // -------------------------
    const float window_x0 = x0 + 0.30f;
    const float window_x1 = x1 - 0.30f;
    const float window_y0 = 0.55f;
    const float window_y1 = 2.75f;

    // Wall segments around the opening (on z = z0).
    scene.objects.push_back(std::make_shared<Quad>(
        Vec3(x0, y0, z0), Vec3(x1 - x0, 0.0f, 0.0f), Vec3(0.0f, window_y0 - y0, 0.0f), wall_mat)); // bottom
    scene.objects.push_back(std::make_shared<Quad>(
        Vec3(x0, window_y1, z0), Vec3(x1 - x0, 0.0f, 0.0f), Vec3(0.0f, y1 - window_y1, 0.0f), wall_mat)); // top
    scene.objects.push_back(std::make_shared<Quad>(
        Vec3(x0, window_y0, z0), Vec3(window_x0 - x0, 0.0f, 0.0f), Vec3(0.0f, window_y1 - window_y0, 0.0f), wall_mat)); // left strip
    scene.objects.push_back(std::make_shared<Quad>(
        Vec3(window_x1, window_y0, z0), Vec3(x1 - window_x1, 0.0f, 0.0f), Vec3(0.0f, window_y1 - window_y0, 0.0f), wall_mat)); // right strip

    // Frame thickness/depth.
    const float frame_t = 0.07f;
    const float frame_d = 0.06f;
    const float frame_z0 = z0;
    const float frame_z1 = z0 + frame_d;

    // Outer frame.
    add_box(Vec3(window_x0, window_y0, frame_z0), Vec3(window_x0 + frame_t, window_y1, frame_z1), metal_dark_mat); // left
    add_box(Vec3(window_x1 - frame_t, window_y0, frame_z0), Vec3(window_x1, window_y1, frame_z1), metal_dark_mat); // right
    add_box(Vec3(window_x0, window_y0, frame_z0), Vec3(window_x1, window_y0 + frame_t, frame_z1), metal_dark_mat); // bottom
    add_box(Vec3(window_x0, window_y1 - frame_t, frame_z0), Vec3(window_x1, window_y1, frame_z1), metal_dark_mat); // top

    // Inner mullions: 5 columns (4 vertical bars), 2 rows (1 horizontal bar).
    const int cols = 5;
    const float w = (window_x1 - window_x0);
    for (int c = 1; c < cols; ++c) {
        const float x = window_x0 + w * (static_cast<float>(c) / static_cast<float>(cols));
        add_box(Vec3(x - 0.5f * frame_t, window_y0, frame_z0),
                Vec3(x + 0.5f * frame_t, window_y1, frame_z1),
                metal_dark_mat);
    }
    const float mid_y = window_y0 + 0.50f * (window_y1 - window_y0);
    add_box(Vec3(window_x0, mid_y - 0.5f * frame_t, frame_z0),
            Vec3(window_x1, mid_y + 0.5f * frame_t, frame_z1),
            metal_dark_mat);

    // Glass quad (slightly inside).
    const float glass_z = z0 + 0.02f;
    scene.objects.push_back(std::make_shared<Quad>(
        Vec3(window_x0, window_y0, glass_z),
        Vec3(window_x1 - window_x0, 0.0f, 0.0f),
        Vec3(0.0f, window_y1 - window_y0, 0.0f),
        glass_mat));

    // -------------------------
    // Furniture
    // -------------------------

    // Bed (right side), oriented along X (faces left wall TV / right wall art).
    // Frame base
    add_box(Vec3(0.85f, 0.0f, -6.70f), Vec3(3.65f, 0.30f, -3.30f), bed_frame_mat);
    // Mattress
    add_box(Vec3(0.95f, 0.30f, -6.60f), Vec3(3.55f, 0.65f, -3.40f), bedding_mat);
    // Headboard (against right wall side, thickness on +X)
    add_box(Vec3(3.65f, 0.0f, -6.75f), Vec3(3.85f, 1.35f, -3.25f), headboard_mat);

    // Nightstands: near / far corners by the headboard (both close to right wall)
    add_box(Vec3(3.40f, 0.0f, -3.25f), Vec3(3.90f, 0.45f, -2.75f), wood_mat); // near (camera side)
    add_box(Vec3(3.40f, 0.0f, -7.25f), Vec3(3.90f, 0.45f, -6.75f), wood_mat); // far (window side)

    // Bedside lamps (warm)
    auto bedside_light_mat = std::make_shared<DiffuseLight>(
        std::make_shared<SolidColor>(Color(6.0f, 5.2f, 4.0f)));
    const float lamp_y = 0.78f;
    const float lamp_s = 0.18f;

    auto add_lamp = [&](float cx, float cz) {
        auto lamp = std::make_shared<Quad>(
            Vec3(cx - 0.5f * lamp_s, lamp_y, cz - 0.5f * lamp_s),
            Vec3(lamp_s, 0.0f, 0.0f),
            Vec3(0.0f, 0.0f, lamp_s),
            bedside_light_mat);
        scene.objects.push_back(lamp);
        scene.lights.add_area_light(lamp);
    };

    add_lamp(3.65f, -3.00f); // near nightstand
    add_lamp(3.65f, -7.00f); // far nightstand

    // TV on left wall + low cabinet (blue), facing into the room (+X).
    add_box(Vec3(x0 + 0.10f, 0.0f, -6.25f), Vec3(x0 + 0.85f, 0.32f, -3.85f), cabinet_blue_mat);

    // TV body (thin box, slightly inside the room)
    add_box(Vec3(x0 + 0.05f, 0.65f, -6.05f), Vec3(x0 + 0.15f, 1.45f, -4.05f), tv_body_mat);

    // TV screen (black for now), ensure normal points +X (into the room)
    scene.objects.push_back(std::make_shared<Quad>(
        Vec3(x0 + 0.151f, 0.70f, -6.00f),
        Vec3(0.0f, 0.72f, 0.0f),
        Vec3(0.0f, 0.0f, 1.90f),
        tv_screen_mat));

    // Floor lamp near left-front corner (realistic size)
    auto lamp_shade_mat = std::make_shared<Lambertian>(
        std::make_shared<SolidColor>(Color(0.86f, 0.86f, 0.84f)));
    auto lamp_pole_mat = metal_dark_mat;

    // Pole
    add_box(Vec3(-3.62f, 0.0f, 0.17f), Vec3(-3.58f, 1.42f, 0.21f), lamp_pole_mat);
    // Shade (small)
    add_box(Vec3(-3.68f, 1.40f, 0.11f), Vec3(-3.52f, 1.52f, 0.27f), lamp_shade_mat);

    // Light inside shade
    auto floorlamp_light_mat = std::make_shared<DiffuseLight>(
        std::make_shared<SolidColor>(Color(7.0f, 6.0f, 4.8f)));
    auto floorlamp = std::make_shared<Quad>(
        Vec3(-3.66f, 1.40f, 0.13f),
        Vec3(0.12f, 0.0f, 0.0f),
        Vec3(0.0f, 0.0f, 0.12f),
        floorlamp_light_mat);
    scene.objects.push_back(floorlamp);
    scene.lights.add_area_light(floorlamp);

    // Ceiling pendant light (center-ish)
    auto pendant_light_mat = std::make_shared<DiffuseLight>(
        std::make_shared<SolidColor>(Color(12.0f, 11.0f, 9.0f)));

    // Pendant rod: connect light to ceiling
    const float rod_r = 0.03f;
    const float rod_cx = 0.0f;
    const float rod_cz = -4.75f;
    add_box(Vec3(rod_cx - rod_r, 2.63f, rod_cz - rod_r),
            Vec3(rod_cx + rod_r, y1,      rod_cz + rod_r),
            metal_dark_mat);

    auto pendant = std::make_shared<Quad>(
        Vec3(-0.30f, 2.62f, -5.05f),
        Vec3(0.60f, 0.0f, 0.0f),
        Vec3(0.0f, 0.0f, 0.60f),
        pendant_light_mat);
    scene.objects.push_back(pendant);
    scene.lights.add_area_light(pendant);

    // Decorative wall art on right wall (normal points -X into the room)
    const float art_x = x1 - 0.01f;
    scene.objects.push_back(std::make_shared<Quad>(
        Vec3(art_x, 1.35f, -3.90f),
        Vec3(0.0f, 1.10f, 0.0f),
        Vec3(0.0f, 0.0f, -1.60f),
        mural_mat));

    // -------------------------
    // Baseboards (subtle)
    // -------------------------
    const float base_h = 0.10f;
    const float base_t = 0.02f;
    add_box(Vec3(x0, 0.0f, z0), Vec3(x0 + base_t, base_h, z1), wood_mat);
    add_box(Vec3(x1 - base_t, 0.0f, z0), Vec3(x1, base_h, z1), wood_mat);
    add_box(Vec3(x0, 0.0f, z1 - base_t), Vec3(x1, base_h, z1), wood_mat);
    add_box(Vec3(x0, 0.0f, z0), Vec3(x1, base_h, z0 + base_t), wood_mat);

    return scene;
}
