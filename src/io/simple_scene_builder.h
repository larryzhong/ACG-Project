#pragma once

#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "core/color.h"
#include "io/gltf_loader.h"
#include "io/obj_loader.h"
#include "scene/material.h"
#include "scene/mesh.h"
#include "scene/moving_sphere.h"
#include "scene/quad.h"
#include "scene/scene.h"
#include "scene/sphere.h"
#include "scene/texture.h"

namespace simple_scene_builder_detail {

inline bool compute_gltf_bounds(const std::vector<GltfMeshInstance>& meshes, AABB& out_bounds) {
    bool has_bounds = false;
    for (const auto& inst : meshes) {
        if (!inst.data) continue;
        const AABB inst_bounds = transform_aabb(inst.data->bounding_box(), inst.transform);
        if (!has_bounds) {
            out_bounds = inst_bounds;
            has_bounds = true;
        } else {
            out_bounds = surrounding_box(out_bounds, inst_bounds);
        }
    }
    return has_bounds;
}

inline bool add_scaled_gltf_to_scene(
    Scene& scene,
    const std::string& path,
    const Vec3& target_pos,
    float target_height,
    std::string* out_err,
    const GltfLoadOptions& options = {})
{
    std::vector<GltfMeshInstance> meshes;
    std::string gltf_err;
    if (!load_gltf_meshes(path, meshes, &gltf_err, options)) {
        if (out_err) *out_err = gltf_err;
        return false;
    }

    AABB bounds;
    if (!compute_gltf_bounds(meshes, bounds)) {
        if (out_err) *out_err = "glTF contained no mesh data.";
        return false;
    }

    const float bounds_h = bounds.max.y - bounds.min.y;
    const float scale = (bounds_h > 1e-6f) ? (target_height / bounds_h) : 1.0f;

    // Pivot at bottom-center.
    const Vec3 pivot(
        (bounds.min.x + bounds.max.x) * 0.5f,
        bounds.min.y,
        (bounds.min.z + bounds.max.z) * 0.5f);

    const Transform world_from_gltf =
        Transform::translate(target_pos) *
        Transform::uniform_scale(scale) *
        Transform::translate(-pivot);

    for (const auto& inst : meshes) {
        if (!inst.data || !inst.material) continue;
        scene.objects.push_back(std::make_shared<Mesh>(inst.data, world_from_gltf * inst.transform, inst.material));
    }

    return true;
}

} // namespace simple_scene_builder_detail

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
// Scene 1b: Principled BSDF Showcase
// ==========================================
inline Scene build_simple_scene_pbr() {
    Scene scene;

    // Cornell-style box to keep lighting controlled and make BRDF differences obvious.
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

    const float x0 = -1.0f, x1 = 1.0f;
    const float y0 = 0.0f, y1 = 2.0f;
    const float z_box_back = -2.0f;
    const float z_studio_back = 4.0f; // Wall behind camera

    // Floor (checker)
    auto floor_even = std::make_shared<SolidColor>(Color(0.8f, 0.8f, 0.8f));
    auto floor_odd = std::make_shared<SolidColor>(Color(0.3f, 0.3f, 0.3f));
    auto floor_checker = std::make_shared<CheckerTexture>(floor_even, floor_odd, 2.0f);
    auto floor_mat = std::make_shared<Lambertian>(floor_checker);

    scene.objects.push_back(std::make_shared<Quad>(
        Vec3(x0, y0, z_box_back),
        Vec3(x1 - x0, 0.0f, 0.0f),
        Vec3(0.0f, 0.0f, z_studio_back - z_box_back),
        floor_mat));

    // Ceiling
    scene.objects.push_back(std::make_shared<Quad>(
        Vec3(x0, y1, z_box_back),
        Vec3(x1 - x0, 0.0f, 0.0f),
        Vec3(0.0f, 0.0f, z_studio_back - z_box_back),
        white));

    // Back wall
    scene.objects.push_back(std::make_shared<Quad>(
        Vec3(x0, y0, z_box_back),
        Vec3(x1 - x0, 0.0f, 0.0f),
        Vec3(0.0f, y1 - y0, 0.0f),
        white));

    // Left wall (red)
    scene.objects.push_back(std::make_shared<Quad>(
        Vec3(x0, y0, z_studio_back),
        Vec3(0.0f, 0.0f, z_box_back - z_studio_back),
        Vec3(0.0f, y1 - y0, 0.0f),
        red));

    // Right wall (green)
    scene.objects.push_back(std::make_shared<Quad>(
        Vec3(x1, y0, z_box_back),
        Vec3(0.0f, 0.0f, z_studio_back - z_box_back),
        Vec3(0.0f, y1 - y0, 0.0f),
        green));

    // Studio back wall (behind camera)
    scene.objects.push_back(std::make_shared<Quad>(
        Vec3(x0, y0, z_studio_back),
        Vec3(x1 - x0, 0.0f, 0.0f),
        Vec3(0.0f, y1 - y0, 0.0f),
        white));

    // Lights
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

    auto studio_light = std::make_shared<Quad>(
        Vec3(lx0, ly, 1.8f),
        Vec3(lx1 - lx0, 0.0f, 0.0f),
        Vec3(0.0f, 0.0f, 0.6f),
        fill_light_mat);
    scene.objects.push_back(studio_light);
    scene.lights.add_area_light(studio_light);

    // --- Objects (Principled BSDF) ---
    // 1) Dielectric-like plastic: metallic=0, medium roughness
    auto plastic_base = std::make_shared<SolidColor>(Color(0.20f, 0.55f, 0.85f));
    auto plastic = std::make_shared<PrincipledBSDF>(
        plastic_base, 0.0f, nullptr, 0.30f, nullptr, nullptr, 1.0f);
    scene.objects.push_back(
        std::make_shared<Sphere>(Vec3(-0.55f, 0.35f, -0.85f), 0.35f, plastic));

    // 2) Metallic: metallic=1, smoother surface for a strong specular lobe
    auto metal_base = std::make_shared<SolidColor>(Color(0.95f, 0.78f, 0.55f)); // warm metal tint
    auto smooth_metal = std::make_shared<PrincipledBSDF>(
        metal_base, 1.0f, nullptr, 0.12f, nullptr, nullptr, 1.0f);
    scene.objects.push_back(
        std::make_shared<Sphere>(Vec3(0.55f, 0.35f, -0.95f), 0.35f, smooth_metal));

    // 3) Textured + normal-mapped Principled: makes shading normal + mipmaps visible
    auto earth_tex = std::make_shared<ImageTexture>("../assets/textures/2k_earth_daymap.jpg");
    auto earth_normal = std::make_shared<NormalMapTexture>("../assets/textures/2k_earth_normal.png");
    auto earth = std::make_shared<PrincipledBSDF>(
        earth_tex, 0.0f, nullptr, 0.65f, nullptr, earth_normal, 1.0f);
    scene.objects.push_back(
        std::make_shared<Sphere>(Vec3(0.0f, 0.55f, -1.45f), 0.30f, earth));

    // 4) Glass (kept as dielectric for refraction; contrasts with Principled metal/plastic)
    auto glass = std::make_shared<Dielectric>(1.5f);
    scene.objects.push_back(
        std::make_shared<Sphere>(Vec3(0.05f, 0.23f, -0.35f), 0.23f, glass));

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
    const float z1 =  10.0f;   // camera side wall

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
    add_box(Vec3(0.85f, 0.0f, -7.90f), Vec3(3.65f, 0.30f, -4.50f), bed_frame_mat);
    // Mattress
    add_box(Vec3(0.95f, 0.30f, -7.80f), Vec3(3.55f, 0.65f, -4.60f), bedding_mat);
    // Rug beside the bed
    auto rug_mat = std::make_shared<Lambertian>(
        std::make_shared<SolidColor>(Color(0.65f, 0.62f, 0.58f)));
    scene.objects.push_back(std::make_shared<Quad>(
        Vec3(-0.5f, 0.003f, -4.0f),
        Vec3(3.5f, 0.0f, 0.0f),
        Vec3(0.0f, 0.0f, 2.5f),
        rug_mat));
    // Pillow (OBJ)
    // Loads ../assets/models/Pillow.obj and auto-scales it to a reasonable size.
    {
        const std::string pillow_path = "../assets/models/Pillow.obj";
        std::vector<ObjMesh> pillow_meshes;
        std::string pillow_error;

        if (load_obj_meshes(pillow_path, pillow_meshes, &pillow_error)) {
            bool bbox_init = false;
            AABB pillow_bbox(Vec3(0.0f), Vec3(0.0f));

            for (const auto& m : pillow_meshes) {
                if (!m.data) continue;
                const AABB b = m.data->bounding_box();
                if (!bbox_init) {
                    pillow_bbox = b;
                    bbox_init = true;
                } else {
                    pillow_bbox = surrounding_box(pillow_bbox, b);
                }
            }

            if (bbox_init) {
                const Vec3 ext = pillow_bbox.max - pillow_bbox.min;
                const float horiz_max = (ext.x > ext.z) ? ext.x : ext.z;
                const float pillow_target_len = 0.85f; // world units
                const float scale = (horiz_max > 1e-12f) ? (pillow_target_len / horiz_max) : 1.0f;

                const float pi = 3.14159265358979323846f;
                const float rot_y = (ext.x > ext.z) ? (0.5f * pi) : 0.0f;

                // Make pillow stand up: rotate the longest axis to +Y.
                Transform upright = Transform::identity();
                if (ext.x >= ext.y && ext.x >= ext.z) {
                    upright = Transform::rotate_z(-0.5f * pi); // X -> Y
                } else if (ext.z >= ext.x && ext.z >= ext.y) {
                    upright = Transform::rotate_x(0.5f * pi);  // Z -> Y
                }

                // Lean slightly toward -X so it "rests" on the bedside table.
                const float lean = -0.20f;

                const Vec3 center_obj = (pillow_bbox.min + pillow_bbox.max) * 0.5f;
                const Transform base = Transform::rotate_z(lean) * upright * Transform::rotate_y(rot_y) *
                                      Transform::uniform_scale(scale);
                const Vec3 center_base = base.apply_point(center_obj);

                // Compute bbox after rotation+scale so we can place it reliably on the mattress.
                const Vec3 p000(pillow_bbox.min.x, pillow_bbox.min.y, pillow_bbox.min.z);
                const Vec3 p001(pillow_bbox.min.x, pillow_bbox.min.y, pillow_bbox.max.z);
                const Vec3 p010(pillow_bbox.min.x, pillow_bbox.max.y, pillow_bbox.min.z);
                const Vec3 p011(pillow_bbox.min.x, pillow_bbox.max.y, pillow_bbox.max.z);
                const Vec3 p100(pillow_bbox.max.x, pillow_bbox.min.y, pillow_bbox.min.z);
                const Vec3 p101(pillow_bbox.max.x, pillow_bbox.min.y, pillow_bbox.max.z);
                const Vec3 p110(pillow_bbox.max.x, pillow_bbox.max.y, pillow_bbox.min.z);
                const Vec3 p111(pillow_bbox.max.x, pillow_bbox.max.y, pillow_bbox.max.z);

                Vec3 tp[8] = {
                    base.apply_point(p000), base.apply_point(p001), base.apply_point(p010), base.apply_point(p011),
                    base.apply_point(p100), base.apply_point(p101), base.apply_point(p110), base.apply_point(p111),
                };
                Vec3 base_min = tp[0];
                for (int i = 1; i < 8; ++i) {
                    base_min.x = std::min(base_min.x, tp[i].x);
                    base_min.y = std::min(base_min.y, tp[i].y);
                    base_min.z = std::min(base_min.z, tp[i].z);
                }

                const float mattress_top_y = 0.65f;
                const float clearance = -0.06f;

                // Place pillow near the near (camera-side) nightstand and lean into it.
                const float desired_cx = 3.45f;
                const float desired_cz = -6.10f;

                const float tx = desired_cx - center_base.x;
                const float tz = desired_cz - center_base.z;
                const float ty = (mattress_top_y + clearance) - base_min.y;

                const Transform pillow_tr = Transform::translate(Vec3(tx, ty, tz)) * base;

                for (const auto& m : pillow_meshes) {
                    if (!m.data) continue;
                    auto obj = std::make_shared<Mesh>(m.data, pillow_tr, m.material);
                    scene.objects.push_back(obj);
                    if (m.emissive) {
                        scene.lights.add_area_light(obj);
                    }
                }
            }
        }
    }
    // Headboard (against right wall side, thickness on +X)
    add_box(Vec3(3.65f, 0.0f, -7.95f), Vec3(3.85f, 1.35f, -4.45f), headboard_mat);

    // Nightstands with legs (make legs readable at the default camera distance).
    const float nightstand_x0 = 3.40f;
    const float nightstand_x1 = 3.90f;
    const float nightstand_top_y = 0.45f;
    const float nightstand_leg_h = 0.24f;
    const float nightstand_leg_t = 0.07f;
    const float nightstand_leg_inset = 0.00f;

    auto leg_mat = std::make_shared<PrincipledBSDF>(
        std::make_shared<SolidColor>(Color(0.60f, 0.58f, 0.55f)),
        /*metallic_factor=*/1.0f,
        /*metallic_tex=*/nullptr,
        /*roughness_factor=*/0.35f,
        /*roughness_tex=*/nullptr,
        /*normal_map=*/nullptr,
        /*normal_strength=*/1.0f);

    auto add_nightstand = [&](float z0_ns, float z1_ns) {
        add_box(Vec3(nightstand_x0, nightstand_leg_h, z0_ns),
                Vec3(nightstand_x1, nightstand_top_y, z1_ns),
                wood_mat);

        const float lx0 = nightstand_x0 + nightstand_leg_inset;
        const float lx1 = lx0 + nightstand_leg_t;
        const float rx1 = nightstand_x1 - nightstand_leg_inset;
        const float rx0 = rx1 - nightstand_leg_t;

        const float bz0 = z0_ns + nightstand_leg_inset;
        const float bz1 = bz0 + nightstand_leg_t;
        const float fz1 = z1_ns - nightstand_leg_inset;
        const float fz0 = fz1 - nightstand_leg_t;

        add_box(Vec3(lx0, 0.0f, bz0), Vec3(lx1, nightstand_leg_h, bz1), leg_mat);
        add_box(Vec3(rx0, 0.0f, bz0), Vec3(rx1, nightstand_leg_h, bz1), leg_mat);
        add_box(Vec3(lx0, 0.0f, fz0), Vec3(lx1, nightstand_leg_h, fz1), leg_mat);
        add_box(Vec3(rx0, 0.0f, fz0), Vec3(rx1, nightstand_leg_h, fz1), leg_mat);
    };

    // Near nightstand (camera side)
    add_nightstand(/*z0_ns=*/-4.45f, /*z1_ns=*/-3.95f);

    // Far nightstand (window side)
    add_nightstand(/*z0_ns=*/-8.45f, /*z1_ns=*/-7.95f);

    // Bedside lamps (warm): shade + glowing bulb (sphere)
    auto bedside_light_mat = std::make_shared<DiffuseLight>(
        std::make_shared<SolidColor>(Color(6.0f, 5.2f, 4.0f)));

    auto bedside_shade_emission = std::make_shared<SolidColor>(Color(2.0f, 1.5f, 0.8f));
    auto bedside_shade_mat = std::make_shared<DiffuseLight>(bedside_shade_emission);

    auto make_open_frustum = [](int segments, float r0, float r1, float h) -> MeshDataPtr {
        if (segments < 3) segments = 3;
        if (r0 < 0.0f) r0 = 0.0f;
        if (r1 < 0.0f) r1 = 0.0f;
        if (h < 0.0f) h = 0.0f;

        std::vector<Vec3> positions;
        std::vector<Vec2> uvs;
        std::vector<std::uint32_t> indices;

        positions.reserve(static_cast<std::size_t>(2 * segments));
        uvs.reserve(static_cast<std::size_t>(2 * segments));
        indices.reserve(static_cast<std::size_t>(segments * 6));

        for (int i = 0; i < segments; ++i) {
            const float t = static_cast<float>(i) / static_cast<float>(segments);
            const float a = 2.0f * kPi * t;
            positions.push_back(Vec3(r0 * std::cos(a), 0.0f, r0 * std::sin(a)));
            uvs.push_back(Vec2(t, 0.0f));
        }
        for (int i = 0; i < segments; ++i) {
            const float t = static_cast<float>(i) / static_cast<float>(segments);
            const float a = 2.0f * kPi * t;
            positions.push_back(Vec3(r1 * std::cos(a), h, r1 * std::sin(a)));
            uvs.push_back(Vec2(t, 1.0f));
        }

        for (int i = 0; i < segments; ++i) {
            const int ni = (i + 1) % segments;
            const std::uint32_t b0 = static_cast<std::uint32_t>(i);
            const std::uint32_t b1 = static_cast<std::uint32_t>(ni);
            const std::uint32_t t0 = static_cast<std::uint32_t>(segments + i);
            const std::uint32_t t1 = static_cast<std::uint32_t>(segments + ni);

            indices.push_back(b0);
            indices.push_back(t0);
            indices.push_back(t1);

            indices.push_back(b0);
            indices.push_back(t1);
            indices.push_back(b1);
        }

        return std::make_shared<MeshData>(std::move(positions), std::vector<Vec3>{}, std::move(uvs), std::move(indices));
    };

    const float bedside_shade_y0 = 0.62f;
    const float bedside_shade_h = 0.30f;
    const float bedside_shade_r0 = 0.18f;
    const float bedside_shade_r1 = 0.11f;
    const int bedside_shade_segments = 28;
    const MeshDataPtr bedside_shade_data =
        make_open_frustum(bedside_shade_segments, bedside_shade_r0, bedside_shade_r1, bedside_shade_h);

    const float bedside_bulb_r = 0.05f;
    const float bedside_bulb_y = bedside_shade_y0 + 0.17f;

    auto add_bedside_lamp = [&](float cx, float cz) {
        const float base_w = 0.06f;
        const float base_h = bedside_shade_y0 - nightstand_top_y;
        add_box(Vec3(cx - base_w, nightstand_top_y, cz - base_w),
                Vec3(cx + base_w, bedside_shade_y0, cz + base_w),
                metal_dark_mat);

        const Transform shade_tr = Transform::translate(Vec3(cx, bedside_shade_y0, cz));
        scene.objects.push_back(std::make_shared<Mesh>(bedside_shade_data, shade_tr, bedside_shade_mat));

        auto bulb = std::make_shared<Sphere>(Vec3(cx, bedside_bulb_y, cz), bedside_bulb_r, bedside_light_mat);
        scene.objects.push_back(bulb);
        scene.lights.add_area_light(bulb);
    };

    add_bedside_lamp(3.65f, -4.20f); // near nightstand
    add_bedside_lamp(3.65f, -8.20f); // far nightstand

    // TV on left wall + low cabinet (blue), facing into the room (+X).
    add_box(Vec3(x0 + 0.10f, 0.0f, -6.85f), Vec3(x0 + 0.85f, 0.32f, -5.55f), cabinet_blue_mat);

    // TV (55 inch, 16:9 ratio: ~1.22m x 0.68m)
    const float tv_w = 1.22f;
    const float tv_h = 0.68f;
    const float tv_cy = 1.10f;
    const float tv_cz = -6.20f;

    // TV body (frame)
    add_box(Vec3(x0 + 0.05f, tv_cy - tv_h/2 - 0.04f, tv_cz - tv_w/2 - 0.04f),
            Vec3(x0 + 0.12f, tv_cy + tv_h/2 + 0.04f, tv_cz + tv_w/2 + 0.04f),
            tv_body_mat);

    // TV screen (16:9, normal points +X)
    scene.objects.push_back(std::make_shared<Quad>(
        Vec3(x0 + 0.121f, tv_cy - tv_h/2, tv_cz - tv_w/2),
        Vec3(0.0f, tv_h, 0.0f),
        Vec3(0.0f, 0.0f, tv_w),
        tv_screen_mat));
    
    // -------------------------
    // Bookshelf (left wall, near camera side of TV cabinet)
    // -------------------------
    // Bookshelf dimensions
    const float shelf_x0 = x0 + 0.10f;      // Against left wall
    const float shelf_x1 = x0 + 0.45f;      // Depth ~0.35m
    const float shelf_z0 = -4.20f;          // Start (window side)
    const float shelf_z1 = -2.60f;          // End (camera side)
    const float shelf_h = 1.50f;            // Total height
    const float shelf_thickness = 0.02f;   // Board thickness
    
    // Bookshelf frame material (dark wood)
    auto shelf_wood_mat = std::make_shared<PrincipledBSDF>(
        std::make_shared<SolidColor>(Color(0.22f, 0.15f, 0.10f)),
        0.0f, nullptr, 0.45f, nullptr, nullptr, 1.0f);

    // Vertical sides
    add_box(Vec3(shelf_x0, 0.0f, shelf_z0), Vec3(shelf_x1, shelf_h, shelf_z0 + shelf_thickness), shelf_wood_mat);  // Left side
    add_box(Vec3(shelf_x0, 0.0f, shelf_z1 - shelf_thickness), Vec3(shelf_x1, shelf_h, shelf_z1), shelf_wood_mat);  // Right side
    
    // Back panel
    add_box(Vec3(shelf_x0, 0.0f, shelf_z0), Vec3(shelf_x0 + shelf_thickness, shelf_h, shelf_z1), shelf_wood_mat);
    
    // Horizontal shelves (4 levels including top and bottom)
    const float shelf_levels[] = {0.0f, 0.38f, 0.76f, 1.14f, shelf_h - shelf_thickness};
    for (float level : shelf_levels) {
        add_box(Vec3(shelf_x0, level, shelf_z0), Vec3(shelf_x1, level + shelf_thickness, shelf_z1), shelf_wood_mat);
    }

    // -------------------------
    // Books on shelves
    // -------------------------
    auto make_book_mat = [](const Color& c) {
        return std::make_shared<PrincipledBSDF>(
            std::make_shared<SolidColor>(c),
            0.0f, nullptr, 0.7f, nullptr, nullptr, 1.0f);
    };

    // Book colors
    auto book_red = make_book_mat(Color(0.6f, 0.15f, 0.12f));
    auto book_blue = make_book_mat(Color(0.12f, 0.18f, 0.45f));
    auto book_green = make_book_mat(Color(0.15f, 0.35f, 0.18f));
    auto book_brown = make_book_mat(Color(0.4f, 0.25f, 0.15f));
    auto book_navy = make_book_mat(Color(0.1f, 0.12f, 0.25f));
    auto book_beige = make_book_mat(Color(0.75f, 0.70f, 0.60f));
    auto book_orange = make_book_mat(Color(0.7f, 0.35f, 0.1f));
    auto book_black = make_book_mat(Color(0.08f, 0.08f, 0.08f));

    // Helper to add a book (standing upright)
    // book_x: depth into shelf, book_z: position along shelf, book_h: height, book_w: width (thickness)
    auto add_book = [&](float z_pos, float level_y, float height, float width, float depth, const MaterialPtr& mat) {
        const float base_y = level_y + shelf_thickness + 0.002f;  // Slightly above shelf
        add_box(Vec3(shelf_x0 + 0.03f, base_y, z_pos),
                Vec3(shelf_x0 + 0.03f + depth, base_y + height, z_pos + width),
                mat);
    };

    // Shelf 1 (bottom, level_y = 0.0)
    add_book(-4.15f, 0.0f, 0.28f, 0.04f, 0.22f, book_brown);
    add_book(-4.10f, 0.0f, 0.30f, 0.03f, 0.20f, book_red);
    add_book(-4.06f, 0.0f, 0.26f, 0.035f, 0.21f, book_navy);
    add_book(-4.02f, 0.0f, 0.29f, 0.025f, 0.19f, book_beige);
    add_book(-3.98f, 0.0f, 0.27f, 0.04f, 0.22f, book_green);
    add_book(-3.55f, 0.0f, 0.32f, 0.05f, 0.24f, book_black);
    add_book(-3.49f, 0.0f, 0.28f, 0.03f, 0.20f, book_blue);
    add_book(-3.45f, 0.0f, 0.30f, 0.04f, 0.22f, book_orange);
    add_book(-3.05f, 0.0f, 0.25f, 0.035f, 0.18f, book_red);
    add_book(-3.00f, 0.0f, 0.28f, 0.04f, 0.21f, book_brown);
    add_book(-2.95f, 0.0f, 0.26f, 0.03f, 0.20f, book_navy);

    // Shelf 2 (level_y = 0.38)
    add_book(-4.15f, 0.38f, 0.26f, 0.035f, 0.20f, book_blue);
    add_book(-4.11f, 0.38f, 0.30f, 0.04f, 0.22f, book_green);
    add_book(-4.06f, 0.38f, 0.28f, 0.03f, 0.19f, book_beige);
    add_book(-4.02f, 0.38f, 0.25f, 0.045f, 0.21f, book_red);
    add_book(-3.60f, 0.38f, 0.29f, 0.03f, 0.20f, book_black);
    add_book(-3.56f, 0.38f, 0.27f, 0.04f, 0.22f, book_orange);
    add_book(-3.51f, 0.38f, 0.31f, 0.035f, 0.21f, book_navy);
    add_book(-3.10f, 0.38f, 0.26f, 0.04f, 0.20f, book_brown);
    add_book(-3.05f, 0.38f, 0.28f, 0.03f, 0.19f, book_blue);
    add_book(-3.01f, 0.38f, 0.24f, 0.04f, 0.22f, book_green);
    add_book(-2.96f, 0.38f, 0.29f, 0.035f, 0.20f, book_beige);

    // Shelf 3 (level_y = 0.76)
    add_book(-4.15f, 0.76f, 0.30f, 0.04f, 0.21f, book_orange);
    add_book(-4.10f, 0.76f, 0.27f, 0.03f, 0.20f, book_navy);
    add_book(-4.06f, 0.76f, 0.29f, 0.035f, 0.22f, book_brown);
    add_book(-3.70f, 0.76f, 0.25f, 0.04f, 0.19f, book_red);
    add_book(-3.65f, 0.76f, 0.28f, 0.03f, 0.21f, book_black);
    add_book(-3.61f, 0.76f, 0.26f, 0.04f, 0.20f, book_beige);
    add_book(-3.20f, 0.76f, 0.31f, 0.035f, 0.22f, book_blue);
    add_book(-3.15f, 0.76f, 0.27f, 0.04f, 0.20f, book_green);
    add_book(-3.10f, 0.76f, 0.29f, 0.03f, 0.21f, book_orange);

    // Shelf 4 (level_y = 1.14)
    add_book(-4.15f, 1.14f, 0.28f, 0.035f, 0.20f, book_green);
    add_book(-4.11f, 1.14f, 0.25f, 0.04f, 0.21f, book_red);
    add_book(-4.06f, 1.14f, 0.27f, 0.03f, 0.19f, book_navy);
    add_book(-3.50f, 1.14f, 0.26f, 0.04f, 0.22f, book_beige);
    add_book(-3.45f, 1.14f, 0.29f, 0.035f, 0.20f, book_brown);
    add_book(-3.40f, 1.14f, 0.24f, 0.03f, 0.21f, book_black);
    add_book(-3.00f, 1.14f, 0.28f, 0.04f, 0.20f, book_blue);
    add_book(-2.95f, 1.14f, 0.26f, 0.035f, 0.22f, book_orange);

    // -------------------------
    // Sofa (facing TV on left wall, against right wall)
    // -------------------------
    auto sofa_fabric_mat = std::make_shared<PrincipledBSDF>(
        std::make_shared<SolidColor>(Color(0.28f, 0.32f, 0.38f)),
        0.0f, nullptr, 0.85f, nullptr, nullptr, 1.0f);

    auto sofa_leg_mat = std::make_shared<PrincipledBSDF>(
        std::make_shared<SolidColor>(Color(0.15f, 0.12f, 0.08f)),
        0.0f, nullptr, 0.5f, nullptr, nullptr, 1.0f);

    const float sofa_w = 1.20f;
    const float sofa_d = 0.70f;
    const float sofa_seat_h = 0.38f;
    const float sofa_back_h = 0.70f;
    const float sofa_arm_w = 0.10f;
    const float sofa_leg_h = 0.10f;

    const float sofa_x1 = x1 - 0.05f;
    const float sofa_x0 = sofa_x1 - sofa_d;
    const float sofa_cz = -2.5f;
    const float sofa_z0 = sofa_cz - sofa_w / 2;
    const float sofa_z1 = sofa_cz + sofa_w / 2;

    add_box(Vec3(sofa_x0, sofa_leg_h, sofa_z0),
            Vec3(sofa_x1, sofa_seat_h, sofa_z1),
            sofa_fabric_mat);

    const float back_thickness = 0.15f;
    add_box(Vec3(sofa_x1 - back_thickness, sofa_seat_h, sofa_z0 + sofa_arm_w),
            Vec3(sofa_x1, sofa_back_h, sofa_z1 - sofa_arm_w),
            sofa_fabric_mat);

    add_box(Vec3(sofa_x0, sofa_seat_h, sofa_z0),
            Vec3(sofa_x1 - back_thickness, sofa_seat_h + 0.18f, sofa_z0 + sofa_arm_w),
            sofa_fabric_mat);

    add_box(Vec3(sofa_x0, sofa_seat_h, sofa_z1 - sofa_arm_w),
            Vec3(sofa_x1 - back_thickness, sofa_seat_h + 0.18f, sofa_z1),
            sofa_fabric_mat);

    const float leg_t = 0.04f;
    const float leg_inset = 0.06f;
    add_box(Vec3(sofa_x0 + leg_inset, 0.0f, sofa_z0 + leg_inset),
            Vec3(sofa_x0 + leg_inset + leg_t, sofa_leg_h, sofa_z0 + leg_inset + leg_t),
            sofa_leg_mat);
    add_box(Vec3(sofa_x0 + leg_inset, 0.0f, sofa_z1 - leg_inset - leg_t),
            Vec3(sofa_x0 + leg_inset + leg_t, sofa_leg_h, sofa_z1 - leg_inset),
            sofa_leg_mat);
    add_box(Vec3(sofa_x1 - leg_inset - leg_t, 0.0f, sofa_z0 + leg_inset),
            Vec3(sofa_x1 - leg_inset, sofa_leg_h, sofa_z0 + leg_inset + leg_t),
            sofa_leg_mat);
    add_box(Vec3(sofa_x1 - leg_inset - leg_t, 0.0f, sofa_z1 - leg_inset - leg_t),
            Vec3(sofa_x1 - leg_inset, sofa_leg_h, sofa_z1 - leg_inset),
            sofa_leg_mat);

    // -------------------------
    // Coffee table in front of sofa
    // -------------------------
    auto table_top_mat = std::make_shared<PrincipledBSDF>(
        std::make_shared<SolidColor>(Color(0.75f, 0.62f, 0.48f)),  // Light oak wood
        0.0f, nullptr, 0.75f, nullptr, nullptr, 1.0f);  // High roughness = no reflection

    auto table_wood_mat = std::make_shared<PrincipledBSDF>(
        std::make_shared<SolidColor>(Color(0.65f, 0.52f, 0.38f)),
        0.0f, nullptr, 0.70f, nullptr, nullptr, 1.0f);

    auto table_leg_mat = std::make_shared<PrincipledBSDF>(
        std::make_shared<SolidColor>(Color(0.60f, 0.48f, 0.35f)),
        0.0f, nullptr, 0.65f, nullptr, nullptr, 1.0f);

    const float table_h = 0.38f;
    const float table_top_thick = 0.03f;
    const float table_w = 0.90f;
    const float table_d = 0.50f;
    const float table_leg_t = 0.04f;
    const float table_leg_inset = 0.05f;

    const float table_cx = sofa_x0 - 0.35f - table_d / 2;
    const float table_cz = sofa_cz;

    const float table_x0 = table_cx - table_d / 2;
    const float table_x1 = table_cx + table_d / 2;
    const float table_z0 = table_cz - table_w / 2;
    const float table_z1 = table_cz + table_w / 2;

    // Wooden tabletop (no reflection)
    add_box(Vec3(table_x0, table_h - table_top_thick, table_z0),
            Vec3(table_x1, table_h, table_z1),
            table_top_mat);

    // Four legs
    const float leg_base_y = table_h - table_top_thick;
    add_box(Vec3(table_x0 + table_leg_inset, 0.0f, table_z0 + table_leg_inset),
            Vec3(table_x0 + table_leg_inset + table_leg_t, leg_base_y, table_z0 + table_leg_inset + table_leg_t),
            table_leg_mat);
    add_box(Vec3(table_x1 - table_leg_inset - table_leg_t, 0.0f, table_z0 + table_leg_inset),
            Vec3(table_x1 - table_leg_inset, leg_base_y, table_z0 + table_leg_inset + table_leg_t),
            table_leg_mat);
    add_box(Vec3(table_x0 + table_leg_inset, 0.0f, table_z1 - table_leg_inset - table_leg_t),
            Vec3(table_x0 + table_leg_inset + table_leg_t, leg_base_y, table_z1 - table_leg_inset),
            table_leg_mat);
    add_box(Vec3(table_x1 - table_leg_inset - table_leg_t, 0.0f, table_z1 - table_leg_inset - table_leg_t),
            Vec3(table_x1 - table_leg_inset, leg_base_y, table_z1 - table_leg_inset),
            table_leg_mat);

    // Lower shelf
    const float tbl_shelf_y = 0.12f;
    const float tbl_shelf_thick = 0.02f;
    add_box(Vec3(table_x0 + 0.06f, tbl_shelf_y, table_z0 + 0.06f),
            Vec3(table_x1 - 0.06f, tbl_shelf_y + tbl_shelf_thick, table_z1 - 0.06f),
            table_wood_mat);

    // Coffee cup on the coffee table
    {
        std::string cup_err;
        const std::string cup_path = "../assets/cup.gltf";

        // Place near the center of the tabletop, slightly above to avoid z-fighting.
        const float tabletop_y = table_h;
        const Vec3 cup_pos(table_cx, tabletop_y + 0.001f, table_cz);

        // Typical cup height in this scene scale. (2x)
        (void)simple_scene_builder_detail::add_scaled_gltf_to_scene(scene, cup_path, cup_pos, 0.14f, &cup_err);
    }
    
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

    // Ceiling pendant light (center-ish): shade + glowing bulb (sphere)
    const float rod_r = 0.03f;
    const float rod_cx = 0.0f;
    const float rod_cz = -4.75f;

    auto pendant_shade_emission = std::make_shared<SolidColor>(Color(4.0f, 3.2f, 2.0f));
    auto pendant_shade_mat = std::make_shared<DiffuseLight>(pendant_shade_emission);

    const float pendant_shade_y0 = 2.05f;
    const float pendant_shade_h  = 0.55f;
    const float pendant_shade_r0 = 0.28f;
    const float pendant_shade_r1 = 0.08f;
    const int pendant_shade_segments = 36;
    const MeshDataPtr pendant_shade_data =
        make_open_frustum(pendant_shade_segments, pendant_shade_r0, pendant_shade_r1, pendant_shade_h);

    const float pendant_shade_y1 = pendant_shade_y0 + pendant_shade_h;

    // Pendant rod: connect shade to ceiling
    add_box(Vec3(rod_cx - rod_r, pendant_shade_y1, rod_cz - rod_r),
            Vec3(rod_cx + rod_r, y1,              rod_cz + rod_r),
            metal_dark_mat);

    const Transform pendant_shade_tr = Transform::translate(Vec3(rod_cx, pendant_shade_y0, rod_cz));
    scene.objects.push_back(std::make_shared<Mesh>(pendant_shade_data, pendant_shade_tr, pendant_shade_mat));

    const float cap_thickness = 0.02f;
    add_box(Vec3(rod_cx - pendant_shade_r1 - 0.02f, pendant_shade_y1, rod_cz - pendant_shade_r1 - 0.02f),
            Vec3(rod_cx + pendant_shade_r1 + 0.02f, pendant_shade_y1 + cap_thickness, rod_cz + pendant_shade_r1 + 0.02f),
            metal_dark_mat);
    auto pendant_light_mat = std::make_shared<DiffuseLight>(
        std::make_shared<SolidColor>(Color(32.0f, 30.0f, 24.0f)));
    const float pendant_bulb_r = 0.07f;
    const float pendant_bulb_y = pendant_shade_y0 + 0.20f;
    auto pendant_bulb = std::make_shared<Sphere>(
        Vec3(rod_cx, pendant_bulb_y, rod_cz),
        pendant_bulb_r,
        pendant_light_mat);
    scene.objects.push_back(pendant_bulb);
    scene.lights.add_area_light(pendant_bulb);

    // Decorative wall art on right wall (normal points -X into the room)
    const float art_x = x1 - 0.01f;
    scene.objects.push_back(std::make_shared<Quad>(
        Vec3(art_x, 1.35f, -6.70f),
        Vec3(0.0f, 0.0f, 1.60f),
        Vec3(0.0f, 1.10f, 0.0f),
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
