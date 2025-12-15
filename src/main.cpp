#include <cstdlib>
#include <iostream>
#include <string>
#include <iomanip>
#include <sstream>

#include "camera/camera.h"
#include "core/color.h"
#include "core/ray.h"
#include "core/rng.h"
#include "core/vec3.h"
#include "io/image_io.h"
#include "io/simple_scene_builder.h"
#include "render/film.h"
#include "render/path_tracer.h"
#include "scene/scene.h"

#define STB_IMAGE_IMPLEMENTATION
#include "io/stb_image.h"

int main(int argc, char** argv) {
    int width = 400;
    int height = 225;
    int samples_per_pixel = 16;
    int max_depth = 5;
    float aperture = 0.0f;
    float focus_dist = 0.0f;
    float shutter_open = 0.0f;
    float shutter_close = 1.0f;
    std::string output = "basic_materials.ppm";
    std::string scene_name = "simple";
    bool turntable_mode = false;
    int turntable_frames = 60;
    float turntable_radius = 0.0f;
    float turntable_height = 0.0f;
    Vec3 turntable_center(0.0f, 0.0f, 0.0f);

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--output" && i + 1 < argc) {
            output = argv[++i];
        } else if (arg == "--width" && i + 1 < argc) {
            width = std::atoi(argv[++i]);
        } else if (arg == "--height" && i + 1 < argc) {
            height = std::atoi(argv[++i]);
        } else if (arg == "--scene" && i + 1 < argc) {
            scene_name = argv[++i];
        } else if (arg == "--spp" && i + 1 < argc) {
            samples_per_pixel = std::atoi(argv[++i]);
        } else if (arg == "--max-depth" && i + 1 < argc) {
            max_depth = std::atoi(argv[++i]);
        } else if (arg == "--aperture" && i + 1 < argc) {
            aperture = static_cast<float>(std::atof(argv[++i]));
        } else if (arg == "--focus-dist" && i + 1 < argc) {
            focus_dist = static_cast<float>(std::atof(argv[++i]));
        } else if (arg == "--shutter-open" && i + 1 < argc) {
            shutter_open = static_cast<float>(std::atof(argv[++i]));
        } else if (arg == "--shutter-close" && i + 1 < argc) {
            shutter_close = static_cast<float>(std::atof(argv[++i]));
        } else if (arg == "--turntable") {
            turntable_mode = true;
        } else if (arg == "--frames" && i + 1 < argc) {
            turntable_frames = std::atoi(argv[++i]);
        } else if (arg == "--turntable-radius" && i + 1 < argc) {
            turntable_radius = static_cast<float>(std::atof(argv[++i]));
        } else if (arg == "--turntable-height" && i + 1 < argc) {
            turntable_height = static_cast<float>(std::atof(argv[++i]));
        } else if (arg == "--turntable-center" && i + 3 < argc) {
            float cx = static_cast<float>(std::atof(argv[++i]));
            float cy = static_cast<float>(std::atof(argv[++i]));
            float cz = static_cast<float>(std::atof(argv[++i]));
            turntable_center = Vec3(cx, cy, cz);
        }
    }

    if (width <= 0 || height <= 0 || samples_per_pixel <= 0 || max_depth <= 0) {
        std::cerr << "Invalid render parameters.\n";
        return 1;
    }

    Film film(width, height);
    Scene scene;
    
    CameraSettings cam_settings;
    cam_settings.aspect_ratio = static_cast<float>(width) / static_cast<float>(height);
    cam_settings.vertical_fov_deg = 40.0f;
    
    if (scene_name == "simple") {
        scene = build_simple_scene_basic();
        cam_settings.look_from = Vec3(0.0f, 1.0f, 2.5f);
        cam_settings.look_at = Vec3(0.0f, 1.0f, -1.0f);
    } 
    else if (scene_name == "dof") {
        scene = build_dof_scene();
        cam_settings.look_from = Vec3(0.0f, 2.0f, 2.5f);
        cam_settings.look_at = Vec3(0.0f, 0.5f, -2.0f); 
        cam_settings.vertical_fov_deg = 35.0f;
    }
    else if (scene_name == "motion") {
        scene = build_motion_blur_scene();
        cam_settings.look_from = Vec3(0.0f, 2.0f, 8.0f);
        cam_settings.look_at = Vec3(0.0f, 1.5f, 0.0f);
        cam_settings.vertical_fov_deg = 30.0f;
    }
    else if (scene_name == "texture") {
        scene = build_texture_scene();
        cam_settings.look_from = Vec3(0.0f, 3.0f, 8.0f);
        cam_settings.look_at = Vec3(0.0f, 2.0f, 0.0f);
        cam_settings.vertical_fov_deg = 35.0f;
    }
    else if (scene_name == "random") {
        scene = build_random_scene();
        cam_settings.look_from = Vec3(13.0f, 2.0f, 3.0f);
        cam_settings.look_at = Vec3(0.0f, 0.0f, 0.0f);
        cam_settings.vertical_fov_deg = 20.0f;
        if (aperture == 0.0f) aperture = 0.1f;
        if (focus_dist < 0.0f) focus_dist = 10.0f;
    }
    else if (scene_name == "solar") {
        scene = build_solar_system_scene();
        cam_settings.look_from = Vec3(0.0f, 0.0f, 8.5f);
        cam_settings.look_at = Vec3(0.0f, 0.0f, 0.0f);
        cam_settings.vertical_fov_deg = 35.0f;
        cam_settings.aperture = 0.05f;
        cam_settings.focus_dist = 8.5f;
    }
    else if (scene_name == "alpha") {
        scene = build_alpha_shadow_scene();
        cam_settings.look_from = Vec3(0.0f, 3.0f, 6.0f);
        cam_settings.look_at = Vec3(0.0f, 1.0f, 0.0f);
        cam_settings.vertical_fov_deg = 40.0f;
    }
    else if (scene_name == "mesh") {
        scene = build_mesh_scene();
        cam_settings.look_from = Vec3(0.0f, 1.6f, 3.5f);
        cam_settings.look_at = Vec3(0.0f, 0.7f, -1.0f);
        cam_settings.vertical_fov_deg = 35.0f;
    }
    else {
        std::cerr << "Unknown scene name: " << scene_name << "\n";
        return 1;
    }

    scene.build_bvh();

    cam_settings.up = Vec3(0.0f, 1.0f, 0.0f);
    cam_settings.aperture = aperture;
    cam_settings.t0 = shutter_open;
    cam_settings.t1 = shutter_close;

    if (focus_dist <= 0.0f) {
        const Vec3 diff = cam_settings.look_from - cam_settings.look_at;
        cam_settings.focus_dist = diff.length();
    } else {
        cam_settings.focus_dist = focus_dist;
    }

    if (!turntable_mode) {
        Camera camera(cam_settings);
        PathTracer integrator(max_depth);

        std::cout << "Rendering scene: " << scene_name << "\n";
        render_image(scene, camera, integrator, film, samples_per_pixel);

        if (output.size() >= 4 && output.substr(output.size() - 4) == ".png") {
            write_png(output, film);
            std::cout << "Wrote PNG image to " << output << "\n";
        } else {
            write_ppm(output, film);
            std::cout << "Wrote PPM image to " << output << "\n";
        }
    } else {
        PathTracer integrator(max_depth);

        if (turntable_radius <= 0.0f) {
            Vec3 diff = cam_settings.look_from - cam_settings.look_at;
            turntable_radius = std::sqrt(diff.x * diff.x + diff.z * diff.z);
        }
        if (turntable_height == 0.0f) {
            turntable_height = cam_settings.look_from.y;
        }
        if (turntable_center.length_squared() == 0.0f) {
            turntable_center = cam_settings.look_at;
        }

        std::cout << "Turntable mode: " << turntable_frames << " frames\n";
        std::cout << "  Center: (" << turntable_center.x << ", " 
                  << turntable_center.y << ", " << turntable_center.z << ")\n";
        std::cout << "  Radius: " << turntable_radius << ", Height: " << turntable_height << "\n";

        std::string base_name = output;
        std::string extension = ".png";
        if (output.size() >= 4) {
            std::string ext = output.substr(output.size() - 4);
            if (ext == ".png" || ext == ".ppm") {
                extension = ext;
                base_name = output.substr(0, output.size() - 4);
            }
        }

        const float angle_step = 2.0f * kPi / static_cast<float>(turntable_frames);

        for (int frame = 0; frame < turntable_frames; ++frame) {
            float angle = frame * angle_step;

            float x = turntable_center.x + turntable_radius * std::sin(angle);
            float z = turntable_center.z + turntable_radius * std::cos(angle);

            cam_settings.look_from = Vec3(x, turntable_height, z);
            cam_settings.look_at = turntable_center;

            Vec3 diff = cam_settings.look_from - cam_settings.look_at;
            cam_settings.focus_dist = diff.length();

            Camera camera(cam_settings);

            Film frame_film(width, height);

            std::cout << "\rRendering frame " << (frame + 1) << "/" << turntable_frames 
                      << " (angle: " << static_cast<int>(angle * 180.0f / kPi) << " deg)" << std::flush;

            render_image(scene, camera, integrator, frame_film, samples_per_pixel);

            std::ostringstream filename;
            filename << base_name << "_" << std::setfill('0') << std::setw(4) << frame << extension;

            if (extension == ".png") {
                write_png(filename.str(), frame_film);
            } else {
                write_ppm(filename.str(), frame_film);
            }
        }

        std::cout << "\nTurntable rendering complete!\n";
        std::cout << "To create GIF, run:\n";
        std::cout << "  ffmpeg -framerate 30 -i " << base_name << "_%04d.png -vf \"fps=30,scale=480:-1:flags=lanczos\" " << base_name << ".gif\n";
    }

    return 0;
}
