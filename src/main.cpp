#include <cstdlib>
#include <iostream>
#include <string>

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
        }
    }

    if (width <= 0 || height <= 0 || samples_per_pixel <= 0 || max_depth <= 0) {
        std::cerr << "Invalid render parameters.\n";
        return 1;
    }

    Film film(width, height);

    Scene scene;
    if (scene_name == "simple") {
        scene = build_simple_scene_basic();
    }

    scene.build_bvh();

    CameraSettings cam_settings;
    cam_settings.aspect_ratio = static_cast<float>(width) / static_cast<float>(height);
    cam_settings.look_from = Vec3(0.0f, 1.0f, -1.05f);
    cam_settings.look_at = Vec3(0.0f, 1.0f, -2.0f);
    cam_settings.up = Vec3(0.0f, 1.0f, 0.0f);
    cam_settings.vertical_fov_deg = 90.0f;
    cam_settings.aperture = aperture;
    if (focus_dist <= 0.0f) {
        const Vec3 diff = cam_settings.look_from - cam_settings.look_at;
        focus_dist = diff.length();
    }
    cam_settings.focus_dist = focus_dist;
    cam_settings.t0 = shutter_open;
    cam_settings.t1 = shutter_close;

    Camera camera(cam_settings);
    PathTracer integrator(max_depth);

    render_image(scene, camera, integrator, film, samples_per_pixel);

    if (output.size() >= 4 && output.substr(output.size() - 4) == ".png") {
        write_png(output, film);
        std::cout << "Wrote PNG image to " << output << "\n";
    } else {
        write_ppm(output, film);
        std::cout << "Wrote PPM image to " << output << "\n";
    }

    return 0;
}
