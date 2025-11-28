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

    CameraSettings cam_settings;
    cam_settings.aspect_ratio = static_cast<float>(width) / static_cast<float>(height);
    cam_settings.look_from = Vec3(0.0f, 1.0f, 2.0f);
    cam_settings.look_at = Vec3(0.0f, 1.0f, -2.0f);
    cam_settings.up = Vec3(0.0f, 1.0f, 0.0f);
    cam_settings.vertical_fov_deg = 40.0f;
    cam_settings.aperture = 0.0f;
    cam_settings.focus_dist = 2.0f;

    Camera camera(cam_settings);
    PathTracer integrator(max_depth);

    render_image(scene, camera, integrator, film, samples_per_pixel);

    write_ppm(output, film);

    std::cout << "Wrote image to " << output << "\n";

    return 0;
}
