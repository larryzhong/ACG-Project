#include <cstdlib>
#include <iostream>
#include <string>
#include <iomanip>
#include <sstream>
#include <cerrno>
#include <limits>

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

namespace {

void print_usage(const char* exe) {
    std::cerr
        << "Usage: " << (exe ? exe : "pathtracer") << " [options]\n"
        << "Options:\n"
        << "  --scene <name>           Built-in scene (simple|dof|motion|texture|random|solar|alpha|mesh)\n"
        << "  --output <path>          Output image path (.png writes PNG, otherwise PPM)\n"
        << "  --width <int>            Image width\n"
        << "  --height <int>           Image height\n"
        << "  --spp <int>              Samples per pixel\n"
        << "  --max-depth <int>        Path tracer max depth\n"
        << "  --env <path>             HDR environment map\n"
        << "  --hide-env-bg            Hide environment background (still lights the scene)\n"
        << "  --aperture <float>       Aperture diameter (0 disables DOF)\n"
        << "  --focus-dist <float>     Focus distance (<=0 uses |look_from-look_at|)\n"
        << "  --shutter-open <float>   Shutter open time\n"
        << "  --shutter-close <float>  Shutter close time\n"
        << "  --look-from x y z        Override camera position\n"
        << "  --look-at x y z          Override camera target\n"
        << "  --up x y z               Override camera up vector\n"
        << "  --turntable              Render turntable animation\n"
        << "  --frames <int>           Turntable frames\n"
        << "  --turntable-radius <f>   Turntable radius (XZ plane)\n"
        << "  --turntable-height <f>   Turntable height (Y)\n"
        << "  --turntable-center x y z Turntable center\n"
        << "  --help                   Show this help\n";
}

bool parse_int_arg(const char* text, int& out) {
    if (!text) {
        return false;
    }
    char* end = nullptr;
    errno = 0;
    const long v = std::strtol(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0') {
        return false;
    }
    if (v < std::numeric_limits<int>::min() || v > std::numeric_limits<int>::max()) {
        return false;
    }
    out = static_cast<int>(v);
    return true;
}

bool parse_float_arg(const char* text, float& out) {
    if (!text) {
        return false;
    }
    char* end = nullptr;
    errno = 0;
    const float v = std::strtof(text, &end);
    if (errno != 0 || end == text || *end != '\0') {
        return false;
    }
    if (!std::isfinite(v)) {
        return false;
    }
    out = v;
    return true;
}

}  // namespace

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
    std::string env_path;
    bool hide_env_bg = false;
    bool turntable_mode = false;
    int turntable_frames = 60;
    float turntable_radius = 0.0f;
    float turntable_height = 0.0f;
    Vec3 turntable_center(0.0f, 0.0f, 0.0f);
    bool turntable_radius_set = false;
    bool turntable_height_set = false;
    bool turntable_center_set = false;
    Vec3 look_from_override(0.0f);
    Vec3 look_at_override(0.0f);
    Vec3 up_override(0.0f, 1.0f, 0.0f);
    bool look_from_set = false;
    bool look_at_set = false;
    bool up_set = false;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto missing_value = [&](int needed) {
            std::cerr << "Missing value for " << arg << " (need " << needed << " argument(s)).\n";
            print_usage(argv[0]);
        };

        if (arg == "--help" || arg == "-h") {
            print_usage(argv[0]);
            return 0;
        } else if (arg == "--output") {
            if (i + 1 >= argc) {
                missing_value(1);
                return 1;
            }
            output = argv[++i];
        } else if (arg == "--width") {
            if (i + 1 >= argc) {
                missing_value(1);
                return 1;
            }
            if (!parse_int_arg(argv[++i], width)) {
                std::cerr << "Invalid --width value.\n";
                return 1;
            }
        } else if (arg == "--height") {
            if (i + 1 >= argc) {
                missing_value(1);
                return 1;
            }
            if (!parse_int_arg(argv[++i], height)) {
                std::cerr << "Invalid --height value.\n";
                return 1;
            }
        } else if (arg == "--scene") {
            if (i + 1 >= argc) {
                missing_value(1);
                return 1;
            }
            scene_name = argv[++i];
        } else if (arg == "--env") {
            if (i + 1 >= argc) {
                missing_value(1);
                return 1;
            }
            env_path = argv[++i];
        } else if (arg == "--hide-env-bg") {
            hide_env_bg = true;
        } else if (arg == "--spp") {
            if (i + 1 >= argc) {
                missing_value(1);
                return 1;
            }
            if (!parse_int_arg(argv[++i], samples_per_pixel)) {
                std::cerr << "Invalid --spp value.\n";
                return 1;
            }
        } else if (arg == "--max-depth") {
            if (i + 1 >= argc) {
                missing_value(1);
                return 1;
            }
            if (!parse_int_arg(argv[++i], max_depth)) {
                std::cerr << "Invalid --max-depth value.\n";
                return 1;
            }
        } else if (arg == "--aperture") {
            if (i + 1 >= argc) {
                missing_value(1);
                return 1;
            }
            if (!parse_float_arg(argv[++i], aperture)) {
                std::cerr << "Invalid --aperture value.\n";
                return 1;
            }
        } else if (arg == "--focus-dist") {
            if (i + 1 >= argc) {
                missing_value(1);
                return 1;
            }
            if (!parse_float_arg(argv[++i], focus_dist)) {
                std::cerr << "Invalid --focus-dist value.\n";
                return 1;
            }
        } else if (arg == "--shutter-open") {
            if (i + 1 >= argc) {
                missing_value(1);
                return 1;
            }
            if (!parse_float_arg(argv[++i], shutter_open)) {
                std::cerr << "Invalid --shutter-open value.\n";
                return 1;
            }
        } else if (arg == "--shutter-close") {
            if (i + 1 >= argc) {
                missing_value(1);
                return 1;
            }
            if (!parse_float_arg(argv[++i], shutter_close)) {
                std::cerr << "Invalid --shutter-close value.\n";
                return 1;
            }
        } else if (arg == "--turntable") {
            turntable_mode = true;
        } else if (arg == "--frames") {
            if (i + 1 >= argc) {
                missing_value(1);
                return 1;
            }
            if (!parse_int_arg(argv[++i], turntable_frames)) {
                std::cerr << "Invalid --frames value.\n";
                return 1;
            }
        } else if (arg == "--turntable-radius") {
            if (i + 1 >= argc) {
                missing_value(1);
                return 1;
            }
            if (!parse_float_arg(argv[++i], turntable_radius)) {
                std::cerr << "Invalid --turntable-radius value.\n";
                return 1;
            }
            turntable_radius_set = true;
        } else if (arg == "--turntable-height") {
            if (i + 1 >= argc) {
                missing_value(1);
                return 1;
            }
            if (!parse_float_arg(argv[++i], turntable_height)) {
                std::cerr << "Invalid --turntable-height value.\n";
                return 1;
            }
            turntable_height_set = true;
        } else if (arg == "--turntable-center") {
            if (i + 3 >= argc) {
                missing_value(3);
                return 1;
            }
            float cx, cy, cz;
            if (!parse_float_arg(argv[++i], cx) ||
                !parse_float_arg(argv[++i], cy) ||
                !parse_float_arg(argv[++i], cz)) {
                std::cerr << "Invalid --turntable-center value(s).\n";
                return 1;
            }
            turntable_center = Vec3(cx, cy, cz);
            turntable_center_set = true;
        } else if (arg == "--look-from") {
            if (i + 3 >= argc) {
                missing_value(3);
                return 1;
            }
            if (!parse_float_arg(argv[++i], look_from_override.x) ||
                !parse_float_arg(argv[++i], look_from_override.y) ||
                !parse_float_arg(argv[++i], look_from_override.z)) {
                std::cerr << "Invalid --look-from value(s).\n";
                return 1;
            }
            look_from_set = true;
        } else if (arg == "--look-at") {
            if (i + 3 >= argc) {
                missing_value(3);
                return 1;
            }
            if (!parse_float_arg(argv[++i], look_at_override.x) ||
                !parse_float_arg(argv[++i], look_at_override.y) ||
                !parse_float_arg(argv[++i], look_at_override.z)) {
                std::cerr << "Invalid --look-at value(s).\n";
                return 1;
            }
            look_at_set = true;
        } else if (arg == "--up") {
            if (i + 3 >= argc) {
                missing_value(3);
                return 1;
            }
            if (!parse_float_arg(argv[++i], up_override.x) ||
                !parse_float_arg(argv[++i], up_override.y) ||
                !parse_float_arg(argv[++i], up_override.z)) {
                std::cerr << "Invalid --up value(s).\n";
                return 1;
            }
            up_set = true;
        } else {
            std::cerr << "Unknown argument: " << arg << "\n";
            print_usage(argv[0]);
            return 1;
        }
    }

    if (width <= 0 || height <= 0 || samples_per_pixel <= 0 || max_depth <= 0) {
        std::cerr << "Invalid render parameters.\n";
        return 1;
    }
    if (aperture < 0.0f) {
        std::cerr << "Invalid --aperture value (must be >= 0).\n";
        return 1;
    }
    if (shutter_close < shutter_open) {
        std::cerr << "Invalid shutter interval (--shutter-close must be >= --shutter-open).\n";
        return 1;
    }
    if (turntable_mode && turntable_frames <= 0) {
        std::cerr << "Invalid --frames value (must be > 0).\n";
        return 1;
    }
    if (turntable_radius_set && turntable_radius <= 0.0f) {
        std::cerr << "Invalid --turntable-radius value (must be > 0).\n";
        return 1;
    }

    Film film(width, height);
    Scene scene;
    
    CameraSettings cam_settings;
    cam_settings.aspect_ratio = static_cast<float>(width) / static_cast<float>(height);
    cam_settings.vertical_fov_deg = 45.0f;
    cam_settings.image_width = width;
    cam_settings.image_height = height;
    
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

    if (look_from_set) {
        cam_settings.look_from = look_from_override;
    }
    if (look_at_set) {
        cam_settings.look_at = look_at_override;
    }
    if (up_set) {
        cam_settings.up = up_override;
    }

    if ((cam_settings.look_from - cam_settings.look_at).length_squared() <= 1e-8f) {
        std::cerr << "Invalid camera: look_from and look_at are too close.\n";
        return 1;
    }
    if (cam_settings.up.length_squared() <= 1e-8f) {
        std::cerr << "Invalid camera: up vector has near-zero length.\n";
        return 1;
    }

    if (!env_path.empty()) {
        scene.environment = std::make_shared<EnvironmentMap>();
        std::string env_error;
        if (!scene.environment->load_hdr(env_path, &env_error)) {
            std::cerr << "Failed to load environment: " << env_path << "\n";
            if (!env_error.empty()) {
                std::cerr << env_error << "\n";
            }
            return 1;
        }
    } else if (hide_env_bg) {
        std::cerr << "Warning: --hide-env-bg set without --env; if the scene has no other lights it will render black.\n";
    }
    scene.hide_environment_background = hide_env_bg;
    if (scene.hide_environment_background && scene.environment && scene.environment->valid()) {
        std::cout << "Environment background hidden (environment still used for lighting).\n";
    }

    scene.build_bvh();

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
        if (!turntable_height_set) {
            turntable_height = cam_settings.look_from.y;
        }
        if (!turntable_center_set) {
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
