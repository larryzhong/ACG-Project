#include <cstdlib>
#include <iostream>
#include <string>
#include <iomanip>
#include <sstream>
#include <cerrno>
#include <limits>
#include <cmath>

#include "camera/camera.h"
#include "core/color.h"
#include "core/ray.h"
#include "core/rng.h"
#include "core/vec3.h"
#include "io/image_io.h"
#include "io/gltf_loader.h"
#include "io/simple_scene_builder.h"
#include "render/film.h"
#include "render/path_tracer.h"
#include "scene/scene.h"

#define STB_IMAGE_IMPLEMENTATION
#include "io/stb_image.h"

namespace {

struct Options {
    int width = 960;
    int height = 540;
    int spp = 64;
    int max_depth = 8;
    bool spp_set = false;
    bool max_depth_set = false;

    float aperture = 0.0f;
    float focus_dist = 0.0f;
    float shutter_open = 0.0f;
    float shutter_close = 1.0f;

    std::string output = "bedroom.png";
    std::string scene_name = "hotel";

    // NOTE:
    // --gltf: for scene "gltf" -> load model into a simple test room
    // --gltf: for scene "hotel" -> OUTSIDE model seen through the window
    std::string gltf_path;

    // New: explicit indoor plant model (optional).
    std::string plant_gltf_path = "../assets/models/potted_plant/potted_plant_02_4k.gltf";

    std::string env_path;
    bool hide_env_bg = false;
    bool hide_env_bg_set = false;

    bool turntable = false;
    int turntable_frames = 60;
    float turntable_radius = 0.0f;
    float turntable_height = 0.0f;
    Vec3 turntable_center = Vec3(0.0f);
    bool turntable_radius_set = false;
    bool turntable_height_set = false;
    bool turntable_center_set = false;

    Vec3 look_from_override = Vec3(0.0f);
    Vec3 look_at_override = Vec3(0.0f);
    Vec3 up_override = Vec3(0.0f, 1.0f, 0.0f);
    bool look_from_set = false;
    bool look_at_set = false;
    bool up_set = false;
};

void print_usage(const char* exe) {
    std::cerr
        << "Usage: " << (exe ? exe : "pathtracer") << " [options]\n"
        << "Options:\n"
        << "  --scene <name>           Built-in scene (simple|dof|motion|texture|random|solar|alpha|mesh|gltf|hotel)\n"
        << "  --gltf <path>            glTF file (scene 'gltf') OR outside model (scene 'hotel')\n"
        << "  --plant-gltf <path>      Indoor plant glTF (scene 'hotel' only)\n"
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
    if (!text) return false;
    char* end = nullptr;
    errno = 0;
    const long v = std::strtol(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0') return false;
    if (v < std::numeric_limits<int>::min() || v > std::numeric_limits<int>::max()) return false;
    out = static_cast<int>(v);
    return true;
}

bool parse_float_arg(const char* text, float& out) {
    if (!text) return false;
    char* end = nullptr;
    errno = 0;
    const float v = std::strtof(text, &end);
    if (errno != 0 || end == text || *end != '\0') return false;
    if (!std::isfinite(v)) return false;
    out = v;
    return true;
}

bool parse_args(int argc, char** argv, Options& opt) {
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto missing_value = [&](int needed) {
            std::cerr << "Missing value for " << arg << " (need " << needed << " argument(s)).\n";
            print_usage(argv[0]);
        };

        if (arg == "--help" || arg == "-h") {
            print_usage(argv[0]);
            return false; // caller should exit 0
        } else if (arg == "--output") {
            if (i + 1 >= argc) { missing_value(1); return false; }
            opt.output = argv[++i];
        } else if (arg == "--width") {
            if (i + 1 >= argc) { missing_value(1); return false; }
            if (!parse_int_arg(argv[++i], opt.width)) { std::cerr << "Invalid --width value.\n"; return false; }
        } else if (arg == "--height") {
            if (i + 1 >= argc) { missing_value(1); return false; }
            if (!parse_int_arg(argv[++i], opt.height)) { std::cerr << "Invalid --height value.\n"; return false; }
        } else if (arg == "--scene") {
            if (i + 1 >= argc) { missing_value(1); return false; }
            opt.scene_name = argv[++i];
        } else if (arg == "--gltf") {
            if (i + 1 >= argc) { missing_value(1); return false; }
            opt.gltf_path = argv[++i];
        } else if (arg == "--plant-gltf") {
            if (i + 1 >= argc) { missing_value(1); return false; }
            opt.plant_gltf_path = argv[++i];
        } else if (arg == "--env") {
            if (i + 1 >= argc) { missing_value(1); return false; }
            opt.env_path = argv[++i];
        } else if (arg == "--hide-env-bg") {
            opt.hide_env_bg = true;
            opt.hide_env_bg_set = true;
        } else if (arg == "--spp") {
            if (i + 1 >= argc) { missing_value(1); return false; }
            if (!parse_int_arg(argv[++i], opt.spp)) { std::cerr << "Invalid --spp value.\n"; return false; }
            opt.spp_set = true;
        } else if (arg == "--max-depth") {
            if (i + 1 >= argc) { missing_value(1); return false; }
            if (!parse_int_arg(argv[++i], opt.max_depth)) { std::cerr << "Invalid --max-depth value.\n"; return false; }
            opt.max_depth_set = true;
        } else if (arg == "--aperture") {
            if (i + 1 >= argc) { missing_value(1); return false; }
            if (!parse_float_arg(argv[++i], opt.aperture)) { std::cerr << "Invalid --aperture value.\n"; return false; }
        } else if (arg == "--focus-dist") {
            if (i + 1 >= argc) { missing_value(1); return false; }
            if (!parse_float_arg(argv[++i], opt.focus_dist)) { std::cerr << "Invalid --focus-dist value.\n"; return false; }
        } else if (arg == "--shutter-open") {
            if (i + 1 >= argc) { missing_value(1); return false; }
            if (!parse_float_arg(argv[++i], opt.shutter_open)) { std::cerr << "Invalid --shutter-open value.\n"; return false; }
        } else if (arg == "--shutter-close") {
            if (i + 1 >= argc) { missing_value(1); return false; }
            if (!parse_float_arg(argv[++i], opt.shutter_close)) { std::cerr << "Invalid --shutter-close value.\n"; return false; }
        } else if (arg == "--turntable") {
            opt.turntable = true;
        } else if (arg == "--frames") {
            if (i + 1 >= argc) { missing_value(1); return false; }
            if (!parse_int_arg(argv[++i], opt.turntable_frames)) { std::cerr << "Invalid --frames value.\n"; return false; }
        } else if (arg == "--turntable-radius") {
            if (i + 1 >= argc) { missing_value(1); return false; }
            if (!parse_float_arg(argv[++i], opt.turntable_radius)) { std::cerr << "Invalid --turntable-radius value.\n"; return false; }
            opt.turntable_radius_set = true;
        } else if (arg == "--turntable-height") {
            if (i + 1 >= argc) { missing_value(1); return false; }
            if (!parse_float_arg(argv[++i], opt.turntable_height)) { std::cerr << "Invalid --turntable-height value.\n"; return false; }
            opt.turntable_height_set = true;
        } else if (arg == "--turntable-center") {
            if (i + 3 >= argc) { missing_value(3); return false; }
            float cx, cy, cz;
            if (!parse_float_arg(argv[++i], cx) ||
                !parse_float_arg(argv[++i], cy) ||
                !parse_float_arg(argv[++i], cz)) {
                std::cerr << "Invalid --turntable-center value(s).\n";
                return false;
            }
            opt.turntable_center = Vec3(cx, cy, cz);
            opt.turntable_center_set = true;
        } else if (arg == "--look-from") {
            if (i + 3 >= argc) { missing_value(3); return false; }
            if (!parse_float_arg(argv[++i], opt.look_from_override.x) ||
                !parse_float_arg(argv[++i], opt.look_from_override.y) ||
                !parse_float_arg(argv[++i], opt.look_from_override.z)) {
                std::cerr << "Invalid --look-from value(s).\n";
                return false;
            }
            opt.look_from_set = true;
        } else if (arg == "--look-at") {
            if (i + 3 >= argc) { missing_value(3); return false; }
            if (!parse_float_arg(argv[++i], opt.look_at_override.x) ||
                !parse_float_arg(argv[++i], opt.look_at_override.y) ||
                !parse_float_arg(argv[++i], opt.look_at_override.z)) {
                std::cerr << "Invalid --look-at value(s).\n";
                return false;
            }
            opt.look_at_set = true;
        } else if (arg == "--up") {
            if (i + 3 >= argc) { missing_value(3); return false; }
            if (!parse_float_arg(argv[++i], opt.up_override.x) ||
                !parse_float_arg(argv[++i], opt.up_override.y) ||
                !parse_float_arg(argv[++i], opt.up_override.z)) {
                std::cerr << "Invalid --up value(s).\n";
                return false;
            }
            opt.up_set = true;
        } else {
            std::cerr << "Unknown argument: " << arg << "\n";
            print_usage(argv[0]);
            return false;
        }
    }
    return true;
}

bool compute_gltf_bounds(const std::vector<GltfMeshInstance>& meshes, AABB& out_bounds) {
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

bool add_scaled_gltf(
    Scene& scene,
    const std::string& path,
    const Vec3& target_pos,
    float target_height,
    std::string* out_err)
{
    std::vector<GltfMeshInstance> meshes;
    std::string gltf_err;
    GltfLoadOptions gltf_options;
    if (!load_gltf_meshes(path, meshes, &gltf_err, gltf_options)) {
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

} // namespace

int main(int argc, char** argv) {
    Options opt;
    if (!parse_args(argc, argv, opt)) {
        // parse_args already printed usage; treat --help as clean exit.
        for (int i = 1; i < argc; ++i) {
            if (std::string(argv[i]) == "--help" || std::string(argv[i]) == "-h") return 0;
        }
        return 1;
    }

    if (opt.width <= 0 || opt.height <= 0 || opt.spp <= 0 || opt.max_depth <= 0) {
        std::cerr << "Invalid render parameters.\n";
        return 1;
    }
    if (opt.aperture < 0.0f) {
        std::cerr << "Invalid --aperture value (must be >= 0).\n";
        return 1;
    }
    if (opt.shutter_close < opt.shutter_open) {
        std::cerr << "Invalid shutter interval (--shutter-close must be >= --shutter-open).\n";
        return 1;
    }
    if (opt.turntable && opt.turntable_frames <= 0) {
        std::cerr << "Invalid --frames value (must be > 0).\n";
        return 1;
    }
    if (opt.turntable_radius_set && opt.turntable_radius <= 0.0f) {
        std::cerr << "Invalid --turntable-radius value (must be > 0).\n";
        return 1;
    }

    Film film(opt.width, opt.height);
    Scene scene;

    CameraSettings cam_settings;
    cam_settings.aspect_ratio = static_cast<float>(opt.width) / static_cast<float>(opt.height);
    cam_settings.vertical_fov_deg = 55.0f;
    cam_settings.image_width = opt.width;
    cam_settings.image_height = opt.height;

    // -------------------------
    // Scene selection
    // -------------------------
    if (opt.scene_name == "simple") {
        scene = build_simple_scene_basic();
        cam_settings.look_from = Vec3(0.0f, 1.0f, 2.5f);
        cam_settings.look_at = Vec3(0.0f, 1.0f, -1.0f);
        cam_settings.vertical_fov_deg = 45.0f;
    } else if (opt.scene_name == "dof") {
        scene = build_dof_scene();
        cam_settings.look_from = Vec3(0.0f, 2.0f, 2.5f);
        cam_settings.look_at = Vec3(0.0f, 0.5f, -2.0f);
        cam_settings.vertical_fov_deg = 35.0f;
    } else if (opt.scene_name == "motion") {
        scene = build_motion_blur_scene();
        cam_settings.look_from = Vec3(0.0f, 2.0f, 8.0f);
        cam_settings.look_at = Vec3(0.0f, 1.5f, 0.0f);
        cam_settings.vertical_fov_deg = 30.0f;
    } else if (opt.scene_name == "texture") {
        scene = build_texture_scene();
        cam_settings.look_from = Vec3(0.0f, 3.0f, 8.0f);
        cam_settings.look_at = Vec3(0.0f, 2.0f, 0.0f);
        cam_settings.vertical_fov_deg = 35.0f;
    } else if (opt.scene_name == "random") {
        scene = build_random_scene();
        cam_settings.look_from = Vec3(13.0f, 2.0f, 3.0f);
        cam_settings.look_at = Vec3(0.0f, 0.0f, 0.0f);
        cam_settings.vertical_fov_deg = 20.0f;
        if (opt.aperture == 0.0f) opt.aperture = 0.1f;
        if (opt.focus_dist < 0.0f) opt.focus_dist = 10.0f;
    } else if (opt.scene_name == "solar") {
        scene = build_solar_system_scene();
        cam_settings.look_from = Vec3(0.0f, 0.0f, 8.5f);
        cam_settings.look_at = Vec3(0.0f, 0.0f, 0.0f);
        cam_settings.vertical_fov_deg = 35.0f;
        cam_settings.aperture = 0.05f;
        cam_settings.focus_dist = 8.5f;
    } else if (opt.scene_name == "alpha") {
        scene = build_alpha_shadow_scene();
        cam_settings.look_from = Vec3(0.0f, 3.0f, 6.0f);
        cam_settings.look_at = Vec3(0.0f, 1.0f, 0.0f);
        cam_settings.vertical_fov_deg = 40.0f;
    } else if (opt.scene_name == "mesh") {
        scene = build_mesh_scene();
        cam_settings.look_from = Vec3(0.0f, 1.6f, 3.5f);
        cam_settings.look_at = Vec3(0.0f, 0.7f, -1.0f);
        cam_settings.vertical_fov_deg = 35.0f;
    } else if (opt.scene_name == "gltf") {
        // Simple test room + ceiling light + glTF.
        auto ground_mat = std::make_shared<Lambertian>(
            std::make_shared<SolidColor>(Color(0.75f, 0.75f, 0.75f)));
        scene.objects.push_back(std::make_shared<Quad>(
            Vec3(-2.0f, 0.0f, -2.0f), Vec3(4.0f, 0.0f, 0.0f), Vec3(0.0f, 0.0f, 4.0f), ground_mat));

        auto back_mat = std::make_shared<Lambertian>(
            std::make_shared<SolidColor>(Color(0.85f, 0.85f, 0.85f)));
        scene.objects.push_back(std::make_shared<Quad>(
            Vec3(-2.0f, 0.0f, -2.0f), Vec3(4.0f, 0.0f, 0.0f), Vec3(0.0f, 3.0f, 0.0f), back_mat));

        auto light_mat = std::make_shared<DiffuseLight>(
            std::make_shared<SolidColor>(Color(25.0f, 25.0f, 25.0f)));
        auto light = std::make_shared<Quad>(
            Vec3(-0.7f, 2.0f, -0.7f), Vec3(1.4f, 0.0f, 0.0f), Vec3(0.0f, 0.0f, 1.4f), light_mat);
        scene.objects.push_back(light);
        scene.lights.add_area_light(light);

        if (opt.gltf_path.empty()) {
            std::cerr << "Scene 'gltf' requires --gltf <path>.\n";
            return 1;
        }
        std::string gltf_err;
        GltfLoadOptions gltf_options;
        if (!append_gltf_to_scene(opt.gltf_path, scene, Transform::identity(), &gltf_err, gltf_options)) {
            std::cerr << "Failed to load glTF: " << opt.gltf_path << "\n";
            if (!gltf_err.empty()) std::cerr << gltf_err << "\n";
            return 1;
        }

        cam_settings.look_from = Vec3(0.0f, 0.55f, 1.6f);
        cam_settings.look_at = Vec3(0.0f, 0.45f, 0.0f);
        cam_settings.vertical_fov_deg = 35.0f;
    } else if (opt.scene_name == "hotel") {
        // Modern bedroom aligned to your reference image.
        scene = build_hotel_room_scene();

        // If user didn't set quality, bump defaults for this interior shot.
        if (!opt.spp_set) opt.spp = 256;
        if (!opt.max_depth_set) opt.max_depth = 10;

        // Lighting: keep HDR for bounce/ambient if present.
        // Background: if we have outside glTF, hiding env bg usually looks better.
        if (!opt.hide_env_bg_set && !opt.gltf_path.empty()) {
            opt.hide_env_bg = true;
        }
        if (opt.env_path.empty()) {
            // Keep your previous default HDR if you have it; otherwise leave empty.
            opt.env_path = "../assets/hdri/venice_sunset_4k.hdr";
        }

        // Place indoor plant (optional).
        if (!opt.plant_gltf_path.empty()) {
            std::string err;
            // Near window-left in the room.
            if (!add_scaled_gltf(scene, opt.plant_gltf_path, Vec3(-2.85f, 0.0f, -7.95f), 1.20f, &err)) {
                std::cerr << "Failed to load plant glTF: " << opt.plant_gltf_path << "\n";
                if (!err.empty()) std::cerr << err << "\n";
                return 1;
            }
        }

        // Place outside model (this is the key requirement).
        if (!opt.gltf_path.empty()) {
            std::string err;
            // Put it outside the window wall (window wall is around z ~= -9 in the builder below).
            if (!add_scaled_gltf(scene, opt.gltf_path, Vec3(0.0f, 0.0f, -13.5f), 8.0f, &err)) {
                std::cerr << "Failed to load outside glTF: " << opt.gltf_path << "\n";
                if (!err.empty()) std::cerr << err << "\n";
                return 1;
            }
        }

        // Camera: straight view facing the window (not diagonal).
        cam_settings.look_from = Vec3(0.0f, 1.30f, 0.80f);
        cam_settings.look_at   = Vec3(0.0f, 1.30f, -5.00f);
        cam_settings.vertical_fov_deg = 55.0f;
    } else {
        std::cerr << "Unknown scene name: " << opt.scene_name << "\n";
        return 1;
    }

    // -------------------------
    // Camera overrides
    // -------------------------
    if (opt.look_from_set) cam_settings.look_from = opt.look_from_override;
    if (opt.look_at_set) cam_settings.look_at = opt.look_at_override;
    if (opt.up_set) cam_settings.up = opt.up_override;

    if ((cam_settings.look_from - cam_settings.look_at).length_squared() <= 1e-8f) {
        std::cerr << "Invalid camera: look_from and look_at are too close.\n";
        return 1;
    }
    if (cam_settings.up.length_squared() <= 1e-8f) {
        std::cerr << "Invalid camera: up vector has near-zero length.\n";
        return 1;
    }

    // -------------------------
    // Environment
    // -------------------------
    if (!opt.env_path.empty()) {
        scene.environment = std::make_shared<EnvironmentMap>();
        std::string env_error;
        if (!scene.environment->load_hdr(opt.env_path, &env_error)) {
            std::cerr << "Failed to load environment: " << opt.env_path << "\n";
            if (!env_error.empty()) std::cerr << env_error << "\n";
            return 1;
        }
    } else if (opt.hide_env_bg) {
        std::cerr << "Warning: --hide-env-bg set without --env; rays escaping outside will render black.\n";
    }
    scene.hide_environment_background = opt.hide_env_bg;

    scene.build_bvh();

    // -------------------------
    // DOF / motion
    // -------------------------
    cam_settings.aperture = opt.aperture;
    cam_settings.t0 = opt.shutter_open;
    cam_settings.t1 = opt.shutter_close;

    if (opt.focus_dist <= 0.0f) {
        cam_settings.focus_dist = (cam_settings.look_from - cam_settings.look_at).length();
    } else {
        cam_settings.focus_dist = opt.focus_dist;
    }

    // -------------------------
    // Render
    // -------------------------
    if (!opt.turntable) {
        Camera camera(cam_settings);
        PathTracer integrator(opt.max_depth);

        std::cout << "Rendering scene: " << opt.scene_name << "\n";
        render_image(scene, camera, integrator, film, opt.spp);

        if (opt.output.size() >= 4 && opt.output.substr(opt.output.size() - 4) == ".png") {
            write_png(opt.output, film);
            std::cout << "Wrote PNG image to " << opt.output << "\n";
        } else {
            write_ppm(opt.output, film);
            std::cout << "Wrote PPM image to " << opt.output << "\n";
        }
    } else {
        PathTracer integrator(opt.max_depth);

        if (opt.turntable_radius <= 0.0f) {
            Vec3 diff = cam_settings.look_from - cam_settings.look_at;
            opt.turntable_radius = std::sqrt(diff.x * diff.x + diff.z * diff.z);
        }
        if (!opt.turntable_height_set) opt.turntable_height = cam_settings.look_from.y;
        if (!opt.turntable_center_set) opt.turntable_center = cam_settings.look_at;

        std::cout << "Turntable mode: " << opt.turntable_frames << " frames\n";
        std::cout << "  Center: (" << opt.turntable_center.x << ", "
                  << opt.turntable_center.y << ", " << opt.turntable_center.z << ")\n";
        std::cout << "  Radius: " << opt.turntable_radius << ", Height: " << opt.turntable_height << "\n";

        std::string base_name = opt.output;
        std::string extension = ".png";
        if (opt.output.size() >= 4) {
            std::string ext = opt.output.substr(opt.output.size() - 4);
            if (ext == ".png" || ext == ".ppm") {
                extension = ext;
                base_name = opt.output.substr(0, opt.output.size() - 4);
            }
        }

        const float angle_step = 2.0f * kPi / static_cast<float>(opt.turntable_frames);

        for (int frame = 0; frame < opt.turntable_frames; ++frame) {
            const float angle = frame * angle_step;

            const float x = opt.turntable_center.x + opt.turntable_radius * std::sin(angle);
            const float z = opt.turntable_center.z + opt.turntable_radius * std::cos(angle);

            cam_settings.look_from = Vec3(x, opt.turntable_height, z);
            cam_settings.look_at = opt.turntable_center;
            cam_settings.focus_dist = (cam_settings.look_from - cam_settings.look_at).length();

            Camera camera(cam_settings);
            Film frame_film(opt.width, opt.height);

            std::cout << "\rRendering frame " << (frame + 1) << "/" << opt.turntable_frames
                      << " (angle: " << static_cast<int>(angle * 180.0f / kPi) << " deg)" << std::flush;

            render_image(scene, camera, integrator, frame_film, opt.spp);

            std::ostringstream filename;
            filename << base_name << "_" << std::setfill('0') << std::setw(4) << frame << extension;

            if (extension == ".png") write_png(filename.str(), frame_film);
            else write_ppm(filename.str(), frame_film);
        }

        std::cout << "\nTurntable rendering complete!\n";
        std::cout << "To create GIF, run:\n";
        std::cout << "  ffmpeg -framerate 30 -i " << base_name
                  << "_%04d.png -vf \"fps=30,scale=480:-1:flags=lanczos\" "
                  << base_name << ".gif\n";
    }

    return 0;
}
