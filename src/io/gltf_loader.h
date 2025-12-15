#pragma once

#include <string>

#include "camera/camera.h"
#include "scene/scene.h"

struct GltfLoadOptions {
    bool use_first_camera = false;
};

bool load_gltf_scene(const std::string& path,
                     Scene& out_scene,
                     CameraSettings* out_camera_settings,
                     std::string* out_error,
                     const GltfLoadOptions& options = {});

