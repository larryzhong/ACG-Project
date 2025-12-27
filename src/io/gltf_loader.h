#pragma once

#include <string>
#include <vector>

#include "math/transform.h"
#include "scene/material.h"
#include "scene/mesh.h"
#include "scene/scene.h"

struct GltfMeshInstance {
    MeshDataPtr data;
    MaterialPtr material;
    Transform transform;
    std::string name;
};

struct GltfLoadOptions {
    // glTF 2.0 uses a top-left texture origin; keep this false unless you know your assets need it.
    bool flip_v = false;
};

bool load_gltf_meshes(const std::string& path,
                      std::vector<GltfMeshInstance>& out_meshes,
                      std::string* out_error,
                      const GltfLoadOptions& options = {});

bool append_gltf_to_scene(const std::string& path,
                          Scene& inout_scene,
                          const Transform& transform,
                          std::string* out_error,
                          const GltfLoadOptions& options = {});

