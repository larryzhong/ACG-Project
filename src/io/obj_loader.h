#pragma once

#include <string>
#include <vector>

#include "math/transform.h"
#include "scene/material.h"
#include "scene/mesh.h"
#include "scene/scene.h"

struct ObjMesh {
    MeshDataPtr data;
    MaterialPtr material;
    bool emissive = false;
    std::string name;
};

struct ObjLoadOptions {
    bool split_by_object = false;
    bool split_by_group = false;
    bool default_smooth_shading = true;
    bool use_illum_for_transmission = true;
};

bool load_obj_meshes(const std::string& path,
                     std::vector<ObjMesh>& out_meshes,
                     std::string* out_error,
                     const ObjLoadOptions& options = {});

bool append_obj_to_scene(const std::string& path,
                         Scene& inout_scene,
                         const Transform& transform,
                         std::string* out_error,
                         const ObjLoadOptions& options = {});

