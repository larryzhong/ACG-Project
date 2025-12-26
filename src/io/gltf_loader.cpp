#include "io/gltf_loader.h"

#include <cmath>
#include <cstring>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <limits>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "scene/material.h"
#include "scene/mesh.h"
#include "scene/texture.h"

#define TINYGLTF_IMPLEMENTATION
#include "io/tinygltf.h"

namespace {

constexpr float kPi = 3.14159265358979323846f;

struct TextureKey {
    int texture_index = -1;
    ImageTexture::ColorSpace color_space = ImageTexture::ColorSpace::sRGB;
    int channel = -1;

    bool operator==(const TextureKey& other) const {
        return texture_index == other.texture_index &&
               color_space == other.color_space &&
               channel == other.channel;
    }
};

struct TextureKeyHash {
    std::size_t operator()(const TextureKey& k) const noexcept {
        std::size_t h = std::hash<int>{}(k.texture_index);
        h ^= (std::hash<int>{}(static_cast<int>(k.color_space)) + 0x9e3779b9 + (h << 6) + (h >> 2));
        h ^= (std::hash<int>{}(k.channel) + 0x9e3779b9 + (h << 6) + (h >> 2));
        return h;
    }
};

static ImageTexture::WrapMode gltf_wrap_to_mode(int wrap) {
    switch (wrap) {
        case TINYGLTF_TEXTURE_WRAP_CLAMP_TO_EDGE:
            return ImageTexture::WrapMode::ClampToEdge;
        case TINYGLTF_TEXTURE_WRAP_MIRRORED_REPEAT:
            return ImageTexture::WrapMode::MirroredRepeat;
        case TINYGLTF_TEXTURE_WRAP_REPEAT:
        default:
            return ImageTexture::WrapMode::Repeat;
    }
}

struct MaterialBuildResult {
    MaterialPtr material;
    bool emissive = false;
};

static float gltf_clamp01(float v) {
    if (!std::isfinite(v)) {
        return 0.0f;
    }
    return clamp_float(v, 0.0f, 1.0f);
}

static int gltf_parse_texture_index(const tinygltf::Value& tex_info) {
    if (!tex_info.IsObject()) {
        return -1;
    }
    if (!tex_info.Has("index")) {
        return -1;
    }
    const tinygltf::Value& idx = tex_info.Get("index");
    if (!idx.IsNumber()) {
        return -1;
    }
    return idx.GetNumberAsInt();
}

static bool gltf_parse_vec3(const tinygltf::Value& arr, Color& out) {
    if (!arr.IsArray() || arr.ArrayLen() < 3) {
        return false;
    }
    const tinygltf::Value& x = arr.Get(0);
    const tinygltf::Value& y = arr.Get(1);
    const tinygltf::Value& z = arr.Get(2);
    if (!x.IsNumber() || !y.IsNumber() || !z.IsNumber()) {
        return false;
    }
    out = Color(static_cast<float>(x.GetNumberAsDouble()),
                static_cast<float>(y.GetNumberAsDouble()),
                static_cast<float>(z.GetNumberAsDouble()));
    return true;
}

static NormalMapPtr load_normal_map(const tinygltf::Model& model,
                                     int texture_index,
                                     const std::filesystem::path& base_dir,
                                     std::unordered_map<int, NormalMapPtr>& cache,
                                    std::string* err) {
    if (texture_index < 0 || texture_index >= static_cast<int>(model.textures.size())) {
        return nullptr;
    }

    if (auto it = cache.find(texture_index); it != cache.end()) {
        return it->second;
    }

    const tinygltf::Texture& tex = model.textures[texture_index];
    ImageTexture::WrapMode wrap_s = ImageTexture::WrapMode::Repeat;
    ImageTexture::WrapMode wrap_t = ImageTexture::WrapMode::Repeat;
    if (tex.sampler >= 0 && tex.sampler < static_cast<int>(model.samplers.size())) {
        const tinygltf::Sampler& s = model.samplers[tex.sampler];
        wrap_s = gltf_wrap_to_mode(s.wrapS);
        wrap_t = gltf_wrap_to_mode(s.wrapT);
    }
        
    if (tex.source < 0 || tex.source >= static_cast<int>(model.images.size())) {
        return nullptr;
    }

    const tinygltf::Image& image = model.images[tex.source];

    NormalMapPtr result;
    // glTF defines (0,0) at the upper-left of the image, so do not flip V here.
    constexpr bool kFlipV = false;
    const bool is_data_uri =
        image.uri.size() >= 5 && image.uri.rfind("data:", 0) == 0;
    if (!image.image.empty() && image.width > 0 && image.height > 0) {
        const int channels = (image.component > 0) ? image.component : 4;
        result = std::make_shared<NormalMapTexture>(image.image, image.width, image.height, channels, kFlipV, wrap_s, wrap_t);
    } else if (is_data_uri) {
        if (err) *err += "Data URI image was not decoded.\n";
        result = nullptr;
    } else if (!image.uri.empty()) {
        std::filesystem::path p = base_dir / std::filesystem::path(image.uri);
        result = std::make_shared<NormalMapTexture>(p.string(), kFlipV, wrap_s, wrap_t);
    } else {
        if (err) *err += "Unsupported image source.\n";
        result = nullptr;
    }

    cache.emplace(texture_index, result);
    return result;
}

static Transform quat_to_transform(const std::vector<double>& q) {
    if (q.size() != 4) {
        return Transform::identity();
    }

    const float x = static_cast<float>(q[0]);
    const float y = static_cast<float>(q[1]);
    const float z = static_cast<float>(q[2]);
    const float w = static_cast<float>(q[3]);

    const float len2 = x * x + y * y + z * z + w * w;
    if (len2 < 1e-20f) {
        return Transform::identity();
    }

    const float inv_len = 1.0f / std::sqrt(len2);
    const float nx = x * inv_len;
    const float ny = y * inv_len;
    const float nz = z * inv_len;
    const float nw = w * inv_len;

    const float xx = nx * nx;
    const float yy = ny * ny;
    const float zz = nz * nz;
    const float xy = nx * ny;
    const float xz = nx * nz;
    const float yz = ny * nz;
    const float wx = nw * nx;
    const float wy = nw * ny;
    const float wz = nw * nz;

    // Column-major Mat3
    const Vec3 c0(1.0f - 2.0f * (yy + zz), 2.0f * (xy + wz), 2.0f * (xz - wy));
    const Vec3 c1(2.0f * (xy - wz), 1.0f - 2.0f * (xx + zz), 2.0f * (yz + wx));
    const Vec3 c2(2.0f * (xz + wy), 2.0f * (yz - wx), 1.0f - 2.0f * (xx + yy));

    return Transform(Mat3(c0, c1, c2), Vec3::zero());
}

static Transform node_local_transform(const tinygltf::Node& node) {
    if (node.matrix.size() == 16) {
        const auto& m = node.matrix;
        // glTF matrices are column-major.
        const Vec3 c0(static_cast<float>(m[0]), static_cast<float>(m[1]), static_cast<float>(m[2]));
        const Vec3 c1(static_cast<float>(m[4]), static_cast<float>(m[5]), static_cast<float>(m[6]));
        const Vec3 c2(static_cast<float>(m[8]), static_cast<float>(m[9]), static_cast<float>(m[10]));
        const Vec3 t(static_cast<float>(m[12]), static_cast<float>(m[13]), static_cast<float>(m[14]));
        return Transform(Mat3(c0, c1, c2), t);
    }

    // glTF node TRS composes as: M = T * R * S.
    Transform s_tr = Transform::identity();
    Transform r_tr = Transform::identity();
    Transform t_tr = Transform::identity();

    if (node.translation.size() == 3) {
        const Vec3 t(static_cast<float>(node.translation[0]),
                     static_cast<float>(node.translation[1]),
                     static_cast<float>(node.translation[2]));
        t_tr = Transform::translate(t);
    }

    if (node.rotation.size() == 4) {
        r_tr = quat_to_transform(node.rotation);
    }

    if (node.scale.size() == 3) {
        const Vec3 s(static_cast<float>(node.scale[0]),
                     static_cast<float>(node.scale[1]),
                     static_cast<float>(node.scale[2]));
        s_tr = Transform::scale(s);
    }

    return t_tr * r_tr * s_tr;
}

static bool read_accessor_vec3(const tinygltf::Model& model,
                               int accessor_index,
                               std::vector<Vec3>& out,
                               std::string* err) {
    if (accessor_index < 0 || accessor_index >= static_cast<int>(model.accessors.size())) {
        if (err) *err += "Missing accessor.\n";
        return false;
    }

    const tinygltf::Accessor& accessor = model.accessors[accessor_index];
    if (accessor.sparse.isSparse) {
        if (err) *err += "Sparse accessors are not supported.\n";
        return false;
    }
    if (accessor.type != TINYGLTF_TYPE_VEC3) {
        if (err) *err += "Expected VEC3 accessor.\n";
        return false;
    }
    if (accessor.componentType != TINYGLTF_COMPONENT_TYPE_FLOAT) {
        if (err) *err += "Expected FLOAT accessor.\n";
        return false;
    }

    if (accessor.bufferView < 0 || accessor.bufferView >= static_cast<int>(model.bufferViews.size())) {
        if (err) *err += "Invalid bufferView index.\n";
        return false;
    }

    const tinygltf::BufferView& view = model.bufferViews[accessor.bufferView];
    if (view.buffer < 0 || view.buffer >= static_cast<int>(model.buffers.size())) {
        if (err) *err += "Invalid buffer index.\n";
        return false;
    }

    const tinygltf::Buffer& buffer = model.buffers[view.buffer];
    const int stride_i = accessor.ByteStride(view);
    if (stride_i < 0) {
        if (err) *err += "Invalid accessor stride.\n";
        return false;
    }
    const std::size_t elem_stride = static_cast<std::size_t>(stride_i);
    const std::size_t elem_size = sizeof(float) * 3;

    const std::size_t base_offset = static_cast<std::size_t>(view.byteOffset + accessor.byteOffset);
    const std::size_t required =
        accessor.count == 0 ? 0 : ((accessor.count - 1) * elem_stride + elem_size);
    if (base_offset < static_cast<std::size_t>(view.byteOffset) ||
        base_offset + required > buffer.data.size() ||
        (view.byteLength > 0 && base_offset + required > static_cast<std::size_t>(view.byteOffset + view.byteLength))) {
        if (err) *err += "Accessor out of bounds.\n";
        return false;
    }

    out.resize(accessor.count);

    for (std::size_t i = 0; i < accessor.count; ++i) {
        float v[3] = {0, 0, 0};
        const unsigned char* ptr = buffer.data.data() + base_offset + i * elem_stride;
        std::memcpy(v, ptr, sizeof(float) * 3);
        out[i] = Vec3(v[0], v[1], v[2]);
    }

    return true;
}

static bool read_accessor_vec2(const tinygltf::Model& model,
                               int accessor_index,
                               std::vector<Vec2>& out,
                               std::string* err) {
    if (accessor_index < 0 || accessor_index >= static_cast<int>(model.accessors.size())) {
        if (err) *err += "Missing accessor.\n";
        return false;
    }

    const tinygltf::Accessor& accessor = model.accessors[accessor_index];
    if (accessor.sparse.isSparse) {
        if (err) *err += "Sparse accessors are not supported.\n";
        return false;
    }
    if (accessor.type != TINYGLTF_TYPE_VEC2) {
        if (err) *err += "Expected VEC2 accessor.\n";
        return false;
    }
    if (accessor.componentType != TINYGLTF_COMPONENT_TYPE_FLOAT) {
        if (err) *err += "Expected FLOAT accessor.\n";
        return false;
    }

    if (accessor.bufferView < 0 || accessor.bufferView >= static_cast<int>(model.bufferViews.size())) {
        if (err) *err += "Invalid bufferView index.\n";
        return false;
    }

    const tinygltf::BufferView& view = model.bufferViews[accessor.bufferView];
    if (view.buffer < 0 || view.buffer >= static_cast<int>(model.buffers.size())) {
        if (err) *err += "Invalid buffer index.\n";
        return false;
    }

    const tinygltf::Buffer& buffer = model.buffers[view.buffer];
    const int stride_i = accessor.ByteStride(view);
    if (stride_i < 0) {
        if (err) *err += "Invalid accessor stride.\n";
        return false;
    }
    const std::size_t elem_stride = static_cast<std::size_t>(stride_i);
    const std::size_t elem_size = sizeof(float) * 2;

    const std::size_t base_offset = static_cast<std::size_t>(view.byteOffset + accessor.byteOffset);
    const std::size_t required =
        accessor.count == 0 ? 0 : ((accessor.count - 1) * elem_stride + elem_size);
    if (base_offset < static_cast<std::size_t>(view.byteOffset) ||
        base_offset + required > buffer.data.size() ||
        (view.byteLength > 0 && base_offset + required > static_cast<std::size_t>(view.byteOffset + view.byteLength))) {
        if (err) *err += "Accessor out of bounds.\n";
        return false;
    }

    out.resize(accessor.count);

    for (std::size_t i = 0; i < accessor.count; ++i) {
        float v[2] = {0, 0};
        const unsigned char* ptr = buffer.data.data() + base_offset + i * elem_stride;
        std::memcpy(v, ptr, sizeof(float) * 2);
        out[i] = Vec2(v[0], v[1]);
    }

    return true;
}

static bool read_accessor_vec4_tangents(const tinygltf::Model& model,
                                        int accessor_index,
                                        std::vector<Vec3>& out_tangent,
                                        std::vector<float>& out_sign,
                                        std::string* err) {
    out_tangent.clear();
    out_sign.clear();

    if (accessor_index < 0 || accessor_index >= static_cast<int>(model.accessors.size())) {
        if (err) *err += "Missing accessor.\n";
        return false;
    }

    const tinygltf::Accessor& accessor = model.accessors[accessor_index];
    if (accessor.sparse.isSparse) {
        if (err) *err += "Sparse accessors are not supported.\n";
        return false;
    }
    if (accessor.type != TINYGLTF_TYPE_VEC4) {
        if (err) *err += "Expected VEC4 accessor.\n";
        return false;
    }
    if (accessor.componentType != TINYGLTF_COMPONENT_TYPE_FLOAT) {
        if (err) *err += "Expected FLOAT accessor.\n";
        return false;
    }

    if (accessor.bufferView < 0 || accessor.bufferView >= static_cast<int>(model.bufferViews.size())) {
        if (err) *err += "Invalid bufferView index.\n";
        return false;
    }

    const tinygltf::BufferView& view = model.bufferViews[accessor.bufferView];
    if (view.buffer < 0 || view.buffer >= static_cast<int>(model.buffers.size())) {
        if (err) *err += "Invalid buffer index.\n";
        return false;
    }

    const tinygltf::Buffer& buffer = model.buffers[view.buffer];
    const int stride_i = accessor.ByteStride(view);
    if (stride_i < 0) {
        if (err) *err += "Invalid accessor stride.\n";
        return false;
    }
    const std::size_t elem_stride = static_cast<std::size_t>(stride_i);
    const std::size_t elem_size = sizeof(float) * 4;

    const std::size_t base_offset = static_cast<std::size_t>(view.byteOffset + accessor.byteOffset);
    const std::size_t required =
        accessor.count == 0 ? 0 : ((accessor.count - 1) * elem_stride + elem_size);
    if (base_offset < static_cast<std::size_t>(view.byteOffset) ||
        base_offset + required > buffer.data.size() ||
        (view.byteLength > 0 && base_offset + required > static_cast<std::size_t>(view.byteOffset + view.byteLength))) {
        if (err) *err += "Accessor out of bounds.\n";
        return false;
    }

    out_tangent.resize(accessor.count);
    out_sign.resize(accessor.count, 1.0f);

    for (std::size_t i = 0; i < accessor.count; ++i) {
        float v[4] = {0, 0, 0, 1};
        const unsigned char* ptr = buffer.data.data() + base_offset + i * elem_stride;
        std::memcpy(v, ptr, sizeof(float) * 4);
        out_tangent[i] = Vec3(v[0], v[1], v[2]);
        out_sign[i] = (v[3] < 0.0f) ? -1.0f : 1.0f;
    }

    return true;
}

static bool read_accessor_indices(const tinygltf::Model& model,
                                  int accessor_index,
                                  std::vector<std::uint32_t>& out,
                                  std::string* err) {
    if (accessor_index < 0 || accessor_index >= static_cast<int>(model.accessors.size())) {
        if (err) *err += "Missing indices accessor.\n";
        return false;
    }

    const tinygltf::Accessor& accessor = model.accessors[accessor_index];
    if (accessor.sparse.isSparse) {
        if (err) *err += "Sparse accessors are not supported.\n";
        return false;
    }
    if (accessor.type != TINYGLTF_TYPE_SCALAR) {
        if (err) *err += "Expected scalar indices accessor.\n";
        return false;
    }

    if (accessor.bufferView < 0 || accessor.bufferView >= static_cast<int>(model.bufferViews.size())) {
        if (err) *err += "Invalid bufferView index.\n";
        return false;
    }

    const tinygltf::BufferView& view = model.bufferViews[accessor.bufferView];
    if (view.buffer < 0 || view.buffer >= static_cast<int>(model.buffers.size())) {
        if (err) *err += "Invalid buffer index.\n";
        return false;
    }

    const tinygltf::Buffer& buffer = model.buffers[view.buffer];
    const std::size_t elem_size = tinygltf::GetComponentSizeInBytes(accessor.componentType);
    const int stride_i = accessor.ByteStride(view);
    if (stride_i < 0) {
        if (err) *err += "Invalid accessor stride.\n";
        return false;
    }
    const std::size_t elem_stride = static_cast<std::size_t>(stride_i);

    const std::size_t base_offset = static_cast<std::size_t>(view.byteOffset + accessor.byteOffset);
    const std::size_t required =
        accessor.count == 0 ? 0 : ((accessor.count - 1) * elem_stride + elem_size);
    if (base_offset < static_cast<std::size_t>(view.byteOffset) ||
        base_offset + required > buffer.data.size() ||
        (view.byteLength > 0 && base_offset + required > static_cast<std::size_t>(view.byteOffset + view.byteLength))) {
        if (err) *err += "Indices accessor out of bounds.\n";
        return false;
    }

    out.resize(accessor.count);

    for (std::size_t i = 0; i < accessor.count; ++i) {
        const unsigned char* ptr = buffer.data.data() + base_offset + i * elem_stride;
        std::uint32_t idx = 0;
        switch (accessor.componentType) {
            case TINYGLTF_COMPONENT_TYPE_UNSIGNED_BYTE: {
                idx = static_cast<std::uint32_t>(ptr[0]);
                break;
            }
            case TINYGLTF_COMPONENT_TYPE_UNSIGNED_SHORT: {
                std::uint16_t v = 0;
                std::memcpy(&v, ptr, sizeof(v));
                idx = static_cast<std::uint32_t>(v);
                break;
            }
            case TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT: {
                std::uint32_t v = 0;
                std::memcpy(&v, ptr, sizeof(v));
                idx = v;
                break;
            }
            default:
                if (err) *err += "Unsupported indices component type.\n";
                return false;
        }
        out[i] = idx;
    }

    return true;
}

static void compute_normals(const std::vector<Vec3>& positions,
                            const std::vector<std::uint32_t>& indices,
                            std::vector<Vec3>& normals_out) {
    normals_out.assign(positions.size(), Vec3(0.0f));

    const std::size_t tri_count = indices.size() / 3;
    for (std::size_t t = 0; t < tri_count; ++t) {
        const std::uint32_t i0 = indices[3 * t + 0];
        const std::uint32_t i1 = indices[3 * t + 1];
        const std::uint32_t i2 = indices[3 * t + 2];
        if (i0 >= positions.size() || i1 >= positions.size() || i2 >= positions.size()) {
            continue;
        }

        const Vec3& p0 = positions[i0];
        const Vec3& p1 = positions[i1];
        const Vec3& p2 = positions[i2];
        Vec3 n = cross(p1 - p0, p2 - p0);
        if (n.length_squared() > 1e-20f) {
            n = normalize(n);
        } else {
            n = Vec3(0.0f, 1.0f, 0.0f);
        }

        normals_out[i0] += n;
        normals_out[i1] += n;
        normals_out[i2] += n;
    }

    for (auto& n : normals_out) {
        if (n.length_squared() > 1e-20f) {
            n = normalize(n);
        } else {
            n = Vec3(0.0f, 1.0f, 0.0f);
        }
    }
}

static TexturePtr load_texture(const tinygltf::Model& model,
                               int texture_index,
                               const std::filesystem::path& base_dir,
                               ImageTexture::ColorSpace color_space,
                               int channel,
                               std::unordered_map<TextureKey, TexturePtr, TextureKeyHash>& cache,
                               std::string* err) {
    if (texture_index < 0 || texture_index >= static_cast<int>(model.textures.size())) {
        return nullptr;
    }

    const TextureKey key{texture_index, color_space, channel};
    if (auto it = cache.find(key); it != cache.end()) {
        return it->second;
    }

    const tinygltf::Texture& tex = model.textures[texture_index];
    ImageTexture::WrapMode wrap_s = ImageTexture::WrapMode::Repeat;
    ImageTexture::WrapMode wrap_t = ImageTexture::WrapMode::Repeat;
    if (tex.sampler >= 0 && tex.sampler < static_cast<int>(model.samplers.size())) {
        const tinygltf::Sampler& s = model.samplers[tex.sampler];
        wrap_s = gltf_wrap_to_mode(s.wrapS);
        wrap_t = gltf_wrap_to_mode(s.wrapT);
    }
    if (tex.source < 0 || tex.source >= static_cast<int>(model.images.size())) {
        return nullptr;
    }

    const tinygltf::Image& image = model.images[tex.source];

    TexturePtr result;
    // glTF defines (0,0) at the upper-left of the image, so do not flip V here.
    constexpr bool kFlipV = false;
    const bool is_data_uri =
        image.uri.size() >= 5 && image.uri.rfind("data:", 0) == 0;
    if (!image.image.empty() && image.width > 0 && image.height > 0) {
        const int channels = (image.component > 0) ? image.component : 4;
        result = std::make_shared<ImageTexture>(
            image.image, image.width, image.height, channels, color_space, channel, kFlipV, wrap_s, wrap_t);
    } else if (is_data_uri) {
        if (err) *err += "Data URI image was not decoded.\n";
        result = nullptr;
    } else if (!image.uri.empty()) {
        std::filesystem::path p = base_dir / std::filesystem::path(image.uri);
        result = std::make_shared<ImageTexture>(p.string(), color_space, channel, kFlipV, wrap_s, wrap_t);
    } else {
        if (err) *err += "Unsupported image source.\n";
        result = nullptr;
    }

    cache.emplace(key, result);
    return result;
}

static MaterialBuildResult build_material(const tinygltf::Model& model,
                                         int material_index,
                                         const std::filesystem::path& base_dir,
                                         std::unordered_map<TextureKey, TexturePtr, TextureKeyHash>& texture_cache,
                                         std::unordered_map<int, NormalMapPtr>& normal_map_cache,
                                         std::string* err) {
    MaterialBuildResult out;

    if (material_index < 0 || material_index >= static_cast<int>(model.materials.size())) {
        out.material = std::make_shared<Lambertian>(std::make_shared<SolidColor>(Color(0.8f)));
        return out;
    }

    const tinygltf::Material& mat = model.materials[material_index];

    float ior = 1.5f;
    if (auto it = mat.extensions.find("KHR_materials_ior"); it != mat.extensions.end()) {
        const tinygltf::Value& ext = it->second;
        if (ext.IsObject() && ext.Has("ior") && ext.Get("ior").IsNumber()) {
            ior = static_cast<float>(ext.Get("ior").GetNumberAsDouble());
        }
    }
    ior = std::max(1.0f, ior);

    float transmission_factor = 0.0f;
    int transmission_tex_index = -1;
    if (auto it = mat.extensions.find("KHR_materials_transmission"); it != mat.extensions.end()) {
        const tinygltf::Value& ext = it->second;
        if (ext.IsObject()) {
            if (ext.Has("transmissionFactor") && ext.Get("transmissionFactor").IsNumber()) {
                transmission_factor = static_cast<float>(ext.Get("transmissionFactor").GetNumberAsDouble());
            }
            if (ext.Has("transmissionTexture")) {
                transmission_tex_index = gltf_parse_texture_index(ext.Get("transmissionTexture"));
            }
        }
    }
    transmission_factor = gltf_clamp01(transmission_factor);

    float thickness_factor = 0.0f;
    int thickness_tex_index = -1;
    Color attenuation_color(1.0f);
    float attenuation_distance = std::numeric_limits<float>::infinity();
    if (auto it = mat.extensions.find("KHR_materials_volume"); it != mat.extensions.end()) {
        const tinygltf::Value& ext = it->second;
        if (ext.IsObject()) {
            if (ext.Has("thicknessFactor") && ext.Get("thicknessFactor").IsNumber()) {
                thickness_factor = static_cast<float>(ext.Get("thicknessFactor").GetNumberAsDouble());
            }
            if (ext.Has("thicknessTexture")) {
                thickness_tex_index = gltf_parse_texture_index(ext.Get("thicknessTexture"));
            }
            if (ext.Has("attenuationColor")) {
                Color c;
                if (gltf_parse_vec3(ext.Get("attenuationColor"), c)) {
                    attenuation_color = c;
                }
            }
            if (ext.Has("attenuationDistance") && ext.Get("attenuationDistance").IsNumber()) {
                attenuation_distance = static_cast<float>(ext.Get("attenuationDistance").GetNumberAsDouble());
            }
        }
    }
    thickness_factor = std::max(0.0f, thickness_factor);
    attenuation_color = Color(gltf_clamp01(attenuation_color.x),
                              gltf_clamp01(attenuation_color.y),
                              gltf_clamp01(attenuation_color.z));
    if (!(attenuation_distance > 0.0f) || !std::isfinite(attenuation_distance)) {
        attenuation_distance = std::numeric_limits<float>::infinity();
    }

    Color base_rgb(1.0f);
    float base_alpha = 1.0f;
    if (mat.pbrMetallicRoughness.baseColorFactor.size() == 4) {
        base_rgb = Color(static_cast<float>(mat.pbrMetallicRoughness.baseColorFactor[0]),
                         static_cast<float>(mat.pbrMetallicRoughness.baseColorFactor[1]),
                         static_cast<float>(mat.pbrMetallicRoughness.baseColorFactor[2]));
        base_alpha = static_cast<float>(mat.pbrMetallicRoughness.baseColorFactor[3]);
    }

    const std::string alpha_mode = mat.alphaMode.empty() ? "OPAQUE" : mat.alphaMode;
    if (alpha_mode == "OPAQUE") {
        base_alpha = 1.0f;
    }

    TexturePtr base_tex;
    if (mat.pbrMetallicRoughness.baseColorTexture.index >= 0) {
        base_tex = load_texture(model,
                                mat.pbrMetallicRoughness.baseColorTexture.index,
                                base_dir,
                                ImageTexture::ColorSpace::sRGB,
                                -1,
                                texture_cache,
                                err);
    } else {
        base_tex = std::make_shared<SolidColor>(base_rgb);
        base_rgb = Color(1.0f);
    }

    if ((base_rgb.x != 1.0f || base_rgb.y != 1.0f || base_rgb.z != 1.0f) || base_alpha != 1.0f) {
        base_tex = std::make_shared<ColorFactorTexture>(base_tex, base_rgb, base_alpha);
    }

    if (alpha_mode == "MASK") {
        const float cutoff =
            (mat.alphaCutoff > 0.0) ? static_cast<float>(mat.alphaCutoff) : 0.5f;
        base_tex = std::make_shared<AlphaCutoffTexture>(base_tex, cutoff);
    } else if (alpha_mode == "OPAQUE") {
        // glTF: ignore baseColor alpha for opaque materials (both factor and texture channel).
        base_tex = std::make_shared<ForceAlphaTexture>(base_tex, 1.0f);
    }

    Color emissive_factor(0.0f);
    if (mat.emissiveFactor.size() == 3) {
        emissive_factor = Color(static_cast<float>(mat.emissiveFactor[0]),
                                static_cast<float>(mat.emissiveFactor[1]),
                                static_cast<float>(mat.emissiveFactor[2]));
    }

    TexturePtr emissive_tex;
    if (mat.emissiveTexture.index >= 0) {
        emissive_tex = load_texture(model,
                                    mat.emissiveTexture.index,
                                    base_dir,
                                    ImageTexture::ColorSpace::sRGB,
                                    -1,
                                    texture_cache,
                                    err);
    }

    if (emissive_tex) {
        emissive_tex = std::make_shared<ColorFactorTexture>(emissive_tex, emissive_factor);
    } else if (emissive_factor.x > 0.0f || emissive_factor.y > 0.0f || emissive_factor.z > 0.0f) {
        emissive_tex = std::make_shared<SolidColor>(emissive_factor);
    }

    const float metallic_factor = static_cast<float>(mat.pbrMetallicRoughness.metallicFactor);
    const float roughness_factor = static_cast<float>(mat.pbrMetallicRoughness.roughnessFactor);

    TexturePtr metallic_tex;
    TexturePtr roughness_tex;
    if (mat.pbrMetallicRoughness.metallicRoughnessTexture.index >= 0) {
        // glTF packs roughness in G, metallic in B.
        roughness_tex = load_texture(model,
                                     mat.pbrMetallicRoughness.metallicRoughnessTexture.index,
                                     base_dir,
                                     ImageTexture::ColorSpace::Linear,
                                     1,
                                     texture_cache,
                                     err);
        metallic_tex = load_texture(model,
                                    mat.pbrMetallicRoughness.metallicRoughnessTexture.index,
                                    base_dir,
                                    ImageTexture::ColorSpace::Linear,
                                    2,
                                    texture_cache,
                                    err);
    }

    TexturePtr occlusion_tex;
    float occlusion_strength = 1.0f;
    if (mat.occlusionTexture.index >= 0) {
        // glTF occlusion is packed in R.
        occlusion_tex = load_texture(model,
                                     mat.occlusionTexture.index,
                                     base_dir,
                                     ImageTexture::ColorSpace::Linear,
                                     0,
                                     texture_cache,
                                     err);
        if (mat.occlusionTexture.strength > 0.0) {
            occlusion_strength = static_cast<float>(mat.occlusionTexture.strength);
        }
    }

    NormalMapPtr normal_map;
    float normal_strength = 1.0f;
    if (mat.normalTexture.index >= 0) {
        normal_map = load_normal_map(model,
                                     mat.normalTexture.index,
                                     base_dir,
                                     normal_map_cache,
                                     err);
        if (mat.normalTexture.scale > 0.0) {
            normal_strength = static_cast<float>(mat.normalTexture.scale);
        }
    }

    TexturePtr transmission_tex;
    if (transmission_tex_index >= 0) {
        transmission_tex = load_texture(model,
                                        transmission_tex_index,
                                        base_dir,
                                        ImageTexture::ColorSpace::Linear,
                                        0,
                                        texture_cache,
                                        err);
    }

    TexturePtr thickness_tex;
    if (thickness_tex_index >= 0) {
        thickness_tex = load_texture(model,
                                     thickness_tex_index,
                                     base_dir,
                                     ImageTexture::ColorSpace::Linear,
                                     0,
                                     texture_cache,
                                     err);
    }

    const bool use_transmission =
        (alpha_mode == "BLEND") || (transmission_factor > 0.0f) || (transmission_tex != nullptr);

    MaterialPtr base_material;
    if (use_transmission) {
        float factor = transmission_factor;
        if (factor <= 0.0f && transmission_tex) {
            factor = 1.0f;
        }
        if (factor <= 0.0f && (alpha_mode == "BLEND")) {
            factor = 1.0f;
        }

        base_material = std::make_shared<DielectricTransmissionBSDF>(ior,
                                                                     base_tex,
                                                                     normal_map,
                                                                     normal_strength,
                                                                     factor,
                                                                     transmission_tex,
                                                                     alpha_mode == "BLEND",
                                                                     thickness_factor,
                                                                     thickness_tex,
                                                                     attenuation_color,
                                                                     attenuation_distance);
    } else {
        base_material = std::make_shared<PrincipledBSDF>(base_tex,
                                                         metallic_factor,
                                                         metallic_tex,
                                                         roughness_factor,
                                                         roughness_tex,
                                                         normal_map,
                                                         normal_strength,
                                                         occlusion_tex,
                                                         occlusion_strength);
    }

    if (emissive_tex) {
        out.material = std::make_shared<EmissiveMaterial>(base_material, emissive_tex, mat.doubleSided);
        out.emissive = true;
        return out;
    }

    out.material = base_material;
    return out;
}

static AABB compute_scene_bounds(const Scene& scene) {
    bool has_any = false;
    AABB bounds;

    for (const auto& obj : scene.objects) {
        if (!obj) continue;
        const AABB b = obj->bounding_box();
        if (!has_any) {
            bounds = b;
            has_any = true;
        } else {
            bounds = surrounding_box(bounds, b);
        }
    }

    if (!has_any) {
        return AABB(Vec3(-1.0f), Vec3(1.0f));
    }
    return bounds;
}

static void auto_camera_from_bounds(const AABB& bounds, CameraSettings& out_cam) {
    const Vec3 center = (bounds.min + bounds.max) * 0.5f;
    const Vec3 extent = bounds.max - bounds.min;
    const float radius = 0.5f * std::sqrt(extent.length_squared());
    const float dist = std::max(1.0f, radius * 2.5f);

    out_cam.look_at = center;
    out_cam.look_from = center + Vec3(0.0f, radius * 0.25f, dist);
    out_cam.up = Vec3(0.0f, 1.0f, 0.0f);
    out_cam.vertical_fov_deg = 40.0f;
}

static bool set_camera_from_node(const tinygltf::Model& model,
                                 int node_index,
                                 const Transform& world_tr,
                                 CameraSettings& out_cam) {
    if (node_index < 0 || node_index >= static_cast<int>(model.nodes.size())) {
        return false;
    }
    const tinygltf::Node& node = model.nodes[node_index];
    if (node.camera < 0 || node.camera >= static_cast<int>(model.cameras.size())) {
        return false;
    }

    const tinygltf::Camera& cam = model.cameras[node.camera];
    if (cam.type != "perspective") {
        return false;
    }

    const Vec3 pos = world_tr.apply_point(Vec3::zero());
    const Vec3 forward = normalize(world_tr.apply_vector(Vec3(0.0f, 0.0f, -1.0f)));
    const Vec3 up = normalize(world_tr.apply_vector(Vec3(0.0f, 1.0f, 0.0f)));

    out_cam.look_from = pos;
    out_cam.look_at = pos + forward;
    out_cam.up = up;

    if (cam.perspective.yfov > 0.0) {
        out_cam.vertical_fov_deg = static_cast<float>(cam.perspective.yfov) * 180.0f / kPi;
    }

    return true;
}

static MeshDataPtr build_mesh_data(const tinygltf::Model& model,
                                  const tinygltf::Primitive& prim,
                                  std::string* err) {
    auto it_pos = prim.attributes.find("POSITION");
    if (it_pos == prim.attributes.end()) {
        if (err) *err += "Primitive missing POSITION.\n";
        return nullptr;
    }

    std::vector<Vec3> positions;
    if (!read_accessor_vec3(model, it_pos->second, positions, err)) {
        return nullptr;
    }

    std::vector<Vec3> normals;
    if (auto it_n = prim.attributes.find("NORMAL"); it_n != prim.attributes.end()) {
        if (!read_accessor_vec3(model, it_n->second, normals, err)) {
            normals.clear();
        }
    }

    std::vector<Vec2> uvs;
    if (auto it_uv = prim.attributes.find("TEXCOORD_0"); it_uv != prim.attributes.end()) {
        if (!read_accessor_vec2(model, it_uv->second, uvs, err)) {
            uvs.clear();
        }
    }

    std::vector<Vec3> tangents;
    std::vector<float> tangent_signs;
    if (auto it_t = prim.attributes.find("TANGENT"); it_t != prim.attributes.end()) {
        if (!read_accessor_vec4_tangents(model, it_t->second, tangents, tangent_signs, err)) {
            tangents.clear();
            tangent_signs.clear();
        }
    }

    if (!tangents.empty() && tangents.size() != positions.size()) {
        if (err) *err += "TANGENT attribute size mismatch.\n";
        tangents.clear();
        tangent_signs.clear();
    }

    std::vector<std::uint32_t> indices;
    if (prim.indices >= 0) {
        if (!read_accessor_indices(model, prim.indices, indices, err)) {
            return nullptr;
        }
    } else {
        indices.resize(positions.size());
        for (std::uint32_t i = 0; i < static_cast<std::uint32_t>(positions.size()); ++i) {
            indices[i] = i;
        }
    }

    if (indices.size() % 3 != 0) {
        indices.resize(indices.size() - (indices.size() % 3));
    }

    if (uvs.size() != positions.size()) {
        uvs.assign(positions.size(), Vec2(0.0f, 0.0f));
    }

    if (normals.size() != positions.size()) {
        compute_normals(positions, indices, normals);
    }

    return std::make_shared<MeshData>(
        std::move(positions),
        std::move(normals),
        std::move(uvs),
        std::move(indices),
        std::move(tangents),
        std::move(tangent_signs));
}

static void traverse_node(const tinygltf::Model& model,
                          int node_index,
                          const Transform& parent,
                          Scene& out_scene,
                          std::unordered_map<int, MaterialBuildResult>& material_cache,
                          std::unordered_map<TextureKey, TexturePtr, TextureKeyHash>& texture_cache,
                          std::unordered_map<int, NormalMapPtr>& normal_map_cache,
                          std::vector<std::pair<int, Transform>>& camera_nodes,
                          const std::filesystem::path& base_dir,
                          std::string* err) {
    if (node_index < 0 || node_index >= static_cast<int>(model.nodes.size())) {
        return;
    }

    const tinygltf::Node& node = model.nodes[node_index];
    const Transform local = node_local_transform(node);
    const Transform world = parent * local;

    if (node.camera >= 0) {
        camera_nodes.emplace_back(node_index, world);
    }

    if (node.mesh >= 0 && node.mesh < static_cast<int>(model.meshes.size())) {
        const tinygltf::Mesh& mesh = model.meshes[node.mesh];

        for (const auto& prim : mesh.primitives) {
            int mode = prim.mode;
            if (mode < 0) {
                mode = TINYGLTF_MODE_TRIANGLES;
            }
            if (mode != TINYGLTF_MODE_TRIANGLES) {
                continue;
            }

            MaterialBuildResult mat_result;
            if (auto it = material_cache.find(prim.material); it != material_cache.end()) {
                mat_result = it->second;
            } else {
                mat_result = build_material(model, prim.material, base_dir, texture_cache, normal_map_cache, err);
                material_cache.emplace(prim.material, mat_result);
            }

            MeshDataPtr data = build_mesh_data(model, prim, err);
            if (!data) {
                continue;
            }

            auto obj = std::make_shared<Mesh>(data, world, mat_result.material);
            out_scene.objects.push_back(obj);
            if (mat_result.emissive) {
                out_scene.lights.add_area_light(obj);
            }
        }
    }

    for (int child : node.children) {
        traverse_node(model,
                      child,
                      world,
                      out_scene,
                      material_cache,
                      texture_cache,
                      normal_map_cache,
                      camera_nodes,
                      base_dir,
                      err);
    }
}

}  // namespace

bool load_gltf_scene(const std::string& path,
                     Scene& out_scene,
                     CameraSettings* out_camera_settings,
                     std::string* out_error,
                     const GltfLoadOptions& options) {
    out_scene = Scene();
    if (out_error) {
        out_error->clear();
    }

    tinygltf::Model model;
    tinygltf::TinyGLTF loader;
    std::string err;
    std::string warn;

    const std::filesystem::path p(path);
    const std::string ext = p.extension().string();
    bool ok = false;
    if (ext == ".glb" || ext == ".GLB") {
        ok = loader.LoadBinaryFromFile(&model, &err, &warn, path);
    } else {
        ok = loader.LoadASCIIFromFile(&model, &err, &warn, path);
    }

    if (!warn.empty()) {
        std::cerr << "glTF warning: " << warn << "\n";
    }
    if (!ok) {
        if (out_error) {
            *out_error = err;
        } else {
            std::cerr << "glTF error: " << err << "\n";
        }
        return false;
    }

    const std::filesystem::path base_dir = p.has_parent_path() ? p.parent_path() : std::filesystem::path(".");

    std::unordered_map<TextureKey, TexturePtr, TextureKeyHash> texture_cache;
    std::unordered_map<int, MaterialBuildResult> material_cache;
    std::unordered_map<int, NormalMapPtr> normal_map_cache;
    std::vector<std::pair<int, Transform>> camera_nodes;

    int scene_index = model.defaultScene >= 0 ? model.defaultScene : 0;
    if (scene_index < 0 || scene_index >= static_cast<int>(model.scenes.size())) {
        scene_index = 0;
    }

    const Transform identity = Transform::identity();
    if (!model.scenes.empty()) {
        const tinygltf::Scene& scene = model.scenes[scene_index];
        for (int node_index : scene.nodes) {
            traverse_node(model,
                          node_index,
                          identity,
                          out_scene,
                          material_cache,
                          texture_cache,
                          normal_map_cache,
                          camera_nodes,
                          base_dir,
                          out_error ? out_error : &err);
        }
    }

    if (out_camera_settings) {
        if (options.use_first_camera) {
            for (const auto& [node_index, world] : camera_nodes) {
                if (set_camera_from_node(model, node_index, world, *out_camera_settings)) {
                    return true;
                }
            }
        }

        const AABB bounds = compute_scene_bounds(out_scene);
        auto_camera_from_bounds(bounds, *out_camera_settings);
    }

    return true;
}
