#include "io/obj_loader.h"

#include <array>
#include <cctype>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

#include "core/color.h"
#include "core/vec2.h"
#include "core/vec3.h"
#include "scene/texture.h"

namespace {

struct MaterialBuildResult {
    MaterialPtr material;
    bool emissive = false;
};

struct TextureKey {
    std::string path;
    ImageTexture::ColorSpace color_space = ImageTexture::ColorSpace::sRGB;
    int channel = -1;
    bool flip_v = true;
    ImageTexture::WrapMode wrap_s = ImageTexture::WrapMode::Repeat;
    ImageTexture::WrapMode wrap_t = ImageTexture::WrapMode::Repeat;

    bool operator==(const TextureKey& other) const {
        return path == other.path &&
               color_space == other.color_space &&
               channel == other.channel &&
               flip_v == other.flip_v &&
               wrap_s == other.wrap_s &&
               wrap_t == other.wrap_t;
    }
};

struct TextureKeyHash {
    std::size_t operator()(const TextureKey& k) const noexcept {
        std::size_t h = std::hash<std::string>{}(k.path);
        h ^= (std::hash<int>{}(static_cast<int>(k.color_space)) + 0x9e3779b9 + (h << 6) + (h >> 2));
        h ^= (std::hash<int>{}(k.channel) + 0x9e3779b9 + (h << 6) + (h >> 2));
        h ^= (std::hash<int>{}(k.flip_v ? 1 : 0) + 0x9e3779b9 + (h << 6) + (h >> 2));
        h ^= (std::hash<int>{}(static_cast<int>(k.wrap_s)) + 0x9e3779b9 + (h << 6) + (h >> 2));
        h ^= (std::hash<int>{}(static_cast<int>(k.wrap_t)) + 0x9e3779b9 + (h << 6) + (h >> 2));
        return h;
    }
};

struct NormalMapKey {
    std::string path;
    bool flip_v = true;
    ImageTexture::WrapMode wrap_s = ImageTexture::WrapMode::Repeat;
    ImageTexture::WrapMode wrap_t = ImageTexture::WrapMode::Repeat;

    bool operator==(const NormalMapKey& other) const {
        return path == other.path &&
               flip_v == other.flip_v &&
               wrap_s == other.wrap_s &&
               wrap_t == other.wrap_t;
    }
};

struct NormalMapKeyHash {
    std::size_t operator()(const NormalMapKey& k) const noexcept {
        std::size_t h = std::hash<std::string>{}(k.path);
        h ^= (std::hash<int>{}(k.flip_v ? 1 : 0) + 0x9e3779b9 + (h << 6) + (h >> 2));
        h ^= (std::hash<int>{}(static_cast<int>(k.wrap_s)) + 0x9e3779b9 + (h << 6) + (h >> 2));
        h ^= (std::hash<int>{}(static_cast<int>(k.wrap_t)) + 0x9e3779b9 + (h << 6) + (h >> 2));
        return h;
    }
};

static std::string_view trim(std::string_view s) {
    while (!s.empty() && std::isspace(static_cast<unsigned char>(s.front()))) {
        s.remove_prefix(1);
    }
    while (!s.empty() && std::isspace(static_cast<unsigned char>(s.back()))) {
        s.remove_suffix(1);
    }
    return s;
}

static bool next_token(std::string_view& in, std::string_view& out) {
    in = trim(in);
    if (in.empty()) {
        return false;
    }
    std::size_t i = 0;
    while (i < in.size() && !std::isspace(static_cast<unsigned char>(in[i]))) {
        ++i;
    }
    out = in.substr(0, i);
    in.remove_prefix(i);
    return true;
}

static bool parse_int(std::string_view s, int& out) {
    s = trim(s);
    if (s.empty()) {
        return false;
    }
    int value = 0;
    const char* begin = s.data();
    const char* end = s.data() + s.size();
    auto [ptr, ec] = std::from_chars(begin, end, value);
    if (ec != std::errc() || ptr != end) {
        return false;
    }
    out = value;
    return true;
}

static bool parse_float(std::string_view s, float& out) {
    s = trim(s);
    if (s.empty()) {
        return false;
    }
    char buf[128];
    if (s.size() >= sizeof(buf)) {
        return false;
    }
    std::memcpy(buf, s.data(), s.size());
    buf[s.size()] = '\0';
    char* end = nullptr;
    const float v = std::strtof(buf, &end);
    if (!end || end != buf + s.size()) {
        return false;
    }
    if (!std::isfinite(v)) {
        return false;
    }
    out = v;
    return true;
}

static bool parse_vec3(std::string_view in, Vec3& out) {
    std::string_view t;
    float x = 0.0f, y = 0.0f, z = 0.0f;
    if (!next_token(in, t) || !parse_float(t, x)) return false;
    if (!next_token(in, t) || !parse_float(t, y)) return false;
    if (!next_token(in, t) || !parse_float(t, z)) return false;
    out = Vec3(x, y, z);
    return true;
}

static bool parse_vec2(std::string_view in, Vec2& out) {
    std::string_view t;
    float x = 0.0f, y = 0.0f;
    if (!next_token(in, t) || !parse_float(t, x)) return false;
    if (!next_token(in, t) || !parse_float(t, y)) return false;
    out = Vec2(x, y);
    return true;
}

static std::string unquote(std::string_view s) {
    s = trim(s);
    if (s.size() >= 2 && ((s.front() == '"' && s.back() == '"') || (s.front() == '\'' && s.back() == '\''))) {
        s = s.substr(1, s.size() - 2);
    }
    return std::string(s);
}

static float cross2(const Vec2& a, const Vec2& b) {
    return a.x * b.y - a.y * b.x;
}

static bool point_in_triangle(const Vec2& p,
                              const Vec2& a,
                              const Vec2& b,
                              const Vec2& c,
                              float winding_sign) {
    constexpr float eps = 1e-8f;
    const float c1 = cross2(b - a, p - a);
    const float c2 = cross2(c - b, p - b);
    const float c3 = cross2(a - c, p - c);
    return winding_sign * c1 >= -eps &&
           winding_sign * c2 >= -eps &&
           winding_sign * c3 >= -eps;
}

static bool triangulate_polygon(const std::vector<Vec3>& positions,
                                const std::vector<int>& pos_indices,
                                std::vector<std::array<int, 3>>& out_tris) {
    out_tris.clear();
    const int n = static_cast<int>(pos_indices.size());
    if (n < 3) {
        return false;
    }
    if (n == 3) {
        out_tris.push_back({0, 1, 2});
        return true;
    }

    Vec3 normal(0.0f);
    for (int i = 0; i < n; ++i) {
        const Vec3& cur = positions[static_cast<std::size_t>(pos_indices[static_cast<std::size_t>(i)])];
        const Vec3& next = positions[static_cast<std::size_t>(pos_indices[static_cast<std::size_t>((i + 1) % n)])];
        normal.x += (cur.y - next.y) * (cur.z + next.z);
        normal.y += (cur.z - next.z) * (cur.x + next.x);
        normal.z += (cur.x - next.x) * (cur.y + next.y);
    }

    const float ax = std::fabs(normal.x);
    const float ay = std::fabs(normal.y);
    const float az = std::fabs(normal.z);
    int axis = 2;
    if (ax >= ay && ax >= az) {
        axis = 0;
    } else if (ay >= ax && ay >= az) {
        axis = 1;
    }

    std::vector<Vec2> poly2;
    poly2.reserve(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        const Vec3& p = positions[static_cast<std::size_t>(pos_indices[static_cast<std::size_t>(i)])];
        if (axis == 0) {
            poly2.emplace_back(p.y, p.z);
        } else if (axis == 1) {
            poly2.emplace_back(p.x, p.z);
        } else {
            poly2.emplace_back(p.x, p.y);
        }
    }

    float area2 = 0.0f;
    for (int i = 0; i < n; ++i) {
        const Vec2& a = poly2[static_cast<std::size_t>(i)];
        const Vec2& b = poly2[static_cast<std::size_t>((i + 1) % n)];
        area2 += a.x * b.y - b.x * a.y;
    }
    if (std::fabs(area2) < 1e-12f) {
        out_tris.reserve(static_cast<std::size_t>(n - 2));
        for (int i = 1; i + 1 < n; ++i) {
            out_tris.push_back({0, i, i + 1});
        }
        return true;
    }
    const float winding_sign = (area2 >= 0.0f) ? 1.0f : -1.0f;

    std::vector<int> v;
    v.reserve(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        v.push_back(i);
    }

    out_tris.reserve(static_cast<std::size_t>(n - 2));
    int guard = 0;
    while (static_cast<int>(v.size()) > 3 && guard++ < n * n) {
        bool clipped = false;
        const int m = static_cast<int>(v.size());
        for (int i = 0; i < m; ++i) {
            const int i_prev = v[static_cast<std::size_t>((i + m - 1) % m)];
            const int i_curr = v[static_cast<std::size_t>(i)];
            const int i_next = v[static_cast<std::size_t>((i + 1) % m)];

            const Vec2& p_prev = poly2[static_cast<std::size_t>(i_prev)];
            const Vec2& p_curr = poly2[static_cast<std::size_t>(i_curr)];
            const Vec2& p_next = poly2[static_cast<std::size_t>(i_next)];

            const float cross = cross2(p_next - p_curr, p_prev - p_curr);
            if (winding_sign * cross <= 1e-12f) {
                continue;
            }

            bool contains = false;
            for (int j = 0; j < m; ++j) {
                const int idx = v[static_cast<std::size_t>(j)];
                if (idx == i_prev || idx == i_curr || idx == i_next) {
                    continue;
                }
                if (point_in_triangle(poly2[static_cast<std::size_t>(idx)], p_prev, p_curr, p_next, winding_sign)) {
                    contains = true;
                    break;
                }
            }
            if (contains) {
                continue;
            }

            out_tris.push_back({i_prev, i_curr, i_next});
            v.erase(v.begin() + i);
            clipped = true;
            break;
        }
        if (!clipped) {
            break;
        }
    }

    if (static_cast<int>(v.size()) == 3) {
        out_tris.push_back({v[0], v[1], v[2]});
        return true;
    }

    out_tris.clear();
    out_tris.reserve(static_cast<std::size_t>(n - 2));
    for (int i = 1; i + 1 < n; ++i) {
        out_tris.push_back({0, i, i + 1});
    }
    return true;
}

struct ObjVertexRef {
    int v = -1;
    int vt = -1;
    int vn = -1;
};

static bool resolve_index(int idx_1_based, int count, int& out_0_based) {
    if (idx_1_based > 0) {
        out_0_based = idx_1_based - 1;
    } else if (idx_1_based < 0) {
        out_0_based = count + idx_1_based;
    } else {
        return false;
    }
    return out_0_based >= 0 && out_0_based < count;
}

static bool parse_vertex_ref(std::string_view tok,
                             int v_count,
                             int vt_count,
                             int vn_count,
                             ObjVertexRef& out) {
    out = ObjVertexRef();
    tok = trim(tok);
    if (tok.empty()) {
        return false;
    }

    const std::size_t s1 = tok.find('/');
    std::string_view a = tok;
    std::string_view b;
    std::string_view c;
    if (s1 != std::string_view::npos) {
        a = tok.substr(0, s1);
        const std::string_view rest = tok.substr(s1 + 1);
        const std::size_t s2 = rest.find('/');
        if (s2 != std::string_view::npos) {
            b = rest.substr(0, s2);
            c = rest.substr(s2 + 1);
        } else {
            b = rest;
        }
    }

    int vi = 0;
    if (!parse_int(a, vi) || !resolve_index(vi, v_count, out.v)) {
        return false;
    }

    if (!b.empty()) {
        int vti = 0;
        if (!parse_int(b, vti) || !resolve_index(vti, vt_count, out.vt)) {
            return false;
        }
    }

    if (!c.empty()) {
        int vni = 0;
        if (!parse_int(c, vni) || !resolve_index(vni, vn_count, out.vn)) {
            return false;
        }
    }

    return true;
}

struct VertexKey {
    int v = -1;
    int vt = -1;
    int vn = -1;
    int smooth = 0;

    bool operator==(const VertexKey& other) const {
        return v == other.v && vt == other.vt && vn == other.vn && smooth == other.smooth;
    }
};

struct VertexKeyHash {
    std::size_t operator()(const VertexKey& k) const noexcept {
        std::size_t h = std::hash<int>{}(k.v);
        h ^= (std::hash<int>{}(k.vt) + 0x9e3779b9 + (h << 6) + (h >> 2));
        h ^= (std::hash<int>{}(k.vn) + 0x9e3779b9 + (h << 6) + (h >> 2));
        h ^= (std::hash<int>{}(k.smooth) + 0x9e3779b9 + (h << 6) + (h >> 2));
        return h;
    }
};

struct MeshBuilder {
    std::string object;
    std::string group;
    std::string material;

    std::vector<Vec3> positions;
    std::vector<Vec3> normals;
    std::vector<Vec2> uvs;
    std::vector<std::uint32_t> indices;
    std::unordered_map<VertexKey, std::uint32_t, VertexKeyHash> vert_map;

    bool any_normals = false;
    bool all_normals = true;

    std::uint32_t add_vertex(const ObjVertexRef& ref,
                             int smooth_key,
                             const std::vector<Vec3>& src_positions,
                             const std::vector<Vec2>& src_uvs,
                             const std::vector<Vec3>& src_normals) {
        VertexKey key;
        key.v = ref.v;
        key.vt = ref.vt;
        key.vn = ref.vn;
        key.smooth = (ref.vn >= 0) ? 0 : smooth_key;

        if (auto it = vert_map.find(key); it != vert_map.end()) {
            return it->second;
        }

        const std::uint32_t out_index = static_cast<std::uint32_t>(positions.size());
        positions.push_back(src_positions[static_cast<std::size_t>(ref.v)]);

        if (ref.vt >= 0 && ref.vt < static_cast<int>(src_uvs.size())) {
            uvs.push_back(src_uvs[static_cast<std::size_t>(ref.vt)]);
        } else {
            uvs.push_back(Vec2::zero());
        }

        if (ref.vn >= 0 && ref.vn < static_cast<int>(src_normals.size())) {
            normals.push_back(src_normals[static_cast<std::size_t>(ref.vn)]);
            any_normals = true;
        } else {
            normals.push_back(Vec3::zero());
            all_normals = false;
        }

        vert_map.emplace(key, out_index);
        return out_index;
    }
};

struct MtlMapInfo {
    std::string filename;
    bool clamp = false;
    float bump_multiplier = 1.0f;
};

struct MtlMaterial {
    std::string name;
    Color kd = Color(0.8f);
    Color ke = Color(0.0f);
    float ns = 0.0f;
    float d = 1.0f;
    bool has_d = false;
    float ni = 1.5f;
    int illum = 2;

    float roughness = -1.0f;  // extension: Pr
    float metallic = -1.0f;   // extension: Pm

    MtlMapInfo map_kd;
    MtlMapInfo map_d;
    MtlMapInfo map_bump;
    MtlMapInfo map_ke;
    MtlMapInfo map_pr;
    MtlMapInfo map_pm;
    MtlMapInfo map_po;
};

static bool is_nonzero(const Color& c) {
    return c.x > 0.0f || c.y > 0.0f || c.z > 0.0f;
}

static bool parse_mtl_map(std::string_view rest, MtlMapInfo& out) {
    out = MtlMapInfo();
    rest = trim(rest);
    if (rest.empty()) {
        return false;
    }

    while (!rest.empty()) {
        std::string_view tok;
        if (!next_token(rest, tok)) {
            break;
        }
        if (tok.empty()) {
            continue;
        }

        if (tok[0] != '-') {
            std::string path(tok);
            rest = trim(rest);
            if (!rest.empty()) {
                path += " ";
                path += std::string(rest);
            }
            out.filename = unquote(path);
            return !out.filename.empty();
        }

        if (tok == "-clamp") {
            std::string_view v;
            if (next_token(rest, v)) {
                out.clamp = (v == "on" || v == "1" || v == "true");
            }
        } else if (tok == "-bm") {
            std::string_view v;
            float bm = 1.0f;
            if (next_token(rest, v) && parse_float(v, bm)) {
                out.bump_multiplier = bm;
            }
        } else if (tok == "-o" || tok == "-s" || tok == "-t") {
            for (int i = 0; i < 3; ++i) {
                std::string_view v;
                if (!next_token(rest, v)) {
                    break;
                }
                float tmp = 0.0f;
                if (!parse_float(v, tmp)) {
                    break;
                }
            }
        } else if (tok == "-mm") {
            std::string_view a, b;
            next_token(rest, a);
            next_token(rest, b);
        } else if (tok == "-texres" || tok == "-cc" || tok == "-blendu" || tok == "-blendv" || tok == "-imfchan" ||
                   tok == "-type" || tok == "-boost") {
            std::string_view v;
            next_token(rest, v);
        } else {
            std::string_view v;
            next_token(rest, v);
        }
    }

    return false;
}

static bool load_mtl_file(const std::filesystem::path& mtl_path,
                          std::unordered_map<std::string, MtlMaterial>& inout_materials,
                          std::string* out_error) {
    std::ifstream f(mtl_path);
    if (!f) {
        if (out_error) {
            *out_error += "Failed to open MTL file: " + mtl_path.string() + "\n";
        }
        return false;
    }

    MtlMaterial* cur = nullptr;

    std::string line;
    while (std::getline(f, line)) {
        std::string_view v(line);
        const std::size_t hash = v.find('#');
        if (hash != std::string_view::npos) {
            v = v.substr(0, hash);
        }
        v = trim(v);
        if (v.empty()) {
            continue;
        }

        std::string_view key;
        if (!next_token(v, key)) {
            continue;
        }

        if (key == "newmtl") {
            const std::string name = std::string(trim(v));
            if (name.empty()) {
                cur = nullptr;
                continue;
            }
            MtlMaterial mat;
            mat.name = name;
            inout_materials[name] = mat;
            cur = &inout_materials[name];
            continue;
        }

        if (!cur) {
            continue;
        }

        if (key == "Kd") {
            Vec3 c;
            if (parse_vec3(v, c)) {
                cur->kd = Color(c.x, c.y, c.z);
            }
        } else if (key == "Ke") {
            Vec3 c;
            if (parse_vec3(v, c)) {
                cur->ke = Color(c.x, c.y, c.z);
            }
        } else if (key == "Ns") {
            float ns = 0.0f;
            std::string_view t;
            if (next_token(v, t) && parse_float(t, ns)) {
                cur->ns = std::max(0.0f, ns);
            }
        } else if (key == "Ni") {
            float ni = 1.5f;
            std::string_view t;
            if (next_token(v, t) && parse_float(t, ni)) {
                if (ni >= 1.0f && std::isfinite(ni)) {
                    cur->ni = ni;
                }
            }
        } else if (key == "d") {
            float d = 1.0f;
            std::string_view t;
            if (next_token(v, t) && parse_float(t, d)) {
                cur->d = clamp_float(d, 0.0f, 1.0f);
                cur->has_d = true;
            }
        } else if (key == "Tr") {
            float tr = 0.0f;
            std::string_view t;
            if (next_token(v, t) && parse_float(t, tr)) {
                if (!cur->has_d) {
                    cur->d = clamp_float(1.0f - tr, 0.0f, 1.0f);
                }
            }
        } else if (key == "illum") {
            int illum = 2;
            std::string_view t;
            if (next_token(v, t) && parse_int(t, illum)) {
                cur->illum = illum;
            }
        } else if (key == "map_Kd") {
            parse_mtl_map(v, cur->map_kd);
        } else if (key == "map_d") {
            parse_mtl_map(v, cur->map_d);
        } else if (key == "map_Ke") {
            parse_mtl_map(v, cur->map_ke);
        } else if (key == "map_Bump" || key == "bump") {
            parse_mtl_map(v, cur->map_bump);
        } else if (key == "Pr") {
            float r = 0.0f;
            std::string_view t;
            if (next_token(v, t) && parse_float(t, r)) {
                cur->roughness = clamp_float(r, 0.0f, 1.0f);
            }
        } else if (key == "Pm") {
            float m = 0.0f;
            std::string_view t;
            if (next_token(v, t) && parse_float(t, m)) {
                cur->metallic = clamp_float(m, 0.0f, 1.0f);
            }
        } else if (key == "map_Pr") {
            parse_mtl_map(v, cur->map_pr);
        } else if (key == "map_Pm") {
            parse_mtl_map(v, cur->map_pm);
        } else if (key == "map_Po" || key == "map_ao") {
            parse_mtl_map(v, cur->map_po);
        }
    }

    return true;
}

static ImageTexture::WrapMode wrap_from_clamp(bool clamp) {
    return clamp ? ImageTexture::WrapMode::ClampToEdge : ImageTexture::WrapMode::Repeat;
}

static TexturePtr load_texture(const std::filesystem::path& base_dir,
                               const MtlMapInfo& map,
                               ImageTexture::ColorSpace color_space,
                               int channel,
                               std::unordered_map<TextureKey, TexturePtr, TextureKeyHash>& cache) {
    if (map.filename.empty()) {
        return nullptr;
    }

    std::filesystem::path p = std::filesystem::path(map.filename);
    if (p.is_relative()) {
        p = base_dir / p;
    }
    p = p.lexically_normal();

    TextureKey key;
    key.path = p.string();
    key.color_space = color_space;
    key.channel = channel;
    key.flip_v = true;
    key.wrap_s = wrap_from_clamp(map.clamp);
    key.wrap_t = wrap_from_clamp(map.clamp);

    if (auto it = cache.find(key); it != cache.end()) {
        return it->second;
    }

    TexturePtr tex = std::make_shared<ImageTexture>(key.path, key.color_space, key.channel, key.flip_v, key.wrap_s, key.wrap_t);
    cache.emplace(key, tex);
    return tex;
}

static NormalMapPtr load_normal_map(const std::filesystem::path& base_dir,
                                    const MtlMapInfo& map,
                                    std::unordered_map<NormalMapKey, NormalMapPtr, NormalMapKeyHash>& cache) {
    if (map.filename.empty()) {
        return nullptr;
    }

    std::filesystem::path p = std::filesystem::path(map.filename);
    if (p.is_relative()) {
        p = base_dir / p;
    }
    p = p.lexically_normal();

    NormalMapKey key;
    key.path = p.string();
    key.flip_v = true;
    key.wrap_s = wrap_from_clamp(map.clamp);
    key.wrap_t = wrap_from_clamp(map.clamp);

    if (auto it = cache.find(key); it != cache.end()) {
        return it->second;
    }

    NormalMapPtr nm = std::make_shared<NormalMapTexture>(key.path, key.flip_v, key.wrap_s, key.wrap_t);
    cache.emplace(key, nm);
    return nm;
}

static float ns_to_roughness(float ns) {
    if (!(ns > 0.0f) || !std::isfinite(ns)) {
        return 1.0f;
    }
    const float alpha = std::sqrt(2.0f / (ns + 2.0f));
    const float r = std::sqrt(std::max(0.0f, alpha));
    return clamp_float(r, 0.0f, 1.0f);
}

static bool illum_is_transmissive(int illum) {
    switch (illum) {
        case 4:
        case 6:
        case 7:
        case 9:
            return true;
        default:
            return false;
    }
}

class AlphaMultiplyTexture : public Texture {
public:
    AlphaMultiplyTexture(const TexturePtr& base, const TexturePtr& mask)
        : base_(base), mask_(mask) {}

    Color value(const HitRecord& hit) const override {
        return base_ ? base_->value(hit) : Color(1.0f);
    }

    float alpha(const HitRecord& hit) const override {
        const float a_base = base_ ? base_->alpha(hit) : 1.0f;
        float a_mask = 1.0f;
        if (mask_) {
            const Color c = mask_->value(hit);
            a_mask = clamp_float(c.x, 0.0f, 1.0f);
        }
        return clamp_float(a_base * a_mask, 0.0f, 1.0f);
    }

private:
    TexturePtr base_;
    TexturePtr mask_;
};

static MaterialBuildResult build_material(const std::filesystem::path& base_dir,
                                         const MtlMaterial* spec,
                                         const ObjLoadOptions& options,
                                         std::unordered_map<TextureKey, TexturePtr, TextureKeyHash>& tex_cache,
                                         std::unordered_map<NormalMapKey, NormalMapPtr, NormalMapKeyHash>& nm_cache) {
    MaterialBuildResult out;

    if (!spec) {
        out.material = std::make_shared<Lambertian>(std::make_shared<SolidColor>(Color(0.8f)));
        return out;
    }

    TexturePtr base_tex;
    if (!spec->map_kd.filename.empty()) {
        base_tex = load_texture(base_dir, spec->map_kd, ImageTexture::ColorSpace::sRGB, -1, tex_cache);
        if ((spec->kd.x != 1.0f || spec->kd.y != 1.0f || spec->kd.z != 1.0f) || spec->d != 1.0f) {
            base_tex = std::make_shared<ColorFactorTexture>(base_tex, spec->kd, spec->d);
        }
    } else {
        base_tex = std::make_shared<SolidColor>(spec->kd);
        if (spec->d != 1.0f) {
            base_tex = std::make_shared<ForceAlphaTexture>(base_tex, spec->d);
        }
    }

    if (!spec->map_d.filename.empty()) {
        TexturePtr mask = load_texture(base_dir, spec->map_d, ImageTexture::ColorSpace::Linear, 0, tex_cache);
        base_tex = std::make_shared<AlphaMultiplyTexture>(base_tex, mask);
    }

    NormalMapPtr normal_map;
    float normal_strength = 1.0f;
    if (!spec->map_bump.filename.empty()) {
        normal_map = load_normal_map(base_dir, spec->map_bump, nm_cache);
        normal_strength = std::max(0.0f, spec->map_bump.bump_multiplier);
    }

    TexturePtr emissive_tex;
    if (!spec->map_ke.filename.empty()) {
        emissive_tex = load_texture(base_dir, spec->map_ke, ImageTexture::ColorSpace::sRGB, -1, tex_cache);
        if (emissive_tex && is_nonzero(spec->ke)) {
            emissive_tex = std::make_shared<ColorFactorTexture>(emissive_tex, spec->ke);
        }
    } else if (is_nonzero(spec->ke)) {
        emissive_tex = std::make_shared<SolidColor>(spec->ke);
    }

    float metallic_factor = 0.0f;
    TexturePtr metallic_tex;
    float roughness_factor = ns_to_roughness(spec->ns);
    TexturePtr roughness_tex;
    TexturePtr occlusion_tex;

    if (spec->metallic >= 0.0f) {
        metallic_factor = clamp_float(spec->metallic, 0.0f, 1.0f);
    }
    if (!spec->map_pm.filename.empty()) {
        metallic_tex = load_texture(base_dir, spec->map_pm, ImageTexture::ColorSpace::Linear, 0, tex_cache);
        if (metallic_factor <= 0.0f) {
            metallic_factor = 1.0f;
        }
    }

    if (spec->roughness >= 0.0f) {
        roughness_factor = clamp_float(spec->roughness, 0.0f, 1.0f);
    }
    if (!spec->map_pr.filename.empty()) {
        roughness_tex = load_texture(base_dir, spec->map_pr, ImageTexture::ColorSpace::Linear, 0, tex_cache);
        if (roughness_factor <= 0.0f) {
            roughness_factor = 1.0f;
        }
    }

    if (!spec->map_po.filename.empty()) {
        occlusion_tex = load_texture(base_dir, spec->map_po, ImageTexture::ColorSpace::Linear, 0, tex_cache);
    }

    const bool treat_as_transmission =
        options.use_illum_for_transmission &&
        illum_is_transmissive(spec->illum) &&
        (spec->d < 1.0f || !spec->map_d.filename.empty());

    MaterialPtr base_material;
    if (treat_as_transmission) {
        const float ior = std::max(1.0f, spec->ni);
        base_material = std::make_shared<DielectricTransmissionBSDF>(ior,
                                                                     base_tex,
                                                                     normal_map,
                                                                     normal_strength,
                                                                     1.0f,
                                                                     nullptr,
                                                                     true);
    } else {
        base_material = std::make_shared<PrincipledBSDF>(base_tex,
                                                         metallic_factor,
                                                         metallic_tex,
                                                         roughness_factor,
                                                         roughness_tex,
                                                         normal_map,
                                                         normal_strength,
                                                         occlusion_tex,
                                                         1.0f);
    }

    if (emissive_tex) {
        out.material = std::make_shared<EmissiveMaterial>(base_material, emissive_tex, false);
        out.emissive = true;
        return out;
    }

    out.material = base_material;
    return out;
}

static std::string make_builder_key(const std::string& object,
                                   const std::string& group,
                                   const std::string& material,
                                   const ObjLoadOptions& options) {
    std::string key;
    if (options.split_by_object) {
        key += object;
        key += '\n';
    }
    if (options.split_by_group) {
        key += group;
        key += '\n';
    }
    key += material;
    return key;
}

static std::string make_mesh_name(const std::filesystem::path& obj_path,
                                 const std::string& object,
                                 const std::string& group,
                                 const std::string& material) {
    std::string name = obj_path.stem().string();
    if (!object.empty()) {
        name += "::";
        name += object;
    }
    if (!group.empty()) {
        name += "/";
        name += group;
    }
    if (!material.empty()) {
        name += " [";
        name += material;
        name += "]";
    }
    return name;
}

}  // namespace

bool load_obj_meshes(const std::string& path,
                     std::vector<ObjMesh>& out_meshes,
                     std::string* out_error,
                     const ObjLoadOptions& options) {
    out_meshes.clear();
    if (out_error) {
        out_error->clear();
    }

    std::ifstream f(path);
    if (!f) {
        if (out_error) {
            *out_error = "Failed to open OBJ file: " + path + "\n";
        }
        return false;
    }

    const std::filesystem::path obj_path(path);
    const std::filesystem::path base_dir =
        obj_path.has_parent_path() ? obj_path.parent_path() : std::filesystem::path(".");

    std::unordered_map<std::string, MtlMaterial> mtl_materials;

    std::vector<Vec3> src_positions;
    std::vector<Vec2> src_uvs;
    std::vector<Vec3> src_normals;
    src_positions.reserve(1024);
    src_uvs.reserve(1024);
    src_normals.reserve(1024);

    std::unordered_map<std::string, MeshBuilder> builders;

    std::string current_object;
    std::string current_group;
    std::string current_material;
    int current_smooth = options.default_smooth_shading ? 1 : 0;
    int face_id = 0;

    auto get_builder = [&](const std::string& object,
                           const std::string& group,
                           const std::string& material) -> MeshBuilder& {
        const std::string key = make_builder_key(object, group, material, options);
        auto it = builders.find(key);
        if (it != builders.end()) {
            return it->second;
        }
        MeshBuilder b;
        b.object = object;
        b.group = group;
        b.material = material;
        auto [ins_it, _] = builders.emplace(key, std::move(b));
        return ins_it->second;
    };

    std::string line;
    int line_number = 0;
    while (std::getline(f, line)) {
        ++line_number;

        while (!line.empty() && line.back() == '\\') {
            line.pop_back();
            std::string cont;
            if (!std::getline(f, cont)) {
                break;
            }
            ++line_number;
            line += cont;
        }

        std::string_view v(line);
        const std::size_t hash = v.find('#');
        if (hash != std::string_view::npos) {
            v = v.substr(0, hash);
        }
        v = trim(v);
        if (v.empty()) {
            continue;
        }

        std::string_view keyword;
        if (!next_token(v, keyword)) {
            continue;
        }

        if (keyword == "v") {
            Vec3 p;
            if (parse_vec3(v, p)) {
                src_positions.push_back(p);
            } else if (out_error) {
                *out_error += "OBJ parse error (v) at line " + std::to_string(line_number) + "\n";
            }
        } else if (keyword == "vt") {
            Vec2 uv;
            if (parse_vec2(v, uv)) {
                src_uvs.push_back(uv);
            } else if (out_error) {
                *out_error += "OBJ parse error (vt) at line " + std::to_string(line_number) + "\n";
            }
        } else if (keyword == "vn") {
            Vec3 n;
            if (parse_vec3(v, n)) {
                src_normals.push_back(n);
            } else if (out_error) {
                *out_error += "OBJ parse error (vn) at line " + std::to_string(line_number) + "\n";
            }
        } else if (keyword == "f") {
            ++face_id;

            std::vector<ObjVertexRef> poly;
            poly.reserve(8);
            std::vector<int> poly_pos_indices;
            poly_pos_indices.reserve(8);

            bool ok = true;
            std::string_view tok;
            while (next_token(v, tok)) {
                ObjVertexRef ref;
                if (!parse_vertex_ref(tok,
                                      static_cast<int>(src_positions.size()),
                                      static_cast<int>(src_uvs.size()),
                                      static_cast<int>(src_normals.size()),
                                      ref)) {
                    ok = false;
                    break;
                }
                poly.push_back(ref);
                poly_pos_indices.push_back(ref.v);
            }

            if (!ok || poly.size() < 3) {
                if (out_error) {
                    *out_error += "OBJ parse error (f) at line " + std::to_string(line_number) + "\n";
                }
                continue;
            }

            std::vector<std::array<int, 3>> tri_indices;
            triangulate_polygon(src_positions, poly_pos_indices, tri_indices);
            if (tri_indices.empty()) {
                continue;
            }

            const int smooth_key = (current_smooth == 0) ? -face_id : current_smooth;
            MeshBuilder& b = get_builder(current_object, current_group, current_material);
            b.indices.reserve(b.indices.size() + tri_indices.size() * 3u);

            for (const auto& tri : tri_indices) {
                const ObjVertexRef& r0 = poly[static_cast<std::size_t>(tri[0])];
                const ObjVertexRef& r1 = poly[static_cast<std::size_t>(tri[1])];
                const ObjVertexRef& r2 = poly[static_cast<std::size_t>(tri[2])];

                const std::uint32_t i0 = b.add_vertex(r0, smooth_key, src_positions, src_uvs, src_normals);
                const std::uint32_t i1 = b.add_vertex(r1, smooth_key, src_positions, src_uvs, src_normals);
                const std::uint32_t i2 = b.add_vertex(r2, smooth_key, src_positions, src_uvs, src_normals);

                if (i0 == i1 || i0 == i2 || i1 == i2) {
                    continue;
                }

                b.indices.push_back(i0);
                b.indices.push_back(i1);
                b.indices.push_back(i2);
            }
        } else if (keyword == "o") {
            current_object = std::string(trim(v));
        } else if (keyword == "g") {
            current_group = std::string(trim(v));
        } else if (keyword == "usemtl") {
            current_material = std::string(trim(v));
        } else if (keyword == "mtllib") {
            std::string_view rest = v;
            std::string_view tok2;
            while (next_token(rest, tok2)) {
                const std::string mtl_name = unquote(tok2);
                if (mtl_name.empty()) {
                    continue;
                }
                std::filesystem::path mtl_path = std::filesystem::path(mtl_name);
                if (mtl_path.is_relative()) {
                    mtl_path = base_dir / mtl_path;
                }
                mtl_path = mtl_path.lexically_normal();
                load_mtl_file(mtl_path, mtl_materials, out_error);
            }
        } else if (keyword == "s") {
            std::string_view tok2;
            if (next_token(v, tok2)) {
                if (tok2 == "off" || tok2 == "0") {
                    current_smooth = 0;
                } else {
                    int g = 1;
                    if (parse_int(tok2, g) && g > 0) {
                        current_smooth = g;
                    } else {
                        current_smooth = options.default_smooth_shading ? 1 : 0;
                    }
                }
            }
        }
    }

    if (builders.empty()) {
        if (out_error) {
            *out_error += "OBJ contains no mesh data.\n";
        }
        return false;
    }

    std::unordered_map<TextureKey, TexturePtr, TextureKeyHash> tex_cache;
    std::unordered_map<NormalMapKey, NormalMapPtr, NormalMapKeyHash> nm_cache;
    std::unordered_map<std::string, MaterialBuildResult> mat_cache;

    auto get_material = [&](const std::string& name) -> MaterialBuildResult {
        if (auto it = mat_cache.find(name); it != mat_cache.end()) {
            return it->second;
        }

        const MtlMaterial* spec = nullptr;
        if (!name.empty()) {
            if (auto it = mtl_materials.find(name); it != mtl_materials.end()) {
                spec = &it->second;
            }
        }

        MaterialBuildResult r = build_material(base_dir, spec, options, tex_cache, nm_cache);
        mat_cache.emplace(name, r);
        return r;
    };

    out_meshes.reserve(builders.size());
    for (auto& [key, b] : builders) {
        (void)key;
        if (b.positions.empty() || b.indices.empty()) {
            continue;
        }

        std::vector<Vec3> normals;
        if (b.any_normals && b.all_normals && b.normals.size() == b.positions.size()) {
            normals = std::move(b.normals);
        }

        std::vector<Vec2> uvs = std::move(b.uvs);
        if (uvs.size() != b.positions.size()) {
            uvs.assign(b.positions.size(), Vec2::zero());
        }

        MeshDataPtr data = std::make_shared<MeshData>(std::move(b.positions),
                                                      std::move(normals),
                                                      std::move(uvs),
                                                      std::move(b.indices));

        MaterialBuildResult mat = get_material(b.material);

        ObjMesh mesh;
        mesh.data = std::move(data);
        mesh.material = mat.material;
        mesh.emissive = mat.emissive;
        mesh.name = make_mesh_name(obj_path, b.object, b.group, b.material);
        out_meshes.push_back(std::move(mesh));
    }

    if (out_meshes.empty()) {
        if (out_error) {
            *out_error += "OBJ contains no triangles.\n";
        }
        return false;
    }

    return true;
}

bool append_obj_to_scene(const std::string& path,
                         Scene& inout_scene,
                         const Transform& transform,
                         std::string* out_error,
                         const ObjLoadOptions& options) {
    std::vector<ObjMesh> meshes;
    if (!load_obj_meshes(path, meshes, out_error, options)) {
        return false;
    }

    for (const auto& m : meshes) {
        if (!m.data) {
            continue;
        }
        auto obj = std::make_shared<Mesh>(m.data, transform, m.material);
        inout_scene.objects.push_back(obj);
        if (m.emissive) {
            inout_scene.lights.add_area_light(obj);
        }
    }

    return true;
}
