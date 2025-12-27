#include "io/gltf_loader.h"

#include <cerrno>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <filesystem>

#include "core/color.h"
#include "scene/mesh.h"
#include "scene/texture.h"

namespace {

struct Json {
    enum class Type {
        Null,
        Bool,
        Number,
        String,
        Array,
        Object
    };

    using Array = std::vector<Json>;
    using Object = std::unordered_map<std::string, Json>;

    Type type = Type::Null;
    bool bool_value = false;
    double number_value = 0.0;
    std::string string_value;
    Array array_value;
    Object object_value;

    static Json make_null() { return Json(); }
    static Json make_bool(bool v) {
        Json j;
        j.type = Type::Bool;
        j.bool_value = v;
        return j;
    }
    static Json make_number(double v) {
        Json j;
        j.type = Type::Number;
        j.number_value = v;
        return j;
    }
    static Json make_string(std::string v) {
        Json j;
        j.type = Type::String;
        j.string_value = std::move(v);
        return j;
    }
    static Json make_array(Array v) {
        Json j;
        j.type = Type::Array;
        j.array_value = std::move(v);
        return j;
    }
    static Json make_object(Object v) {
        Json j;
        j.type = Type::Object;
        j.object_value = std::move(v);
        return j;
    }

    bool is_null() const { return type == Type::Null; }
    bool is_bool() const { return type == Type::Bool; }
    bool is_number() const { return type == Type::Number; }
    bool is_string() const { return type == Type::String; }
    bool is_array() const { return type == Type::Array; }
    bool is_object() const { return type == Type::Object; }

    const Json* get(const std::string& key) const {
        if (!is_object()) {
            return nullptr;
        }
        auto it = object_value.find(key);
        if (it == object_value.end()) {
            return nullptr;
        }
        return &it->second;
    }

    const Json* at(std::size_t idx) const {
        if (!is_array()) {
            return nullptr;
        }
        if (idx >= array_value.size()) {
            return nullptr;
        }
        return &array_value[idx];
    }
};

class JsonParser {
public:
    explicit JsonParser(const std::string& src) : src_(src) {}

    bool parse(Json& out, std::string* out_error) {
        if (out_error) {
            out_error->clear();
        }
        pos_ = 0;
        line_ = 1;
        col_ = 1;

        skip_ws();
        if (!parse_value(out)) {
            set_error(out_error);
            return false;
        }
        skip_ws();
        if (pos_ != src_.size()) {
            error_ = "Unexpected trailing characters";
            set_error(out_error);
            return false;
        }
        return true;
    }

private:
    const std::string& src_;
    std::size_t pos_ = 0;
    int line_ = 1;
    int col_ = 1;
    std::string error_;

    void set_error(std::string* out_error) const {
        if (!out_error) {
            return;
        }
        std::ostringstream oss;
        oss << "JSON parse error at line " << line_ << ", col " << col_ << ": " << error_;
        *out_error = oss.str();
    }

    char peek() const {
        if (pos_ >= src_.size()) {
            return '\0';
        }
        return src_[pos_];
    }

    char getch() {
        if (pos_ >= src_.size()) {
            return '\0';
        }
        const char c = src_[pos_++];
        if (c == '\n') {
            ++line_;
            col_ = 1;
        } else {
            ++col_;
        }
        return c;
    }

    void skip_ws() {
        while (pos_ < src_.size()) {
            const char c = peek();
            if (c == ' ' || c == '\t' || c == '\r' || c == '\n') {
                getch();
                continue;
            }
            break;
        }
    }

    bool consume(char expected) {
        if (peek() != expected) {
            return false;
        }
        getch();
        return true;
    }

    bool parse_value(Json& out) {
        const char c = peek();
        if (c == '{') {
            return parse_object(out);
        }
        if (c == '[') {
            return parse_array(out);
        }
        if (c == '"') {
            std::string s;
            if (!parse_string(s)) {
                return false;
            }
            out = Json::make_string(std::move(s));
            return true;
        }
        if (c == 't' || c == 'f') {
            return parse_bool(out);
        }
        if (c == 'n') {
            return parse_null(out);
        }
        if (c == '-' || std::isdigit(static_cast<unsigned char>(c))) {
            return parse_number(out);
        }
        error_ = "Unexpected character";
        return false;
    }

    bool parse_null(Json& out) {
        if (src_.compare(pos_, 4, "null") == 0) {
            pos_ += 4;
            col_ += 4;
            out = Json::make_null();
            return true;
        }
        error_ = "Invalid token (expected null)";
        return false;
    }

    bool parse_bool(Json& out) {
        if (src_.compare(pos_, 4, "true") == 0) {
            pos_ += 4;
            col_ += 4;
            out = Json::make_bool(true);
            return true;
        }
        if (src_.compare(pos_, 5, "false") == 0) {
            pos_ += 5;
            col_ += 5;
            out = Json::make_bool(false);
            return true;
        }
        error_ = "Invalid token (expected true/false)";
        return false;
    }

    bool parse_number(Json& out) {
        const std::size_t start = pos_;
        const char* begin = src_.c_str() + static_cast<std::ptrdiff_t>(start);
        char* end = nullptr;
        errno = 0;
        const double v = std::strtod(begin, &end);
        if (end == begin) {
            error_ = "Invalid number";
            return false;
        }
        if (errno != 0 || !std::isfinite(v)) {
            error_ = "Invalid number (overflow)";
            return false;
        }
        const std::size_t consumed = static_cast<std::size_t>(end - begin);
        pos_ += consumed;
        col_ += static_cast<int>(consumed);
        out = Json::make_number(v);
        return true;
    }

    static int hex_value(char c) {
        if (c >= '0' && c <= '9') return c - '0';
        if (c >= 'a' && c <= 'f') return 10 + (c - 'a');
        if (c >= 'A' && c <= 'F') return 10 + (c - 'A');
        return -1;
    }

    bool parse_string(std::string& out) {
        out.clear();
        if (!consume('"')) {
            error_ = "Expected string opening quote";
            return false;
        }
        while (pos_ < src_.size()) {
            const char c = getch();
            if (c == '"') {
                return true;
            }
            if (c == '\0') {
                break;
            }
            if (c == '\\') {
                const char e = getch();
                switch (e) {
                    case '"': out.push_back('"'); break;
                    case '\\': out.push_back('\\'); break;
                    case '/': out.push_back('/'); break;
                    case 'b': out.push_back('\b'); break;
                    case 'f': out.push_back('\f'); break;
                    case 'n': out.push_back('\n'); break;
                    case 'r': out.push_back('\r'); break;
                    case 't': out.push_back('\t'); break;
                    case 'u': {
                        int code = 0;
                        for (int i = 0; i < 4; ++i) {
                            const char h = getch();
                            const int v = hex_value(h);
                            if (v < 0) {
                                error_ = "Invalid \\\\u escape";
                                return false;
                            }
                            code = (code << 4) | v;
                        }
                        if (code <= 0x7f) {
                            out.push_back(static_cast<char>(code));
                        } else {
                            out.push_back('?');
                        }
                        break;
                    }
                    default:
                        error_ = "Invalid escape sequence";
                        return false;
                }
                continue;
            }
            out.push_back(c);
        }
        error_ = "Unterminated string";
        return false;
    }

    bool parse_array(Json& out) {
        if (!consume('[')) {
            error_ = "Expected '['";
            return false;
        }
        skip_ws();
        Json::Array arr;
        if (consume(']')) {
            out = Json::make_array(std::move(arr));
            return true;
        }
        while (true) {
            skip_ws();
            Json v;
            if (!parse_value(v)) {
                return false;
            }
            arr.push_back(std::move(v));
            skip_ws();
            if (consume(']')) {
                break;
            }
            if (!consume(',')) {
                error_ = "Expected ',' or ']'";
                return false;
            }
        }
        out = Json::make_array(std::move(arr));
        return true;
    }

    bool parse_object(Json& out) {
        if (!consume('{')) {
            error_ = "Expected '{'";
            return false;
        }
        skip_ws();
        Json::Object obj;
        if (consume('}')) {
            out = Json::make_object(std::move(obj));
            return true;
        }
        while (true) {
            skip_ws();
            if (peek() != '"') {
                error_ = "Expected object key string";
                return false;
            }
            std::string key;
            if (!parse_string(key)) {
                return false;
            }
            skip_ws();
            if (!consume(':')) {
                error_ = "Expected ':'";
                return false;
            }
            skip_ws();
            Json v;
            if (!parse_value(v)) {
                return false;
            }
            obj.emplace(std::move(key), std::move(v));
            skip_ws();
            if (consume('}')) {
                break;
            }
            if (!consume(',')) {
                error_ = "Expected ',' or '}'";
                return false;
            }
        }
        out = Json::make_object(std::move(obj));
        return true;
    }
};

static bool read_file_text(const std::string& path, std::string& out, std::string* out_error) {
    out.clear();
    std::ifstream f(path, std::ios::in | std::ios::binary);
    if (!f) {
        if (out_error) {
            *out_error = "Failed to open file: " + path;
        }
        return false;
    }
    std::ostringstream ss;
    ss << f.rdbuf();
    out = ss.str();
    return true;
}

static bool read_file_bytes(const std::string& path,
                            std::vector<std::uint8_t>& out,
                            std::string* out_error) {
    out.clear();
    std::ifstream f(path, std::ios::in | std::ios::binary);
    if (!f) {
        if (out_error) {
            *out_error = "Failed to open file: " + path;
        }
        return false;
    }
    f.seekg(0, std::ios::end);
    const std::streamoff size = f.tellg();
    if (size < 0) {
        if (out_error) {
            *out_error = "Failed to get file size: " + path;
        }
        return false;
    }
    f.seekg(0, std::ios::beg);
    out.resize(static_cast<std::size_t>(size));
    if (!out.empty()) {
        f.read(reinterpret_cast<char*>(out.data()), size);
        if (!f) {
            if (out_error) {
                *out_error = "Failed to read file: " + path;
            }
            out.clear();
            return false;
        }
    }
    return true;
}

static int json_int(const Json* j, int def) {
    if (!j || !j->is_number()) {
        return def;
    }
    const double v = j->number_value;
    if (!std::isfinite(v)) {
        return def;
    }
    if (v < static_cast<double>(std::numeric_limits<int>::min()) ||
        v > static_cast<double>(std::numeric_limits<int>::max())) {
        return def;
    }
    return static_cast<int>(v);
}

static float json_float(const Json* j, float def) {
    if (!j || !j->is_number()) {
        return def;
    }
    const double v = j->number_value;
    if (!std::isfinite(v)) {
        return def;
    }
    if (v < -static_cast<double>(std::numeric_limits<float>::max()) ||
        v > static_cast<double>(std::numeric_limits<float>::max())) {
        return def;
    }
    return static_cast<float>(v);
}

static std::string json_string(const Json* j, const std::string& def) {
    if (!j || !j->is_string()) {
        return def;
    }
    return j->string_value;
}

static bool json_bool(const Json* j, bool def) {
    if (!j || !j->is_bool()) {
        return def;
    }
    return j->bool_value;
}

static bool decode_vec(const Json* j, int n, std::vector<float>& out) {
    out.clear();
    if (!j || !j->is_array()) {
        return false;
    }
    if (static_cast<int>(j->array_value.size()) != n) {
        return false;
    }
    out.reserve(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        const Json& e = j->array_value[static_cast<std::size_t>(i)];
        if (!e.is_number()) {
            return false;
        }
        out.push_back(static_cast<float>(e.number_value));
    }
    return true;
}

static ImageTexture::WrapMode gltf_wrap_mode(int wrap) {
    // REPEAT=10497, CLAMP_TO_EDGE=33071, MIRRORED_REPEAT=33648
    switch (wrap) {
        case 33071:
            return ImageTexture::WrapMode::ClampToEdge;
        case 33648:
            return ImageTexture::WrapMode::MirroredRepeat;
        case 10497:
        default:
            return ImageTexture::WrapMode::Repeat;
    }
}

struct BufferView {
    int buffer = 0;
    std::size_t byte_offset = 0;
    std::size_t byte_length = 0;
    std::size_t byte_stride = 0;
};

struct Accessor {
    int buffer_view = -1;
    std::size_t byte_offset = 0;
    int component_type = 0;
    std::size_t count = 0;
    std::string type;
    bool normalized = false;
};

struct ImageDef {
    std::string uri;
};

struct SamplerDef {
    int wrap_s = 10497;
    int wrap_t = 10497;
};

struct TextureDef {
    int sampler = -1;
    int source = -1;
};

struct MaterialDef {
    MaterialPtr material;
};

struct PrimitiveDef {
    MeshDataPtr data;
    int material = -1;
    std::string name;
};

struct MeshDef {
    std::vector<PrimitiveDef> primitives;
};

struct NodeDef {
    std::string name;
    int mesh = -1;
    std::vector<int> children;
    Transform local = Transform::identity();
};

struct SceneDef {
    std::vector<int> nodes;
};

struct GltfDoc {
    std::filesystem::path base_dir;
    std::vector<std::vector<std::uint8_t>> buffers;
    std::vector<BufferView> buffer_views;
    std::vector<Accessor> accessors;
    std::vector<ImageDef> images;
    std::vector<SamplerDef> samplers;
    std::vector<TextureDef> textures;
    std::vector<MaterialDef> materials;
    std::vector<MeshDef> meshes;
    std::vector<NodeDef> nodes;
    std::vector<SceneDef> scenes;
    int default_scene = 0;
};

static int type_components(const std::string& type) {
    if (type == "SCALAR") return 1;
    if (type == "VEC2") return 2;
    if (type == "VEC3") return 3;
    if (type == "VEC4") return 4;
    if (type == "MAT2") return 4;
    if (type == "MAT3") return 9;
    if (type == "MAT4") return 16;
    return 0;
}

static std::size_t component_size(int component_type) {
    switch (component_type) {
        case 5120: return 1; // BYTE
        case 5121: return 1; // UNSIGNED_BYTE
        case 5122: return 2; // SHORT
        case 5123: return 2; // UNSIGNED_SHORT
        case 5125: return 4; // UNSIGNED_INT
        case 5126: return 4; // FLOAT
        default: return 0;
    }
}

static bool read_accessor_bytes(const GltfDoc& doc,
                                const Accessor& acc,
                                const std::uint8_t*& out_ptr,
                                std::size_t& out_stride,
                                std::size_t& out_elem_size,
                                std::string* out_error) {
    if (out_error) {
        out_error->clear();
    }

    if (acc.buffer_view < 0 || acc.buffer_view >= static_cast<int>(doc.buffer_views.size())) {
        if (out_error) {
            *out_error = "Invalid accessor bufferView index.";
        }
        return false;
    }
    const BufferView& view = doc.buffer_views[static_cast<std::size_t>(acc.buffer_view)];
    if (view.buffer < 0 || view.buffer >= static_cast<int>(doc.buffers.size())) {
        if (out_error) {
            *out_error = "Invalid bufferView buffer index.";
        }
        return false;
    }
    const std::vector<std::uint8_t>& buf = doc.buffers[static_cast<std::size_t>(view.buffer)];

    const int comps = type_components(acc.type);
    if (comps <= 0) {
        if (out_error) {
            *out_error = "Unsupported accessor type: " + acc.type;
        }
        return false;
    }
    const std::size_t csize = component_size(acc.component_type);
    if (csize == 0) {
        if (out_error) {
            *out_error = "Unsupported accessor componentType: " + std::to_string(acc.component_type);
        }
        return false;
    }
    const std::size_t elem_size = csize * static_cast<std::size_t>(comps);
    const std::size_t stride = (view.byte_stride != 0) ? view.byte_stride : elem_size;

    const std::size_t start = view.byte_offset + acc.byte_offset;
    if (start > buf.size()) {
        if (out_error) {
            *out_error = "Accessor starts past end of buffer.";
        }
        return false;
    }
    if (acc.count == 0) {
        out_ptr = nullptr;
        out_stride = stride;
        out_elem_size = elem_size;
        return true;
    }

    const std::size_t last = start + stride * (acc.count - 1) + elem_size;
    const std::size_t view_end = view.byte_offset + view.byte_length;

    if (last > buf.size()) {
        if (out_error) {
            *out_error = "Accessor reads past end of buffer.";
        }
        return false;
    }
    if (view.byte_length != 0 && last > view_end) {
        if (out_error) {
            *out_error = "Accessor reads past end of bufferView.";
        }
        return false;
    }

    out_ptr = buf.data() + start;
    out_stride = stride;
    out_elem_size = elem_size;
    return true;
}

static float decode_normalized(int component_type, std::int64_t v) {
    switch (component_type) {
        case 5120: { // BYTE
            const float denom = 127.0f;
            return clamp_float(static_cast<float>(v) / denom, -1.0f, 1.0f);
        }
        case 5121: { // UNSIGNED_BYTE
            const float denom = 255.0f;
            return clamp_float(static_cast<float>(v) / denom, 0.0f, 1.0f);
        }
        case 5122: { // SHORT
            const float denom = 32767.0f;
            return clamp_float(static_cast<float>(v) / denom, -1.0f, 1.0f);
        }
        case 5123: { // UNSIGNED_SHORT
            const float denom = 65535.0f;
            return clamp_float(static_cast<float>(v) / denom, 0.0f, 1.0f);
        }
        default:
            return 0.0f;
    }
}

static bool read_accessor_floats(const GltfDoc& doc,
                                 int accessor_index,
                                 int expected_components,
                                 std::vector<float>& out,
                                 std::string* out_error) {
    out.clear();
    if (out_error) {
        out_error->clear();
    }
    if (accessor_index < 0 || accessor_index >= static_cast<int>(doc.accessors.size())) {
        if (out_error) {
            *out_error = "Invalid accessor index.";
        }
        return false;
    }
    const Accessor& acc = doc.accessors[static_cast<std::size_t>(accessor_index)];
    const int comps = type_components(acc.type);
    if (comps != expected_components) {
        if (out_error) {
            *out_error = "Accessor has unexpected type: " + acc.type;
        }
        return false;
    }

    const std::uint8_t* ptr = nullptr;
    std::size_t stride = 0;
    std::size_t elem_size = 0;
    if (!read_accessor_bytes(doc, acc, ptr, stride, elem_size, out_error)) {
        return false;
    }

    out.resize(acc.count * static_cast<std::size_t>(comps));
    if (acc.count == 0) {
        return true;
    }

    const std::size_t csize = component_size(acc.component_type);
    if (csize == 0) {
        if (out_error) {
            *out_error = "Unsupported accessor componentType.";
        }
        return false;
    }

    for (std::size_t i = 0; i < acc.count; ++i) {
        const std::uint8_t* e = ptr + i * stride;
        for (int c = 0; c < comps; ++c) {
            float v = 0.0f;
            const std::uint8_t* src = e + static_cast<std::size_t>(c) * csize;
            if (acc.component_type == 5126) {
                float f = 0.0f;
                std::memcpy(&f, src, sizeof(float));
                v = f;
            } else {
                std::int64_t iv = 0;
                if (acc.component_type == 5120) {
                    std::int8_t t = 0;
                    std::memcpy(&t, src, 1);
                    iv = t;
                } else if (acc.component_type == 5121) {
                    std::uint8_t t = 0;
                    std::memcpy(&t, src, 1);
                    iv = t;
                } else if (acc.component_type == 5122) {
                    std::int16_t t = 0;
                    std::memcpy(&t, src, 2);
                    iv = t;
                } else if (acc.component_type == 5123) {
                    std::uint16_t t = 0;
                    std::memcpy(&t, src, 2);
                    iv = t;
                } else if (acc.component_type == 5125) {
                    std::uint32_t t = 0;
                    std::memcpy(&t, src, 4);
                    iv = t;
                } else {
                    if (out_error) {
                        *out_error = "Unsupported accessor componentType for float conversion.";
                    }
                    return false;
                }
                if (acc.normalized && acc.component_type != 5125) {
                    v = decode_normalized(acc.component_type, iv);
                } else {
                    v = static_cast<float>(iv);
                }
            }
            out[i * static_cast<std::size_t>(comps) + static_cast<std::size_t>(c)] = v;
        }
    }

    (void)elem_size;
    return true;
}

static bool read_accessor_indices(const GltfDoc& doc,
                                  int accessor_index,
                                  std::vector<std::uint32_t>& out,
                                  std::string* out_error) {
    out.clear();
    if (out_error) {
        out_error->clear();
    }
    if (accessor_index < 0 || accessor_index >= static_cast<int>(doc.accessors.size())) {
        if (out_error) {
            *out_error = "Invalid indices accessor index.";
        }
        return false;
    }

    const Accessor& acc = doc.accessors[static_cast<std::size_t>(accessor_index)];
    if (acc.type != "SCALAR") {
        if (out_error) {
            *out_error = "Indices accessor must be SCALAR.";
        }
        return false;
    }
    if (acc.component_type != 5121 && acc.component_type != 5123 && acc.component_type != 5125) {
        if (out_error) {
            *out_error = "Unsupported index componentType: " + std::to_string(acc.component_type);
        }
        return false;
    }

    const std::uint8_t* ptr = nullptr;
    std::size_t stride = 0;
    std::size_t elem_size = 0;
    if (!read_accessor_bytes(doc, acc, ptr, stride, elem_size, out_error)) {
        return false;
    }

    out.resize(acc.count);
    for (std::size_t i = 0; i < acc.count; ++i) {
        const std::uint8_t* src = ptr + i * stride;
        std::uint32_t v = 0;
        if (acc.component_type == 5121) {
            std::uint8_t t = 0;
            std::memcpy(&t, src, 1);
            v = t;
        } else if (acc.component_type == 5123) {
            std::uint16_t t = 0;
            std::memcpy(&t, src, 2);
            v = t;
        } else if (acc.component_type == 5125) {
            std::uint32_t t = 0;
            std::memcpy(&t, src, 4);
            v = t;
        }
        out[i] = v;
    }

    (void)elem_size;
    return true;
}

static Transform quat_transform(float x, float y, float z, float w) {
    const float xx = x * x;
    const float yy = y * y;
    const float zz = z * z;
    const float xy = x * y;
    const float xz = x * z;
    const float yz = y * z;
    const float wx = w * x;
    const float wy = w * y;
    const float wz = w * z;

    const float m00 = 1.0f - 2.0f * (yy + zz);
    const float m01 = 2.0f * (xy + wz);
    const float m02 = 2.0f * (xz - wy);

    const float m10 = 2.0f * (xy - wz);
    const float m11 = 1.0f - 2.0f * (xx + zz);
    const float m12 = 2.0f * (yz + wx);

    const float m20 = 2.0f * (xz + wy);
    const float m21 = 2.0f * (yz - wx);
    const float m22 = 1.0f - 2.0f * (xx + yy);

    Mat3 linear(
        Vec3(m00, m10, m20),
        Vec3(m01, m11, m21),
        Vec3(m02, m12, m22));
    return Transform(linear, Vec3::zero(), transpose(linear));
}

static Transform node_local_transform(const Json& node_obj) {
    const Json* matrix = node_obj.get("matrix");
    if (matrix && matrix->is_array() && matrix->array_value.size() == 16) {
        float m[16];
        for (std::size_t i = 0; i < 16; ++i) {
            const Json& e = matrix->array_value[i];
            m[i] = e.is_number() ? static_cast<float>(e.number_value) : 0.0f;
        }
        Mat3 linear(
            Vec3(m[0], m[1], m[2]),
            Vec3(m[4], m[5], m[6]),
            Vec3(m[8], m[9], m[10]));
        Vec3 t(m[12], m[13], m[14]);
        return Transform(linear, t);
    }

    Vec3 t(0.0f);
    Vec3 s(1.0f, 1.0f, 1.0f);
    float rx = 0.0f, ry = 0.0f, rz = 0.0f, rw = 1.0f;

    std::vector<float> v;
    if (decode_vec(node_obj.get("translation"), 3, v)) {
        t = Vec3(v[0], v[1], v[2]);
    }
    if (decode_vec(node_obj.get("scale"), 3, v)) {
        s = Vec3(v[0], v[1], v[2]);
    }
    if (decode_vec(node_obj.get("rotation"), 4, v)) {
        rx = v[0];
        ry = v[1];
        rz = v[2];
        rw = v[3];
    }

    const Transform T = Transform::translate(t);
    const Transform R = quat_transform(rx, ry, rz, rw);
    const Transform S = Transform::scale(s);
    return T * R * S;
}

static bool parse_gltf_doc(const std::string& path,
                           const GltfLoadOptions& /*options*/,
                           GltfDoc& out,
                           Json& out_root,
                           std::string* out_error) {
    if (out_error) {
        out_error->clear();
    }
    out = GltfDoc();
    out_root = Json();

    std::string text;
    if (!read_file_text(path, text, out_error)) {
        return false;
    }

    Json root;
    JsonParser parser(text);
    std::string parse_err;
    if (!parser.parse(root, &parse_err)) {
        if (out_error) {
            *out_error = parse_err;
        }
        return false;
    }
    if (!root.is_object()) {
        if (out_error) {
            *out_error = "glTF root must be a JSON object.";
        }
        return false;
    }

    out_root = root;
    out.base_dir = std::filesystem::path(path).parent_path();
    out.default_scene = json_int(root.get("scene"), 0);

    const Json* buffers = root.get("buffers");
    if (buffers && buffers->is_array()) {
        out.buffers.resize(buffers->array_value.size());
        for (std::size_t i = 0; i < buffers->array_value.size(); ++i) {
            const Json* b = buffers->at(i);
            if (!b || !b->is_object()) {
                continue;
            }
            const std::string uri = json_string(b->get("uri"), "");
            if (uri.empty()) {
                if (out_error) {
                    *out_error = "Only external .bin buffers are supported (missing uri).";
                }
                return false;
            }
            if (uri.rfind("data:", 0) == 0) {
                if (out_error) {
                    *out_error = "data: URIs are not supported for buffers.";
                }
                return false;
            }
            const std::filesystem::path buf_path = out.base_dir / std::filesystem::path(uri);
            std::vector<std::uint8_t> bytes;
            std::string io_err;
            if (!read_file_bytes(buf_path.string(), bytes, &io_err)) {
                if (out_error) {
                    *out_error = io_err;
                }
                return false;
            }
            out.buffers[i] = std::move(bytes);
        }
    }

    const Json* views = root.get("bufferViews");
    if (views && views->is_array()) {
        out.buffer_views.resize(views->array_value.size());
        for (std::size_t i = 0; i < views->array_value.size(); ++i) {
            const Json* v = views->at(i);
            if (!v || !v->is_object()) {
                continue;
            }
            BufferView bv;
            bv.buffer = json_int(v->get("buffer"), 0);
            bv.byte_offset = static_cast<std::size_t>(json_int(v->get("byteOffset"), 0));
            bv.byte_length = static_cast<std::size_t>(json_int(v->get("byteLength"), 0));
            bv.byte_stride = static_cast<std::size_t>(json_int(v->get("byteStride"), 0));
            out.buffer_views[i] = bv;
        }
    }

    const Json* accessors = root.get("accessors");
    if (accessors && accessors->is_array()) {
        out.accessors.resize(accessors->array_value.size());
        for (std::size_t i = 0; i < accessors->array_value.size(); ++i) {
            const Json* a = accessors->at(i);
            if (!a || !a->is_object()) {
                continue;
            }
            Accessor acc;
            acc.buffer_view = json_int(a->get("bufferView"), -1);
            acc.byte_offset = static_cast<std::size_t>(json_int(a->get("byteOffset"), 0));
            acc.component_type = json_int(a->get("componentType"), 0);
            acc.count = static_cast<std::size_t>(json_int(a->get("count"), 0));
            acc.type = json_string(a->get("type"), "");
            acc.normalized = json_bool(a->get("normalized"), false);
            out.accessors[i] = std::move(acc);
        }
    }

    const Json* images = root.get("images");
    if (images && images->is_array()) {
        out.images.resize(images->array_value.size());
        for (std::size_t i = 0; i < images->array_value.size(); ++i) {
            const Json* im = images->at(i);
            if (!im || !im->is_object()) {
                continue;
            }
            ImageDef idef;
            idef.uri = json_string(im->get("uri"), "");
            out.images[i] = std::move(idef);
        }
    }

    const Json* samplers = root.get("samplers");
    if (samplers && samplers->is_array()) {
        out.samplers.resize(samplers->array_value.size());
        for (std::size_t i = 0; i < samplers->array_value.size(); ++i) {
            const Json* s = samplers->at(i);
            if (!s || !s->is_object()) {
                continue;
            }
            SamplerDef sd;
            sd.wrap_s = json_int(s->get("wrapS"), 10497);
            sd.wrap_t = json_int(s->get("wrapT"), 10497);
            out.samplers[i] = sd;
        }
    }

    const Json* textures = root.get("textures");
    if (textures && textures->is_array()) {
        out.textures.resize(textures->array_value.size());
        for (std::size_t i = 0; i < textures->array_value.size(); ++i) {
            const Json* t = textures->at(i);
            if (!t || !t->is_object()) {
                continue;
            }
            TextureDef td;
            td.sampler = json_int(t->get("sampler"), -1);
            td.source = json_int(t->get("source"), -1);
            out.textures[i] = td;
        }
    }

    return true;
}

static bool build_gltf_materials(const Json& root,
                                 const GltfLoadOptions& options,
                                 GltfDoc& doc,
                                 std::string* out_error) {
    if (out_error) {
        out_error->clear();
    }

    std::unordered_map<std::string, TexturePtr> image_cache;
    std::unordered_map<std::string, NormalMapPtr> normal_cache;

    auto resolve_texture = [&](int texture_index,
                               std::filesystem::path& out_path,
                               ImageTexture::WrapMode& out_wrap_s,
                               ImageTexture::WrapMode& out_wrap_t) -> bool {
        out_path.clear();
        out_wrap_s = ImageTexture::WrapMode::Repeat;
        out_wrap_t = ImageTexture::WrapMode::Repeat;

        if (texture_index < 0 || texture_index >= static_cast<int>(doc.textures.size())) {
            return false;
        }
        const TextureDef& tex = doc.textures[static_cast<std::size_t>(texture_index)];
        if (tex.source < 0 || tex.source >= static_cast<int>(doc.images.size())) {
            return false;
        }

        const std::string uri = doc.images[static_cast<std::size_t>(tex.source)].uri;
        if (uri.empty() || uri.rfind("data:", 0) == 0) {
            return false;
        }

        out_path = doc.base_dir / std::filesystem::path(uri);
        if (tex.sampler >= 0 && tex.sampler < static_cast<int>(doc.samplers.size())) {
            const SamplerDef& s = doc.samplers[static_cast<std::size_t>(tex.sampler)];
            out_wrap_s = gltf_wrap_mode(s.wrap_s);
            out_wrap_t = gltf_wrap_mode(s.wrap_t);
        }
        return true;
    };

    auto load_image_texture = [&](int texture_index,
                                  ImageTexture::ColorSpace color_space,
                                  int channel) -> TexturePtr {
        std::filesystem::path tex_path;
        ImageTexture::WrapMode wrap_s = ImageTexture::WrapMode::Repeat;
        ImageTexture::WrapMode wrap_t = ImageTexture::WrapMode::Repeat;
        if (!resolve_texture(texture_index, tex_path, wrap_s, wrap_t)) {
            return nullptr;
        }

        std::ostringstream key;
        key << tex_path.string()
            << "|cs=" << static_cast<int>(color_space)
            << "|ch=" << channel
            << "|flip=" << (options.flip_v ? 1 : 0)
            << "|ws=" << static_cast<int>(wrap_s)
            << "|wt=" << static_cast<int>(wrap_t);

        const std::string k = key.str();
        auto it = image_cache.find(k);
        if (it != image_cache.end()) {
            return it->second;
        }

        auto tex = std::make_shared<ImageTexture>(
            tex_path.string(), color_space, channel, options.flip_v, wrap_s, wrap_t);
        image_cache.emplace(k, tex);
        return tex;
    };

    auto load_normal_map = [&](int texture_index) -> NormalMapPtr {
        std::filesystem::path tex_path;
        ImageTexture::WrapMode wrap_s = ImageTexture::WrapMode::Repeat;
        ImageTexture::WrapMode wrap_t = ImageTexture::WrapMode::Repeat;
        if (!resolve_texture(texture_index, tex_path, wrap_s, wrap_t)) {
            return nullptr;
        }

        std::ostringstream key;
        key << tex_path.string()
            << "|flip=" << (options.flip_v ? 1 : 0)
            << "|ws=" << static_cast<int>(wrap_s)
            << "|wt=" << static_cast<int>(wrap_t);

        const std::string k = key.str();
        auto it = normal_cache.find(k);
        if (it != normal_cache.end()) {
            return it->second;
        }

        auto nm = std::make_shared<NormalMapTexture>(tex_path.string(), options.flip_v, wrap_s, wrap_t);
        normal_cache.emplace(k, nm);
        return nm;
    };

    const Json* mats = root.get("materials");
    if (!mats || !mats->is_array()) {
        doc.materials.clear();
        return true;
    }

    doc.materials.resize(mats->array_value.size());
    for (std::size_t i = 0; i < mats->array_value.size(); ++i) {
        const Json* m = mats->at(i);
        if (!m || !m->is_object()) {
            auto fallback = std::make_shared<Lambertian>(
                std::make_shared<SolidColor>(Color(0.8f, 0.8f, 0.8f)));
            doc.materials[i].material = fallback;
            continue;
        }

        std::string alpha_mode = json_string(m->get("alphaMode"), "OPAQUE");
        const float alpha_cutoff = json_float(m->get("alphaCutoff"), 0.5f);

        float bc_r = 1.0f;
        float bc_g = 1.0f;
        float bc_b = 1.0f;
        float bc_a = 1.0f;
        int base_color_tex = -1;
        int metallic_roughness_tex = -1;
        float metallic_factor = 1.0f;
        float roughness_factor = 1.0f;

        const Json* pbr = m->get("pbrMetallicRoughness");
        if (pbr && pbr->is_object()) {
            std::vector<float> bc;
            if (decode_vec(pbr->get("baseColorFactor"), 4, bc)) {
                bc_r = bc[0];
                bc_g = bc[1];
                bc_b = bc[2];
                bc_a = bc[3];
            }

            const Json* bct = pbr->get("baseColorTexture");
            if (bct && bct->is_object()) {
                base_color_tex = json_int(bct->get("index"), -1);
            }

            metallic_factor = json_float(pbr->get("metallicFactor"), 1.0f);
            roughness_factor = json_float(pbr->get("roughnessFactor"), 1.0f);

            const Json* mrt = pbr->get("metallicRoughnessTexture");
            if (mrt && mrt->is_object()) {
                metallic_roughness_tex = json_int(mrt->get("index"), -1);
            }
        }

        TexturePtr base_tex;
        if (base_color_tex >= 0) {
            base_tex = load_image_texture(base_color_tex, ImageTexture::ColorSpace::sRGB, -1);
        }
        if (!base_tex) {
            base_tex = std::make_shared<SolidColor>(Color(1.0f, 1.0f, 1.0f));
        }

        TexturePtr base = std::make_shared<ColorFactorTexture>(
            base_tex, Color(bc_r, bc_g, bc_b), bc_a);

        if (alpha_mode == "MASK") {
            base = std::make_shared<AlphaCutoffTexture>(base, alpha_cutoff);
        } else if (alpha_mode == "OPAQUE") {
            base = std::make_shared<ForceAlphaTexture>(base, 1.0f);
        }

        TexturePtr metallic_tex;
        TexturePtr roughness_tex;
        if (metallic_roughness_tex >= 0) {
            roughness_tex = load_image_texture(
                metallic_roughness_tex, ImageTexture::ColorSpace::Linear, 1);
            if (metallic_factor > 0.0f) {
                metallic_tex = load_image_texture(
                    metallic_roughness_tex, ImageTexture::ColorSpace::Linear, 2);
            }
        }

        NormalMapPtr normal_map;
        float normal_strength = 1.0f;
        const Json* nt = m->get("normalTexture");
        if (nt && nt->is_object()) {
            const int idx = json_int(nt->get("index"), -1);
            normal_strength = json_float(nt->get("scale"), 1.0f);
            normal_map = load_normal_map(idx);
        }

        TexturePtr occlusion_tex;
        float occlusion_strength = 1.0f;
        const Json* ot = m->get("occlusionTexture");
        if (ot && ot->is_object()) {
            const int idx = json_int(ot->get("index"), -1);
            occlusion_strength = json_float(ot->get("strength"), 1.0f);
            occlusion_tex = load_image_texture(idx, ImageTexture::ColorSpace::Linear, 0);
        }

        doc.materials[i].material = std::make_shared<PrincipledBSDF>(
            base,
            metallic_factor,
            metallic_tex,
            roughness_factor,
            roughness_tex,
            normal_map,
            normal_strength,
            occlusion_tex,
            occlusion_strength);
    }

    return true;
}

static bool build_gltf_meshes(const Json& root, GltfDoc& doc, std::string* out_error) {
    if (out_error) {
        out_error->clear();
    }

    const Json* meshes = root.get("meshes");
    if (!meshes || !meshes->is_array()) {
        doc.meshes.clear();
        return true;
    }

    doc.meshes.resize(meshes->array_value.size());
    for (std::size_t mi = 0; mi < meshes->array_value.size(); ++mi) {
        const Json* mesh_obj = meshes->at(mi);
        if (!mesh_obj || !mesh_obj->is_object()) {
            continue;
        }

        const std::string mesh_name = json_string(mesh_obj->get("name"), "");
        const Json* prims = mesh_obj->get("primitives");
        if (!prims || !prims->is_array()) {
            continue;
        }

        MeshDef mdef;
        mdef.primitives.reserve(prims->array_value.size());

        for (std::size_t pi = 0; pi < prims->array_value.size(); ++pi) {
            const Json* prim = prims->at(pi);
            if (!prim || !prim->is_object()) {
                continue;
            }

            const int mode = json_int(prim->get("mode"), 4);
            if (mode != 4) {
                if (out_error) {
                    *out_error = "Only TRIANGLES primitives are supported (mode=4).";
                }
                return false;
            }

            const Json* attrs = prim->get("attributes");
            if (!attrs || !attrs->is_object()) {
                if (out_error) {
                    *out_error = "Primitive is missing attributes.";
                }
                return false;
            }

            const int pos_acc = json_int(attrs->get("POSITION"), -1);
            if (pos_acc < 0) {
                if (out_error) {
                    *out_error = "Primitive is missing POSITION attribute.";
                }
                return false;
            }

            const int nrm_acc = json_int(attrs->get("NORMAL"), -1);
            const int uv0_acc = json_int(attrs->get("TEXCOORD_0"), -1);
            const int idx_acc = json_int(prim->get("indices"), -1);
            const int mat_idx = json_int(prim->get("material"), -1);

            std::vector<float> pos_f;
            std::string err;
            if (!read_accessor_floats(doc, pos_acc, 3, pos_f, &err)) {
                if (out_error) {
                    *out_error = "POSITION: " + err;
                }
                return false;
            }
            std::vector<Vec3> positions;
            positions.reserve(pos_f.size() / 3u);
            for (std::size_t i = 0; i + 2 < pos_f.size(); i += 3) {
                positions.emplace_back(pos_f[i + 0], pos_f[i + 1], pos_f[i + 2]);
            }

            std::vector<Vec3> normals;
            if (nrm_acc >= 0) {
                std::vector<float> nrm_f;
                if (!read_accessor_floats(doc, nrm_acc, 3, nrm_f, &err)) {
                    if (out_error) {
                        *out_error = "NORMAL: " + err;
                    }
                    return false;
                }
                normals.reserve(nrm_f.size() / 3u);
                for (std::size_t i = 0; i + 2 < nrm_f.size(); i += 3) {
                    normals.emplace_back(nrm_f[i + 0], nrm_f[i + 1], nrm_f[i + 2]);
                }
            }

            std::vector<Vec2> uvs;
            if (uv0_acc >= 0) {
                std::vector<float> uv_f;
                if (!read_accessor_floats(doc, uv0_acc, 2, uv_f, &err)) {
                    if (out_error) {
                        *out_error = "TEXCOORD_0: " + err;
                    }
                    return false;
                }
                uvs.reserve(uv_f.size() / 2u);
                for (std::size_t i = 0; i + 1 < uv_f.size(); i += 2) {
                    uvs.emplace_back(uv_f[i + 0], uv_f[i + 1]);
                }
            }

            std::vector<std::uint32_t> indices;
            if (idx_acc >= 0) {
                if (!read_accessor_indices(doc, idx_acc, indices, &err)) {
                    if (out_error) {
                        *out_error = "indices: " + err;
                    }
                    return false;
                }
            } else {
                indices.reserve(positions.size());
                for (std::uint32_t v = 0; v < static_cast<std::uint32_t>(positions.size()); ++v) {
                    indices.push_back(v);
                }
            }

            if (indices.size() % 3u != 0u) {
                if (out_error) {
                    *out_error = "Primitive indices are not a multiple of 3 (expected triangles).";
                }
                return false;
            }

            const std::size_t vert_count = positions.size();
            for (std::size_t i = 0; i < indices.size(); ++i) {
                if (indices[i] >= vert_count) {
                    if (out_error) {
                        *out_error = "Primitive has out-of-range vertex indices.";
                    }
                    return false;
                }
            }

            PrimitiveDef pdef;
            pdef.material = mat_idx;
            if (!mesh_name.empty()) {
                if (prims->array_value.size() > 1) {
                    pdef.name = mesh_name + "_prim" + std::to_string(pi);
                } else {
                    pdef.name = mesh_name;
                }
            } else {
                pdef.name = "mesh" + std::to_string(mi) + "_prim" + std::to_string(pi);
            }

            pdef.data = std::make_shared<MeshData>(
                std::move(positions),
                std::move(normals),
                std::move(uvs),
                std::move(indices));

            mdef.primitives.push_back(std::move(pdef));
        }

        doc.meshes[mi] = std::move(mdef);
    }

    return true;
}

static bool build_gltf_nodes_and_scenes(const Json& root, GltfDoc& doc, std::string* out_error) {
    if (out_error) {
        out_error->clear();
    }

    const Json* nodes = root.get("nodes");
    if (nodes && nodes->is_array()) {
        doc.nodes.resize(nodes->array_value.size());
        for (std::size_t i = 0; i < nodes->array_value.size(); ++i) {
            const Json* n = nodes->at(i);
            if (!n || !n->is_object()) {
                continue;
            }
            NodeDef nd;
            nd.name = json_string(n->get("name"), "");
            nd.mesh = json_int(n->get("mesh"), -1);
            nd.local = node_local_transform(*n);

            const Json* children = n->get("children");
            if (children && children->is_array()) {
                nd.children.reserve(children->array_value.size());
                for (std::size_t ci = 0; ci < children->array_value.size(); ++ci) {
                    const Json& c = children->array_value[ci];
                    if (!c.is_number()) {
                        continue;
                    }
                    nd.children.push_back(static_cast<int>(c.number_value));
                }
            }
            doc.nodes[i] = std::move(nd);
        }
    } else {
        doc.nodes.clear();
    }

    const Json* scenes = root.get("scenes");
    if (scenes && scenes->is_array()) {
        doc.scenes.resize(scenes->array_value.size());
        for (std::size_t si = 0; si < scenes->array_value.size(); ++si) {
            const Json* s = scenes->at(si);
            if (!s || !s->is_object()) {
                continue;
            }
            SceneDef sd;
            const Json* snodes = s->get("nodes");
            if (snodes && snodes->is_array()) {
                sd.nodes.reserve(snodes->array_value.size());
                for (std::size_t ni = 0; ni < snodes->array_value.size(); ++ni) {
                    const Json& idx = snodes->array_value[ni];
                    if (!idx.is_number()) {
                        continue;
                    }
                    sd.nodes.push_back(static_cast<int>(idx.number_value));
                }
            }
            doc.scenes[si] = std::move(sd);
        }
    } else {
        doc.scenes.clear();
    }

    return true;
}

static void traverse_gltf_node(const GltfDoc& doc,
                               int node_index,
                               const Transform& parent,
                               int depth,
                               const MaterialPtr& default_material,
                               std::vector<GltfMeshInstance>& out_meshes) {
    if (depth > 256) {
        return;
    }
    if (node_index < 0 || node_index >= static_cast<int>(doc.nodes.size())) {
        return;
    }
    const NodeDef& node = doc.nodes[static_cast<std::size_t>(node_index)];
    const Transform world = parent * node.local;

    if (node.mesh >= 0 && node.mesh < static_cast<int>(doc.meshes.size())) {
        const MeshDef& mesh = doc.meshes[static_cast<std::size_t>(node.mesh)];
        for (const auto& prim : mesh.primitives) {
            if (!prim.data) {
                continue;
            }

            MaterialPtr mat = default_material;
            if (prim.material >= 0 && prim.material < static_cast<int>(doc.materials.size())) {
                const MaterialPtr& candidate = doc.materials[static_cast<std::size_t>(prim.material)].material;
                if (candidate) {
                    mat = candidate;
                }
            }

            std::string name;
            if (!node.name.empty()) {
                name = node.name;
            } else {
                name = prim.name;
            }

            out_meshes.push_back({prim.data, mat, world, std::move(name)});
        }
    }

    for (int child : node.children) {
        traverse_gltf_node(doc, child, world, depth + 1, default_material, out_meshes);
    }
}

static bool build_gltf_mesh_instances(const GltfDoc& doc,
                                      std::vector<GltfMeshInstance>& out_meshes,
                                      std::string* out_error) {
    out_meshes.clear();
    if (out_error) {
        out_error->clear();
    }

    const MaterialPtr default_mat = std::make_shared<Lambertian>(
        std::make_shared<SolidColor>(Color(0.8f, 0.8f, 0.8f)));

    std::vector<int> roots;
    int scene_index = doc.default_scene;
    if (!doc.scenes.empty()) {
        if (scene_index < 0 || scene_index >= static_cast<int>(doc.scenes.size())) {
            scene_index = 0;
        }
        roots = doc.scenes[static_cast<std::size_t>(scene_index)].nodes;
    } else {
        for (int i = 0; i < static_cast<int>(doc.nodes.size()); ++i) {
            roots.push_back(i);
        }
    }

    const Transform identity = Transform::identity();
    for (int root : roots) {
        traverse_gltf_node(doc, root, identity, 0, default_mat, out_meshes);
    }

    return true;
}

}  // namespace

bool load_gltf_meshes(const std::string& path,
                      std::vector<GltfMeshInstance>& out_meshes,
                      std::string* out_error,
                      const GltfLoadOptions& options) {
    out_meshes.clear();
    if (out_error) {
        out_error->clear();
    }

    GltfDoc doc;
    Json root;
    std::string err;
    if (!parse_gltf_doc(path, options, doc, root, &err)) {
        if (out_error) {
            *out_error = err;
        }
        return false;
    }
    if (!build_gltf_materials(root, options, doc, &err)) {
        if (out_error) {
            *out_error = err;
        }
        return false;
    }
    if (!build_gltf_meshes(root, doc, &err)) {
        if (out_error) {
            *out_error = err;
        }
        return false;
    }
    if (!build_gltf_nodes_and_scenes(root, doc, &err)) {
        if (out_error) {
            *out_error = err;
        }
        return false;
    }
    if (!build_gltf_mesh_instances(doc, out_meshes, &err)) {
        if (out_error) {
            *out_error = err;
        }
        return false;
    }

    return true;
}

bool append_gltf_to_scene(const std::string& path,
                          Scene& inout_scene,
                          const Transform& transform,
                          std::string* out_error,
                          const GltfLoadOptions& options) {
    std::vector<GltfMeshInstance> meshes;
    std::string err;
    if (!load_gltf_meshes(path, meshes, &err, options)) {
        if (out_error) {
            *out_error = err;
        }
        return false;
    }

    for (const auto& inst : meshes) {
        if (!inst.data || !inst.material) {
            continue;
        }
        inout_scene.objects.push_back(
            std::make_shared<Mesh>(inst.data, transform * inst.transform, inst.material));
    }

    return true;
}
