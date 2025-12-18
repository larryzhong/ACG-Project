#pragma once

#include <cstdint>
#include <memory>

#include "core/ray.h"
#include "core/vec3.h"
#include "math/aabb.h"

class Material;
class Hittable;

struct HitRecord {
    Vec3 point;
    Vec3 normal;
    Vec3 tangent;
    float tangent_sign = 1.0f;
    float t;
    float u;
    float v;
    const Material* material;
    const Hittable* object = nullptr;
    bool front_face;
    Vec3 dpdu = Vec3(0.0f);
    Vec3 dpdv = Vec3(0.0f);
    float ray_footprint = 0.0f; // world-space radius at hit, for texture LOD
    std::uint32_t primitive_id = 0xffffffffu;

    void set_face_normal(const Ray& r, const Vec3& outward_normal) {
        front_face = dot(r.direction, outward_normal) < 0.0f;
        normal = front_face ? outward_normal : -outward_normal;
    }
};

class Hittable {
public:
    virtual ~Hittable() = default;

    virtual bool hit(const Ray& r, float t_min, float t_max, HitRecord& rec) const = 0;

    virtual AABB bounding_box() const = 0;
};

using HittablePtr = std::shared_ptr<Hittable>;
