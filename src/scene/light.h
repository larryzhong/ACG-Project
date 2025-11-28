#pragma once

#include <memory>
#include <vector>

#include "scene/hittable.h"

enum class LightType {
    Area,
    Point
};

struct Light {
    LightType type;
    HittablePtr shape;

    Light() = default;

    Light(LightType t, const HittablePtr& s)
        : type(t), shape(s) {}
};

class LightCollection {
public:
    void add_area_light(const HittablePtr& shape) {
        add_light(LightType::Area, shape);
    }

    void add_point_light(const HittablePtr& shape) {
        add_light(LightType::Point, shape);
    }

    const std::vector<Light>& lights() const {
        return lights_;
    }

    bool empty() const {
        return lights_.empty();
    }

    std::size_t size() const {
        return lights_.size();
    }

private:
    void add_light(LightType type, const HittablePtr& shape) {
        if (shape) {
            lights_.push_back(Light(type, shape));
        }
    }

    std::vector<Light> lights_;
};

