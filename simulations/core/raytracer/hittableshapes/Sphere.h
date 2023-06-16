#pragma once

#include "HittableShape.h"

namespace RT {
    class Sphere : public HittableShape {
    public:
        Sphere() = default;
        Sphere(Point3 center, double radius, std::shared_ptr<Material> material): m_center(center), m_radius(radius),m_material(std::move(material)) {};

        virtual bool hit(
            const Ray& ray, double tMin, double tMax, HitRecord& record) const override; 

        Point3 m_center;
        double m_radius;
        std::shared_ptr<Material> m_material;
    };
}