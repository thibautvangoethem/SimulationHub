#pragma once

#include "HittableShape.h"

namespace RT {
    class Sphere : public HittableShape {
    public:
        Sphere() = default;
        Sphere(Point3 center, double radius): m_center(center), m_radius(radius) {};

        virtual bool hit(
            const Ray& ray, double tMin, double tMax, HitRecord& record) const override;

    public:
        Point3 m_center;
        double m_radius;
    };
}