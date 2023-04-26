#include "Sphere.h"

using namespace RT;



bool Sphere::hit(const Ray& ray, double tMin, double tMax, HitRecord& record) const {
    const Vec3 oc = ray.origin() - m_center;
    const auto a = ray.direction().length_squared();
    const auto half_b = dot(oc, ray.direction());
    const auto c = oc.length_squared() - m_radius * m_radius;

    const auto discriminant = half_b * half_b - a * c;
    if (discriminant < 0) return false;
    const auto sqrtd = sqrt(discriminant);

    // Find the nearest root that lies in the acceptable range.
    auto root = (-half_b - sqrtd) / a;
    if (root < tMin || tMax < root) {
        root = (-half_b + sqrtd) / a;
        if (root < tMin || tMax < root)
            return false;
    }

    record.t = root;
    record.p = ray.at(record.t);
    const Vec3 outward_normal = (record.p - m_center) / m_radius;
    record.setFaceNormal(ray, outward_normal);

    return true;
}