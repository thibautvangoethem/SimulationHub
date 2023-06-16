#include "Sphere.h"

using namespace RT;



bool Sphere::hit(const Ray& ray, double tMin, double tMax, HitRecord& record) const {
    const Vec3 oc = ray.origin() - m_center;
    const auto a = ray.direction().squaredNorm();
    const auto half_b = oc.dot(ray.direction());
    const auto c = oc.squaredNorm() - m_radius * m_radius;
    const auto discriminant = (half_b * half_b) - (a * c);


    //Believe it or not but adding the likely/unlickely here causes a 3% performance uplift
    if (discriminant < 0) [[likely]]
        return false;
    else [[unlikely]] {
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
        record.material = m_material;
        const Vec3 outward_normal = (record.p - m_center) / m_radius;
        record.setFaceNormal(ray, outward_normal);

        return true;
    }
}