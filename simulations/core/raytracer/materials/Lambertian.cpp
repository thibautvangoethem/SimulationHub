#include "Lambertian.h"

#include "../hittableshapes/HittableShape.h"

using namespace RT;

bool Lambertian::scatter(const Ray& rIn, const HitRecord& record, Color& attenuation, Ray& scattered) const
{
    auto scatterDirection = record.normal + Vec3::randomInUnitSPhere();
    if (scatterDirection.nearZero())
        scatterDirection = record.normal;

    scattered = Ray(record.p, scatterDirection);
    attenuation = m_albedo;
    return true;
}