#include "Lambertian.h"

#include "../hittableshapes/HittableShape.h"

using namespace RT;

bool Lambertian::scatter(const Ray& rIn, const HitRecord& record, Color& attenuation, Ray& scattered) const
{
    //Auto doesnt work here due to templates in eigen
    Vec3 scatterDirection = record.normal + randomInUnitSPhere();
    if (nearZero(scatterDirection))
        scatterDirection = record.normal;

    scattered = Ray(record.p, scatterDirection);
    attenuation = m_albedo;
    return true;
}