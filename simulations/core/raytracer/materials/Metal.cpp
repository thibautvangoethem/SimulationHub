#include "Metal.h"

#include "../hittableshapes/HittableShape.h"

using namespace RT;

bool Metal::scatter(const Ray& rIn, const HitRecord& record, Color& attenuation, Ray& scattered) const
{
	Vec3 reflected = reflect(unit_vector(rIn.direction()), record.normal);
	scattered = Ray{ record.p, reflected+m_fuzz* Vec3::randomInUnitSPhere() };
	attenuation = m_albedo;
	return (dot(scattered.direction(), record.normal) > 0);
}