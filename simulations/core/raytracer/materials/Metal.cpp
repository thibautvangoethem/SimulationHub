#include "Metal.h"

#include "../hittableshapes/HittableShape.h"

using namespace RT;

bool Metal::scatter(const Ray& rIn, const HitRecord& record, Color& attenuation, Ray& scattered) const
{
	Vec3 reflected = reflect(unit_vector(rIn.direction()), record.normal);
	scattered = Ray{ record.p, reflected+m_fuzz* randomInUnitSPhere() };
	attenuation = m_albedo;
	return (scattered.direction().dot( record.normal) > 0);
}