#include "Dielectric.h"

#include "../hittableshapes/HittableShape.h"

using namespace RT;

bool Dielectric::scatter(const Ray& rIn, const HitRecord& record, Color& attenuation, Ray& scattered) const
{
	attenuation = Color(1.0, 1.0, 1.0);
	const double refactRatio = record.frontFace ? (1.0 / m_refraction) : m_refraction;

	const Vec3 unitDirection = unit_vector(rIn.direction());

	double ct = fmin(-unitDirection.dot(record.normal), 1.0);
	double st = sqrt(1.0 - ct * ct);

	const bool cannotRefract = refactRatio * st > 1.0;

	Vec3 direction;
	static std::uniform_real_distribution<double> distribution(0.0, 1.0);
	static std::mt19937 generator;

	if(cannotRefract || reflectance(ct,refactRatio)> distribution(generator))
	{
		direction = reflect(unitDirection, record.normal);
	}else
	{
		direction = refract(unitDirection, record.normal, refactRatio);
	}

	scattered = Ray{ record.p,direction };
	return true;
}

double Dielectric::reflectance(double cosine, double refIdx)
{
	auto r0 = (1 - refIdx) / (1 + refIdx);
	r0 = r0 * r0;
	return r0 + (1 - r0) * pow((1 - cosine), 5);
}
