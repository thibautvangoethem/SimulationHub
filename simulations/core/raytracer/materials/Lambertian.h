#pragma once
#include "Material.h"
#include "../utility/Vec3.h"

namespace RT
{
	class Lambertian : public Material
	{
	public:
		explicit Lambertian(const Color& col): m_albedo(col)
		{
		}

		bool scatter(const Ray& rIn, const HitRecord& record, Color& attenuation, Ray& scattered) const final;

		Color m_albedo;
	};
}
