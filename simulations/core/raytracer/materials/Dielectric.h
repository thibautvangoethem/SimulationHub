#pragma once
#include "Material.h"
#include "../utility/Vec3.h"

namespace RT
{
	class Dielectric : public Material
	{
	public:
		explicit Dielectric(double refraction) : m_refraction(refraction)
		{
		}

		bool scatter(const Ray& rIn, const HitRecord& record, Color& attenuation, Ray& scattered) const final;

	private:
		static inline double reflectance(double cosine, double refIdx);
		double m_refraction;
	};
}
