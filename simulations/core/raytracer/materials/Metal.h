#pragma once
#include "Material.h"
#include "../utility/vec3Util.h"

namespace RT
{
	class Metal : public Material
	{
	public:
		explicit Metal(const Color& col,double fuzz) : m_albedo(col),m_fuzz(fuzz)
		{
		}

		bool scatter(const Ray& rIn, const HitRecord& record, Color& attenuation, Ray& scattered) const final;

		Color m_albedo;
		double m_fuzz;
	};
}
