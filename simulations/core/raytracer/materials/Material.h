#pragma once
#include "../utility/Ray.h"
#include "../utility/Vec3.h"

namespace RT
{
	struct HitRecord;
}

namespace RT
{

	class Material
	{
	public:
		virtual ~Material() = default;
		virtual bool scatter(const Ray& rIn, const HitRecord& record, Color& attenuation, Ray& scattered) const = 0;
	};
}