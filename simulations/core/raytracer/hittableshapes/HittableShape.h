#pragma once
#include "../utility/Ray.h"

namespace RT
{
	struct HitRecord {
		Point3 p;
		Vec3 normal;
		double t;

		bool frontFace;

		inline void setFaceNormal(const Ray& ray, const Vec3& outwardNormal) {
			frontFace = dot(ray.direction(), outwardNormal) < 0;
			normal = frontFace ? outwardNormal : -outwardNormal;
		}
	};

	class HittableShape
	{
	public:
		virtual ~HittableShape() = default;
		virtual bool hit(const Ray& ray, double tMin, double tMax, HitRecord& record) const = 0;
	};
}