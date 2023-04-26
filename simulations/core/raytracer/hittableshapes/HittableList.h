#pragma once
#include "HittableShape.h"
#include <vector>
#include <memory>

namespace RT
{
	class hittableList : public HittableShape
	{
    public:
        hittableList() = default;
        hittableList(std::shared_ptr<HittableShape> object) { add(object); }

        void clear() { m_objects.clear(); }
        void add(std::shared_ptr<HittableShape> object) { m_objects.push_back(object); }

        virtual bool hit(
            const Ray& ray, double tMin, double tMax, HitRecord& record) const override;

        std::vector<std::shared_ptr<HittableShape>> m_objects;
	};
}