#include "HittableList.h"

using namespace RT;

bool hittableList::hit(
    const Ray& ray, double tMin, double tMax, HitRecord& record) const
{
    HitRecord tempRec;
    bool hitAnything = false;
    auto closestSoFar = tMax;

    for (const auto& object : m_objects) {
        if (object->hit(ray, tMin, closestSoFar, tempRec)) {
            hitAnything = true;
            closestSoFar = tempRec.t;
            record = tempRec;
        }
    }

    return hitAnything;
}