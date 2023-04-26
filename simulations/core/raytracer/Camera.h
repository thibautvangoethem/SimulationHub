#pragma once

#include "utility/vec3Util.h"
#include "utility/Ray.h"

namespace RT {
    class Camera
    {
    public:
        Camera();

        Ray getRay(const double u,const double v) const;

    private:
        Point3 m_origin;
        Point3 m_lowerLeftCorner;
        Vec3 m_horizontal;
        Vec3 m_vertical;
    };
}