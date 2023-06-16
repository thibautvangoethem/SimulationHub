#pragma once

#include "utility/vec3Util.h"
#include "utility/Ray.h"

namespace RT {
    class Camera
    {
    public:
        Camera(Point3 from, Point3 at, Point3 viewUp,double verticalFov,double aspectRatio,double aperture, double focusDistance);

        Ray getRay(const double u,const double v) const;

    private:
        Point3 m_origin;
        Point3 m_lowerLeftCorner;
        Vec3 m_horizontal;
        Vec3 m_vertical;
        Vec3 m_u;
        Vec3 m_v;
        Vec3 m_w;
        double m_lensRadius;
    };
}