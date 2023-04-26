#include "Camera.h"
using namespace RT;

Camera::Camera()
{
    //I am simply drawing square images here, this can be changed to fit the 16/9 aspect ratio of a default screen
    auto aspect_ratio = 1.0;
    auto viewport_height = 2.0;
    auto viewport_width = aspect_ratio * viewport_height;;
    auto focal_length = 1.0;

    m_origin = Point3(0, 0, 0);
    m_horizontal = Vec3(viewport_width, 0.0, 0.0);
    m_vertical = Vec3(0.0, viewport_height, 0.0);
    m_lowerLeftCorner= m_origin - m_horizontal / 2 - m_vertical / 2 - Vec3(0, 0, focal_length);
}

Ray Camera::getRay(const double u, const double v) const
{
    return Ray(m_origin, m_lowerLeftCorner+ u * m_horizontal + v * m_vertical - m_origin);
}