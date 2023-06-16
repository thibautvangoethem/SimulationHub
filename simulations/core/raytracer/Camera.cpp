#include "Camera.h"
using namespace RT;

Camera::Camera(Point3 from, Point3 at, Point3 viewUp, double verticalFov, double aspectRatio, double aperture, double focusDistance)
{
    //I am simply drawing square images here, this can be changed to fit the 16/9 aspect ratio of a default screen
    const auto theta = degrees_to_radians(verticalFov);
    const auto h = tan(theta / 2);

    const auto viewport_height = 2.0*h;
    const auto viewport_width = aspectRatio * viewport_height;

    m_w = unit_vector(from - at);
    m_u = unit_vector(viewUp.cross(m_w));
    m_v = m_w.cross(m_u);

    m_origin = from;
    m_horizontal = focusDistance * viewport_width * m_u;
    m_vertical = focusDistance* viewport_height * m_v;
    m_lowerLeftCorner= m_origin - m_horizontal / 2 - m_vertical / 2 - focusDistance*m_w;

    m_lensRadius = aperture / 2;
}

Ray Camera::getRay(const double u, const double v) const
{
    const Vec3 rd = m_lensRadius * randomInUnitDisk();
    const Vec3 offset = m_u * rd.x() + m_v * rd.y();
    return Ray(m_origin+offset, m_lowerLeftCorner+ u * m_horizontal + v * m_vertical - m_origin-offset);
}