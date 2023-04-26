#define _USE_MATH_DEFINES

#pragma once
#include "Vec3.h"

#include <corecrt_math_defines.h>
#include <limits>

namespace RT
{
    const double infinity = std::numeric_limits<double>::infinity();

    inline std::ostream& operator<<(std::ostream& out, const Vec3& v) {
        return out << v.m_values[0] << ' ' << v.m_values[1] << ' ' << v.m_values[2];
    }

    inline Vec3 operator+(const Vec3& u, const Vec3& v) {
        return Vec3(u.m_values[0] + v.m_values[0], u.m_values[1] + v.m_values[1], u.m_values[2] + v.m_values[2]);
    }

    inline Vec3 operator-(const Vec3& u, const Vec3& v) {
        return Vec3(u.m_values[0] - v.m_values[0], u.m_values[1] - v.m_values[1], u.m_values[2] - v.m_values[2]);
    }

    inline Vec3 operator*(const Vec3& u, const Vec3& v) {
        return Vec3(u.m_values[0] * v.m_values[0], u.m_values[1] * v.m_values[1], u.m_values[2] * v.m_values[2]);
    }

    inline Vec3 operator*(double t, const Vec3& v) {
        return Vec3(t * v.m_values[0], t * v.m_values[1], t * v.m_values[2]);
    }

    inline Vec3 operator*(const Vec3& v, double t) {
        return t * v;
    }

    inline Vec3 operator/(Vec3 v, double t) {
        return (1 / t) * v;
    }

    inline double dot(const Vec3& u, const Vec3& v) {
        return u.m_values[0] * v.m_values[0]
            + u.m_values[1] * v.m_values[1]
            + u.m_values[2] * v.m_values[2];
    }

    inline Vec3 cross(const Vec3& u, const Vec3& v) {
        return Vec3(u.m_values[1] * v.m_values[2] - u.m_values[2] * v.m_values[1],
            u.m_values[2] * v.m_values[0] - u.m_values[0] * v.m_values[2],
            u.m_values[0] * v.m_values[1] - u.m_values[1] * v.m_values[0]);
    }

    inline Vec3 unit_vector(Vec3 v) {
        return v / v.length();
    }

    inline double degrees_to_radians(const double degrees) {
        return degrees * M_PI / 180.0;
    }
}