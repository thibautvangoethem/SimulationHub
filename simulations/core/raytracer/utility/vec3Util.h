#pragma once

#define _USE_MATH_DEFINES
#include "Vec3.h"

#include <corecrt_math_defines.h>
#include <limits>

namespace RT
{
    const double infinity = std::numeric_limits<double>::infinity();

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

    inline Vec3 reflect(const Vec3& v,const Vec3& n )
    {
        return v - 2 * dot(v, n) * n;
    }

    inline Vec3 refract(const Vec3& uv, const Vec3& n,double frac)
    {
        const auto sub1 = fmin(dot(-uv, n), 1.0);
        const Vec3 perp = frac * (uv + sub1 * n);
        const Vec3 parallel = -sqrt(fabs(1.0 - perp.length_squared())) * n;
        return perp + parallel;
    }

    inline Vec3 randomInUnitDisk() {
        static std::uniform_real_distribution<double> distribution(-1.0, 1.0);
        static std::mt19937 generator;
        while (true) {
            auto p = Vec3(distribution(generator), distribution(generator), 0);
            if (p.length_squared() >= 1) continue;
            return p;
        }
    }
}