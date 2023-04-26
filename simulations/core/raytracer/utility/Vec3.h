#pragma once

#include <cmath>
#include <iostream>
#include <array>

using std::sqrt;
//copied from https://raytracing.github.io/books/RayTracingInOneWeekend.html, modified towards my preference
namespace RT {
    class Vec3 {
    public:
        Vec3() : m_values{ 0,0,0 } {}
        Vec3(const double e0, const double e1,const double e2) : m_values{ e0, e1, e2 } {}

        double x() const { return m_values[0]; }
        double y() const { return m_values[1]; }
        double z() const { return m_values[2]; }

        Vec3 operator-() const { return Vec3(-m_values[0], -m_values[1], -m_values[2]); }
        double operator[](const int i) const { return m_values[i]; }
        double& operator[](const int i) { return m_values[i]; }

        Vec3& operator+=(const Vec3& v) {
            m_values[0] += v.m_values[0];
            m_values[1] += v.m_values[1];
            m_values[2] += v.m_values[2];   
            return *this;
        }

        Vec3& operator*=(const double t) {
            m_values[0] *= t;
            m_values[1] *= t;
            m_values[2] *= t;
            return *this;
        }

        Vec3& operator/=(const double t) {
            return *this *= 1 / t;
        }

        double length() const;
        double length_squared() const;
    
        std::array<double, 3> m_values;
    };

    // Type aliases for vec3
    using Point3 = Vec3;   // 3D point
    using Color = Vec3;    // RGB color
}