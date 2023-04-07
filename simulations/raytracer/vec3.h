#pragma once

#include <cmath>
#include <iostream>
#include <array>

using std::sqrt;
//copied from https://raytracing.github.io/books/RayTracingInOneWeekend.html, modified towards my preference
namespace SIM {
    class vec3 {
    public:
        vec3() : m_values{ 0,0,0 } {}
        vec3(const double e0, const double e1,const double e2) : m_values{ e0, e1, e2 } {}

        double x() const { return m_values[0]; }
        double y() const { return m_values[1]; }
        double z() const { return m_values[2]; }

        vec3 operator-() const { return vec3(-m_values[0], -m_values[1], -m_values[2]); }
        double operator[](const int i) const { return m_values[i]; }
        double& operator[](const int i) { return m_values[i]; }

        vec3& operator+=(const vec3& v) {
            m_values[0] += v.m_values[0];
            m_values[1] += v.m_values[1];
            m_values[2] += v.m_values[2];   
            return *this;
        }

        vec3& operator*=(const double t) {
            m_values[0] *= t;
            m_values[1] *= t;
            m_values[2] *= t;
            return *this;
        }

        vec3& operator/=(const double t) {
            return *this *= 1 / t;
        }

        double length() const;
        double length_squared() const;
    private:
        std::array<double, 3> m_values;
    };

    // Type aliases for vec3
    using point3 = vec3;   // 3D point
    using color = vec3;    // RGB color
}