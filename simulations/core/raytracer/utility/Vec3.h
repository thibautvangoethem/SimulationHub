#pragma once

#include <cmath>
#include <array>
#include <random>

using std::sqrt;
//copied from https://raytracing.github.io/books/RayTracingInOneWeekend.html, modified towards my preference
namespace RT {
    class Vec3 {
    public:
        Vec3() : m_values{ 0,0,0 } {}
        Vec3(const double e0, const double e1,const double e2) : m_values{ e0, e1, e2 } {}

        [[nodiscard]] double x() const { return m_values[0]; }
        [[nodiscard]] double y() const { return m_values[1]; }
        [[nodiscard]] double z() const { return m_values[2]; }

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

        [[nodiscard]] bool nearZero() const;

        inline static Vec3 getRandomVec()
        {
            static std::uniform_real_distribution<double> distribution(0.0, 1.0);
            static std::mt19937 generator;
            return Vec3(distribution(generator), distribution(generator), distribution(generator));
        }

        inline static Vec3 random(double min, double max) {
            //TODO asses performance impact of creating the distributione aech time, if too slow replace with a 0-1 distribution and multiply, this would be faster but less accurate
            static std::uniform_real_distribution<double> distribution(0.0, 1.0);
            static std::mt19937 generator;
            return Vec3(min + (max - min) * distribution(generator), min + (max - min) * distribution(generator), min + (max - min) * distribution(generator));
        }


        static Vec3 randomInUnitSPhere();

        double length() const;
        double length_squared() const;

		
		//not private for ease of use and code readabilty
        std::array<double, 3> m_values;
    };

    // Type aliases for vec3
    using Point3 = Vec3;   // 3D point
    using Color = Vec3;    // RGB color
}