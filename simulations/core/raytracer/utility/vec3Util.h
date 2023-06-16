#pragma once

//#define USE_MATH_DEFINES

#include <random>
#include <corecrt_math_defines.h>
#include <limits>

#include <Eigen/Dense>

namespace RT
{
    //syntactic sugar for everyone
    using Point3 = Eigen::Vector3d;   
    using Color = Eigen::Vector3d;
    using Vec3 = Eigen::Vector3d;

    const double infinity = std::numeric_limits<double>::infinity();

    inline Vec3 unit_vector(Vec3 v) {
        return v / v.norm();
    }

    inline double degrees_to_radians(const double degrees) {
        return degrees * M_PI / 180.0;
    }

    inline Vec3 reflect(const Vec3& v,const Vec3& n )
    {
        return v - 2 * v.dot(n) * n;
    }

    inline Vec3 refract(const Vec3& uv, const Vec3& n,double frac)
    {
        const auto sub1 = fmin(-uv.dot(n), 1.0);
        const Vec3 perp = frac * (uv + sub1 * n);
        const Vec3 parallel = -sqrt(fabs(1.0 - perp.squaredNorm())) * n;
        return perp + parallel;
    }

    inline Vec3 randomInUnitDisk() {
        static std::uniform_real_distribution<double> distribution(-1.0, 1.0);
        static std::mt19937 generator;
        while (true) {
            auto p = Vec3(distribution(generator), distribution(generator), 0);
            
            if (p.squaredNorm() >= 1) continue;
            return p;
        }
    }

    inline static Vec3 getRandomVec()
    {
        static std::uniform_real_distribution<double> distribution(0.0, 1.0);
        static std::mt19937 generator;
        return Vec3(distribution(generator), distribution(generator), distribution(generator));
    }

    inline static Vec3 randomVector(double min, double max) {
        //TODO asses performance impact of creating the distributione aech time, if too slow replace with a 0-1 distribution and multiply, this would be faster but less accurate
        static std::uniform_real_distribution<double> distribution(0.0, 1.0);
        static std::mt19937 generator;
        return Vec3(min + (max - min) * distribution(generator), min + (max - min) * distribution(generator), min + (max - min) * distribution(generator));
    }

    inline static bool nearZero(const Vec3& vec)
        {
            const auto s = 1e-8;
            return (fabs(vec[0]) < s) && (fabs(vec[1]) < s) && (fabs(vec[2]) < s);
        }

    inline static Vec3 randomInUnitSPhere()
{
    while (true) {
        auto randomVec = randomVector(-1, 1);
        if (randomVec.squaredNorm() >= 1) continue;
        return randomVec;
    }
}
        
}