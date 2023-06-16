#include "Vec3.h"

#include <random>
using namespace RT;

double Vec3::length_squared() const
{
	return m_values[0] * m_values[0] + m_values[1] * m_values[1] + m_values[2] * m_values[2];
}

double Vec3::length() const
{
    return sqrt(length_squared());
}

Vec3 Vec3::randomInUnitSPhere()
{
    while (true) {
        auto randomVec = Vec3::random(-1, 1);
        if (randomVec.length_squared() >= 1) continue;
        return randomVec;
    }
}

bool Vec3::nearZero() const
{
    const auto s = 1e-8;
    return (fabs(m_values[0]) < s) && (fabs(m_values[1]) < s) && (fabs(m_values[2]) < s);
}