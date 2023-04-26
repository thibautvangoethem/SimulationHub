#include "Vec3.h"
using namespace RT;


double Vec3::length_squared() const
{
	return m_values[0] * m_values[0] + m_values[1] * m_values[1] + m_values[2] * m_values[2];
}

double Vec3::length() const
{
    return sqrt(length_squared());
}