#include "vec3.h"
using namespace SIM;


double vec3::length_squared() const
{
	return m_values[0] * m_values[0] + m_values[1] * m_values[1] + m_values[2] * m_values[2];
}

double vec3::length() const
{
    return sqrt(length_squared());
}