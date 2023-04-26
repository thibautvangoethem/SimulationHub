#include "RayTracerImpl.h"

using namespace FSIM;

void RayTracerImpl::updateImage(sf::Image& simImage)
{
	
	for (auto i = 0; i < m_size; ++i)
	{
		for (auto j = 0; j < m_size; ++j)
		{
			const auto& val = m_currentState[j][i];

			auto col= sf::Color{ val.r,val.g,val.b };
			for (auto is = 0; is < m_scale; ++is)
			{
				for (auto js = 0; js < m_scale; ++js)
				{
					simImage.setPixel(i * m_scale + is, j * m_scale + js, col);
				}
			}
		}
	}
}