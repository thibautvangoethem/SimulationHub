#include "RayTracerImpl.h"

using namespace FSIM;

void RayTracerImpl::updateImage(sf::Image& simImage)
{
	
	//for (auto i = 0; i < m_size; ++i)
	//{
	//	for (auto j = 0; j < m_size; ++j)
	//	{
	//		//The idea here is that higher density means that the cell has more white, and higher speed means that the cell has more red
	//		int val = m_density[i * m_size + j] * 32;
	//		if (val > 255)val = 255;
	//		if (val < 0)val = 0;
	//		auto valCasted = static_cast<uint8_t>(val);
	//		//TODO asses the performance impact of this calculation
	//		auto combinedSpeedVal = 1 - sqrt(m_velY[i * m_size + j] * m_velY[i * m_size + j] + m_velX[i * m_size + j] * m_velX[i * m_size + j]);

	//		valCasted = std::clamp(valCasted, static_cast<uint8_t>(0), static_cast<uint8_t>(255));
	//		auto red = valCasted * combinedSpeedVal;
	//		red = std::clamp(red, static_cast<float>(0), static_cast <float>(255));


	//		auto col= sf::Color{ static_cast<uint8_t>(red),valCasted,valCasted };
	//		for (int is = 0; is < m_scale; ++is)
	//		{
	//			for (int js = 0; js < m_scale; ++js)
	//			{
	//				simImage.setPixel(i * m_scale + is, j * m_scale + js, col);
	//			}
	//		}
	//	}
	//}

	//if (m_storedClick)
	//{
	//	auto i = m_storedClick->first;
	//	auto j = m_storedClick->second;
	//	for (int is = 0; is < m_scale; ++is)
	//	{
	//		for (int js = 0; js < m_scale; ++js)
	//		{
	//			simImage.setPixel(i * m_scale + is, j * m_scale + js, sf::Color{0,255,0});
	//		}
	//	}
	//}
}