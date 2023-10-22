#include "SlimeSimImpl.h"

using namespace FSIM;

void SlimeSimImpl::updateImage(sf::Image& simImage)
{
	static int size = getSettings()->getSize();
	for (auto i = 0; i < size; ++i)
	{
		for (auto j = 0; j < size; ++j)		{
			const auto& pix = m_field[j+i*size];
			const auto pixcol = sf::Color{ uint8_t(pix*255), uint8_t(pix * 255), uint8_t(pix * 255)};
			for (auto is = 0; is < m_scale; ++is)
			{
				for (auto js = 0; js < m_scale; ++js)
				{
					simImage.setPixel(i * m_scale + is, j * m_scale + js, pixcol);
				}
			}
		}
	}
}
