#include "FiresimImpl.h"

using namespace FSIM;

void FireSimImpl::updateImage(sf::Image& simImage)
{
	static int size = getSettings()->getSize();
	for (auto i = 0; i < size; ++i)
	{
		for (auto j = 0; j < size; ++j)
		{
			const auto& pix = m_currentState[j][i];
			const auto pixcol = sf::Color{ pix.r, pix.g, pix.b };
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