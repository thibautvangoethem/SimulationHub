#pragma once

#include "SFMLSimInterface.h"
#include "../core/slimesim/SlimeSim.h"

namespace FSIM
{
	class SlimeSimImpl :public SFMLSimInterface, SSIM::SlimeSim
	{

	public:
		SlimeSimImpl(std::shared_ptr<SIM::SimulationSettings> settings, int scale) :SSIM::SlimeSim(settings), m_scale(scale) {};
		[[nodiscard]] std::shared_ptr<SIM::SimulationSettings> getInternalSettings() const final { return getSettings(); }
		void click(const bool isLeftClick, const int xpos, const int ypos) final { handleClick(isLeftClick, xpos, ypos); }
		void advanceSim(const double timestep) final { advance(timestep); }
		void updateImage(sf::Image& simImage) final;

	private:
		int m_scale = 1;
	};
}