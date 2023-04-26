#pragma once
#include <SFML/Graphics.hpp>

#pragma once

#include "../core/baseinterface/SimulationSettings.h"

#include <vector>

namespace FSIM {
	class SFMLSimInterface
	{

	public:
		virtual void updateImage(sf::Image& simImage) = 0;

		[[nodiscard]] virtual std::shared_ptr<SIM::SimulationSettings> getInternalSettings() const = 0;
		virtual void click(const bool isLeftClick, const int xpos, const int ypos) = 0;
		virtual void advanceSim(const double timestep) = 0;

	};
}