#pragma once

#include "SimulationSettings.h"

#include<vector>
#include <memory>

 
namespace SIM {
	struct colour {
		uint8_t r;
		uint8_t g;
		uint8_t b;
	};

	class Simulation {
	public:
		Simulation(std::shared_ptr<SimulationSettings> settings) :m_settings{std::move(settings)} {}
		virtual void advance(const double timestep) = 0;
		virtual std::vector<std::vector<colour>>&  getCurrentState()=0;
		virtual void handleClick(const bool isLeftClick,const int xpos,const int ypos)=0;

		std::shared_ptr<SimulationSettings> getSettings() const { return m_settings; };
	protected:
		const std::shared_ptr<SimulationSettings> m_settings;
	};
}