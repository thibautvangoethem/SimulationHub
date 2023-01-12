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
		Simulation(std::unique_ptr<SimulationSettings> settings) :m_settings{std::move(settings)} {}
		virtual void advance(double timestep) = 0;
		virtual std::vector<std::vector<colour>>&  getCurrentState()=0;
		virtual void handleClick(bool isLeftClick,int xpos,int ypos)=0;
	protected:
		const std::unique_ptr<SimulationSettings> m_settings;
	};
}