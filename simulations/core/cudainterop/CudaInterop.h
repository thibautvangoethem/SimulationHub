#pragma once

#include "../baseinterface/Simulation.h"

#include <vector>
#include <memory>
#include <random>
#include <corecrt_math_defines.h>

namespace SIM {
	class SimulationSetting;
}
//Slime simulation; cpu only reimplementation of https://github.com/SebLague/Slime-Simulation 
namespace CSIM {

	class CudaInterop : public SIM::Simulation {

	public:
		CudaInterop(std::shared_ptr<SIM::SimulationSettings> settings);

		void advance(const double timestep) final;

		void handleClick(const bool isLeftClick, const int xpos, const int ypos) final;

	protected:
		
		const int m_size;
		//TODO put in settings
		const double speed = 100.0;
		const double evaporate = 0.99;
		const double steerStrength = 8;
		const int  sensorSize = 16;
		const double sensorDirOffset = 8.0;
	};
}