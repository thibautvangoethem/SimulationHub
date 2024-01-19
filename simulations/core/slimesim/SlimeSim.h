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
namespace SSIM {
	struct SlimeAgent {
		float posx;
		float posy;
		float angle;

		SlimeAgent(float posx, float posy, float angle) :posx{posx}, posy{posy}, angle{angle} {};
	};

	class SlimeSim : public SIM::Simulation {

	public:
		SlimeSim(std::shared_ptr<SIM::SimulationSettings> settings);

		void advance(const double timestep) final;

		void handleClick(const bool isLeftClick, const int xpos, const int ypos) final;

	protected:
		void advanceAgents(const double timestep) noexcept;
		void blurDiffuse(const double timestep);
		void slimeSteerUpdate(SlimeAgent &agent, const double timestep);
		float sense(SlimeAgent& agent,float angle);

		const int m_size;
		std::vector < std::unique_ptr<SlimeAgent>> m_agents;
		std::vector<float> m_field;

		//TODO put in settings
		const double speed = 100.0;
		const double evaporate= 0.95;
		const double steerStrength = 2;
		const int  sensorSize = 16;
		const double sensorDirOffset = 8.0;
	};
}