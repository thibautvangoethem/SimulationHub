#pragma once
#include "../baseinterface/Simulation.h"
#include "vec3.h"
#include <utility> 
namespace SIM {
	class SimulationSetting;
}


namespace SIM {
	//This raytracer was made following the tutorial on https://raytracing.github.io/books/RayTracingInOneWeekend.html
	class RayTracer :public Simulation {
	public:
		RayTracer(std::shared_ptr<SIM::SimulationSettings> settings);
		void advance(const double timestep) final;
		std::vector<std::vector<SIM::colour>>& getCurrentState() final;
		void handleClick(const bool isLeftClick, const int xpos, const int ypos) final;
	private:
		std::vector<std::vector<SIM::colour>> m_currentState;
		int m_size;
	};
}