#pragma once
#include "../baseinterface/Simulation.h"
#include <utility> 
namespace SIM {
	class SimulationSetting;
}

namespace SIM {
	struct Flow {
		std::pair<int, int> location;
		std::pair<float, float> velocity;
	};

	class FluidSim :public Simulation {
		static constexpr SIM::colour EMPTY = SIM::colour{ 0,0,0 };
	public:
		FluidSim(std::unique_ptr<SIM::SimulationSettings> settings);
		void advance(double timestep) final;
		std::vector<std::vector<SIM::colour>>& getCurrentState() final;
		void handleClick(bool isLeftClick, int xpos, int ypos) final;
	private:
		void addDensity(int x, int y, float amount);
		void addVelocity(int x, int y, float amount_x, float amount_y);
		void applySources();

		void setBoundary(int b, std::vector<std::vector<float>>& array);
		void linSolve(int b, std::vector<std::vector<float>>& array, std::vector<std::vector<float>>& prevArray, float a, float c);
		void diffuse(int b, std::vector<std::vector<float>>& array, std::vector<std::vector<float>>& prevArray, float diff, float dt);
		void project(std::vector<std::vector<float>>& velocityX, std::vector<std::vector<float>>& velocityY, std::vector<std::vector<float>>& clearVector, std::vector<std::vector<float>>& targetVectory);
		void advect(int b, std::vector<std::vector<float>>& d, std::vector<std::vector<float>>& prevD, std::vector<std::vector<float>>& velocityX, std::vector<std::vector<float>>& velocitY, float dt);


		typedef std::vector<std::vector<float>> grid;

		int m_size;

		std::vector<Flow> m_flows = std::vector<Flow>();
		std::vector<std::pair<int, int>> m_paintSources = std::vector<std::pair<int, int>>();

		grid m_s;
		grid m_density;

		grid m_velX;
		grid m_velY;

		grid m_oldVelX;
		grid m_oldVelY;

		std::vector<std::vector<SIM::colour>> m_currentState;

		std::unique_ptr<std::pair<int, int>> m_storedClick=nullptr;
	};
}