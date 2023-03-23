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
		void addDensity(const int x, const int y, const float amount);
		void addVelocity(const int x, const int y, const float amount_x, const float amount_y);
		void applySources();

		void setBoundary(const int b, std::vector<float>& array);
		void linSolve(const int b, std::vector<float>& array, const std::vector<float>& prevArray, const float a, const float c);
		void diffuse(const int b, std::vector<float>& array, const std::vector<float>& prevArray, const float diff, const float dt);
		void project(std::vector<float>& velocityX, std::vector<float>& velocityY, std::vector<float>& clearVector, std::vector<float>& targetVectory);
		void advect(const int b, std::vector<float>& d, const std::vector<float>& prevD, const std::vector<float>& velocityX, const std::vector<float>& velocityY, const float dt);


		//Functions that use avx2 for acceleration
#ifdef AVX
		void linSolveAvx(const int b, std::vector<float>& array, const std::vector<float>& prevArray, const float a, const float c);
		void advectAvx(const int b, std::vector<float>& d, const std::vector<float>& prevD, const std::vector<float>& velocityX, const std::vector<float>& velocityY, const float dt);
		void projectAvx(std::vector<float>& velocityX, std::vector<float>& velocityY, std::vector<float>& clearVector, std::vector<float>& targetVectory);
#endif


		

		int m_size;

		std::vector<Flow> m_flows = std::vector<Flow>();
		std::vector<std::pair<int, int>> m_paintSources = std::vector<std::pair<int, int>>();

		std::vector<float> m_s;
		std::vector<float> m_density;

		std::vector<float> m_velX;
		std::vector<float> m_velY;

		std::vector<float> m_oldVelX;
		std::vector<float> m_oldVelY;

		std::vector<std::vector<SIM::colour>> m_currentState;

		std::unique_ptr<std::pair<int, int>> m_storedClick=nullptr;
	};
}