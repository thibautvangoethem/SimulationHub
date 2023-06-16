#pragma once
#include "../baseinterface/Simulation.h"
#include <utility>
#include <vector>
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
		/*#### Sim::Simulation functions ####*/

		/*
		 * This constructor will as example also generate a flow at the middle of the simulation fuild
		 */
		FluidSim(std::shared_ptr<SIM::SimulationSettings> settings);

		/*
		 *Advances the simulation by the given time. Very high timesteps (>0.5) cause a high inaccuracy in the simulation
		 */
		void advance(const double timestep) final;
		
		/*
		 * Handles a click for this simulation.
		 * left click -> add a new flow(velocity source) between 2 clicks
		 * right click -> add new density source
		 */
		void handleClick(const bool isLeftClick, const int xpos, const int ypos) final;
	protected:
		void addDensity(const int x, const int y, const float amount);
		void addVelocity(const int x, const int y, const float amountX, const float amountY);
		/*
		 * Applies the current density/velocity sources to the simulation
		 */
		void applySources();

		/*
		 * Handles the boundary conditions, so fluid can not exit the simulation, and also gets rebound a little bit.
		 */
		void setBoundary(const int b, std::vector<float>& array);
		void linSolve(const int b, std::vector<float>& array, const std::vector<float>& prevArray, const float a, const float c);
		void diffuse(const int b, std::vector<float>& array, const std::vector<float>& prevArray, const float diff, const float dt);
		void project(std::vector<float>& velocityX, std::vector<float>& velocityY, std::vector<float>& clearVector, std::vector<float>& targetVectory);
		void advect(const int b, std::vector<float>& d, const std::vector<float>& prevD, const std::vector<float>& velocityX, const std::vector<float>& velocityY, const float dt);

		//stored here so it doesnt need to come from the settings each time
		const int m_size;

		std::vector<Flow> m_flows = std::vector<Flow>();
		std::vector<std::pair<int, int>> m_paintSources = std::vector<std::pair<int, int>>();

		std::vector<float> m_tempDensity;
		std::vector<float> m_density;

		std::vector<float> m_velX;
		std::vector<float> m_velY;

		std::vector<float> m_tempVelX;
		std::vector<float> m_tempVelY;

		std::vector<std::vector<SIM::colour>> m_currentState;

		std::unique_ptr<std::pair<int, int>> m_storedClick=nullptr;
	};
}