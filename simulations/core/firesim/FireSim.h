#pragma once
#include "../baseinterface/Simulation.h"
#include "RandomNumberQueue.h"
namespace SIM {
	class SimulationSetting;
}

namespace SIM {
	class FireSim: public SIM::Simulation {
		static constexpr SIM::colour EMPTY = SIM::colour{ 0,0,0 };
	public:
		/*#### Sim::Simulation functions ####*/
		FireSim(std::shared_ptr<SIM::SimulationSettings> settings);
		void advance(const double timestep) final;
		//not used here
		void handleClick(const bool isLeftClick,const int xpos,const int ypos) final{};
	protected:
		
		void advancePixel(int x, int y);

		
		std::vector<std::vector<SIM::colour>> m_currentState;
		std::vector<std::vector<int>> m_currentStateNonMapped;
		std::vector<SIM::colour> m_colorbuffer;

		RandomNumberQueue m_random1{ 1000000,2};
		RandomNumberQueue m_random2{ 1000000,2};
	};
}