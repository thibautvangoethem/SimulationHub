#pragma once

namespace SIM{
	class SimulationSettings {
	public:
		SimulationSettings() = default;
		explicit SimulationSettings(int size, bool timeBased) :m_size{ size }, m_timeBased{ timeBased } {
		};
		~SimulationSettings() = default;
		/*SimulationSettings(SimulationSettings& settings) = delete;
		SimulationSettings(SimulationSettings&& settings) = delete;
		SimulationSettings& operator=(SimulationSettings& settings) = delete;
		SimulationSettings& operator=(SimulationSettings&& settings) = delete;*/

		void setSize(int size) { m_size = size; };
		int getSize() { return m_size; };
		void setTimeBased(bool timeBased) { m_timeBased = timeBased; };
		int isTimeBased() { return m_timeBased; };
		void setDiffusion(double diff) { m_diffusion = diff; };
		double getDiffusion() { return m_diffusion; };
		void setViscosity(double vis) { m_viscosity = vis; };
		double getViscosity() { return m_viscosity; };
	private:

		int m_size;
		bool m_timeBased;
		double m_diffusion = 0;
		double m_viscosity = 0;
	};
}