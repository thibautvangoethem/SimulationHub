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

		void setSize(int size) { m_size = size; }
		int getSize() const{ return m_size; }
		void setTimeBased(bool timeBased) { m_timeBased = timeBased; }
		int isTimeBased() const{ return m_timeBased; }
		void setunlockedTime(bool unlockedTime) { m_unlockedTime = unlockedTime; }
		int isUnlockedTime() const{ return m_unlockedTime; }
		void setDiffusion(double diff) { m_diffusion = diff; }
		double getDiffusion() const{ return m_diffusion; }
		void setViscosity(double vis) { m_viscosity = vis; }
		double getViscosity() const{ return m_viscosity; }
		void setConstantTickrate(double newRate) { m_constantTickRate = newRate; }
		double getConstantTickrate() const { return m_constantTickRate; }
	private:

		int m_size;
		bool m_timeBased;
		bool m_unlockedTime;
		double m_diffusion = 0;
		double m_viscosity = 0;
		double m_constantTickRate = 0.013;
	};
}