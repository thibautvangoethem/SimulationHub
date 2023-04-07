#include "RayTracer.h"

#include <iostream>

using namespace SIM;

RayTracer::RayTracer(std::shared_ptr<SIM::SimulationSettings> settings) : Simulation(std::move(settings))
{
    m_size = this->m_settings->getSize();
	m_currentState = std::vector<std::vector<SIM::colour> >(
		m_size,
		std::vector<SIM::colour>(m_size));
}

void RayTracer::advance(const double timestep)
{
	//todo
}

std::vector<std::vector<SIM::colour>>& RayTracer::getCurrentState()
{
    for (int j = 0; j<m_size; ++j) {
        for (int i = 0; i < m_size; ++i) {
            auto r = double(i) / (m_size - 1);
            auto g = double(j) / (m_size - 1);
            auto b = 0.25;

            auto ir = static_cast<uint8_t>(255.999 * r);
            auto ig = static_cast<uint8_t>(255.999 * g);
            auto ib = static_cast<uint8_t>(255.999 * b);

            m_currentState[m_size-j-1][i] = colour{ ir,ig,ib };
        }
    }
    return m_currentState;
}

void RayTracer::handleClick(const bool isLeftClick, const int xpos, const int ypos)
{
	//todo
}