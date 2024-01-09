#include "CudaInterop.h"

CSIM::CudaInterop::CudaInterop(std::shared_ptr<SIM::SimulationSettings> settings):SIM::Simulation(std::move(settings)), m_size(this->m_settings->getSize())
{
}

void CSIM::CudaInterop::advance(const double timestep)
{
}

void CSIM::CudaInterop::handleClick(const bool isLeftClick, const int xpos, const int ypos)
{
}
