#include "SlimeSim.h"

#include <algorithm>

SSIM::SlimeSim::SlimeSim(std::shared_ptr<SIM::SimulationSettings> settings):Simulation(std::move(settings)), m_size(this->m_settings->getSize())
{
	m_field = std::vector<float>(m_size * m_size, 0);
	m_agents =std::vector<std::unique_ptr<SlimeAgent>>();
	
	for (int i = 0; i < 4096; i++) {
		//m_agents.push_back(std::make_unique<SlimeAgent>(m_size / 2, m_size / 2, (static_cast <float> (rand()) / static_cast <float> (RAND_MAX))*2*M_PI));
		m_agents.push_back(std::make_unique<SlimeAgent>(m_size* ((rand()*1.0)/ RAND_MAX), m_size * ((rand() * 1.0) / RAND_MAX), (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) * 2 * M_PI));
	}
}

void SSIM::SlimeSim::advance(const double timestep)
{
	advanceAgents(timestep);
	blurDiffuse(timestep);
}

void SSIM::SlimeSim::handleClick(const bool isLeftClick, const int xpos, const int ypos)
{
}

void SSIM::SlimeSim::advanceAgents(const double timestep) noexcept
{
	for (auto& agent : m_agents) {

		slimeSteerUpdate(*agent, timestep);

		const auto dirx = cos(agent->angle);
		auto newx = agent->posx + dirx * speed * timestep;
		
		if (newx < 0.0 || newx >= m_size * 1.0) {
			newx = std::clamp(newx, 0.0, (m_size - 1) * 1.0);
			agent->angle = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) * 2 * M_PI;
		}

		const auto diry = sin(agent->angle);
		auto newy = agent->posy + diry * speed * timestep;
		
		if (newy < 0.0 || newy >= m_size * 1.0) {
			newy = std::clamp(newy, 0.0, (m_size-1) * 1.0);
			agent->angle = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) * 2 * M_PI;
		}
		agent->posy = newy;
		agent->posx = newx;

		m_field[int(newx)+int(newy)*m_size]=1.0;

	}
}

void SSIM::SlimeSim::blurDiffuse(const double timestep)
{
	
	for (int i = 0; i < m_size; i++) {
		for (int j = 0; j < m_size; j++) {
			m_field[i + j * m_size] = std::max(0.0, m_field[i + j * m_size] - (evaporate * timestep));
		}
	}
	for (int i = 0; i < m_size; i++) {
		for (int j = 0; j < m_size; j++) {
			float sum = 0.0;
			for (int dx = -1; dx <= 1; dx++) {
				for (int dy = -1; dy <= 1; dy++) {
					int sx = i + dx;
					int sy = j + dy;
					if (sx >= 0 && sx < m_size && sy >= 0 && sy < m_size) {
						sum += m_field[sx + sy * m_size];
					}
				}
			}
			m_field[i + j * m_size] = std::lerp(m_field[i + j * m_size], sum / 9.0, 10 * timestep);

		}
	}

}

void SSIM::SlimeSim::slimeSteerUpdate(SlimeAgent& agent, const double timestep)
{
	float  wForward = sense(agent, 0);
	float  wLeft = sense(agent, M_PI/4);
	float  wRight = sense(agent, -(M_PI /4));

	//TODO make random
	float randomStrength = 1;


	if (wForward > wLeft && wForward > wRight) {
		return;
	}
	
	if (wLeft > wRight) {
		agent.angle += randomStrength * steerStrength *timestep;
	}
	else {
		agent.angle -= randomStrength * steerStrength * timestep;
	}
	
}

float SSIM::SlimeSim::sense(SlimeAgent& agent, float angle)
{
	float sAngle = agent.angle + angle;
	float sDirx = cos(sAngle);
	int sPosx = agent.posx + sDirx * sensorDirOffset;
	float sDiry = sin(sAngle);
	int sPosy = agent.posy + sDiry * sensorDirOffset;

	float sum = 0.0;

	for (int i = -sensorSize; i <= sensorSize; i++) {
		for (int j = -sensorSize; j <= sensorSize; j++) {
			int checkPosx = sPosx + i;
			int checkPosy = sPosy + j;

			if (checkPosx >= 0 && checkPosx < m_size && checkPosy >= 0 && checkPosy < m_size) {
				sum += m_field[checkPosx + checkPosy * m_size];
			}
		}
	}
	return sum;
}
