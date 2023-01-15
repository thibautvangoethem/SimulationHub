#include "FluidSim.h"
#include <cmath>
#include <iostream>

using namespace SIM;

FluidSim::FluidSim(std::unique_ptr<SIM::SimulationSettings> settings) : Simulation(std::move(settings)) {
	m_size = this->m_settings->getSize();
	m_s = grid(
		m_size,
		std::vector<float>(m_size));
	m_density = grid(
		m_size,
		std::vector<float>(m_size));
	m_velX = grid(
		m_size,
		std::vector<float>(m_size));
	m_velY = grid(
		m_size,
		std::vector<float>(m_size));
	m_oldVelX = grid(
		m_size,
		std::vector<float>(m_size));
	m_oldVelY = grid(
		m_size,
		std::vector<float>(m_size));

	m_currentState = std::vector<std::vector<SIM::colour> >(
		m_size,
		std::vector<SIM::colour>(m_size));


	//temp
	m_paintSources.emplace_back(m_size / 2, m_size / 2);
	auto temp = Flow(std::make_pair(m_size / 2, m_size / 2), std::make_pair(1,1));
	m_flows.emplace_back(temp);
}
void FluidSim::advance(double timestep) {
	applySources();

	diffuse(1, m_oldVelX, m_velX, m_settings->getViscosity(), timestep);
	diffuse(2, m_oldVelY, m_velY, m_settings->getViscosity(), timestep);

	project(m_oldVelX, m_oldVelY, m_velX, m_velY);

	advect(1, m_velX, m_oldVelX, m_oldVelX, m_oldVelY, timestep);
	advect(2, m_velY, m_oldVelY, m_oldVelX, m_oldVelY, timestep);

	project(m_velX, m_velY, m_oldVelX, m_oldVelY);

	diffuse(0, m_s, m_density, m_settings->getDiffusion(), timestep);
	advect(0, m_density, m_s, m_velX, m_velY, timestep);


}
std::vector<std::vector<SIM::colour>>& FluidSim::getCurrentState() {
	for (int i = 0; i < m_size; ++i)
	{
		for (int j = 0; j < m_size; ++j)
		{

			int val = m_s[i][j]* 52;
			if (val > 255)val = 255;
			if (val < 0)val = 0;
			auto valCasted = (uint8_t)(val);
			//TODO asses the performance impact of this calculation
			float combinedSpeedVal =1-sqrt(m_velY[i][j]* m_velY[i][j] + m_velX[i][j]* m_velX[i][j]);
			
			m_currentState[j][i] = colour{ valCasted,(uint8_t)(valCasted* combinedSpeedVal),(uint8_t)(valCasted* combinedSpeedVal) };
		}
	}
	return m_currentState;
}

void FluidSim::addDensity(int x, int y, float amount) {
	m_density[x][y] += amount;
}

void FluidSim::addVelocity(int x, int y, float amount_x, float amount_y) {
	m_velX[x][y] += amount_x;
	m_velY[x][y] += amount_y;
}

void FluidSim::applySources() {
	for (auto& flow : m_flows) {
		m_velX[flow.location.first][flow.location.second] = flow.velocity.first;
		m_velY[flow.location.first][flow.location.second] = flow.velocity.second;
	}

	for (auto& source : m_paintSources) {
			addDensity(source.first, source.second, 2);
	}

}


void FluidSim::setBoundary(int b, std::vector<std::vector<float>>& array) {
	for (unsigned int i = 1; i < m_size - 1; i++)
	{
		array[i][0] = (b == 2 ? -array[i][1] : array[i][1]);
		array[i][m_size - 1] = (b == 2 ? -array[i][m_size - 2] : array[i][m_size - 2]);
	}
	for (unsigned int i = 1; i < m_size - 1; i++)
	{
		array[0][i] = (b == 1 ? -array[1][i] : array[1][i]);
		array[m_size - 1][i] = (b == 1 ? -array[m_size - 2][i] : array[m_size - 2][i]);
	}
	array[0][0] = 0.5 * array[1][0] + array[0][1];
	array[0][m_size - 1] = 0.5 * (array[1][m_size - 1] + array[0][m_size - 2]);
	array[m_size - 1][0] = 0.5 * (array[m_size - 2][0] + array[m_size - 1][1]);
	array[m_size - 1][m_size - 1] = 0.5 * (array[m_size - 2][m_size - 1] + array[m_size - 1][m_size - 2]);
}

void FluidSim::linSolve(int b, std::vector<std::vector<float>>& array, std::vector<std::vector<float>>& prevArray, float a, float c) {
	float cRecip = 1.0 / c;
	for (unsigned int i = 1; i < m_size - 1; i++)
	{
		for (unsigned int j = 1; j < m_size - 1; j++)
		{
			array[i][j] = (prevArray[i][j] + a * (array[i + 1][j] + array[i - 1][j] + array[i][j + 1] + array[i][j - 1])) * cRecip;
		}
	}
	this->setBoundary(b, array);
}
void FluidSim::diffuse(int b, std::vector<std::vector<float>>& array, std::vector<std::vector<float>>& prevArray, float diff, float dt) {
	float a = dt * diff * (m_size - 2) * (m_size - 2);
	linSolve(b, array, prevArray, a, 1 + 6 * a);

}
void FluidSim::project(std::vector<std::vector<float>>& velocityX, std::vector<std::vector<float>>& velocityY, std::vector<std::vector<float>>& clearVector, std::vector<std::vector<float>>& targetVectory) {
	for (unsigned int j = 1; j < m_size - 1; j++)
	{
		for (unsigned int i = 1; i < m_size - 1; i++)
		{
			targetVectory[j][i] = - 0.5 * (velocityX[j + 1][i] - velocityX[j - 1][i] + velocityY[j][i + 1] - velocityY[j][i - 1]) / m_size;
			clearVector[j][i] = 0;
		}
	}
	this->setBoundary(0, targetVectory);
	this->setBoundary(0, clearVector);
	this->linSolve(0, clearVector, targetVectory, 1, 6);

	for (unsigned int i = 1; i < m_size - 1; i++)
	{
		for (unsigned int j = 1; j < m_size - 1; j++)
		{
			velocityX[i][j] = velocityX[i][j] - 0.5 * (clearVector[i + 1][j] - clearVector[i - 1][j]) * m_size;
			velocityY[i][j] = velocityY[i][j] - 0.5 * (clearVector[i][j + 1] - clearVector[i][j - 1]) * m_size;
		}
	}
	this->setBoundary(1, velocityX);
	this->setBoundary(2, velocityY);
}
void FluidSim::advect(int b, std::vector<std::vector<float>>& d, std::vector<std::vector<float>>& prevD, std::vector<std::vector<float>>& velocityX, std::vector<std::vector<float>>& velocitY, float dt) {
	float i0=0, i1=0, j0=0, j1=0;
	float dtx = dt * (m_size - 2);
	float dty = dt * (m_size - 2);

	float s0=0, s1=0, t0=0, t1=0;
	float tmp1=0, tmp2=0, x=0, y=0;

	float Nfloat = m_size - 2;
	for (unsigned int i = 1; i < m_size - 1; i++)
	{
		for (unsigned int j = 1; j < m_size - 1; j++)
		{
			tmp1 = dtx * velocityX[i][j];
			tmp2 = dty * velocitY[i][j];
			x = i - tmp1;
			y = j - tmp2;
			if (x < 0.5)x = 0.5;
			if (x > Nfloat + 0.5)x = Nfloat + 0.5;
			i0 = std::floor(x);
			i1 = i0 + 1;
			if (y < 0.5)x = 0.5;
			if (y > Nfloat + 0.5)y = Nfloat + 0.5;
			j0 = std::floor(y);
			j1 = j0 + 1;

			s1 = x - i0;
			s0 = 1 - s1;
			t1 = y - j0;
			t0 = 1 - t1;

			//std::round is extremely slow for some reason, this should be nearly equivalent within the needed accuracy
			int i0i = std::floor(i0+0.5);
			int i1i = std::floor(i1+0.5);
			int j0i = std::floor(j0+0.5);
			int j1i = std::floor(j1+0.5);
			//This is by far the slowest line in the entire program
			d[i][j] = s0 * (t0 * prevD[i0i][j0i] + t1 * prevD[i0i][j1i]) + s1 * (t0 * prevD[i1i][j0i] + t1 * prevD[i1i][j1i]);
			
		}
	}
	this->setBoundary(b, d);
}

void FluidSim::handleClick(bool isLeftClick, int xpos, int ypos)
{
	if (xpos >= m_size || ypos >= m_size)return;
	if(!isLeftClick)
	{
		m_paintSources.emplace_back(xpos, ypos);
	}else
	{
		if(m_storedClick==nullptr)
		{
			m_storedClick = std::make_unique<std::pair<int, int>>(xpos, ypos);
		}else
		{
			float xvel= -(m_storedClick->first- xpos) / (m_size/ 10);
			float yvel = -(m_storedClick->second - ypos) / (m_size / 10);
			m_flows.emplace_back(std::make_pair(m_storedClick->first, m_storedClick->second), std::make_pair(xvel,yvel));
			
			m_storedClick = nullptr;
		}
	}

}
