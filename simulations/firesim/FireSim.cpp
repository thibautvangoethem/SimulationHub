/**
 * @file
 * an implementation of a fire effect almost identitcal to the one found here http://fabiensanglard.net/doom_fire_psx/
 * @author Thibaut Van Goethem
 */
#include "FireSim.h"

#include "../baseinterface/SimulationSettings.h"

#include <iostream>
 //#include <SFML/Graphics.hpp>
#include <iostream>
#include <random>

using namespace FSIM;

FireSim::FireSim(std::shared_ptr<SIM::SimulationSettings> settings) : SIM::Simulation(std::move(settings))
{

	m_currentState = std::vector<std::vector<SIM::colour> >(
		this->m_settings->getSize(),
		std::vector<SIM::colour>(this->m_settings->getSize()));

	m_colorbuffer.emplace_back(7, 7, 7);
	m_colorbuffer.emplace_back(31, 7, 7);
	m_colorbuffer.emplace_back(47, 15, 7);
	m_colorbuffer.emplace_back(71, 15, 7);
	m_colorbuffer.emplace_back(87, 23, 7);
	m_colorbuffer.emplace_back(103, 31, 7);
	m_colorbuffer.emplace_back(119, 31, 7);
	m_colorbuffer.emplace_back(143, 39, 7);
	m_colorbuffer.emplace_back(159, 47, 7);
	m_colorbuffer.emplace_back(175, 63, 7);
	m_colorbuffer.emplace_back(191, 71, 7);
	m_colorbuffer.emplace_back(199, 71, 7);
	m_colorbuffer.emplace_back(223, 79, 7);
	m_colorbuffer.emplace_back(223, 87, 7);
	m_colorbuffer.emplace_back(223, 87, 7);
	m_colorbuffer.emplace_back(215, 95, 7);
	m_colorbuffer.emplace_back(215, 95, 7);
	m_colorbuffer.emplace_back(215, 103, 15);
	m_colorbuffer.emplace_back(207, 111, 15);
	m_colorbuffer.emplace_back(207, 119, 15);
	m_colorbuffer.emplace_back(207, 127, 15);
	m_colorbuffer.emplace_back(207, 135, 23);
	m_colorbuffer.emplace_back(199, 135, 23);
	m_colorbuffer.emplace_back(199, 143, 23);
	m_colorbuffer.emplace_back(199, 151, 31);
	m_colorbuffer.emplace_back(191, 159, 31);
	m_colorbuffer.emplace_back(191, 159, 31);
	m_colorbuffer.emplace_back(191, 167, 39);
	m_colorbuffer.emplace_back(191, 167, 39);
	m_colorbuffer.emplace_back(191, 175, 47);
	m_colorbuffer.emplace_back(183, 175, 47);
	m_colorbuffer.emplace_back(183, 183, 47);
	m_colorbuffer.emplace_back(183, 183, 55);
	m_colorbuffer.emplace_back(207, 207, 111);
	m_colorbuffer.emplace_back(223, 223, 159);
	m_colorbuffer.emplace_back(239, 239, 199);
	m_colorbuffer.emplace_back(255, 255, 255);

	m_currentStateNonMapped = std::vector<std::vector<int> >(
		this->m_settings->getSize(),
		std::vector<int>(this->m_settings->getSize()));
	m_currentStateNonMapped.push_back(std::vector<int>(this->m_settings->getSize(), m_colorbuffer.size()-1));
	m_currentStateNonMapped.push_back(std::vector<int>(this->m_settings->getSize(), m_colorbuffer.size()-1));
}

void FireSim::advance(double timestep) {
	for (int i = 0; i < m_settings->getSize(); ++i) {
		for (int j = 0; j < m_settings->getSize()+1; ++j) {
			advancePixel(i, j);
		}
	}
}

std::vector<std::vector<SIM::colour>>& FireSim::getCurrentState() {
	return m_currentState;
}


/**
 * this function calculates the spread of the fire for a single pixel
 * @param pixelArray the array of fire pixels
 * @param x
 * @param y
 */
void FireSim::advancePixel(int x, int y) {
    //2 random variables are used one to randomise the amount of decay and one to randomize the movement of the fire going up
    /*int rand = std::rand() % 3;
    int dst = static_cast<int>(round((std::rand() % 2)));*/
	int rand = m_random1.getNext();
	int dst=m_random2.getNext();

    int targetx = x + dst;
    if (targetx >m_settings->getSize() - 1)targetx -= (m_settings->getSize());
    if (y <= 0) return;
	m_currentStateNonMapped[y - 1][x] = m_currentStateNonMapped[y][targetx] - rand;
	if (m_currentStateNonMapped[y - 1][x] > 0) {
		m_currentState[y - 1][x] = m_colorbuffer[m_currentStateNonMapped[y - 1][x]];
	}
	else {
		m_currentState[y - 1][x] = EMPTY;
	}
}