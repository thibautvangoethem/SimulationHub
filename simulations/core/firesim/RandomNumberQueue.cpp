#include "RandomNumberQueue.h"

#include <chrono>
#include <random>

using namespace SIM;

RandomNumberQueue::RandomNumberQueue(int size, int max) :m_size{size} {
	std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());
	for (unsigned int i = 0; i < size; i++) {
		m_queue.push_back(rng() % max);
	}
}

int RandomNumberQueue::getNext() noexcept {
	m_currentIndex++;
	m_currentIndex %= m_size;
	return m_queue[m_currentIndex];

}