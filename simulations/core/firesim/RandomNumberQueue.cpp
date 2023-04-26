#include "RandomNumberQueue.h"

#include <random>

using namespace SIM;

RandomNumberQueue::RandomNumberQueue(int size, int max) :m_size{size} {
	for (unsigned int i = 0; i < size; i++) {
		m_queue.push_back(std::rand() % max);
	}
}

int RandomNumberQueue::getNext() noexcept {
	m_currentIndex++;
	m_currentIndex %= m_size;
	return m_queue[m_currentIndex];

}