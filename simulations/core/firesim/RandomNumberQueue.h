#pragma once

#include <vector>
//A class that caches random numbers to use as a faster pseudorandom function
namespace SIM {
	class RandomNumberQueue {
	public:
		RandomNumberQueue(int size, int max);
		int getNext() noexcept;
	private:
		int m_currentIndex{0};
		int m_size;
		std::vector<int> m_queue;

	};
}