#include "FluidSim.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <immintrin.h>
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
	if(m_storedClick)
	{
		m_currentState[m_storedClick->second][m_storedClick->first]= colour{ 0,255,0 };
	}
	return m_currentState;
}

void FluidSim::addDensity(const int x,const int y,const float amount) {
	m_density[x][y] += amount;
}

void FluidSim::addVelocity(const int x,const int y,const float amount_x,const float amount_y) {
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


void FluidSim::setBoundary(const int b, std::vector<std::vector<float>>& array) {
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

/*
 * float cRecip = 1.0 / c;
	__m256 a_vec = _mm256_set1_ps(a); // Create a vector with all elements set to 'a'
	__m256 cRecip_vec = _mm256_set1_ps(cRecip); // Create a vector with all elements set to 'cRecip'
	for (unsigned int i = 1; i < m_size - 1; i++)
	{
		for (unsigned int j = 1; j < m_size - 1; j += 8) // Process 8 elements at a time
		{
			// Load 8 elements from the arrays into 256-bit vectors
			__m256 prevArray_vec = _mm256_loadu_ps(&prevArray[i][j]);
			__m256 array_up_vec = _mm256_loadu_ps(&array[i - 1][j]);
			__m256 array_down_vec = _mm256_loadu_ps(&array[i + 1][j]);
			__m256 array_left_vec = _mm256_loadu_ps(&array[i][j - 1]);
			__m256 array_right_vec = _mm256_loadu_ps(&array[i][j + 1]);
			
			// Add the 5 vectors together and multiply by 'a_vec'
			__m256 sum_vec = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(_mm256_add_ps(prevArray_vec, array_up_vec), array_down_vec), array_left_vec), array_right_vec);
			__m256 prod_vec = _mm256_mul_ps(a_vec, sum_vec);
			
			// Multiply by 'cRecip_vec' and store the result back into the array
			__m256 result_vec = _mm256_mul_ps(cRecip_vec, prod_vec);
			_mm256_storeu_ps(&array[i][j], result_vec);
		}
	}
 */
#ifdef AVX
void FluidSim::linSolveAvx(const int b, std::vector<std::vector<float>>& array, const std::vector<std::vector<float>>& prevArray, const float a, const float c)
{
	const float cRecip = 1.0 / c;
	const auto a_vec = _mm256_set1_ps(a); 
	const auto cRecip_vec = _mm256_set1_ps(cRecip); 
	for (unsigned int i = 1; i < m_size - 1; i++)
	{
		for (unsigned int j = 1; j < m_size - 1; j += 8) // Process 8 elements at a time due to avx2 taking 8 floats per operation
		{
			const auto prevArray_vec = _mm256_loadu_ps(&prevArray[i][j]);
			const auto array_up_vec = _mm256_loadu_ps(&array[i - 1][j]);
			const auto array_down_vec = _mm256_loadu_ps(&array[i + 1][j]);
			const auto array_left_vec = _mm256_loadu_ps(&array[i][j - 1]);
			const auto array_right_vec = _mm256_loadu_ps(&array[i][j + 1]);

			const auto sum_vec = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(array_up_vec, array_down_vec), array_left_vec), array_right_vec);
			const auto prod_vec = _mm256_mul_ps(a_vec, sum_vec);
			const auto sum2_vec = _mm256_add_ps(prevArray_vec, prod_vec);

			const auto result_vec = _mm256_mul_ps(cRecip_vec, sum2_vec);
			_mm256_storeu_ps(&array[i][j], result_vec);
		}
	}
}

void FluidSim::advectAvx(const int b, std::vector<std::vector<float>>& d, const std::vector<std::vector<float>>& prevD, const std::vector<std::vector<float>>& velocityX, const std::vector<std::vector<float>>& velocityY, const float dt)
{
	std::vector<int> resulti0i(8);
	std::vector<int> resulti1i(8);
	std::vector<int> resultj0i(8);
	std::vector<int> resultj1i(8);
	std::vector<float> prevDFirst(8, 0);
	std::vector<float> prevDSecond(8, 0);
	std::vector<float> prevDThird(8, 0);
	std::vector<float> prevDFourth(8, 0);
	// Initialize variables used for spatial interpolation
	const auto dtx_vec = _mm256_set1_ps(dt * (m_size - 2));
	const auto dty_vec = _mm256_set1_ps(dt * (m_size - 2));

	const auto one_vec = _mm256_set1_epi32(1);
	const auto one_vecps = _mm256_set1_ps(1.0f);
	const auto half_vecps = _mm256_set1_ps(0.5f);

	// Compute the size of the fluid grid
	float Nfloat = m_size - 2;
	auto lowerBound_vec = _mm256_set1_ps(0.5);
	auto higherBound_vec = _mm256_set1_ps(Nfloat + 0.5);
	auto indexAditor_vec= _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);

	// Iterate over each grid cell and advect its fluid property
	static int s = m_size - 1;
	for (unsigned int i = 1; i < s; i++)
	{
		for (unsigned int j = 1; j < m_size - 1; j+=8) // Process 8 elements at a time due to avx2 taking 8 floats per operation
		{
			auto velx_vec = _mm256_loadu_ps(&velocityX[i][j]);
			auto vely_vec = _mm256_loadu_ps(&velocityY[i][j]);
			// Compute the position at which to interpolate the density field
			auto tmp1_vec = _mm256_mul_ps(dtx_vec, velx_vec);
			auto tmp2_vec = _mm256_mul_ps(dty_vec, vely_vec);

			auto i_vec = _mm256_set1_ps(i);
			auto j_vec = _mm256_cvtepi32_ps(_mm256_add_epi32(indexAditor_vec, _mm256_set1_epi32(j)));
			auto x_vec = _mm256_fnmadd_ps(dtx_vec, velx_vec, i_vec); /*-(dtx*velx)+i*/
			auto y_vec = _mm256_fnmadd_ps(dty_vec, vely_vec, j_vec);/*-(dty*vely)+j*/

			
			// Clamp the position within the fluid grid
			auto clampedX_vec=_mm256_min_ps(_mm256_max_ps(x_vec, lowerBound_vec), higherBound_vec);
			
			auto i0_vec = _mm256_cvttps_epi32(clampedX_vec);
			auto i1_vec = _mm256_add_epi32(i0_vec, one_vec);

			auto clampedY_vec = _mm256_min_ps(_mm256_max_ps(y_vec, lowerBound_vec), higherBound_vec);
			auto j0_vec = _mm256_cvttps_epi32(clampedY_vec);
			auto j1_vec = _mm256_add_epi32(j0_vec, one_vec);

			auto i0_vecps = _mm256_cvtepi32_ps(i0_vec);
			auto i1_vecps = _mm256_cvtepi32_ps(i1_vec);
			auto j0_vecps = _mm256_cvtepi32_ps(j0_vec);
			auto j1_vecps = _mm256_cvtepi32_ps(j1_vec);
			// Compute the interpolation coefficients
			
			auto s1_vec = _mm256_sub_ps(clampedX_vec, i0_vecps);
			auto s0_vec = _mm256_sub_ps(one_vecps, s1_vec);
			auto t1_vec = _mm256_sub_ps(clampedY_vec, j0_vecps);
			auto t0_vec = _mm256_sub_ps(one_vecps, t1_vec);


			// Perform bilinear interpolation of the density field at the interpolated position
			// and store the result in the current density field

			auto i0i_vec = _mm256_cvttps_epi32(_mm256_add_ps(i0_vecps,half_vecps));
			auto i1i_vec = _mm256_cvttps_epi32(_mm256_add_ps(i1_vecps, half_vecps));
			auto j0i_vec = _mm256_cvttps_epi32(_mm256_add_ps(j0_vecps, half_vecps));
			auto j1i_vec = _mm256_cvttps_epi32(_mm256_add_ps(j1_vecps, half_vecps));

			//Annoyingly have to go out of avx mode as I used nested vectors which cant be accessed in a flattened way  for use with gather avx function
			//TODO change this so _mm256_i32gather_epi32 can be used
			
			_mm256_storeu_epi32(&resulti0i[0], i0i_vec);
			_mm256_storeu_epi32(&resulti1i[0], i1i_vec);
			_mm256_storeu_epi32(&resultj0i[0], j0i_vec);
			_mm256_storeu_epi32(&resultj1i[0], j1i_vec);
			for(int prevDIterator=0; prevDIterator <8; prevDIterator++)
			{
				prevDFirst[prevDIterator] = prevD[resulti0i[prevDIterator]][resultj0i[prevDIterator]];
				prevDSecond[prevDIterator] = prevD[resulti0i[prevDIterator]][resultj1i[prevDIterator]];
				prevDThird[prevDIterator] = prevD[resulti1i[prevDIterator]][resultj0i[prevDIterator]];
				prevDFourth[prevDIterator] = prevD[resulti1i[prevDIterator]][resultj1i[prevDIterator]];
			}
			auto prevDFirst_vec = _mm256_loadu_ps(&prevDFirst[0]);
			auto prevDSecond_vec = _mm256_loadu_ps(&prevDSecond[0]);
			auto prevDThird_vec = _mm256_loadu_ps(&prevDThird[0]);
			auto prevDFourth_vec = _mm256_loadu_ps(&prevDFourth[0]);

			auto part1 = _mm256_fmadd_ps(t0_vec, prevDFirst_vec, _mm256_mul_ps(t1_vec, prevDSecond_vec));
			auto part2 = _mm256_fmadd_ps(t0_vec, prevDThird_vec, _mm256_mul_ps(t1_vec, prevDFourth_vec));
			auto result_vec = _mm256_fmadd_ps(s0_vec, part1, _mm256_mul_ps(s1_vec, part2));
			_mm256_storeu_ps(&d[i][j], result_vec);
			
		}
	}
}

#endif

void FluidSim::linSolve(const int b, std::vector<std::vector<float>>& array, const std::vector<std::vector<float>>& prevArray,const float a,const float c) {
	
#ifdef AVX
	linSolveAvx(b, array, prevArray, a, c);
#else
	float cRecip = 1.0 / c;
	
	for (unsigned int i = 1; i < m_size - 1; i++)
	{
		for (unsigned int j = 1; j < m_size - 1; j++)
		{
			array[i][j] = (prevArray[i][j] + a * (array[i + 1][j] + array[i - 1][j] + array[i][j + 1] + array[i][j - 1])) * cRecip;
		}
	}
#endif

	this->setBoundary(b, array);
}
void FluidSim::diffuse(const int b, std::vector<std::vector<float>>& array, const std::vector<std::vector<float>>& prevArray, const float diff, const float dt) {
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
	__m256i mask = _mm256_setr_epi32(-20, -72, -48, -9, -100, 3, 5, 8);
	this->setBoundary(1, velocityX);
	this->setBoundary(2, velocityY);
}

// This function performs advection of a fluid property (represented by the density field 'd') 
// using the velocity field (represented by 'velocityX' and 'velocityY') over a timestep 'dt'
// for a given boundary 'b'. It uses the previous state of the density field ('prevD') as input.
void FluidSim::advect(const int b, std::vector<std::vector<float>>& d, const std::vector<std::vector<float>>& prevD, const std::vector<std::vector<float>>& velocityX,const std::vector<std::vector<float>>& velocityY, const float dt) {
#ifdef AVX
	advectAvx(b, d, prevD, velocityX,velocityY, dt);
#else
	// Initialize variables used for spatial interpolation
	float i0 = 0, i1 = 0, j0 = 0, j1 = 0;
	float dtx = dt * (m_size - 2);
	float dty = dt * (m_size - 2);

	float s0 = 0, s1 = 0, t0 = 0, t1 = 0;
	float tmp1 = 0, tmp2 = 0, x = 0, y = 0;

	// Compute the size of the fluid grid
	float Nfloat = m_size - 2;

	// Iterate over each grid cell and advect its fluid property
	static int s = m_size - 1;
	for (unsigned int i = 1; i < s; i++)
	{
		for (unsigned int j = 1; j < m_size - 1; j++)
		{
			// Compute the position at which to interpolate the density field
			tmp1 = dtx * velocityX[i][j];
			tmp2 = dty * velocityY[i][j];
			x = i - tmp1;
			y = j - tmp2;

			// Clamp the position within the fluid grid
			std::clamp(x, static_cast<float>(0.5), static_cast<float>(Nfloat + 0.5));
			i0 = static_cast<int>(x);
			i1 = i0 + 1;
			std::clamp(y, static_cast<float>(0.5), static_cast<float>(Nfloat + 0.5));
			j0 = static_cast<int>(y);
			j1 = j0 + 1;

			// Compute the interpolation coefficients
			s1 = x - i0;
			s0 = 1 - s1;
			t1 = y - j0;
			t0 = 1 - t1;

			// Perform bilinear interpolation of the density field at the interpolated position
			// and store the result in the current density field
			int i0i = static_cast<int>(i0 + 0.5);
			int i1i = static_cast<int>(i1 + 0.5);
			int j0i = static_cast<int>(j0 + 0.5);
			int j1i = static_cast<int>(j1 + 0.5);
			auto te= s0 * (t0 * prevD[i0i][j0i] + t1 * prevD[i0i][j1i]) + s1 * (t0 * prevD[i1i][j0i] + t1 * prevD[i1i][j1i]);
			d[i][j] = te;
		}
	}
#endif


	// Apply the specified boundary conditions to the current density field
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
