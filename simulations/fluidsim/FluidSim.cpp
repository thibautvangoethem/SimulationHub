#include "FluidSim.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <immintrin.h>
using namespace SIM;

FluidSim::FluidSim(std::shared_ptr<SIM::SimulationSettings> settings) : Simulation(std::move(settings)) {
	m_size = this->m_settings->getSize();
	m_s = std::vector<float>(m_size * m_size, 0);
	m_density = std::vector<float>(m_size * m_size, 0);
	m_velX = std::vector<float>(m_size * m_size, 0);
	m_velY = std::vector<float>(m_size * m_size, 0);
	m_oldVelX = std::vector<float>(m_size * m_size, 0);
	m_oldVelY = std::vector<float>(m_size * m_size, 0);

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

			int val = m_s[(i)*m_size + (j)]* 52;
			if (val > 255)val = 255;
			if (val < 0)val = 0;
			auto valCasted = (uint8_t)(val);
			//TODO asses the performance impact of this calculation
			float combinedSpeedVal =1-sqrt(m_velY[(i)*m_size + (j)]* m_velY[(i)*m_size + (j)] + m_velX[(i)*m_size + (j)]* m_velX[(i)*m_size + (j)]);
			std::clamp(valCasted, (uint8_t)0, (uint8_t)255);
			auto g = (uint8_t)(valCasted * combinedSpeedVal);
			std::clamp(g, (uint8_t)0, (uint8_t)255);
			
			
			m_currentState[(j)][(i)] = colour{ valCasted,g,g };
		}
	}
	if(m_storedClick)
	{
		m_currentState[(m_storedClick->second)][(m_storedClick->first)]= colour{ 0,255,0 };
	}
	return m_currentState;
}

void FluidSim::addDensity(const int x,const int y,const float amount) {
	m_density[(x)*m_size + (y)] += amount;
}

void FluidSim::addVelocity(const int x,const int y,const float amount_x,const float amount_y) {
	m_velX[(x)*m_size + (y)] += amount_x;
	m_velY[(x)*m_size + (y)] += amount_y;
}

void FluidSim::applySources() {
	for (auto& flow : m_flows) {
		m_velX[(flow.location.first)*m_size + (flow.location.second)] = flow.velocity.first;
		m_velY[(flow.location.first)*m_size + (flow.location.second)] = flow.velocity.second;
	}

	for (auto& source : m_paintSources) {
			addDensity(source.first, source.second, 2);
	}

}


void FluidSim::setBoundary(const int b, std::vector<float>& array) {
	for (unsigned int i = 1; i < m_size - 1; i++)
	{
		array[(i)*m_size + (0)] = (b == 2 ? -array[(i)*m_size + (1)] : array[(i)*m_size + (1)]);
		array[(i)*m_size + (m_size - 1)] = (b == 2 ? -array[(i)*m_size + (m_size - 2)] : array[(i)*m_size + (m_size - 2)]);
	}
	for (unsigned int i = 1; i < m_size - 1; i++)
	{
		array[(0)*m_size + (i)] = (b == 1 ? -array[(1)*m_size + (i)] : array[(1)*m_size + (i)]);
		array[(m_size - 1)*m_size + (i)] = (b == 1 ? -array[(m_size - 2)*m_size + (i)] : array[(m_size - 2)*m_size + (i)]);
	}
	array[(0)*m_size + (0)] = 0.5 * array[(1)*m_size + (0)] + array[(0)*m_size + (1)];
	array[(0)*m_size + ((m_size - 1))] = 0.5 * (array[((1))*m_size + ((m_size - 1))] + array[(0)*m_size + (m_size - 2)]);
	array[(m_size - 1)*m_size + (0)] = 0.5 * (array[(m_size - 2)*m_size + (0)] + array[(m_size - 1)*m_size + (1)]);
	array[(m_size - 1)*m_size + (m_size - 1)] = 0.5 * (array[(m_size - 2)*m_size + (m_size - 1)] + array[(m_size - 1)*m_size + (m_size - 2)]);
}

void FluidSim::linSolve(const int b, std::vector<float>& array, const std::vector<float>& prevArray,const float a,const float c) {
	
#ifdef AVX
	linSolveAvx(b, array, prevArray, a, c);
#else
	float cRecip = 1.0 / c;
	
	for (unsigned int i = 1; i < m_size - 1; i++)
	{
		for (unsigned int j = 1; j < m_size - 1; j++)
		{
			array[(i)*m_size + (j)] = (prevArray[(i)*m_size + (j)] + a * (array[(i + 1)*m_size + (j)] + array[(i - 1)*m_size + (j)] + array[(i)*m_size + (j + 1)] + array[(i)*m_size + (j - 1)])) * cRecip;
		}
	}
#endif

	this->setBoundary(b, array);
}
void FluidSim::diffuse(const int b, std::vector<float>& array, const std::vector<float>& prevArray, const float diff, const float dt) {
	float a = dt * diff * (m_size - 2) * (m_size - 2);
	linSolve(b, array, prevArray, a, 1 + 6 * a);

}
void FluidSim::project(std::vector<float>& velocityX, std::vector<float>& velocityY, std::vector<float>& clearVector, std::vector<float>& targetVectory) {
#ifdef AVX
	projectAvx(velocityX, velocityY, clearVector, targetVectory);
#else
	for (unsigned int j = 1; j < m_size - 1; j++)
	{
		for (unsigned int i = 1; i < m_size - 1; i++)
		{
			targetVectory[(j)*m_size + (i)] = - 0.5 * (velocityX[(j + 1)*m_size + (i)] - velocityX[(j - 1)*m_size + (i)] + velocityY[(j)*m_size + (i + 1)] - velocityY[(j)*m_size + (i - 1)]) / m_size;
		}
	}
	std::fill(clearVector.begin(), clearVector.end(), 0);
	this->setBoundary(0, targetVectory);
	this->setBoundary(0, clearVector);
	this->linSolve(0, clearVector, targetVectory, 1, 6);

	for (unsigned int i = 1; i < m_size - 1; i++)
	{
		for (unsigned int j = 1; j < m_size - 1; j++)
		{
			velocityX[(i)*m_size + (j)] = velocityX[(i)*m_size + (j)] - 0.5 * (clearVector[(i + 1) * m_size + (j)] - clearVector[(i - 1) * m_size + (j)]) * m_size;
			velocityY[(i)*m_size + (j)] = velocityY[(i)*m_size + (j)] - 0.5 * (clearVector[(i)*m_size + (j + 1)] - clearVector[(i)*m_size + (j - 1)]) * m_size;
		}
	}
#endif

	this->setBoundary(1, velocityX);
	this->setBoundary(2, velocityY);
}

// This function performs advection of a fluid property (represented by the density field 'd') 
// using the velocity field (represented by 'velocityX' and 'velocityY') over a timestep 'dt'
// for a given boundary 'b'. It uses the previous state of the density field ('prevD') as input.
void FluidSim::advect(const int b, std::vector<float>& d, const std::vector<float>& prevD, const std::vector<float>& velocityX,const std::vector<float>& velocityY, const float dt) {
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
			tmp1 = dtx * velocityX[(i)*m_size + (j)];
			tmp2 = dty * velocityY[(i)*m_size + (j)];
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
			auto te= s0 * (t0 * prevD[(i0i)*m_size + (j0i)] + t1 * prevD[(i0i)*m_size + (j1i)]) + s1 * (t0 * prevD[(i1i)*m_size + (j0i)] + t1 * prevD[(i1i)*m_size + (j1i)]);
			d[(i)*m_size + (j)] = te;
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
			float xvel= static_cast<float>(- (m_storedClick->first - xpos)) / (m_size / 10);
			float yvel = static_cast<float> (- (m_storedClick->second - ypos)) / (m_size / 10);
			m_flows.emplace_back(std::make_pair(m_storedClick->first, m_storedClick->second), std::make_pair(xvel,yvel));
			
			m_storedClick = nullptr;
		}
	}

}

#ifdef AVX
void FluidSim::linSolveAvx(const int b, std::vector<float>& array, const std::vector<float>& prevArray, const float a, const float c)
{
	const float cRecip = 1.0 / c;
	const auto a_vec = _mm256_set1_ps(a);
	const auto cRecip_vec = _mm256_set1_ps(cRecip);
	for (unsigned int i = 1; i < m_size - 1; i++)
	{
		for (unsigned int j = 1; j < m_size - 1; j += 8) // Process 8 elements at a time due to avx2 taking 8 floats per operation
		{
			const auto prevArray_vec = _mm256_loadu_ps(&prevArray[(i)*m_size + (j)]);
			const auto array_up_vec = _mm256_loadu_ps(&array[(i - 1) * m_size + (j)]);
			const auto array_down_vec = _mm256_loadu_ps(&array[(i + 1) * m_size + (j)]);
			const auto array_left_vec = _mm256_loadu_ps(&array[(i)*m_size + (j - 1)]);
			const auto array_right_vec = _mm256_loadu_ps(&array[(i)*m_size + (j + 1)]);

			const auto sum_vec = _mm256_add_ps(_mm256_add_ps(_mm256_add_ps(array_up_vec, array_down_vec), array_left_vec), array_right_vec);
			const auto prod_vec = _mm256_mul_ps(a_vec, sum_vec);
			const auto sum2_vec = _mm256_add_ps(prevArray_vec, prod_vec);

			const auto result_vec = _mm256_mul_ps(cRecip_vec, sum2_vec);
			_mm256_storeu_ps(&array[(i)*m_size + (j)], result_vec);
		}
	}
}

void FluidSim::advectAvx(const int b, std::vector<float>& d, const std::vector<float>& prevD, const std::vector<float>& velocityX, const std::vector<float>& velocityY, const float dt)
{
	// Initialize variables used for spatial interpolation
	const auto dtx_vec = _mm256_set1_ps(dt * (m_size - 2));
	const auto dty_vec = _mm256_set1_ps(dt * (m_size - 2));

	const auto one_vec = _mm256_set1_epi32(1);
	const auto mSize_vec = _mm256_set1_epi32(m_size);
	const auto one_vecps = _mm256_set1_ps(1.0f);
	const auto half_vecps = _mm256_set1_ps(0.5f);

	// Compute the size of the fluid grid
	const float Nfloat = m_size - 2;
	const auto lowerBound_vec = _mm256_set1_ps(0.5);
	const auto higherBound_vec = _mm256_set1_ps(Nfloat + 0.5);
	const auto indexAdditor_vec = _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);

	// Iterate over each grid cell and advect its fluid property
	static int s = m_size - 1;
	for (unsigned int i = 1; i < s; i++)
	{
		for (unsigned int j = 1; j < m_size - 1; j += 8) // Process 8 elements at a time due to avx2 taking 8 floats per operation
		{
			const auto velx_vec = _mm256_loadu_ps(&velocityX[(i)*m_size + (j)]);
			const auto vely_vec = _mm256_loadu_ps(&velocityY[(i)*m_size + (j)]);
			// Compute the position at which to interpolate the density field
			const auto tmp1_vec = _mm256_mul_ps(dtx_vec, velx_vec);
			const auto tmp2_vec = _mm256_mul_ps(dty_vec, vely_vec);

			const auto i_vec = _mm256_set1_ps(i);
			const auto j_vec = _mm256_cvtepi32_ps(_mm256_add_epi32(indexAdditor_vec, _mm256_set1_epi32(j)));
			const auto x_vec = _mm256_fnmadd_ps(dtx_vec, velx_vec, i_vec); /*-(dtx*velx)+i*/
			const auto y_vec = _mm256_fnmadd_ps(dty_vec, vely_vec, j_vec);/*-(dty*vely)+j*/


			// Clamp the position within the fluid grid
			const auto clampedX_vec = _mm256_min_ps(_mm256_max_ps(x_vec, lowerBound_vec), higherBound_vec);

			const auto i0_vec = _mm256_cvttps_epi32(clampedX_vec);
			const auto i1_vec = _mm256_add_epi32(i0_vec, one_vec);

			const auto clampedY_vec = _mm256_min_ps(_mm256_max_ps(y_vec, lowerBound_vec), higherBound_vec);
			const auto j0_vec = _mm256_cvttps_epi32(clampedY_vec);
			const auto j1_vec = _mm256_add_epi32(j0_vec, one_vec);

			const auto i0_vecps = _mm256_cvtepi32_ps(i0_vec);
			const auto i1_vecps = _mm256_cvtepi32_ps(i1_vec);
			const auto j0_vecps = _mm256_cvtepi32_ps(j0_vec);
			const auto j1_vecps = _mm256_cvtepi32_ps(j1_vec);
			// Compute the interpolation coefficients

			const auto s1_vec = _mm256_sub_ps(clampedX_vec, i0_vecps);
			const auto s0_vec = _mm256_sub_ps(one_vecps, s1_vec);
			const auto t1_vec = _mm256_sub_ps(clampedY_vec, j0_vecps);
			const auto t0_vec = _mm256_sub_ps(one_vecps, t1_vec);


			// Perform bilinear interpolation of the density field at the interpolated position
			// and store the result in the current density field

			const auto i0i_vec = _mm256_cvttps_epi32(_mm256_add_ps(i0_vecps, half_vecps));
			const auto i1i_vec = _mm256_cvttps_epi32(_mm256_add_ps(i1_vecps, half_vecps));
			const auto j0i_vec = _mm256_cvttps_epi32(_mm256_add_ps(j0_vecps, half_vecps));
			const auto j1i_vec = _mm256_cvttps_epi32(_mm256_add_ps(j1_vecps, half_vecps));

			//no fnmadd for epi32 :(
			const auto prevDFirstIndex_vex = _mm256_add_epi32(_mm256_mullo_epi32(i0i_vec, mSize_vec), j0i_vec);
			const auto prevDSecondIndex_vex = _mm256_add_epi32(_mm256_mullo_epi32(i0i_vec, mSize_vec), j1i_vec);
			const auto prevDThirdIndex_vex = _mm256_add_epi32(_mm256_mullo_epi32(i1i_vec, mSize_vec), j0i_vec);
			const auto prevDFourthIndex_vex = _mm256_add_epi32(_mm256_mullo_epi32(i1i_vec, mSize_vec), j1i_vec);

			const auto prevDFirst_vec = _mm256_i32gather_ps(&prevD[0], prevDFirstIndex_vex, 4);
			const auto prevDSecond_vec = _mm256_i32gather_ps(&prevD[0], prevDSecondIndex_vex, 4);
			const auto prevDThird_vec = _mm256_i32gather_ps(&prevD[0], prevDThirdIndex_vex, 4);
			const auto prevDFourth_vec = _mm256_i32gather_ps(&prevD[0], prevDFourthIndex_vex, 4);


			const auto part1 = _mm256_fmadd_ps(t0_vec, prevDFirst_vec, _mm256_mul_ps(t1_vec, prevDSecond_vec));
			const auto part2 = _mm256_fmadd_ps(t0_vec, prevDThird_vec, _mm256_mul_ps(t1_vec, prevDFourth_vec));
			const auto result_vec = _mm256_fmadd_ps(s0_vec, part1, _mm256_mul_ps(s1_vec, part2));
			_mm256_storeu_ps(&d[(i)*m_size + (j)], result_vec);

		}
	}
}

void FluidSim::projectAvx(std::vector<float>& velocityX, std::vector<float>& velocityY, std::vector<float>& clearVector, std::vector<float>& targetVectory)
{
	/*const auto half_vecps = _mm256_set1_ps(-0.5f);
	const auto size_vecps = _mm256_set1_ps(m_size);*/
	for (unsigned int j = 1; j < m_size - 1; j++)
	{
		for (unsigned int i = 1; i < m_size - 1; i+=1)
		{
			targetVectory[(j)*m_size + (i)] = -0.5 * (velocityX[(j + 1) * m_size + (i)] - velocityX[(j - 1) * m_size + (i)] + velocityY[(j)*m_size + (i + 1)] - velocityY[(j)*m_size + (i - 1)]) / m_size;
			//you cant use avx here, as each update is reliant on the surrounding pixels, so doing it in groups of 8 breaks it :(
			// Idea: calculate it 8 at a time, but never 2 adjacent ones in the same calculation?
			// Not needed that much, as most time is spent elsewhere
			/*auto first_vec = _mm256_loadu_ps(&velocityX[(j + 1) * m_size + (i)]);
			auto second_vec = _mm256_loadu_ps(&velocityX[(j - 1) * m_size + (i)]);
			auto third_vec = _mm256_loadu_ps(&velocityY[(j)*m_size + (i + 1)]);
			auto four_vec = _mm256_loadu_ps(&velocityY[(j)*m_size + (i - 1)]);
			auto result_vec = _mm256_div_ps(_mm256_mul_ps(half_vecps, _mm256_add_ps(_mm256_sub_ps(first_vec, second_vec), _mm256_sub_ps(third_vec, four_vec))), size_vecps);
			_mm256_storeu_ps(&targetVectory[(j)*m_size + (i)], result_vec);*/
		}
	}
	std::fill(clearVector.begin(), clearVector.end(), 0);
	this->setBoundary(0, targetVectory);
	this->setBoundary(0, clearVector);
	this->linSolve(0, clearVector, targetVectory, 1, 6);

	for (unsigned int i = 1; i < m_size - 1; i++) {
		for (unsigned int j = 1; j < m_size - 1; j += 8) {
			const auto one_vec = _mm256_loadu_ps(&clearVector[(i + 1) * m_size + (j)]);
			const auto two_vec = _mm256_loadu_ps(&clearVector[(i - 1) * m_size + (j)]);
			const auto three_vec = _mm256_loadu_ps(&clearVector[(i)*m_size + (j + 1)]);
			const auto four_vec = _mm256_loadu_ps(&clearVector[(i)*m_size + (j - 1)]);
			const auto velX_vec = _mm256_loadu_ps(&velocityX[(i)*m_size + (j)]);
			const auto vely_vec = _mm256_loadu_ps(&velocityY[(i)*m_size + (j)]);
			const auto factor_vec= _mm256_set1_ps(0.5 * m_size);

			const auto resultOne_vec= _mm256_sub_ps(velX_vec, _mm256_mul_ps(_mm256_sub_ps(one_vec, two_vec), factor_vec));
			const auto resultTwo_vec = _mm256_sub_ps(vely_vec, _mm256_mul_ps(_mm256_sub_ps(three_vec, four_vec), factor_vec));

			_mm256_storeu_ps(&velocityX[(i)*m_size + (j)], resultOne_vec);
			_mm256_storeu_ps(&velocityY[(i)*m_size + (j)], resultTwo_vec);
		}
	}
}

#endif