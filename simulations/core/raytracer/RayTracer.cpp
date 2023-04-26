#include "RayTracer.h"

#include <iostream>
#include <algorithm>
#include <chrono>
#include <random>
using namespace RT;


RayTracer::RayTracer(std::shared_ptr<SIM::SimulationSettings> settings) : SIM::Simulation(std::move(settings))
{
    m_size = this->m_settings->getSize();
	m_currentState = std::vector<std::vector<SIM::colour> >(
		m_size,
		std::vector<SIM::colour>(m_size));

    
    m_world.add(std::make_shared<Sphere>(Point3(0, 0, -1), 0.5));
    m_world.add(std::make_shared<Sphere>(Point3(0, -100.5, -1), 100));
    m_samplesPerPixel = 50;
}

void RayTracer::advance(const double timestep)
{
    static std::uniform_real_distribution<double> distribution(0.0, 1.0);
    static std::mt19937 generator;
    for (int j = 0; j < m_size; ++j) {
        for (int i = 0; i < m_size; ++i) {
            Color pixelCol(0, 0, 0);
            for (int s = 0; s < m_samplesPerPixel; ++s) {
                auto u = (i + distribution(generator)) / (m_size- 1);
                auto v = (j + distribution(generator)) / (m_size - 1);
                Ray r = m_camera.getRay(u, v);
                pixelCol += rayHit(r, m_world);
            }

            auto scale = 1.0 / m_samplesPerPixel;

            auto ir = static_cast<uint8_t>(std::clamp(pixelCol.x()*scale,0.0,0.999)*256);
            auto ig = static_cast<uint8_t>(std::clamp(pixelCol.y() * scale , 0.0, 0.999) * 256);
            auto ib = static_cast<uint8_t>(std::clamp(pixelCol.z()*scale, 0.0, 0.999) * 256);

            m_currentState[m_size - j - 1][i] = SIM::colour{ ir,ig,ib};
        }
    }
}
//
//std::vector<std::vector<SIM::colour>>& RayTracer::getCurrentState()
//{
//    return m_currentState;
//}

void RayTracer::handleClick(const bool isLeftClick, const int xpos, const int ypos)
{
	//todo
}

Color RayTracer::rayHit(const Ray& ray,const HittableShape& world) const
{
    HitRecord record;
    if (world.hit(ray, 0, infinity, record)) {
        return 0.5 * (record.normal + Color(1, 1, 1));
    }
    Vec3 unitDirection = unit_vector(ray.direction());
    auto t = 0.5 * (unitDirection.y() + 1.0);
    return (1.0 - t) * Color(1.0, 1.0, 1.0) + t * Color(0.5, 0.7, 1.0);
}

