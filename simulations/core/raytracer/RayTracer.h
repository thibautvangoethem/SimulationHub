#pragma once
#include "../baseinterface/Simulation.h"
#include "../utility/vec3Util.h"
#include "utility/Ray.h"
#include "hittableshapes/HittableShape.h"
#include "hittableshapes/HittableList.h"
#include "hittableshapes/Sphere.h"
#include "Camera.h"

#include <utility> 
namespace SIM {
	class SimulationSetting;
}


namespace RT {
	//This raytracer was made following the tutorial on https://raytracing.github.io/books/RayTracingInOneWeekend.html
	class RayTracer :public SIM::Simulation {
	public:
		RayTracer(std::shared_ptr<SIM::SimulationSettings> settings);
		void advance(const double timestep) final;
		void handleClick(const bool isLeftClick, const int xpos, const int ypos) final;
	protected:
		Color rayHit(const Ray& ray,const HittableShape& world, int depth) const;
		std::vector<std::vector<SIM::colour>> m_currentState;
		int m_size;
		hittableList m_world;
		int m_samplesPerPixel;
		std::unique_ptr<Camera> m_camera;
		int m_maxDepth = 20;
		bool m_done = false;
	};
}