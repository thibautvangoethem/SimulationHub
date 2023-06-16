#include "RayTracer.h"

#include <iostream>
#include <algorithm>
#include <chrono>
#include <random>

#include "materials/Lambertian.h"
#include "materials/Metal.h"
#include "materials/Dielectric.h"

using namespace RT;

RayTracer::RayTracer(std::shared_ptr<SIM::SimulationSettings> settings) : SIM::Simulation(std::move(settings))
{
    m_size = this->m_settings->getSize();
    m_currentState = std::vector<std::vector<SIM::colour> >(
        m_size,
        std::vector<SIM::colour>(m_size));

    Point3 lookfrom(13, 2, 3);
    Point3 lookat(0, 0, 0);
    Vec3 vup(0, 1, 0);
    auto dist_to_focus = 10.0;
    auto aperture = 0.1;

    m_camera = std::make_unique<Camera>( lookfrom, lookat, vup, 20, 1, aperture, dist_to_focus);
    //auto groundMat = std::make_shared<Lambertian>(Color{ 0.8,0.8,0.0 });
    //auto centerMat = std::make_shared<Lambertian>(Color{ 0.1,0.2,0.5 });
    ///*auto leftMat = std::make_shared<Metal>(Color{ 0.8,0.8,0.8 },0.3);*/
    ////auto centerMat = std::make_shared<Dielectric>(1.5);
    //auto leftMat = std::make_shared<Dielectric>(1.5);
    //auto rightMat = std::make_shared<Metal>(Color{ 0.8,0.6,0.2 },0.0);
    //
    //m_world.add(std::make_shared<Sphere>(Point3{ 0.0, -100.5, -1.0 }, 100.0, groundMat));
    //m_world.add(std::make_shared<Sphere>(Point3{ 0.0, 0.0, -1.0 }, 0.5, centerMat));
    //m_world.add(std::make_shared<Sphere>(Point3{ -1.0, 0.0, -1.0 }, 0.5, leftMat));
    //m_world.add(std::make_shared<Sphere>(Point3{ -1.0, 0.0, -1.0 }, -0.4, leftMat));
    //m_world.add(std::make_shared<Sphere>(Point3{ 1.0, 0.0, -1.0 }, 0.5, rightMat));

   /* auto R = cos(M_PI / 4);
    auto material_left = std::make_shared<Lambertian>(Color{ 0, 0, 1 });
    auto material_right = std::make_shared<Lambertian>(Color{1, 0, 0});


    m_world.add(std::make_shared<Sphere>(Point3(-R, 0, -1), R, material_left));
    m_world.add(std::make_shared<Sphere>(Point3(R, 0, -1), R, material_right));*/


   /* auto material_ground = std::make_shared<Lambertian>(Color(0.8, 0.8, 0.0));
    auto material_center = std::make_shared<Lambertian>(Color(0.1, 0.2, 0.5));
    auto material_left = std::make_shared<Dielectric>(1.5);
    auto material_right = std::make_shared<Metal>(Color(0.8, 0.6, 0.2), 0.0);

    m_world.add(std::make_shared<Sphere>(Point3(0.0, -100.5, -1.0), 100.0, material_ground));
    m_world.add(std::make_shared<Sphere>(Point3(0.0, 0.0, -1.0), 0.5, material_center));
    m_world.add(std::make_shared<Sphere>(Point3(-1.0, 0.0, -1.0), 0.5, material_left));
    m_world.add(std::make_shared<Sphere>(Point3(-1.0, 0.0, -1.0), -0.45, material_left));
    m_world.add(std::make_shared<Sphere>(Point3(1.0, 0.0, -1.0), 0.5, material_right));*/

    auto ground_material = std::make_shared<Lambertian>(Color(0.5, 0.5, 0.5));
    m_world.add(std::make_shared<Sphere>(Point3(0, -1000, 0), 1000, ground_material));

    static std::uniform_real_distribution<double> distribution(0.0, 1.0);
    static std::uniform_real_distribution<double> distrFuzz(0, 0.5);
    static std::mt19937 generator;
    for (int a = -11; a < 11; a++) {
        for (int b = -11; b < 11; b++) {
            auto choose_mat = distribution(generator);
            Point3 center(a + 0.9 * distribution(generator), 0.2, b + 0.9 * distribution(generator));

            if ((center - Point3(4, 0.2, 0)).norm() > 0.9) {
                std::shared_ptr<Material> sphere_material;

                if (choose_mat < 0.8) {
                    // diffuse
                    auto albedo = randomVector(0,1).cwiseProduct(randomVector(0,1));
                    sphere_material = std::make_shared<Lambertian>(albedo);
                    m_world.add(std::make_shared<Sphere>(center, 0.2, sphere_material));
                }
                else if (choose_mat < 0.95) {
                    // metal
                    auto albedo = randomVector(0.5, 1);
                    auto fuzz = distrFuzz(generator);
                    sphere_material = std::make_shared<Metal>(albedo, fuzz);
                    m_world.add(std::make_shared<Sphere>(center, 0.2, sphere_material));
                }
                else {
                    // glass
                    sphere_material = std::make_shared<Dielectric>(1.5);
                    m_world.add(std::make_shared<Sphere>(center, 0.2, sphere_material));
                }
            }
        }
    }

    auto material1 = std::make_shared<Dielectric>(1.5);
    m_world.add(std::make_shared<Sphere>(Point3(0, 1, 0), 1.0, material1));

    auto material2 = std::make_shared<Lambertian>(Color(0.4, 0.2, 0.1));
    m_world.add(std::make_shared<Sphere>(Point3(-4, 1, 0), 1.0, material2));

    auto material3 = std::make_shared<Metal>(Color(0.7, 0.6, 0.5), 0.0);
    m_world.add(std::make_shared<Sphere>(Point3(4, 1, 0), 1.0, material3));
    
    m_samplesPerPixel = 50;
}

void RayTracer::advance(const double timestep)
{
    if (m_done) return;
    static std::uniform_real_distribution<double> distribution(0.0, 1.0);
    static std::mt19937 generator;
    for (auto j = 0; j < m_size; ++j) {
        for (auto  i = 0; i < m_size; ++i) {
            Color pixelCol(0, 0, 0);
            for (int s = 0; s < m_samplesPerPixel; ++s) {
                auto u = (i + distribution(generator)) / (m_size- 1);
                auto v = (j + distribution(generator)) / (m_size - 1);
                Ray r = m_camera->getRay(u, v);
                pixelCol += rayHit(r, m_world,m_maxDepth);
            }

            auto scale = 1.0 / m_samplesPerPixel;
            auto r = std::sqrt(pixelCol.x() * scale);
            auto g = std::sqrt(pixelCol.y() * scale);
            auto b = std::sqrt(pixelCol.z() * scale);

            auto ir = static_cast<uint8_t>(std::clamp(r,0.0,0.999)*256);
            auto ig = static_cast<uint8_t>(std::clamp(g , 0.0, 0.999) * 256);
            auto ib = static_cast<uint8_t>(std::clamp(b, 0.0, 0.999) * 256);

            m_currentState[m_size - j - 1][i] = SIM::colour{ ir,ig,ib};
        }
    }
    m_done = true;
}

void RayTracer::handleClick(const bool isLeftClick, const int xpos, const int ypos)
{
	//todo
}

Color RayTracer::rayHit(const Ray& ray,const HittableShape& world, int depth) const
{
    // Gaurd against infinite recursion by simply cutting of the recursion
    if (depth <= 0)
        return Color(0, 0, 0);
    HitRecord record;
    if (world.hit(ray, 0.001, infinity, record)) {
        Ray scattered;
        Color attenuation;
        if(record.material->scatter(ray,record,attenuation,scattered))
            return attenuation.cwiseProduct(rayHit(scattered, world, depth - 1));
        return Color{ 0,0,0 };
    }
    Vec3 unitDirection = unit_vector(ray.direction());
    auto t = 0.5 * (unitDirection.y() + 1.0);
    return (1.0 - t) * Color { 1.0, 1.0, 1.0 } + t * Color{ 0.5, 0.7, 1.0 };
}

