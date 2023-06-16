# CMakeList.txt : CMake project for Fluidsim, include source and define
# project specific logic here.
#	
SET(USEAVX true)
cmake_minimum_required (VERSION 3.8)
project ("SimulationHub")
set(LIB RAYTRACER)

# Add source to this project's executable.
SET(SRC
    "RayTracer.h"
    "RayTracer.cpp"
    "utility/vec3.h"
    "utility/vec3.cpp"
    "utility/vec3Util.h"
    "utility/ray.h" 
    "hittableshapes/HittableShape.h"
    "hittableshapes/Sphere.h" "hittableshapes/Sphere.cpp" "hittableshapes/HittableList.h" "hittableshapes/HittableList.cpp"
"Camera.h" "Camera.cpp" "materials/Material.h" "materials/Lambertian.cpp" "materials/Lambertian.h" "materials/Metal.h" "materials/Metal.cpp" "materials/Dielectric.h" "materials/Dielectric.cpp")


ADD_LIBRARY(${LIB}  ${SRC})

TARGET_INCLUDE_DIRECTORIES(${LIB} PUBLIC include)

#target_link_libraries(main PRIVATE sfml-main)

# TODO: Add tests and install targets if needed
if(USEAVX)
	target_compile_definitions(${LIB} PUBLIC AVX)
endif()