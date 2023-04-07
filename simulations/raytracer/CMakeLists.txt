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
    "vec3.h"
    "vec3.cpp"
    "vec3Util.h"
     )


ADD_LIBRARY(${LIB}  ${SRC})

TARGET_INCLUDE_DIRECTORIES(${LIB} PUBLIC include)

#target_link_libraries(main PRIVATE sfml-main)

# TODO: Add tests and install targets if needed
if(USEAVX)
	target_compile_definitions(${LIB} PUBLIC AVX)
endif()