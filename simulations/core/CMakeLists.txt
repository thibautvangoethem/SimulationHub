# CMakeList.txt : CMake project for Fluidsim, include source and define
# project specific logic here.
#	
SET(USEAVX true)
cmake_minimum_required (VERSION 3.8)
project ("SimulationHub")
set(LIB CORE)

add_subdirectory("baseinterface")
add_subdirectory("firesim")
add_subdirectory("fluidsim")
add_subdirectory("raytracer")
add_subdirectory("slimesim")

SET(SRC
    "_.cpp")


ADD_LIBRARY(${LIB}  ${SRC} )

TARGET_INCLUDE_DIRECTORIES(${LIB} PUBLIC include)

