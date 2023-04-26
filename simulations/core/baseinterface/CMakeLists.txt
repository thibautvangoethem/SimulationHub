# CMakeList.txt : CMake project for Fluidsim, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)
project ("SimulationHub")


set(LIB BASE)

# Add source to this project's executable.
SET(SRC 
	"Simulation.h"
	"SimulationSettings.h" 
	"_.cpp")


ADD_LIBRARY(${LIB} ${SRC})

TARGET_INCLUDE_DIRECTORIES(${LIB} PUBLIC include)

#this needs to be done as there are no cpp files in this directory
set_target_properties(BASE PROPERTIES LINKER_LANGUAGE CXX)