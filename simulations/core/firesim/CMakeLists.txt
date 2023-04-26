# CMakeList.txt : CMake project for Fluidsim, include source and define
# project specific logic here.
#	
cmake_minimum_required (VERSION 3.8)
project ("SimulationHub")
set(LIB FIRESIM)

# Add source to this project's executable.
SET(SRC
      "FireSim.cpp"
      "FireSim.h"
      "RandomNumberQueue.h"
      "RandomNumberQueue.cpp")


ADD_LIBRARY(${LIB}  ${SRC})

TARGET_INCLUDE_DIRECTORIES(${LIB} PUBLIC include)

#target_link_libraries(main PRIVATE sfml-main)

# TODO: Add tests and install targets if needed.
