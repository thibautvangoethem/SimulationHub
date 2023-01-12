# CMakeList.txt : CMake project for Fluidsim, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)
project ("SimulationHub")

add_subdirectory ("baseinterface")

add_subdirectory ("firesim")

add_subdirectory ("fluidsim")
# Add source to this project's executable.
SET(SRC "main.cpp")


add_executable (simulationHub ${SRC})

TARGET_LINK_LIBRARIES(simulationHub LINK_PUBLIC BASE FIRESIM)
TARGET_LINK_LIBRARIES(simulationHub LINK_PUBLIC BASE FLUIDSIM)

find_package(SFML COMPONENTS system window graphics CONFIG REQUIRED)
target_link_libraries(simulationHub PRIVATE sfml-system sfml-network sfml-graphics sfml-window)

# TODO: Add tests and install targets if needed.
