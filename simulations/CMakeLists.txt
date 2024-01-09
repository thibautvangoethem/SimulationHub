# CMakeList.txt : CMake project for Fluidsim, include source and define
# project specific logic here.
#
SET(USEAVX true)

add_subdirectory("core")
add_subdirectory("sfmlFrontend")
# Add source to this project's executable.
SET(SRC "main.cpp")


add_executable(simulationHub ${SRC})

TARGET_LINK_LIBRARIES(simulationHub LINK_PUBLIC BASE)
TARGET_LINK_LIBRARIES(simulationHub LINK_PUBLIC FIRESIM)
TARGET_LINK_LIBRARIES(simulationHub LINK_PUBLIC FLUIDSIM)
TARGET_LINK_LIBRARIES(simulationHub LINK_PUBLIC RAYTRACER)
TARGET_LINK_LIBRARIES(simulationHub LINK_PUBLIC SLIMESIM)
TARGET_LINK_LIBRARIES(simulationHub LINK_PUBLIC CUDAINTEROP)

TARGET_LINK_LIBRARIES(simulationHub LINK_PUBLIC SFMLFRONTEND)

find_package(SFML COMPONENTS system window graphics CONFIG REQUIRED)

target_link_libraries(simulationHub PRIVATE sfml-system sfml-network sfml-graphics sfml-window)


if(USEAVX)
	message("AVX2 will be used, linking in executable")
	SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /arch:AVX2")
	SET_TARGET_PROPERTIES(simulationHub PROPERTIES COMPILE_FLAGS "/arch:AVX2")
	target_compile_definitions(simulationHub  PUBLIC AVX)
else()  
	message("AVX2 not used")
endif()