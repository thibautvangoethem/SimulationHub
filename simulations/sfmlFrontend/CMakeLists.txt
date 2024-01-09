set(LIB SFMLFRONTEND)

SET(SRC
      "SFMLSimInterface.h"
      "FireSimImpl.h" 
      "FireSimImpl.cpp"
      "FluidSimImpl.h"
      "FluidSimImpl.cpp"
      "RayTracerImpl.h"
      "RayTracerImpl.cpp"
      "SlimeSimImpl.h"
      "SlimeSimImpl.cpp")


ADD_LIBRARY(${LIB} STATIC  ${SRC})

TARGET_LINK_LIBRARIES(${LIB} LINK_PUBLIC BASE FIRESIM)
TARGET_LINK_LIBRARIES(${LIB} LINK_PUBLIC BASE FLUIDSIM)
TARGET_LINK_LIBRARIES(${LIB} LINK_PUBLIC BASE RAYTRACER)
TARGET_LINK_LIBRARIES(${LIB} LINK_PUBLIC BASE SLIMESIM)

find_package(SFML COMPONENTS system window graphics CONFIG REQUIRED)
target_link_libraries(${LIB} PRIVATE sfml-system sfml-network sfml-graphics sfml-window)
