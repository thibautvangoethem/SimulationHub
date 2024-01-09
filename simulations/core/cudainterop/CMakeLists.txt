cmake_minimum_required (VERSION 3.8)
project ("SimulationHub")
set(LIB CUDAINTEROP)

SET(SRC
      "CudaInterop.cpp"
      "CudaInterop.h")


ADD_LIBRARY(${LIB} STATIC ${SRC})

TARGET_INCLUDE_DIRECTORIES(${LIB} PUBLIC include)