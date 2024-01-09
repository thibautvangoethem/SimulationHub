SET(USEAVX true)
set(LIB FLUIDSIM)

SET(SRC
    "FluidSim.h"
    "FluidSim.cpp"
     )


ADD_LIBRARY(${LIB} STATIC  ${SRC})

TARGET_INCLUDE_DIRECTORIES(${LIB} PUBLIC include)

if(USEAVX)
	target_compile_definitions(${LIB} PUBLIC AVX)
endif()