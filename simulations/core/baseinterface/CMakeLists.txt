set(LIB BASE)

SET(SRC 
	"Simulation.h"
	"SimulationSettings.h" 
	"_.cpp")


ADD_LIBRARY(${LIB} STATIC ${SRC})

TARGET_INCLUDE_DIRECTORIES(${LIB} PUBLIC include)

set_target_properties(BASE PROPERTIES LINKER_LANGUAGE CXX)