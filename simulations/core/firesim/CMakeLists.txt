set(LIB FIRESIM)

SET(SRC
      "FireSim.cpp"
      "FireSim.h"
      "RandomNumberQueue.h"
      "RandomNumberQueue.cpp")


ADD_LIBRARY(${LIB} STATIC ${SRC})

TARGET_INCLUDE_DIRECTORIES(${LIB} PUBLIC include)

