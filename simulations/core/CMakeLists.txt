SET(USEAVX true)
set(LIB CORE)

add_subdirectory("baseinterface")
add_subdirectory("firesim")
add_subdirectory("fluidsim")
add_subdirectory("raytracer")
add_subdirectory("slimesim")

SET(SRC
    "_.cpp")


ADD_LIBRARY(${LIB} STATIC ${SRC} )

TARGET_INCLUDE_DIRECTORIES(${LIB} PUBLIC include)

