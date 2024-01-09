SET(USEAVX true)
set(LIB RAYTRACER)

# Add source to this project's executable.
SET(SRC
    "RayTracer.h"
    "RayTracer.cpp"
    "utility/vec3Util.h"
    "utility/ray.h" 
    "hittableshapes/HittableShape.h"
    "hittableshapes/Sphere.h" "hittableshapes/Sphere.cpp" "hittableshapes/HittableList.h" "hittableshapes/HittableList.cpp"
"Camera.h" "Camera.cpp" "materials/Material.h" "materials/Lambertian.cpp" "materials/Lambertian.h" "materials/Metal.h" "materials/Metal.cpp" "materials/Dielectric.h" "materials/Dielectric.cpp")


ADD_LIBRARY(${LIB} STATIC ${SRC})

TARGET_INCLUDE_DIRECTORIES(${LIB} PUBLIC include)

find_package(Eigen3 CONFIG REQUIRED)
target_link_libraries(${LIB} PRIVATE Eigen3::Eigen)


#target_link_libraries(main PRIVATE sfml-main)

# TODO: Add tests and install targets if needed
if(USEAVX)
	target_compile_definitions(${LIB} PUBLIC AVX)
endif()