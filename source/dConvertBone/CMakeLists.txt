project(DoodleConvertBone)


# ËÑË÷eng¿â
find_package(Eigen3)
find_path(EIGEN_INCLUDE_DIR
    NAMES 
        "Eigen/Dense"
        "Eigen/Sparse"
        "Eigen/StdVector"
)
include_directories("${EIGEN_INCLUDE_DIR}")
# ËÑË÷openmp¿â
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# ¸ü¸Ä±àÒëÉèÖÃ
if(MSVC)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /ZI")
set(LINK_FLAGS "${LINK_FLAGS} /force:multiple")
endif()
#string(REGEX REPLACE "/ZI" ""
#         CMAKE_CXX_FLAGS_RELEASE
#         "${CMAKE_CXX_FLAGS_RELEASE}")


include($ENV{DEVKIT_LOCATION}/cmake/pluginEntry.cmake)

include($ENV{DEVKIT_LOCATION}/cmake/pluginEntry.cmake)

# set SOURCE_FILES
set(SOURCE_FILES
   doodleConvert.cpp
   doodleConvert.h
   doodleConvertPlug.cpp
   doodleConvertPlug.h
   doodleWeight.cpp
   doodleWeight.h
   plugMain.cpp
)

# set linking libraries
set(LIBRARIES
     OpenMaya
     Foundation
     OpenMayaAnim
)
# Build plugin
build_plugin()
# add_library(dConvertBone doodleConvert.cpp doodleConvert.h)
install(TARGETS DoodleConvertBone
    DESTINATION "C:/Users/teXiao/Documents/maya/plug-ins")