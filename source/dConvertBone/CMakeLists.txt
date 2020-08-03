# ËÑË÷eng¿â
find_package(Eigen3)
find_path(EIGEN_INCLUDE_DIR
    NAMES 
        "Eigen/Dense"
        "Eigen/Sparse"
        "Eigen/StdVector"
)
include_directories("${EIGEN_INCLUDE_DIR}")

include($ENV{DEVKIT_LOCATION}/cmake/pluginEntry.cmake)
add_library(dConvertBone doodleConvert.cpp doodleConvert.h)