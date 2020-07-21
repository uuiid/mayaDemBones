#-
# ==========================================================================
# Copyright (c) 2018 Autodesk, Inc.
# All rights reserved.
# 
# These coded instructions, statements, and computer programs contain
# unpublished proprietary information written by Autodesk, Inc., and are
# protected by Federal copyright law. They may not be disclosed to third
# parties or copied or duplicated in any form, in whole or in part, without
# the prior written consent of Autodesk, Inc.
# ==========================================================================
#+


cmake_minimum_required(VERSION 2.8)
project("doodle_bone")
# include the project setting file
include($ENV{DEVKIT_LOCATION}/cmake/pluginEntry.cmake)
include_directories("${CMAKE_CURRENT_SOURCE_DIR}")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")
# specify project name

find_package(Eigen3)
find_path(EIGEN_INCLUDE_DIR
    NAMES 
        "Eigen/Dense"
        "Eigen/Sparse"
        "Eigen/StdVector"
)
include_directories("${EIGEN_INCLUDE_DIR}")
# specify project name
set(PROJECT_NAME DoodleConvertBone)



# set SOURCE_FILES
set(SOURCE_FILES
   DoodleConvertBone.cpp
   DoodleConvertBone.h
)

# set linking libraries
set(LIBRARIES
     OpenMaya
     Foundation
     OpenMayaAnim
)





# Build plugin
build_plugin()

install(TARGETS DoodleConvertBone
    DESTINATION "C:/Users/teXiao/Documents/maya/plug-ins")