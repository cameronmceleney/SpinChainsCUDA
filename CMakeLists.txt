cmake_minimum_required(VERSION 3.25 FATAL_ERROR)
project(SpinChainsCUDA LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CUDA_STANDARD 23)

add_executable(SpinChainsCUDA "main.cu"
               "example_kernel.cu"
               "Numerical_Methods_Class.cu" "Numerical_Methods_Class.cuh"
               "SpinChainEigenSolverClass.cpp" "SpinChainEigenSolverClass.h"
               "GlobalVariables.cpp" "GlobalVariables.h"
               "CommonLibs.h"
               "linspace.cpp" "linspace.h"
               "progressbar.hpp"
               "Numerical_Methods_Class_backup.h" "Numerical_Methods_Class_backup.cpp")

set_target_properties(SpinChainsCUDA
        PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)


if(APPLE)
    # Help the static cuda runtime find the driver (libcuda.dyllib) at runtime.
    set_property(TARGET <target_name> PROPERTY BUILD_RPATH ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})
endif()
