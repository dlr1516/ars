cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(ars-cuda LANGUAGES CXX CUDA)
#add_definitions(-std=c++0x)

add_compile_options(-std=c++14)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CUDA_FLAGS_DEBUG "-g -G")
#add_compile_options(-lineinfo)

set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -gencode arch=compute_75,code=sm_75)



include_directories(include)
include_directories(src)

include(CTest)


#list(APPEND ars_LINK_LIBS boost_iostreams boost_system boost_filesystem)
find_library(CUDART_LIBRARY cudart ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})


# Eigen library v. 3
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})


# PCL library 
#message(STATUS "pcl dir: ${PCL_DIR}")
find_package(PCL 1.8 REQUIRED COMPONENTS search KdTree)
if(${PCL_FOUND})
  include_directories(${PCL_INCLUDE_DIRS})
  link_directories(${PCL_LIBRARY_DIRS})
  add_definitions(${PCL_DEFINITIONS})
endif()


find_package(rofl REQUIRED)
if(${rofl_FOUND})
    include_directories(${rofl_INCLUDE_DIRS})
    link_directories(${rofl_LIBRARY_DIRS})
endif()



# Library files
add_library(cuars STATIC
    src/ars2d.cpp
    src/ars2d.cu
    
    src/definitions.cpp
    src/utils.cpp
    src/functions.cpp


    src/GaussianMixtureEstimator.cpp
    src/DisjointSet.cpp
    src/MortonSort.cpp


    src/ArsKernelIsotropic2d.cpp
    src/ArsKernelAnisotropic2d.cpp
    src/BBOptimizer1d.cpp
    src/Profiler.cpp
)

add_library(cuars_mpeg7 STATIC
    src/ars2d.cpp
    src/ars2d.cu

    src/definitions.cpp
    src/utils.cpp
    src/functions.cpp

    src/ArsKernelIsotropic2d.cpp
    src/ArsKernelAnisotropic2d.cpp

    src/BBOptimizer1d.cpp
    src/Profiler.cpp

    src/mpeg7RW.cpp
    src/mpeg7_io.cpp
)



####################################################################################
### Isotropic cuArs

add_executable(gpu_testIsotropicArs app/gpu_testIsotropicArs.cu )
target_link_libraries(gpu_testIsotropicArs ${ars_LINK_LIBS} cuars)

#target_compile_features(gpu_testIsotropicArs PUBLIC cxx_std_11)
set_target_properties(gpu_testIsotropicArs PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
set_target_properties(gpu_testIsotropicArs PROPERTIES CUDA_ARCHITECTURES "52")

####################################################################################
 
add_executable(gpu_testArsIsoMpeg7 app/gpu_testArsIsoMpeg7.cu)
target_link_libraries(gpu_testArsIsoMpeg7 ${ars_LINK_LIBS} cuars_mpeg7 rofl_common ${PCL_LIBRARIES} -lstdc++fs)

target_compile_features(gpu_testArsIsoMpeg7 PUBLIC cxx_std_11)
set_target_properties(gpu_testArsIsoMpeg7 PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
set_target_properties(gpu_testArsIsoMpeg7 PROPERTIES CUDA_ARCHITECTURES "52")

####################################################################################
 
add_executable(gpu_fullArsIsoMpeg7 app/gpu_fullArsIsoMpeg7.cu)
target_link_libraries(gpu_fullArsIsoMpeg7 ${ars_LINK_LIBS} cuars_mpeg7 rofl_common ${PCL_LIBRARIES} -lstdc++fs)

target_compile_features(gpu_fullArsIsoMpeg7 PUBLIC cxx_std_11)
set_target_properties(gpu_fullArsIsoMpeg7 PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
set_target_properties(gpu_fullArsIsoMpeg7 PROPERTIES CUDA_ARCHITECTURES "52")

#####################################################################################
### Anisotropic cuArs

#add_executable(gpu_testAnisotropicArs app/gpu_testAnisotropicArs.cu )
#target_link_libraries(gpu_testAnisotropicArs ${ars_LINK_LIBS} cuars rofl_common -lboost_iostreams -lboost_system -lboost_filesystem)
#
##target_compile_features(gpu_testAnisotropicArs PUBLIC cxx_std_11)
#set_target_properties(gpu_testAnisotropicArs PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
#set_target_properties(gpu_testAnisotropicArs PROPERTIES CUDA_ARCHITECTURES "52")




