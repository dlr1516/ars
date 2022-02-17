cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(ars-cuda LANGUAGES CXX CUDA)
#add_definitions(-std=c++0x)

set(CMAKE_BUILD_TYPE Debug)

set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -gencode arch=compute_75,code=sm_75)

add_definitions(-lineinfo)


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
find_package(PCL 1.8 REQUIRED)
if(${PCL_FOUND})
  include_directories(${PCL_INCLUDE_DIRS})
  link_directories(${PCL_LIBRARY_DIRS})
  add_definitions(${PCL_DEFINITIONS})
  add_definitions(-DENABLE_PCL_REGISTRATION)
endif()


find_package(rofl REQUIRED)
if(${rofl_FOUND})
    include_directories(${rofl_INCLUDE_DIRS})
    link_directories(${rofl_LIBRARY_DIRS})
endif()



# Library files
add_library(cuars_base STATIC
    src/ars2d.cpp
#    src/definitions.cpp #since using double2 and double4, those definitions are sorta obsolete
    src/functions.cpp

    src/ArsKernelIsotropic2d.cpp
    src/BBOptimizer1d.cpp
    src/Profiler.cpp
)

add_library(cuars_mpeg7 STATIC
    src/ars2d.cpp
#    src/definitions.cpp #since using double2 and double4, those definitions are sorta obsolete
    src/functions.cpp
    

    src/ArsKernelIsotropic2d.cpp
    src/BBOptimizer1d.cpp
    src/Profiler.cpp

    src/mpeg7RW.cpp
)

#add_executable(testFunctions app/testFunctions.cpp )
#target_link_libraries(testFunctions ${ars_LINK_LIBS} ars)
#
#add_executable(testIsotropicArs app/testIsotropicArs.cpp )
#target_link_libraries(testIsotropicArs ${ars_LINK_LIBS} ars)



add_executable(gpu_testIsotropicArs app/gpu_testIsotropicArs.cu )
target_link_libraries(gpu_testIsotropicArs ${ars_LINK_LIBS} cuars_base)

target_compile_features(gpu_testIsotropicArs PUBLIC cxx_std_11)
set_target_properties(gpu_testIsotropicArs PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
set_target_properties(gpu_testIsotropicArs PROPERTIES CUDA_ARCHITECTURES "52")

####################################################################################
 
add_executable(gpu_testArsIsoMpeg7 app/gpu_testArsIsoMpeg7.cu)
target_link_libraries(gpu_testArsIsoMpeg7 ${ars_LINK_LIBS} cuars_mpeg7 rofl_common ${PCL_LIBRARIES})

target_compile_features(gpu_testArsIsoMpeg7 PUBLIC cxx_std_11)
set_target_properties(gpu_testArsIsoMpeg7 PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
set_target_properties(gpu_testArsIsoMpeg7 PROPERTIES CUDA_ARCHITECTURES "52")

#add_executable(testAnisotropicArs app/testAnisotropicArs.cpp)
#target_link_libraries(testAnisotropicArs ${ars_LINK_LIBS} ars)
#
#add_executable(testCorrelationArs app/testCorrelationArs.cpp)
#target_link_libraries(testCorrelationArs ${ars_LINK_LIBS} ars)
#
#add_executable(testAnisotropicKernel app/testAnisotropicKernel.cpp)
#target_link_libraries(testAnisotropicKernel ${ars_LINK_LIBS} ars)
#
#add_executable(testGaussianMixtureEstimator app/testGaussianMixtureEstimator.cpp)
#target_link_libraries(testGaussianMixtureEstimator ${ars_LINK_LIBS} ars)

#add_executable(testMortonOrderedPoints app/testMortonOrderedPoints.cpp)
#target_link_libraries(testMortonOrderedPoints ${ars_LINK_LIBS} ars)

# Option "make install": copy binaries
#INSTALL(TARGETS ars 
#  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
#  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
#  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
#)

# Option "make install": copy headers
#FILE(GLOB headers "${CMAKE_CURRENT_SOURCE_DIR}/include/ars/*.h")
#FILE(GLOB headers_3rdparty "${CMAKE_CURRENT_SOURCE_DIR}/include/ars/thirdparty/*.h")
#INSTALL(FILES ${headers} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/ars/)
#INSTALL(FILES ${headers_3rdparty} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/ars/thirdparty/)

# Option "make install": copy cmake script
#FILE(GLOB cmake_script "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/arsConfig.cmake")
#message(STATUS "cmake_script " ${cmake_script})
#INSTALL(FILES ${cmake_script} DESTINATION ${CMAKE_INSTALL_PREFIX}/share/ars/)

