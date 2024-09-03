cmake_minimum_required(VERSION 3.10)
project(ars)
#add_definitions(-std=c++0x)
set(CMAKE_CXX_FLAGS='-std=c++17')
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_BUILD_TYPE Debug)

# Set the output directory for the build executables and libraries
SET(ars_RUNTIME_OUTPUT_DIRECTORY ${ars_SOURCE_DIR}/bin CACHE PATH "Target for the binaries")
SET(ars_LIBRARY_OUTPUT_DIRECTORY ${ars_SOURCE_DIR}/lib CACHE PATH "Target for the libraries")
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${ars_LIBRARY_OUTPUT_DIRECTORY})
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${ars_LIBRARY_OUTPUT_DIRECTORY})
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${ars_RUNTIME_OUTPUT_DIRECTORY})
SET(ars_ENABLE_PCL 0)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")

include_directories(include)
include_directories(src)

list(APPEND ars_LINK_LIBS boost_iostreams boost_system boost_filesystem)

# Eigen library v. 3
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# rofl
find_package(rofl REQUIRED)
include_directories(${rofl_INCLUDE_DIR})

if(${rofl_FOUND})
  list(APPEND ${rofl_LIBRARIES})
endif()

# Library files
add_library(ars
  src/ars2d.cpp
  src/functions.cpp
  src/utils.cpp
  src/ArsKernelAnisotropic2d.cpp
  src/ArsKernelIsotropic2d.cpp
  src/BBOptimizer1d.cpp
  src/BBTranslation.cpp
  src/ConsensusTranslationEstimator.cpp
  src/DisjointSet.cpp
  src/GaussianMixtureEstimator.cpp
  src/HistogramCircularCorrelation.cpp
  src/HoughSpectrum.cpp

  # src/MortonOrderedPoints2d.cpp
  src/MortonSort.cpp
  src/ParamMap.cpp
  src/Profiler.cpp
  src/ProcrustesUmeyama.cpp
  src/TranslationRefiner.cpp
)

add_executable(testFunctions app/testFunctions.cpp)
target_link_libraries(testFunctions ${ars_LINK_LIBS} ars)

add_executable(testIsotropicArs app/testIsotropicArs.cpp)
target_link_libraries(testIsotropicArs ${ars_LINK_LIBS} ars)

add_executable(testAnisotropicArs app/testAnisotropicArs.cpp)
target_link_libraries(testAnisotropicArs ${ars_LINK_LIBS} ars)

add_executable(testCorrelationArs app/testCorrelationArs.cpp)
target_link_libraries(testCorrelationArs ${ars_LINK_LIBS} ars)

add_executable(testAnisotropicKernel app/testAnisotropicKernel.cpp)
target_link_libraries(testAnisotropicKernel ${ars_LINK_LIBS} ars)

add_executable(testGaussianMixtureEstimator app/testGaussianMixtureEstimator.cpp)
target_link_libraries(testGaussianMixtureEstimator ${ars_LINK_LIBS} ars)

add_executable(testConsensusTranslationEstimator app/testConsensusTranslationEstimator.cpp)
target_link_libraries(testConsensusTranslationEstimator ${ars_LINK_LIBS} ars)

add_executable(testProcrustes app/testProcrustes.cpp)
target_link_libraries(testProcrustes ${ars_LINK_LIBS} ars)

add_executable(testProcrustesUmeyama2d app/testProcrustesUmeyama2d.cpp)
target_link_libraries(testProcrustesUmeyama2d ${ars_LINK_LIBS} ars)

add_executable(testProcrustesUmeyama3d app/testProcrustesUmeyama3d.cpp)
target_link_libraries(testProcrustesUmeyama3d ${ars_LINK_LIBS} ars)

add_executable(testBBTranslation app/testBBTranslation.cpp)
target_link_libraries(testBBTranslation ${ars_LINK_LIBS} ars)

add_executable(testTranslationRefiner app/testTranslationRefiner.cpp)
target_link_libraries(testTranslationRefiner ${ars_LINK_LIBS} ars)

# add_executable(testMortonOrderedPoints app/testMortonOrderedPoints.cpp)
# target_link_libraries(testMortonOrderedPoints ${ars_LINK_LIBS} ars)

# Option "make install": copy binaries
INSTALL(TARGETS ars
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)

# Option "make install": copy headers
FILE(GLOB headers "${CMAKE_CURRENT_SOURCE_DIR}/include/ars/*.h")
FILE(GLOB headers_3rdparty "${CMAKE_CURRENT_SOURCE_DIR}/include/ars/thirdparty/*.h")
INSTALL(FILES ${headers} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/ars/)
INSTALL(FILES ${headers_3rdparty} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/ars/thirdparty/)

# Option "make install": copy cmake script
FILE(GLOB cmake_script "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/arsConfig.cmake")
message(STATUS "cmake_script " ${cmake_script})
INSTALL(FILES ${cmake_script} DESTINATION ${CMAKE_INSTALL_PREFIX}/share/ars/)

# uninstall target
if(NOT TARGET uninstall)
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_modules/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

  add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_modules/cmake_uninstall.cmake)
endif()
