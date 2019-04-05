cmake_minimum_required(VERSION 2.8)
project(ars)
add_definitions(-std=c++0x)

set(CMAKE_BUILD_TYPE Release)

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

# Library files
add_library(ars
  src/ars2d.cpp
  src/functions.cpp
  src/BBOptimizer1d.cpp
  src/HistogramCircularCorrelation.cpp
  src/HoughSpectrum.cpp
  src/Profiler.cpp
)

add_executable(testFunctions test/testFunctions.cpp )
target_link_libraries(testFunctions ${ars_LINK_LIBS} ars)

add_executable(testIsotropicARS test/testIsotropicARS.cpp )
target_link_libraries(testIsotropicARS ${ars_LINK_LIBS} ars)

add_executable(testAnisotropicArs test/testAnisotropicArs.cpp)
target_link_libraries(testAnisotropicArs ${ars_LINK_LIBS} ars)

add_executable(testCorrelationArs test/testCorrelationArs.cpp)
target_link_libraries(testCorrelationArs ${ars_LINK_LIBS} ars)

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
