# ARS - Angular Radon Spectrum
#### Copyright (C) 2017 Dario Lodi Rizzini.


OVERVIEW
-------------------------------------------------

Library **ars** implements the Angular Radon Spectrum method 
for estimation of rotation. 
It has been kept to a minimal design. 

If you use this library, please cite the following paper: 

D. Lodi Rizzini. 
Angular Radon Spectrum for Rotation Estimation. 
Pattern Recognition, Volume 84, Dec. 2018, Pages 182-196, 
DOI [10.1016/j.patcog.2018.07.017][https://doi.org/10.1016/j.patcog.2018.07.017].

or the most relevant associated publications by visiting: 
http://rimlab.ce.unipr.it/


DEPENDENCIES
-------------------------------------------------

The software depends on the following external libraries

- Boost (submodule lexical_cast)
- Eigen 3.0 

Other dependencies are placed in directory thirdparty. 
Some examples require the external application "gnuplot" to display 
results. 


HOW TO COMPILE
-------------------------------------------------

Let ${ars_ROOT} be the install directory of your local copy 
of library ars. 
The following standard commands are required to compile it:

-  cd ${ars_ROOT}
-  mkdir build
-  cd build
-  cmake ..
-  make

You can also install the library into a system directory. 
To change the install directory you must set cmake environment
variable ${CMAKE_INSTALL_PREFIX} (e.g. using command "ccmake .."
before calling "cmake .."). 
Its default value on UNIX-like/Linux systems is "/usr/local".
After compiling library ars, run the command:

-  sudo make install

The command "sudo" is required only if ${CMAKE_INSTALL_PREFIX} 
is a system diretory managed by administrator user root.
Such command copies:
- header files of ${ars_ROOT}/include/ars to
   ${CMAKE_INSTALL_PREFIX}/include/ars/
- library files ${ars_ROOT}/lib/libars.a to
   ${CMAKE_INSTALL_PREFIX}/lib/
- cmake script ${ars_ROOT}/cmake_modules/arsConfig.cmake to
   ${CMAKE_INSTALL_PREFIX}/share/ars/


HOW TO USE LIBRARY ars IN YOUR PROJECT
-------------------------------------------------

If library ars has been installed in system directory "/usr/local",
then it is straighforward to use it in your projects.
You need to add the following lines to your project as in this example:


> CMAKE_MINIMUM_REQUIRED(VERSION 2.8)  
> PROJECT(foobar)  
> 
> find_package(ars REQUIRED)  
> message(STATUS "ars_FOUND ${ars_FOUND}")  
> message(STATUS "ars_INCLUDE_DIRS ${ars_INCLUDE_DIRS}")  
> message(STATUS "ars_LIBRARY_DIRS ${ars_LIBRARY_DIRS}")  
> message(STATUS "ars_LIBRARIES ${ars_LIBRARIES}")  
>  
> if(${ars_FOUND})   
>   include_directories(${ars_INCLUDE_DIRS})  
>   link_directories(${ars_LIBRARY_DIRS})  
> endif()  
> 
> add_executable(foobar foobar.cpp)  
> target_link_libraries(foobar ${ars_LIBRARIES})  

The above example uses the variables defined in arsConfig.cmake:

-  ars_FOUND - system has ars module
-  ars_INCLUDE_DIRS - the ars include directories
-  ars_LIBRARY_DIRS - the ars library directories
-  ars_LIBRARIES - link these to use ars


