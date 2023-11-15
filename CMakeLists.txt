CMAKE_MINIMUM_REQUIRED (VERSION 3.1...3.27)
PROJECT (NanoSim CXX)

###################
## Compiler flags
set (CMAKE_CXX_STANDARD 20)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_EXTENSIONS OFF)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MARCH_NATIVE)
  message(STATUS Using -march=native)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()


#########################################
### Find the SUNDIALS library
find_package(SUNDIALS REQUIRED)
find_package(Matplot++ REQUIRED)
find_package(Boost REQUIRED)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif()


###################
## Optimization library






###################
## Find include files in this project
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/lib ${CMAKE_CURRENT_SOURCE_DIR}/lib/eigen)

###################
## Define interface library to tells executable how to compile properly
add_library(libNanoSim INTERFACE)
target_include_directories(libNanoSim
    INTERFACE 
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
find_package(OpenMP)
if (OpenMP_CXX_FOUND)
  message(STATUS Found OpenMP)
  target_link_libraries(libNanoSim
      INTERFACE
      SUNDIALS::cvode
      OpenMP::OpenMP_CXX
      Matplot++::matplot
      -lnlopt -lm)
else ()
  target_link_libraries(libNanoSim
    INTERFACE
    SUNDIALS::cvode)
endif()




###################
## Set up the test suite machinery
enable_testing()
add_subdirectory(tests)

###################
## Add the tutorial programs
add_subdirectory(tutorial)
add_subdirectory(examples)