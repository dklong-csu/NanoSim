CMAKE_MINIMUM_REQUIRED (VERSION 3.1)
PROJECT (NanoSim CXX)

###################
## Compiler flags
set (CMAKE_CXX_STANDARD 20)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_EXTENSIONS OFF)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MARCH_NATIVE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()


#########################################
### Find the SUNDIALS library
FIND_PACKAGE (
        _sundials
        NAMES SUNDIALS
        HINTS ${SUNDIALS_DIR} 
        REQUIRED
)
IF ("${_sundials}" STREQUAL "_sundials-NOTFOUND")
  MESSAGE(FATAL_ERROR
          "The SUNDIALS library was not found. You have to specify a path "
          "to that library by setting the SUNDIALS_DIR environment variable, "
          "or by passing '-DSUNDIALS_DIR=...' as an argument to 'cmake'.")
ELSE()
  MESSAGE(STATUS "Found SUNDIALS")
ENDIF()


FIND_PATH(_sundials_include_dir
        NAMES cvode/cvode.h
        HINTS ${SUNDIALS_DIR}/include)
IF ("${_sundials_include_dir}" STREQUAL "_sundials_include_dir-NOTFOUND")
  MESSAGE(FATAL_ERROR
          "The SUNDIALS headers were not found. You have to specify a path "
          "to that library by setting the SUNDIALS_DIR environment variable, "
          "or by passing '-DSUNDIALS_DIR=...' as an argument to 'cmake'.")
ELSE()
  MESSAGE(STATUS "Found SUNDIALS headers at ${_sundials_include_dir}")
  INCLUDE_DIRECTORIES(${_sundials_include_dir})
ENDIF()


FIND_LIBRARY(_sundials_cvode_lib
        NAMES sundials_cvode
        HINTS ${SUNDIALS_DIR}/lib)
IF ("${_sundials_cvode_lib}" STREQUAL "_sundials_cvode_lib-NOTFOUND")
  MESSAGE(FATAL_ERROR
          "The SUNDIALS CVODE library was not found. You have to specify a path "
          "to that library by setting the SUNDIALS_DIR environment variable, "
          "or by passing '-DSUNDIALS_DIR=...' as an argument to 'cmake'.")
ELSE()
  MESSAGE(STATUS "Found SUNDIALS CVODE library at ${_sundials_cvode_lib}")
  #INCLUDE_DIRECTORIES(${_sundials_cvode_lib})
  LINK_DIRECTORIES(${SUNDIALS_DIR}/lib)
ENDIF()

#########################################
### Find the EIGEN library
if(NOT DEFINED EIGEN_DIR)
  message(FATAL_ERROR
    "No path to the Eigen software is provided."
    "You have to specify a path to that library by "
    "passing '-DEIGEN_DIR=...' as an argument to 'cmake'."
    "The directory specified should access header files via "
    "'Eigen/Dense', 'Eigen/Sparse', etc.")
endif()
INCLUDE_DIRECTORIES(${EIGEN_DIR})



# #########################################
# ### Add optimlib for optimization
if(NOT DEFINED OPTIMLIB_DIR)
  message(FATAL_ERROR
    "Provide the path to the folder containing 'optim.hpp' with '-DOPTIMLIB_DIR=/path/to/optimlib/')")
endif()
include_directories(${OPTIMLIB_DIR})

##################
### MUQ for uncertainty quantification
FIND_PACKAGE (
        _muq
        NAMES MUQ
        HINTS ${MUQ_DIR} 
        REQUIRED
)
IF ("${_muq}" STREQUAL "_muq-NOTFOUND")
  MESSAGE(FATAL_ERROR
          "The MUQ library was not found. You have to specify a path "
          "to that library by setting the MUQ_DIR environment variable, "
          "or by passing '-DMUQ_DIR=...' as an argument to 'cmake'.")
ELSE()
  MESSAGE(STATUS "Found MUQ")
  include_directories(${MUQ_INCLUDE_DIRS})
  MESSAGE(STATUS "Include MUQ directories: ${MUQ_INCLUDE_DIRS}")
ENDIF()


MESSAGE(STATUS "MUQ libraries: ${MUQ_LIBRARIES}")
MESSAGE(STATUS "MUQ link libraries: ${MUQ_LINK_LIBRARIES}")


#########################################
### Add Boost library
find_package(Boost REQUIRED)
if (Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  message(STATUS "Found Boost include directories at ${Boost_INCLUDE_DIRS}.")
  message(STATUS "Found Boost libraries at ${Boost_LIBRARY_DIRS}.")
  message(STATUS "Including libraries: ${Boost_LIBRARIES}.")

else ()
  message(FATAL_ERROR
    "Boost was not found. Install Boost, ensure it can be found in your PATH variable, and try again.")
endif ()

###################
## Find include files in this project
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/lib)

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
      ${_sundials_cvode_lib}
      ${Boost_LIBRARIES}
      ${MUQ_LIBRARIES}
      ${MUQ_LINK_LIBRARIES}
      OpenMP::OpenMP_CXX)
else ()
  target_link_libraries(libNanoSim
    INTERFACE
    ${_sundials_cvode_lib}
    ${Boost_LIBRARIES}
    ${MUQ_LIBRARIES}
    ${MUQ_LINK_LIBRARIES})
endif()




###################
## Set up the test suite machinery
enable_testing()
add_subdirectory(tests)

###################
## Add the tutorial programs
add_subdirectory(tutorial)
add_subdirectory(code_gallery)