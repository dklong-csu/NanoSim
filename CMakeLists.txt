CMAKE_MINIMUM_REQUIRED (VERSION 3.1)
PROJECT (NanoSim CXX)

###################
## Compiler flags
set (CMAKE_CXX_STANDARD 20)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_EXTENSIONS OFF)


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
  MESSAGE(STATUS "Found SUNDIALS CVODE library at ${_sundials_CVode_lib}")
  INCLUDE_DIRECTORIES(${_sundials_include_dir})
  LINK_DIRECTORIES(${SUNDIALS_DIR}/lib)
ENDIF()


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
target_link_libraries(libNanoSim
    INTERFACE
    ${_sundials_cvode_lib})

###################
## Set up the test suite machinery
enable_testing()
add_subdirectory(tests)