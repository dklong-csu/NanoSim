CMAKE_MINIMUM_REQUIRED (VERSION 3.1)
MESSAGE(STATUS "Setting up tests")

FIND_PROGRAM(NUMDIFF_EXECUTABLE
  NAMES numdiff
  HINTS ${NUMDIFF_DIR}
  PATH_SUFFIXES bin
)
mark_as_advanced(NUMDIFF_EXECUTABLE)

IF( "${NUMDIFF_EXECUTABLE}" MATCHES "NUMDIFF_EXECUTABLE-NOTFOUND")
  MESSAGE(FATAL_ERROR
    "Could not find 'numdiff', which is required for running the testsuite.\n"
    "Please specify NUMDIFF_DIR to a location containing the binary, or\n"
    "install 'numdiff' via a package manage.")
ENDIF()

# Create a target that we can use for 'make check'.
ADD_CUSTOM_TARGET(check COMMAND ${CMAKE_CTEST_COMMAND})

# Loop over all .cc files in this directory and make tests out of them.
FILE(GLOB _testfiles "*cpp")
FOREACH(_testfile ${_testfiles})
  STRING(REPLACE ".cpp" "" _testname ${_testfile})
  STRING(REPLACE "${CMAKE_CURRENT_SOURCE_DIR}/" "" _testname ${_testname})
  MESSAGE(STATUS "  ${_testname}")

  # For each test, add the corresponding executable and the commands
  # necessary to run the test. Compile the executable with -Werror
  # to ensure that there are no warnings in either the tests, or the
  # header files these tests #include.
  ADD_EXECUTABLE(${_testname}
                 EXCLUDE_FROM_ALL
                 ${_testfile})
  IF(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    TARGET_COMPILE_OPTIONS(${_testname} PRIVATE "/WX")
  ELSE()
    TARGET_COMPILE_OPTIONS(${_testname} PRIVATE "-Werror")
  ENDIF()
  TARGET_LINK_LIBRARIES (${_testname} libNanoSim)


  # Then specify what it means to run a test:
  #
  # Step 1 is to:
  # - remove the previous .result file just in case the executable aborts before
  #   it produces any output (which leaves the previous output file in place)
  # - execute the test and write the output to a .result file
  ADD_CUSTOM_COMMAND(
    OUTPUT ${_testname}.result
    COMMAND rm -f ${_testname}.result
    COMMAND ${_testname} > ${_testname}.result
    DEPENDS ${_testname}
    COMMENT "Running test <${_testname}>...")

  # Step 2 is to:
  # - remove the .ok file
  # - compare the .result file to the .output file
  # - if this succeeds, write a .ok file
  ADD_CUSTOM_COMMAND(
    OUTPUT ${_testname}.ok
    COMMAND rm -f ${_testname}.ok
    COMMAND ${NUMDIFF_EXECUTABLE} -a 1e-6 -r 1e-8 -s "' \\t\\n\\r=,:;<>[](){}^'" ${CMAKE_CURRENT_SOURCE_DIR}/${_testname}.output ${_testname}.result
    COMMAND touch ${_testname}.ok
    DEPENDS ${_testname}.result ${_testname}.output
    COMMENT "Comparing test <${_testname}>...")

  # Now make success of the test be dependent on there being an .ok file:
  ADD_CUSTOM_TARGET(run_${_testname}
                    DEPENDS ${_testname}.ok)

  # Define running the test as executing the cmake program
  # with the given target.
  ADD_TEST(NAME ${_testname} 
           COMMAND "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target run_${_testname})

  # Now also make sure that the test is part of the 'make check' command.
  ADD_DEPENDENCIES(check ${_testname})
ENDFOREACH()