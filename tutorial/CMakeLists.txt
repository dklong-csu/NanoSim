CMAKE_MINIMUM_REQUIRED (VERSION 3.1)

add_executable(step01
  step01-solveRxnODE.cpp)
target_link_libraries(step01
  libNanoSim)

add_executable(step02
  step02-parametrizedRxn.cpp)
target_link_libraries(step02
  libNanoSim)

add_executable(step03
  step03-costFunction.cpp)
target_link_libraries(step03
  libNanoSim)

add_executable(step04
  step04-optimization.cpp)
target_link_libraries(step04
  libNanoSim)

add_executable(step05
  step05-solveRxnFast.cpp)
target_link_libraries(step05
  libNanoSim)