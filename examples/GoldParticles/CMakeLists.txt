CMAKE_MINIMUM_REQUIRED (VERSION 3.1...3.27)

add_executable(GoldParticleSimulation
GoldParticleSimulation.cpp)
target_link_libraries(GoldParticleSimulation
  libNanoSim)