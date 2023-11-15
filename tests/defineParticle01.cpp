/*
Tests that defining a particle creates the appropriate chemical species
*/

#include "particleSystem.h"
#include <iostream>

int main(){
  NanoSim::particleSystem<double> my_rxns;

  std::function<realtype(const int)> atoms2diameter 
    = [](const int atoms){ return 0.3 * std::cbrt(1.0*atoms);};

  my_rxns.defineParticle(1,3, atoms2diameter);

  my_rxns.finalizeReactions();
  my_rxns.printChemicalSpecies();
}