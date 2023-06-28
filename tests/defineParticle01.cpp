/*
Tests that defining a particle creates the appropriate chemical species
*/

#include "particleSystem.h"
#include <iostream>

int main(){
  NanoSim::particleSystem<double> my_rxns;


  my_rxns.defineParticle(1,3,2);

  my_rxns.finalizeReactions();
  my_rxns.printChemicalSpecies();
}