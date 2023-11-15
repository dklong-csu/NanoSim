/*
Tests that adding a single reaction constructs the intended reaction and ODE system
*/

#include "particleSystem.h"
#include "linearAlgebraEigen.h"
#include "testingUtilities.h"
#include <iostream>

using EigenMatrix = Eigen::Matrix<realtype, Eigen::Dynamic, Eigen::Dynamic>;

int main(){
  NanoSim::particleSystem<double> my_rxns;
  NanoSim::eigenLinearAlgebraOperations<realtype, EigenMatrix> lin_alg;

  std::function<realtype(const int)> atoms2diameter 
    = [](const int atoms){ return 0.3 * std::cbrt(1.0*atoms);};

  my_rxns.defineParticle(1,6, atoms2diameter);

  my_rxns.addReaction({{2,"A"}, {3,"B"}},
  {{4,"C"}},
  5.0);

  my_rxns.addNucleation({{1,"A"}},
    {{1,"C"}},
    2.0);

  const std::function<realtype(const unsigned int)> growth_kernel = [](const unsigned int i){
    return 5.0 * i;
  };
  my_rxns.addGrowth({{1,"A"}},{},growth_kernel);

  const std::function<realtype(const unsigned int, const unsigned int)> agglom_kernel 
  = [](const unsigned int i, const unsigned int j){
    if (i > 3 || j > 3){
      return 0.0;
    } else {
      return 5.0 * i * j;
    }
  };
  my_rxns.addAgglomeration({},{},agglom_kernel);

  my_rxns.finalizeReactions();

  sunindextype length = 9;

  N_Vector x = lin_alg.createNewVector(length);
  N_Vector x_dot = lin_alg.createNewVector(length);

  for (unsigned int i=0; i<length;++i){
    lin_alg.vectorInsert(x,i+2.0,i);
  }

  const auto rhs_fcn = my_rxns.composeRHSfunction();
  
  void * user_data = static_cast<void *>(&lin_alg);
  const int err = rhs_fcn(0.0, x, x_dot, user_data);

  const auto jac_fcn = my_rxns.composeJacobianfunction();
 
  N_Vector tmp1 = lin_alg.createNewVector(length);
  N_Vector tmp2 = lin_alg.createNewVector(length);
  N_Vector tmp3 = lin_alg.createNewVector(length);

  SUNMatrix J = lin_alg.createNewMatrix(length,length);

  const int errJ = jac_fcn(0.0, x, x_dot, J, user_data, tmp1, tmp2, tmp3);

  NanoSim::Testing::printParticleSystemOutput<realtype>(my_rxns, x, x_dot, J, lin_alg);

  N_VDestroy(x);
  N_VDestroy(x_dot);
  N_VDestroy(tmp1);
  N_VDestroy(tmp2);
  N_VDestroy(tmp3);
  SUNMatDestroy(J);
}