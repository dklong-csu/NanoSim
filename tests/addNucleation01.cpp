/*
Tests that adding a nucleation event creates the proper reactions
*/

#include "particleSystem.h"
#include "linearAlgebraEigen.h"
#include "testingUtilities.h"
#include <iostream>

using EigenMatrix = Eigen::Matrix<realtype, Eigen::Dynamic, Eigen::Dynamic>;

int main(){
  NanoSim::particleSystem<double> my_rxns;
  NanoSim::eigenLinearAlgebraOperations<realtype, EigenMatrix> lin_alg;

  /* 
    Nucleation is
    2A ->[2] B_2 + C
  */
 std::function<realtype(const int)> atoms2diameter 
    = [](const int atoms){ return 0.3 * std::cbrt(1.0*atoms);};
  my_rxns.defineParticle(2,3, atoms2diameter);
  my_rxns.addNucleation({{2,"A"}},
    {{1,"C"}},
    2.0);

  my_rxns.finalizeReactions();

  sunindextype length = 4;

  N_Vector x = lin_alg.createNewVector(length);
  N_Vector x_dot = lin_alg.createNewVector(length);

  for (sunindextype i=0; i<length; ++i){
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