/*
Tests that adding a nucleation event creates the proper reactions
*/

#include "particleSystem.h"
#include "linearAlgebraSUNDense.h"
#include "testingUtilities.h"

#include <sundials/sundials_nvector.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <iostream>

int main(){
  NanoSim::particleSystem<double> my_rxns;

  /* 
    Nucleation is
    2A ->[2] B_2 + C
    Other reaction is
    A + D ->[3] E
  */
  my_rxns.addReaction({{1,"A"},{1,"D"}},
    {{1,"E"}},
    3.0);
  my_rxns.defineParticle(2,2,2);
  my_rxns.addNucleation({{2,"A"}},
    {{1,"C"}},
    2.0);

  my_rxns.finalizeReactions();

  N_Vector x, x_dot, tmp1, tmp2, tmp3;
  sunindextype length = 5;

  x = N_VNew_Serial(length);
  x_dot = N_VNew_Serial(length);

  NanoSim::sunDenseLinearAlgebraOperations<realtype> lin_alg;
  for (sunindextype i=0; i<length; ++i){
    lin_alg.vectorInsert(x,i+2.0,i);
  }

  const auto rhs_fcn = my_rxns.composeRHSfunction();
  void * user_data = static_cast<void *>(&lin_alg);
  const int err = rhs_fcn(0.0, x, x_dot, user_data);

  const auto jac_fcn = my_rxns.composeJacobianfunction();
 
  tmp1 = N_VNew_Serial(length);
  tmp2 = N_VNew_Serial(length);
  tmp3 = N_VNew_Serial(length);


  SUNMatrix J = SUNDenseMatrix(length,length);

  const int errJ = jac_fcn(0.0, x, x_dot, J, user_data, tmp1, tmp2, tmp3);

  NanoSim::Testing::printParticleSystemOutput<realtype>(my_rxns, x, x_dot, J, lin_alg);

  N_VDestroy(x);
  N_VDestroy(x_dot);
  N_VDestroy(tmp1);
  N_VDestroy(tmp2);
  N_VDestroy(tmp3);
  SUNMatDestroy(J);
}