/*
Tests that adding a single reaction constructs the intended reaction and ODE system
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
    Do reaction:
    2A + 3B ->[5.0] 4C  
  */
 my_rxns.addReaction({{2,"A"}, {3,"B"}},
  {{4,"C"}},
  5.0);

  my_rxns.finalizeReactions();

  N_Vector x, x_dot, tmp1, tmp2, tmp3;
  sunindextype length = 3;

  x = N_VNew_Serial(length);
  x_dot = N_VNew_Serial(length);

  NanoSim::sunDenseLinearAlgebraOperations<realtype> lin_alg;
  lin_alg.vectorInsert(x, 1.0, 0);
  lin_alg.vectorInsert(x, 2.0, 1);
  lin_alg.vectorInsert(x, 3.0, 2);

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