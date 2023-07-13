/*
Tests that having missing reactants or products results in intended behavior
*/

#include "particleSystem.h"
#include "linearAlgebraEigen.h"
#include "testingUtilities.h"
#include <iostream>

using EigenMatrix = Eigen::Matrix<realtype, Eigen::Dynamic, Eigen::Dynamic>;

int main(){
  {
    NanoSim::particleSystem<double> my_rxns;
    NanoSim::eigenLinearAlgebraOperations<realtype, EigenMatrix> lin_alg;

    /* 
      Do reaction:
      ->[5.0] 2A  
    */
    my_rxns.addReaction({},
      {{2,"A"}},
      5.0);

    my_rxns.finalizeReactions();

    sunindextype length = 1;

    N_Vector x = lin_alg.createNewVector(length);
    N_Vector x_dot = lin_alg.createNewVector(length);

    lin_alg.vectorInsert(x,1.0,0);

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

  {
    NanoSim::particleSystem<double> my_rxns;
    NanoSim::eigenLinearAlgebraOperations<realtype, EigenMatrix> lin_alg;

    /* 
      Do reaction:
      2A ->[5.0]   
    */
    my_rxns.addReaction({{2,"A"}},
      {},
      5.0);

    my_rxns.finalizeReactions();

    sunindextype length = 1;

    N_Vector x = lin_alg.createNewVector(length);
    N_Vector x_dot = lin_alg.createNewVector(length);

    lin_alg.vectorInsert(x,1.0,0);

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
}