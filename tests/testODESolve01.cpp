// Tests that the ODE solver is working properly

#include "particleSystem.h"
#include "linearAlgebraEigen.h"
#include "linearSolverEigen.h"
#include "odeSolver.h"

#include <iostream>
#include <iomanip>

using EigenMatrix = Eigen::Matrix<realtype, Eigen::Dynamic, Eigen::Dynamic>;
using SolverType = Eigen::PartialPivLU< EigenMatrix >;
int main(){
  NanoSim::particleSystem<double> my_rxns;
  NanoSim::eigenLinearAlgebraOperations<realtype, EigenMatrix> lin_alg;

  /* 
    Do reaction:
    A->[1.0] B
    Initial condition: A=1, B=0
    Known exact solution:
    A(t) = exp(-t)
    B(t) = 1 - exp(-t)  
  */
 my_rxns.addReaction({{1,"A"}},
  {{1,"B"}},
  1.0);

  my_rxns.finalizeReactions();

  const sunindextype length = my_rxns.getNumberOfSpecies();

  N_Vector ic = lin_alg.createNewVector(length);
  lin_alg.vectorInsert(ic, 1.0, 0);
  lin_alg.vectorInsert(ic, 0.0, 1);


  SUNMatrix template_matrix = lin_alg.createNewMatrix(length,length);

  SUNLinearSolver lin_solve = NanoSim::createLinearSolverEigenDense<realtype, SolverType>();

  NanoSim::cvodeOptions<realtype> opts("SUNDIALS_errors.txt",
    1e-8,
    1e-14,
    1000,
    CV_ADAMS,
    SUNFALSE);
  auto cvode_mem = NanoSim::prepareODESolver<realtype>(ic, 
    template_matrix, 
    lin_solve, 
    opts);

  std::pair< NanoSim::particleSystem<realtype>*, NanoSim::abstractLinearAlgebraOperations<realtype>* >
  data_pair = {&my_rxns, &lin_alg};
  void * user_data = static_cast<void *>(&data_pair);
  auto flag = CVodeSetUserData(cvode_mem, user_data);

  for (unsigned int i=0;i<5;++i){
    N_Vector sol = ic->ops->nvclone(ic);
    realtype t;
    auto err = CVode(cvode_mem, 0.2 + i*0.2, sol, &t, CV_NORMAL);

    for (sunindextype i=0; i<length;++i){
      std::cout << std::setprecision(20) << lin_alg.vectorGetValue(sol, i) << "\n";
    }

    sol->ops->nvdestroy(sol);
  }

  CVodeFree(&cvode_mem);
  ic->ops->nvdestroy(ic);
  SUNLinSolFree(lin_solve);
  SUNMatDestroy(template_matrix);

}