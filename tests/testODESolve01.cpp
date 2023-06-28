// Tests that the ODE solver is working properly

#include "particleSystem.h"
#include "linearAlgebraSUNDense.h"
#include "odeSolver.h"

#include <sundials/sundials_nvector.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>

#include <iostream>
#include <iomanip>

int main(){
  NanoSim::particleSystem<double> my_rxns;

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

  N_Vector ic;
  NanoSim::sunDenseLinearAlgebraOperations<realtype> lin_alg;
  const sunindextype length = my_rxns.getNumberOfSpecies();
  ic = N_VNew_Serial(length);
  lin_alg.vectorInsert(ic, 1.0, 0);
  lin_alg.vectorInsert(ic, 0.0, 1);
  // for (sunindextype i=0;i<length;++i){
  //   lin_alg.vectorInsert(ic,1.0*(i+1),i);
  // }

  SUNMatrix template_matrix = SUNDenseMatrix(length,length);

  SUNLinearSolver lin_solve = SUNLinSol_Dense(ic, template_matrix);

  NanoSim::cvodeOptions<realtype> opts("SUNDIALS_errors.txt",
    1e-8,
    1e-14,
    1000,
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