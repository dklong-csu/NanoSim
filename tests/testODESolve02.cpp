// Tests that the ODE solver is working properly
// Reference solution is from matlab code:
/*
y0 = [1;0;0];
tspan = [0, 500];
options = odeset('RelTol',1e-10,'AbsTol',1e-15);
tic
sol = ode15s(@rhsfcn,tspan,y0, options);
toc

y1 = deval(sol,100);
fprintf("%.20f\n",y1)
y2 = deval(sol,200);
fprintf("%.20f\n",y2)
y3 = deval(sol,300);
fprintf("%.20f\n",y3)
y4 = deval(sol,400);
fprintf("%.20f\n",y4)
y5 = deval(sol,500);
fprintf("%.20f\n",y5)

function dy = rhsfcn(t,Y)
    dy = Y;
    x = Y(1);
    y = Y(2);
    z = Y(3);
    dy(1) = -0.04*x + 10^4 * y*z;
    dy(2) = 0.04 * x - 10^4 * y *z - 3e7 * y * y;
    dy(3) = 3e7 * y * y;
end
*/

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


 my_rxns.addReaction({{1,"X"}},
  {{1,"Y"}},
  0.04);
  my_rxns.addReaction({{2,"Y"}},
    {{2,"Z"}},
    1.5e7);
  my_rxns.addReaction({{1,"Y"},{1,"Z"}},
    {{1,"X"},{1,"Z"}},
    1.0e4);

  my_rxns.finalizeReactions();

  N_Vector ic;
  NanoSim::sunDenseLinearAlgebraOperations<realtype> lin_alg;
  const sunindextype length = my_rxns.getNumberOfSpecies();
  ic = N_VNew_Serial(length);
  lin_alg.vectorInsert(ic, 1.0, 0);
  lin_alg.vectorInsert(ic, 0.0, 1);
  lin_alg.vectorInsert(ic, 0.0, 2);

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
    auto err = CVode(cvode_mem, 100.0 + i*100.0, sol, &t, CV_NORMAL);

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