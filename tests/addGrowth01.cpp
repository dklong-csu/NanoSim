/*
Tests that adding a single reaction constructs the intended reaction and ODE system
*/

#include "particleSystem.h"
#include "operationsSUNDenseMatrix.h"

#include <sundials/sundials_nvector.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <iostream>

int main(){
  NanoSim::particleSystem<double> my_rxns;

  my_rxns.defineParticle(1,3,3);
  const std::function<realtype(const unsigned int)> growth_kernel = [](const unsigned int i){
    return 5.0 * i;
  };
  my_rxns.addGrowth({{1,"A"}},{},growth_kernel);

  my_rxns.finalizeReactions();
  my_rxns.printChemicalReactions();

  N_Vector x, x_dot, tmp1, tmp2, tmp3;
  sunindextype length = 4;

  x = N_VNew_Serial(length);
  x_dot = N_VNew_Serial(length);

  realtype *x_data = N_VGetArrayPointer(x);
  for (unsigned int i=0; i<length;++i){
    x_data[i] = i + 2.0;
  }

  realtype *x_dot_data = N_VGetArrayPointer(x_dot);

  const auto rhs_fcn = my_rxns.composeRHSfunction();
  
  const int err = rhs_fcn(0.0, x, x_dot, nullptr);

  // x values should be untouched by application of rhs function
  std::cout << "x = \n";
  for (unsigned int i=0;i<length;++i){
    std::cout << x_data[i] << "\n";
  }
  
 
  std::cout << "dx = \n";
  for (unsigned int i=0;i<length;++i){
    std::cout << x_dot_data[i] << "\n";
  }

  const auto jac_fcn = my_rxns.composeJacobianfunction();
 
  tmp1 = N_VNew_Serial(length);
  tmp2 = N_VNew_Serial(length);
  tmp3 = N_VNew_Serial(length);


  SUNMatrix J = SUNDenseMatrix(length,length);


  std::function<int(SUNMatrix,sunindextype,sunindextype,realtype)> matrixInsertAdd = [](SUNMatrix A, sunindextype row, sunindextype col, realtype value){
    return NanoSim::matrixInsertAdd(A, row, col, value);
  };
  
  void * user_data = static_cast<void*>(&matrixInsertAdd);
  const int errJ = jac_fcn(0.0, x, x_dot, J, user_data, tmp1, tmp2, tmp3);
  realtype* J_data = SUNDenseMatrix_Data(J);

  std::cout << "J =\n";
  for(unsigned int row=0; row<length; ++row){
    for(unsigned int col=0; col<length; ++col){
      std::cout << J_data[col*length + row] << " ";
      if (col == length-1){
        std::cout << "\n";
      }
    }
  }

  N_VDestroy(x);
  N_VDestroy(x_dot);
  N_VDestroy(tmp1);
  N_VDestroy(tmp2);
  N_VDestroy(tmp3);
  SUNMatDestroy(J);
}