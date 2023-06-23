/*
Tests that adding a multiple reactions constructions the correct system
*/

#include "particleSystem.h"
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
    2C ->[7.0] -> 3A + 4B  
  */
 my_rxns.addReaction({{2,"A"}, {3,"B"}},
  {{4,"C"}},
  5.0);
  my_rxns.addReaction({{2,"C"}},
    {{3,"A"},{4,"B"}},
    7.0);

  my_rxns.finalizeReactions();
  my_rxns.printChemicalReactions();

  N_Vector x, x_dot, tmp1, tmp2, tmp3;
  sunindextype length = 3;

  x = N_VNew_Serial(length);
  x_dot = N_VNew_Serial(length);

  realtype *x_data = N_VGetArrayPointer(x);
  x_data[0] = 2.0;
  x_data[1] = 3.0;
  x_data[2] = 4.0;

  realtype *x_dot_data = N_VGetArrayPointer(x_dot);

  const auto rhs_fcn = my_rxns.composeRHSfunction();
  
  const int err = rhs_fcn(0.0, x, x_dot, nullptr);

  // should be [2, 3, 4] because untouched by rhs function
  std::cout << "x = \n" << x_data[0] << "\n" 
                       << x_data[1] << "\n" 
                       << x_data[2] << "\n";
  

  // should be:
  // -2*5*(2^2)*(3^3) + 3*7*(4^2) = -744
  // -3*5*(2^2)*(3^3) + 4*7*(4^2) = -1172
  //  4*5*(2^2)*(3^3) - 2*7*(4^2) = 1936
  std::cout << "dx = \n" << x_dot_data[0] << "\n" 
                           << x_dot_data[1] << "\n" 
                           << x_dot_data[2] << "\n";

  const auto jac_fcn = my_rxns.composeJacobianfunction();
 
  tmp1 = N_VNew_Serial(length);
  tmp2 = N_VNew_Serial(length);
  tmp3 = N_VNew_Serial(length);


  SUNMatrix J = SUNDenseMatrix(3,3);


  std::function<int(SUNMatrix,sunindextype,sunindextype,realtype)> matrixInsertAdd = [](SUNMatrix A, sunindextype row, sunindextype col, realtype value){
    // Specific to SUNMatrixDense shipped with SUNDIALS
    realtype* A_data = SUNDenseMatrix_Data(A);
    // Organized as array indexed with [col*n_rows + row] with row and col counting from 0
    const sunindextype n_rows = SUNDenseMatrix_Rows(A);
    A_data[col*n_rows + row] = value;
    return 0;
  };
  void * user_data = static_cast<void*>(&matrixInsertAdd);
  const int errJ = jac_fcn(0.0, x, x_dot, J, user_data, tmp1, tmp2, tmp3);
  realtype* J_data = SUNDenseMatrix_Data(J);
  
  // Should be
  /*
    -2*5*2*(2)*(3^3) = -1080    -3*5*2*(2^2)*(3^2) = -1080      2*7*3*(4) = 168
    -2*5*3*(2)*(3^3) = -1620    -3*5*3*(2^2)*(3^2) = -1620      2*7*4*(4) = 224
     2*5*4*(2)*(3^3) =  2160     3*5*4*(2^2)*(3^2) =  2160     -2*7*2*(4) = -112
  */
  std::cout << "J =\n";
  for(unsigned int row=0; row<3; ++row){
    for(unsigned int col=0; col<3; ++col){
      std::cout << J_data[col*3 + row] << " ";
      if (col == 2){
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