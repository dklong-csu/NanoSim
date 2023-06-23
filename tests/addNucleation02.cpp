/*
Tests that adding a nucleation event creates the proper reactions
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
  // Order should be [A, D, E, C, B_2]
  my_rxns.printChemicalReactions();

  N_Vector x, x_dot, tmp1, tmp2, tmp3;
  sunindextype length = 5;

  x = N_VNew_Serial(length);
  x_dot = N_VNew_Serial(length);

  realtype *x_data = N_VGetArrayPointer(x);
  x_data[0] = 2.0;
  x_data[1] = 3.0;
  x_data[2] = 4.0;
  x_data[3] = 5.0;
  x_data[4] = 6.0;

  realtype *x_dot_data = N_VGetArrayPointer(x_dot);

  const auto rhs_fcn = my_rxns.composeRHSfunction();
  
  const int err = rhs_fcn(0.0, x, x_dot, nullptr);

  // should be [2, 3, 4, 5, 6] because untouched by rhs function
  std::cout << "x = \n" << x_data[0] << "\n" 
                       << x_data[1] << "\n" 
                       << x_data[2] << "\n"
                       << x_data[3] << "\n"
                       << x_data[4] << "\n";
  

  // should be:
  // dA/dt = -2*2*(A=2^2) -3*(A=2)*(D=3) = -34
  // dD/dt = -3*(A=2)*(D=3) = -18
  // dE/dt = 3*(A=2)*(D=3) =  18
  // dC/dt = 1*2*(A=2^2) =   8
  // dB_2/dt = 1*2*(A=2^2) =   8
  std::cout << "dx = \n" << x_dot_data[0] << "\n" 
                           << x_dot_data[1] << "\n" 
                           << x_dot_data[2] << "\n"
                           << x_dot_data[3] << "\n"
                           << x_dot_data[4] << "\n";

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
  
  // Should be
  /*
    -2*2*2*2 -3*3 = -25     -3*2 = -6   0   0   0
    -3*3 = -9               -3*2 = -6   0   0   0
    3*3 = 9                 3*2 = 6     0   0   0
    2*2*2 = 8               0           0   0   0
    2*2*2 = 8               0           0   0   0
  */
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