/*
Tests that having missing reactants or products results in intended behavior
*/

#include "particleSystem.h"
#include "operationsSUNDenseMatrix.h"

#include <sundials/sundials_nvector.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <iostream>

int main(){
  {
    NanoSim::particleSystem<double> my_rxns;

    /* 
      Do reaction:
      ->[5.0] 2A  
    */
  my_rxns.addReaction({},
    {{2,"A"}},
    5.0);

    my_rxns.finalizeReactions();
    my_rxns.printChemicalReactions();

    N_Vector x, x_dot, tmp1, tmp2, tmp3;
    sunindextype length = 1;

    x = N_VNew_Serial(length);
    x_dot = N_VNew_Serial(length);

    realtype *x_data = N_VGetArrayPointer(x);
    x_data[0] = 1.0;

    realtype *x_dot_data = N_VGetArrayPointer(x_dot);

    const auto rhs_fcn = my_rxns.composeRHSfunction();
    
    const int err = rhs_fcn(0.0, x, x_dot, nullptr);

    // should be [1] because untouched by rhs function
    std::cout << "x = \n" << x_data[0] << "\n";
    

    // should be:
    // 2*5 = 10 
    std::cout << "dx = \n" << x_dot_data[0] << "\n";

    const auto jac_fcn = my_rxns.composeJacobianfunction();
  
    tmp1 = N_VNew_Serial(length);
    tmp2 = N_VNew_Serial(length);
    tmp3 = N_VNew_Serial(length);


    SUNMatrix J = SUNDenseMatrix(1,1);


    std::function<int(SUNMatrix,sunindextype,sunindextype,realtype)> matrixInsertAdd = [](SUNMatrix A, sunindextype row, sunindextype col, realtype value){
      return NanoSim::matrixInsertAdd(A, row, col, value);
    };

    void * user_data = static_cast<void*>(&matrixInsertAdd);
    const int errJ = jac_fcn(0.0, x, x_dot, J, user_data, tmp1, tmp2, tmp3);
    realtype* J_data = SUNDenseMatrix_Data(J);

    // Should be: 0
    std::cout << "J =\n";
    for(unsigned int row=0; row<1; ++row){
      for(unsigned int col=0; col<1; ++col){
        std::cout << J_data[col*1 + row] << " ";
        if (col == 0){
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

  {
    NanoSim::particleSystem<double> my_rxns;

    /* 
      Do reaction:
      2A ->[5.0]   
    */
  my_rxns.addReaction({{2,"A"}},
    {},
    5.0);

    my_rxns.finalizeReactions();
    my_rxns.printChemicalReactions();

    N_Vector x, x_dot, tmp1, tmp2, tmp3;
    sunindextype length = 1;

    x = N_VNew_Serial(length);
    x_dot = N_VNew_Serial(length);

    realtype *x_data = N_VGetArrayPointer(x);
    x_data[0] = 1.0;

    realtype *x_dot_data = N_VGetArrayPointer(x_dot);

    const auto rhs_fcn = my_rxns.composeRHSfunction();
    
    const int err = rhs_fcn(0.0, x, x_dot, nullptr);

    // should be [1] because untouched by rhs function
    std::cout << "x = \n" << x_data[0] << "\n";
    

    // should be:
    // -2*(1^2)*5 = -10
    std::cout << "dx = \n" << x_dot_data[0] << "\n";

    const auto jac_fcn = my_rxns.composeJacobianfunction();
  
    tmp1 = N_VNew_Serial(length);
    tmp2 = N_VNew_Serial(length);
    tmp3 = N_VNew_Serial(length);


    SUNMatrix J = SUNDenseMatrix(1,1);


    std::function<int(SUNMatrix,sunindextype,sunindextype,realtype)> matrixInsertAdd = [](SUNMatrix A, sunindextype row, sunindextype col, realtype value){
      return NanoSim::matrixInsertAdd(A, row, col, value);
    };

    void * user_data = static_cast<void*>(&matrixInsertAdd);
    const int errJ = jac_fcn(0.0, x, x_dot, J, user_data, tmp1, tmp2, tmp3);
    realtype* J_data = SUNDenseMatrix_Data(J);

    // Should be: -4 * 1 * 5 = -20
    std::cout << "J =\n";
    for(unsigned int row=0; row<1; ++row){
      for(unsigned int col=0; col<1; ++col){
        std::cout << J_data[col*1 + row] << " ";
        if (col == 0){
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
}