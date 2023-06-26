#include "operationsSUNDenseMatrix.h"

#include <sundials/sundials_types.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <iostream>

int main(){
  const unsigned int length = 3;
  SUNMatrix A = SUNDenseMatrix(length,length);
  realtype* A_data = SUNDenseMatrix_Data(A);

  // Start matrix at zero
  A->ops->zero(A);
  std::cout << "A =\n";
  for(unsigned int row=0; row<length; ++row){
    for(unsigned int col=0; col<length; ++col){
      std::cout << A_data[col*length + row] << " ";
      if (col == length-1){
        std::cout << "\n";
      }
    }
  }

  // Insert some elements
  for (unsigned int i=0;i<length;++i){
    for (unsigned int j=0;j<length;++j){
      NanoSim::matrixInsertAdd(A,i,j, 1.0*(i+1)*(j+1));
    }
  }
  std::cout << "A =\n";
  for(unsigned int row=0; row<length; ++row){
    for(unsigned int col=0; col<length; ++col){
      std::cout << A_data[col*length + row] << " ";
      if (col == length-1){
        std::cout << "\n";
      }
    }
  }
}