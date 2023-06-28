#ifndef NANOSIM_OPERATIONSSUNDENSEMATRIX_H
#define NANOSIM_OPERATIONSSUNDENSEMATRIX_H

#include <sundials/sundials_types.h>
#include <sunmatrix/sunmatrix_dense.h>

namespace NanoSim {
  template<typename Real>
  int
  matrixInsertAdd(SUNMatrix A, sunindextype row, sunindextype col, Real value){
    Real* A_data = SUNDenseMatrix_Data(A);
    const sunindextype n_rows = SUNDenseMatrix_Rows(A);
    A_data[col*n_rows + row] += value;
    return 0;
  }
}

#endif