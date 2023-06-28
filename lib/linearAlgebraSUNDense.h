#ifndef NANOSIM_LINEARALGEBRASUNDENSE_H
#define NANOSIM_LINEARALGEBRASUNDENSE_H

#include "linearAlgebra.h"
#include <sunmatrix/sunmatrix_dense.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>


namespace NanoSim{
  template <typename Real>
  class sunDenseLinearAlgebraOperations : public abstractLinearAlgebraOperations<Real> {
    public:
    void vectorInsertAdd(N_Vector x, Real value, sunindextype i);

    Real vectorGetValue(N_Vector x, sunindextype i);

    void vectorInsert(N_Vector x, Real value, sunindextype i);

    N_Vector createNewVector(sunindextype N);

    SUNMatrix createNewMatrix(sunindextype N, sunindextype M);

    void matrixInsertAdd(SUNMatrix A, Real value, sunindextype row, sunindextype col);

    Real matrixGetValue(SUNMatrix A, sunindextype row, sunindextype col);
  };

  template<typename Real>
  void
  sunDenseLinearAlgebraOperations<Real>::vectorInsertAdd(N_Vector x, Real value, sunindextype i){
    Real *x_data = x->ops->nvgetarraypointer(x);
    x_data[i] += value;
  }



  template<typename Real>
  Real
  sunDenseLinearAlgebraOperations<Real>::vectorGetValue(N_Vector x, sunindextype i){
    const Real *x_data = x->ops->nvgetarraypointer(x);
    return x_data[i];
  }



  template<typename Real>
  void
  sunDenseLinearAlgebraOperations<Real>::vectorInsert(N_Vector x, Real value, sunindextype i){
    Real *x_data = x->ops->nvgetarraypointer(x);
    x_data[i] = value;
  }



  template<typename Real>
  N_Vector
  sunDenseLinearAlgebraOperations<Real>::createNewVector(sunindextype N){
    return N_VNew_Serial(N);
  }



  template<typename Real>
  SUNMatrix
  sunDenseLinearAlgebraOperations<Real>::createNewMatrix(sunindextype N, sunindextype M){
    return SUNDenseMatrix(N,M);
  }



  template<typename Real>
  void
  sunDenseLinearAlgebraOperations<Real>::matrixInsertAdd(SUNMatrix A, Real value, sunindextype row, sunindextype col){
    Real* A_data = SUNDenseMatrix_Data(A);
    const sunindextype n_rows = SUNDenseMatrix_Rows(A);
    A_data[col*n_rows + row] += value;
  }



  template<typename Real>
  Real
  sunDenseLinearAlgebraOperations<Real>::matrixGetValue(SUNMatrix A, sunindextype row, sunindextype col){
    const Real* A_data = SUNDenseMatrix_Data(A);
    const sunindextype n_rows = SUNDenseMatrix_Rows(A);
    return A_data[col*n_rows + row];
  }
}
#endif