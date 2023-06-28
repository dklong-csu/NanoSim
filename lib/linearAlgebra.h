#ifndef NANOSIM_LINEARALGEBRA_H
#define NANOSIM_LINEARALGEBRA_H

#include <sundials/sundials_nvector.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_types.h>


namespace NanoSim{
  template <typename Real>
  class abstractLinearAlgebraOperations {
    public:
    virtual void vectorInsertAdd(N_Vector x, Real value, sunindextype i) = 0;

    virtual Real vectorGetValue(N_Vector x, sunindextype i) = 0;

    virtual void vectorInsert(N_Vector x, Real value, sunindextype i) = 0;

    virtual N_Vector createNewVector(sunindextype N) = 0;

    virtual SUNMatrix createNewMatrix(sunindextype N, sunindextype M) = 0;

    virtual void matrixInsertAdd(SUNMatrix A, Real value, sunindextype row, sunindextype col) = 0;

    virtual Real matrixGetValue(SUNMatrix A, sunindextype row, sunindextype col) = 0;
  };
}
#endif