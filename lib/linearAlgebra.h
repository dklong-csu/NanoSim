#ifndef NANOSIM_LINEARALGEBRA_H
#define NANOSIM_LINEARALGEBRA_H

#include <sundials/sundials_nvector.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_types.h>

#include <cassert>


namespace NanoSim{
  template <typename Real>
  class abstractLinearAlgebraOperations {
    public:
    virtual void vectorInsertAdd(N_Vector x, Real value, sunindextype i) = 0;

    virtual Real vectorGetValue(N_Vector x, sunindextype i) = 0;

    virtual void vectorInsert(N_Vector x, Real value, sunindextype i) = 0;

    void stackVectors(N_Vector top, N_Vector bottom, N_Vector combo){
      // Make sure vector dimensions are fine
      const sunindextype dim_top = N_VGetLength(top);
      const sunindextype dim_bottom = N_VGetLength(bottom);
      const sunindextype dim_combo = N_VGetLength(combo);
      assert(dim_top + dim_bottom == dim_combo);

      for (sunindextype i=0; i<dim_top; ++i){
        vectorInsert(combo, vectorGetValue(top, i), i);
      }

      for (sunindextype i=0; i<dim_bottom; ++i){
        vectorInsert(combo, vectorGetValue(bottom, i), dim_top + i);
      }
    }

    void getSubvector(N_Vector full_vec, N_Vector sub_vec, sunindextype idx0){
      // Make sure vector dimensions are fine
      const sunindextype dim_full = N_VGetLength(full_vec);
      const sunindextype dim_sub  = N_VGetLength(sub_vec);
      assert(dim_full >= dim_sub);
      assert(dim_sub + idx0 <= dim_full);

      for (sunindextype i=0; i<dim_sub; ++i){
        vectorInsert(sub_vec, vectorGetValue(full_vec, idx0+i), i);
      }
    }

    virtual N_Vector createNewVector(sunindextype N) = 0;

    virtual SUNMatrix createNewMatrix(sunindextype N, sunindextype M) = 0;

    virtual void matrixInsertAdd(SUNMatrix A, Real value, sunindextype row, sunindextype col) = 0;

    virtual Real matrixGetValue(SUNMatrix A, sunindextype row, sunindextype col) = 0;

    virtual void matrixInsert(SUNMatrix A, Real value, sunindextype row, sunindextype col) = 0;

    void matrixInsertColumn(SUNMatrix A, N_Vector x, sunindextype c){
      const sunindextype n_rows = N_VGetLength(x);
      for (sunindextype r=0; r<n_rows; ++r){
        matrixInsert(A, vectorGetValue(x,r), r, c);
      }
    }

    void matrixGetColumn(SUNMatrix A, N_Vector x, sunindextype c){
      const sunindextype n_rows = N_VGetLength(x);
      for (sunindextype r=0; r<n_rows; ++r){
        vectorInsert(x, matrixGetValue(A, r, c), r);
      }
    }

    void getSubmatrix(SUNMatrix full_mat, SUNMatrix sub_mat, sunindextype row0, sunindextype col0, sunindextype row1, sunindextype col1){
      assert(row0 <= row1);
      assert(col0 <= col1);
      const sunindextype R = row1 - row0 + 1;
      const sunindextype C = col1 - col0 + 1;

      for (sunindextype r=0; r<R; ++r){
        for (sunindextype c=0; c<C; ++c){
          matrixInsert(sub_mat, matrixGetValue(full_mat, r+row0, c+col0), r, c);
        }
      }
    }
  };
}
#endif