#ifndef NANOSIM_SUNMATRIXEIGENDENSE_H
#define NANOSIM_SUNMATRIXEIGENDENSE_H

#include <sundials/sundials_types.h>
#include <sundials/sundials_matrix.h>

namespace NanoSim{
  SUNMatrix_ID
  sunmat_eigen_getid(SUNMatrix A)
  {
    return SUNMATRIX_CUSTOM;
  }



  void
  sunmat_eigen_copy_ops(SUNMatrix A, SUNMatrix B)
  {
    B->ops->getid = A->ops->getid;
    B->ops->clone = A->ops->clone;
    B->ops->destroy = A->ops->destroy;
    B->ops->zero = A->ops->zero;
    B->ops->copy = A->ops->copy;
    B->ops->scaleaddi = A->ops->scaleaddi;
  }



  template<typename EigenMatrix>
  SUNMatrix
  sunmat_eigen_clone(SUNMatrix A)
  {
    sundials::Context cntxt;
    SUNMatrix B = SUNMatNewEmpty(cntxt);
    sunmat_eigen_copy_ops(A, B);

    auto A_ptr = static_cast< EigenMatrix *>(A->content);
    auto n_cols = A_ptr->cols();
    auto n_rows = A_ptr->rows();
    auto cloned = new EigenMatrix(n_rows,n_cols);
    B->content = (void*)cloned;
    return B;
  }



  template<typename EigenMatrix>
  void
  sunmat_eigen_destroy(SUNMatrix A)
  {
    if (A->content != nullptr)
    {
      auto content = static_cast< EigenMatrix *>(A->content);
      delete content;
      A->content = nullptr;
    }

    SUNMatFreeEmpty(A);
  }



  template<typename EigenMatrix>
  int
  sunmat_eigen_zero(SUNMatrix A)
  {
    auto A_ptr = static_cast< EigenMatrix*>(A->content);
    A_ptr->setZero();
    return SUNMAT_SUCCESS;
  }



  template<typename EigenMatrix>
  int
  sunmat_eigen_copy(SUNMatrix A, SUNMatrix B)
  {
    auto A_ptr = static_cast< EigenMatrix*>(A->content);
    if (B->content == nullptr)
    {
      auto rows = A_ptr->rows();
      auto cols = A_ptr->cols();
      EigenMatrix* mat = new EigenMatrix(rows, cols);
      B->content = (void*)mat;
    }
    auto B_ptr = static_cast<EigenMatrix*>(B->content);

    *B_ptr = *A_ptr;

    return SUNMAT_SUCCESS;
  }



  template<typename Real, typename EigenMatrix>
  int
  sunmat_eigen_scaleaddI(Real c, SUNMatrix A)
  {
    auto A_ptr = static_cast<EigenMatrix*>(A->content);

    *A_ptr = c * (*A_ptr);
    auto n = A_ptr->cols();
    for (unsigned int i=0; i<n; ++i)
    {
      A_ptr->coeffRef(i,i) += 1.;
    }

    return SUNMAT_SUCCESS;
  }


  template<typename Real, typename EigenMatrix>
  SUNMatrix
  SUNMatNewEigenMatrix(sunindextype rows, sunindextype cols){
    sundials::Context cntxt;
    SUNMatrix A = SUNMatNewEmpty(cntxt);

    EigenMatrix* mat = new EigenMatrix(rows, cols);
    mat->setZero();
    A->content = (void*)mat;

    A->ops->getid = sunmat_eigen_getid;
    A->ops->clone = sunmat_eigen_clone<EigenMatrix>;
    A->ops->destroy = sunmat_eigen_destroy<EigenMatrix>;
    A->ops->zero = sunmat_eigen_zero<EigenMatrix>;
    A->ops->copy = sunmat_eigen_copy<EigenMatrix>;
    A->ops->scaleaddi = sunmat_eigen_scaleaddI<Real, EigenMatrix>; 

    return A;
  }
}

#endif