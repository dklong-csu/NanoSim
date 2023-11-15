#ifndef NANOSIM_LINEARALGEBRAEIGEN_H
#define NANOSIM_LINEARALGEBRAEIGEN_H

#include "linearAlgebra.h"
#include "sunmatrixEigen.h"
#include <nvectorEigen.h>


namespace NanoSim{
  template <typename Real, typename EigenMatrix>
  class eigenLinearAlgebraOperations : public abstractLinearAlgebraOperations<Real> {
    public:
    void vectorInsertAdd(N_Vector x, Real value, sunindextype i);

    Real vectorGetValue(N_Vector x, sunindextype i);

    void vectorInsert(N_Vector x, Real value, sunindextype i);

    N_Vector createNewVector(sunindextype N);

    SUNMatrix createNewMatrix(sunindextype N, sunindextype M);

    void matrixInsertAdd(SUNMatrix A, Real value, sunindextype row, sunindextype col);

    Real matrixGetValue(SUNMatrix A, sunindextype row, sunindextype col);

    void matrixInsert(SUNMatrix A, Real value, sunindextype row, sunindextype col);
  };



  template<typename Real, typename EigenMatrix>
  void
  eigenLinearAlgebraOperations<Real, EigenMatrix>::vectorInsertAdd(N_Vector x, Real value, sunindextype i){
    Eigen::Vector<Real,Eigen::Dynamic>* x_data = static_cast< Eigen::Vector<Real,Eigen::Dynamic>* >(x->content);
    x_data->coeffRef(i) += value;
  }



  template<typename Real, typename EigenMatrix>
  Real
  eigenLinearAlgebraOperations<Real, EigenMatrix>::vectorGetValue(N_Vector x, sunindextype i){
    Eigen::Vector<Real,Eigen::Dynamic>* x_data = static_cast< Eigen::Vector<Real,Eigen::Dynamic>* >(x->content);
    return x_data->coeff(i);
  }



  template<typename Real, typename EigenMatrix>
  void
  eigenLinearAlgebraOperations<Real, EigenMatrix>::vectorInsert(N_Vector x, Real value, sunindextype i){
    Eigen::Vector<Real,Eigen::Dynamic>* x_data 
      = static_cast< Eigen::Vector<Real,Eigen::Dynamic>* >(x->content);
    x_data->coeffRef(i) = value;
  }



  template<typename Real, typename EigenMatrix>
  N_Vector
  eigenLinearAlgebraOperations<Real, EigenMatrix>::createNewVector(sunindextype N){
    return N_VNew_Eigen<Real>(N);
  }



  template<typename Real, typename EigenMatrix>
  SUNMatrix
  eigenLinearAlgebraOperations<Real, EigenMatrix>::createNewMatrix(sunindextype N, sunindextype M){
    return SUNMatNewEigenMatrix<Real, EigenMatrix>(N,M); 
  }



  template<typename Real, typename EigenMatrix>
  void
  eigenLinearAlgebraOperations<Real, EigenMatrix>::matrixInsertAdd(SUNMatrix A, Real value, sunindextype row, sunindextype col){
    EigenMatrix* A_data = static_cast< EigenMatrix* >(A->content);
    A_data->coeffRef(row, col) += value;
  }



  template<typename Real, typename EigenMatrix>
  Real
  eigenLinearAlgebraOperations<Real, EigenMatrix>::matrixGetValue(SUNMatrix A, sunindextype row, sunindextype col){
    EigenMatrix* A_data = static_cast< EigenMatrix* >(A->content);
    return A_data->coeff(row, col);
  }



  template<typename Real, typename EigenMatrix>
  void
  eigenLinearAlgebraOperations<Real, EigenMatrix>::matrixInsert(SUNMatrix A, Real value, sunindextype row, sunindextype col){
    EigenMatrix* A_data = static_cast< EigenMatrix* >(A->content);
    A_data->coeffRef(row, col) = value;
  }
}
#endif