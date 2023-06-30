#ifndef NANOSIM_LINEARALGEBRAEIGENDENSE_H
#define NANOSIM_LINEARALGEBRAEIGENDENSE_H

#include "linearAlgebra.h"
#include <Eigen/Dense>


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

    void matrixInsert(SUNMatrix A, Real value, sunindextype row, sunindextype col);
  };

  template<typename Real>
  void
  sunDenseLinearAlgebraOperations<Real>::vectorInsertAdd(N_Vector x, Real value, sunindextype i){
    Eigen::Vector<Real,Eigen::Dynamic>* x_data = static_cast< Eigen::Vector<Real,Eigen::Dynamic>* >(x->content);
    x_data->operator()(i) += value;
  }



  template<typename Real>
  Real
  sunDenseLinearAlgebraOperations<Real>::vectorGetValue(N_Vector x, sunindextype i){
    Eigen::Vector<Real,Eigen::Dynamic>* x_data = static_cast< Eigen::Vector<Real,Eigen::Dynamic>* >(x->content);
    return x_data->operator()(i);
  }



  template<typename Real>
  void
  sunDenseLinearAlgebraOperations<Real>::vectorInsert(N_Vector x, Real value, sunindextype i){
    Eigen::Vector<Real,Eigen::Dynamic>* x_data 
      = static_cast< Eigen::Vector<Real,Eigen::Dynamic>* >(x->content);
    x_data->operator()(i) = value;
  }



  template<typename Real>
  N_Vector
  sunDenseLinearAlgebraOperations<Real>::createNewVector(sunindextype N){
    return N_VNew_Eigen(N); // TODO
  }



  template<typename Real>
  SUNMatrix
  sunDenseLinearAlgebraOperations<Real>::createNewMatrix(sunindextype N, sunindextype M){
    return SUNDenseEigenMatrix(N,M); // TODO
  }



  template<typename Real>
  void
  sunDenseLinearAlgebraOperations<Real>::matrixInsertAdd(SUNMatrix A, Real value, sunindextype row, sunindextype col){
    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>* A_data
      = static_cast< Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>* >(A->content);
    A_data->operator()(row, col) += value;
  }



  template<typename Real>
  Real
  sunDenseLinearAlgebraOperations<Real>::matrixGetValue(SUNMatrix A, sunindextype row, sunindextype col){
    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>* A_data
      = static_cast< Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>* >(A->content);
    return A_data->operator()(row, col);
  }



  template<typename Real>
  void
  sunDenseLinearAlgebraOperations<Real>::matrixInsert(SUNMatrix A, Real value, sunindextype row, sunindextype col){
    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>* A_data
      = static_cast< Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>* >(A->content);
    A_data->operator()(row, col) = value;
  }
}
#endif