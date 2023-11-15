#include "sunmatrixEigen.h"
#include "Eigen/Dense"
#include <iostream>

using Matrix = Eigen::Matrix<realtype,Eigen::Dynamic, Eigen::Dynamic>;

int main ()
{
  SUNMatrix A = NanoSim::SUNMatNewEigenMatrix<realtype, Matrix>(2,2);
  auto A_matrix = static_cast<Matrix*>(A->content);
  A_matrix->operator()(0,0) = 1;
  A_matrix->operator()(1,0) = 2;

  {
    auto mat = *static_cast<Matrix*>(A->content);

    auto rows = mat.rows();
    auto cols = mat.cols();

    for (unsigned int i=0; i<rows; ++i)
      for (unsigned int j=0; j<cols; ++j)
      {
        std::cout << mat.coeff(i,j);
        if (j < cols-1)
          std::cout << ", ";
        else
          std::cout << std::endl;
      }
  }



  SUNMatrix B = NanoSim::SUNMatNewEigenMatrix<realtype, Matrix>(2,2);

  

  // test the matrix with initialized memory
  {
    auto ierr = A->ops->copy(A, B);
    std::cout << ierr << std::endl;

    auto mat = *static_cast<Matrix*>(B->content);

    auto rows = mat.rows();
    auto cols = mat.cols();

    for (unsigned int i=0; i<rows; ++i)
      for (unsigned int j=0; j<cols; ++j)
      {
        std::cout << mat.coeff(i,j);
        if (j < cols-1)
          std::cout << ", ";
        else
          std::cout << std::endl;
      }
  }

  A->ops->destroy(A);
  B->ops->destroy(B);
  
}