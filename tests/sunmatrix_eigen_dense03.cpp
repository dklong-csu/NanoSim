// A = 1,0;2,0. c = 3.
// result =
// 4, 0
// 6, 1

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

  // make sure the matrix has the correct initial coefficients
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

  realtype c = 3;
  std::cout << c << std::endl;


  A->ops->scaleaddi(c, A);
  // make sure the matrix is now scaled and has the identity added
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

  A->ops->destroy(A);
}