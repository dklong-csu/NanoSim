#include "nvector_eigen_dense.h"
#include "Eigen/Dense"
#include <iostream>

using Vector = Eigen::Matrix<realtype, Eigen::Dynamic, 1>;

int main ()
{
  // constants to scale by
  realtype c [2] = {2,3};


  // vectors to be scaled
  N_Vector y = NanoSim::N_VNew_Eigen<realtype>(2);
  auto y_vec = static_cast<Vector*>(y->content);
  *y_vec << 1,2;


  N_Vector z = NanoSim::N_VNew_Eigen<realtype>(2);
  auto z_vec = static_cast<Vector*>(z->content);
  *z_vec << 3,4;


  N_Vector X [2] = {y,z};


  // vectors to store results
  N_Vector a = NanoSim::N_VNew_Eigen<realtype>(2);

  N_Vector b = NanoSim::N_VNew_Eigen<realtype>(2);

  N_Vector C [2] = {a,b};

  // test
  auto result = y->ops->nvscalevectorarray(2, c, X, C);
  std::cout << result << std::endl;
  std::cout << *(static_cast<Vector*>(a->content)) << std::endl;
  std::cout << *(static_cast<Vector*>(b->content)) << std::endl;

  result = y->ops->nvscalevectorarray(0, c, X, C);
  std::cout << result << std::endl;

  y->ops->nvdestroy(y);
  z->ops->nvdestroy(z);
  a->ops->nvdestroy(a);
  b->ops->nvdestroy(b);
}