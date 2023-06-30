#include "nvector_eigen_dense.h"
#include "Eigen/Dense"
#include <iostream>

using Vector = Eigen::Matrix<realtype, Eigen::Dynamic, 1>;

int main ()
{
  N_Vector x = NanoSim::N_VNew_Eigen<realtype>(2);
  auto x_vec = static_cast<Vector*>(x->content);
  *x_vec << 2, 3;

  N_Vector y = NanoSim::N_VNew_Eigen<realtype>(2);
  auto y_vec = static_cast<Vector*>(y->content);
  *y_vec << 4, 5;

  N_Vector z = NanoSim::N_VNew_Eigen<realtype>(2);

  // z = x .* y = (2*4, 3*5) = (8, 15)
  x->ops->nvprod(x,y,z);

  auto z_vec = static_cast<Vector*>(z->content);
  std::cout << *z_vec << std::endl;

  x->ops->nvdestroy(x);
  y->ops->nvdestroy(y);
  z->ops->nvdestroy(z);
}