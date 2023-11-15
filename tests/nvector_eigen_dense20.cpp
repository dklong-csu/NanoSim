#include "nvectorEigen.h"
#include "Eigen/Dense"
#include <iostream>

using Vector = Eigen::Matrix<realtype, Eigen::Dynamic, 1>;

int main ()
{
  N_Vector x = NanoSim::N_VNew_Eigen<realtype>(2);
  auto x_vec = static_cast<Vector*>(x->content);
  *x_vec << 2, 3;

  realtype c = 5;

  N_Vector y = NanoSim::N_VNew_Eigen<realtype>(2);

  x->ops->nvscale(c, x, y);
  auto y_vec = static_cast<Vector*>(y->content);
  std::cout << *y_vec << std::endl;

  x->ops->nvdestroy(x);
  y->ops->nvdestroy(y);
}