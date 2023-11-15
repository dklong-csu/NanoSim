#include "nvectorEigen.h"
#include "Eigen/Dense"
#include <iostream>

using Vector = Eigen::Matrix<realtype, Eigen::Dynamic, 1>;

int main ()
{
  realtype c = 1;

  N_Vector x = NanoSim::N_VNew_Eigen<realtype>(3);
  auto x_vec = static_cast<Vector*>(x->content);
  *x_vec << 0, 1, 2;

  N_Vector y = NanoSim::N_VNew_Eigen<realtype>(3);

  x->ops->nvcompare(c,x,y);
  std::cout << *static_cast<Vector*>(y->content) << std::endl;

  x->ops->nvdestroy(x);
  y->ops->nvdestroy(y);
}