#include "nvector_eigen_dense.h"
#include "Eigen/Dense"
#include <iostream>

using Vector = Eigen::Matrix<realtype, Eigen::Dynamic, 1>;

int main ()
{
  N_Vector x = NanoSim::N_VNew_Eigen<realtype>(2);
  auto x_vec = static_cast<Vector*>(x->content);
  *x_vec << 2, 4;


  N_Vector y = NanoSim::N_VNew_Eigen<realtype>(2);

  x->ops->nvinv(x,y);

  auto y_vec = static_cast<Vector*>(y->content);
  std::cout << *y_vec << std::endl;

  x->ops->nvdestroy(x);
  y->ops->nvdestroy(y);
}