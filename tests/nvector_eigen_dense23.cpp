#include "nvector_eigen_dense.h"
#include "Eigen/Dense"
#include <iostream>

using Vector = Eigen::Matrix<realtype, Eigen::Dynamic, 1>;

int main ()
{
  N_Vector x = NanoSim::N_VNew_Eigen<realtype>(4);
  auto x_vec = static_cast<Vector*>(x->content);
  *x_vec << 1, 2, 4, 8;

  N_Vector w = NanoSim::N_VNew_Eigen<realtype>(4);
  auto w_vec = static_cast<Vector*>(w->content);
  *w_vec << 2, 1, .5, .25;


  // Result should be 2
  auto result = x->ops->nvwrmsnorm(x,w);
  std::cout << result << std::endl;

  x->ops->nvdestroy(x);
  w->ops->nvdestroy(w);
}