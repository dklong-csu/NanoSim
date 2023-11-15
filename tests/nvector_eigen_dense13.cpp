#include "nvectorEigen.h"
#include "Eigen/Dense"
#include <iostream>

using Vector = Eigen::Matrix<realtype, Eigen::Dynamic, 1>;

int main ()
{
  N_Vector x = NanoSim::N_VNew_Eigen<realtype>(2);
  auto x_vec = static_cast<Vector*>(x->content);
  *x_vec << 2, 4;


  N_Vector y = NanoSim::N_VNew_Eigen<realtype>(2);

  // pass the test
  auto result = x->ops->nvinvtest(x,y);
  std::cout << result << std::endl;
  std::cout << *static_cast<Vector*>(y->content) << std::endl;

  // fail the test
  *x_vec << 0, 4;
  result = x->ops->nvinvtest(x,y);
  std::cout << result << std::endl;

  x->ops->nvdestroy(x);
  y->ops->nvdestroy(y);
}