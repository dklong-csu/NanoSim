#include "nvectorEigen.h"
#include "Eigen/Dense"
#include <iostream>

using Vector = Eigen::Matrix<realtype, Eigen::Dynamic, 1>;

int main ()
{
  realtype c[] {2, 3, 4};


  N_Vector x = NanoSim::N_VNew_Eigen<realtype>(2);
  auto x_vec = static_cast<Vector*>(x->content);
  *x_vec << 1, 2;


  N_Vector y = NanoSim::N_VNew_Eigen<realtype>(2);
  auto y_vec = static_cast<Vector*>(y->content);
  *y_vec << 2, 3;


  N_Vector z = NanoSim::N_VNew_Eigen<realtype>(2);
  auto z_vec = static_cast<Vector*>(z->content);
  *z_vec << 3, 4;

  N_Vector X [3] = {x, y, z};


  N_Vector r = NanoSim::N_VNew_Eigen<realtype>(2);


  auto result = x->ops->nvlinearcombination(3, c, X, r);
  std::cout << result << std::endl;
  std::cout << *(static_cast<Vector*>(r->content)) << std::endl;

  result = x->ops->nvlinearcombination(0, c, X, r);
  std::cout << result << std::endl;

  x->ops->nvdestroy(x);
  y->ops->nvdestroy(y);
  z->ops->nvdestroy(z);
  r->ops->nvdestroy(r);
}