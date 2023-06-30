#include "nvector_eigen_dense.h"
#include "Eigen/Dense"
#include <iostream>

using Vector = Eigen::Matrix<realtype, Eigen::Dynamic, 1>;

int main ()
{
  N_Vector v = NanoSim::N_VNew_Eigen<realtype>(2);

  N_Vector x = v->ops->nvclone(v);
  // nvclone is only supposed to allocate storage, not modify the values. So the output should be a vector of zeros
  // since Eigen will default those values to zero.
  std::cout << *static_cast<Vector*>(x->content) << std::endl;
  v->ops->nvdestroy(v);
  x->ops->nvdestroy(x);
}