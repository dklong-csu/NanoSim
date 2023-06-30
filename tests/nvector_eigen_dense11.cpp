#include "nvector_eigen_dense.h"
#include "Eigen/Dense"
#include <iostream>

using Vector = Eigen::Matrix<realtype, Eigen::Dynamic, 1>;

int main ()
{
  N_Vector v = NanoSim::N_VNew_Eigen<realtype>(2);

  sunindextype n = v->ops->nvgetlength(v);
  std::cout << n << std::endl;

  v->ops->nvdestroy(v);
}