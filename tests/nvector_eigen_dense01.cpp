#include "nvector_eigen_dense.h"
#include "Eigen/Dense"
#include <iostream>

using Vector = Eigen::Matrix<realtype, Eigen::Dynamic, 1>;

int main ()
{
  N_Vector v = NanoSim::N_VNew_Eigen<realtype>(2);
  auto v_vec = static_cast<Vector*>(v->content);
  *v_vec << -1, 2;

  N_Vector x = v->ops->nvclone(v);
  v->ops->nvabs(v, x);

  auto vec = *(static_cast<Vector *>(x->content));
  std::cout << vec << std::endl;
  v->ops->nvdestroy(v);
  x->ops->nvdestroy(x);
}