#include "nvectorEigen.h"
#include "Eigen/Dense"
#include <iostream>

using Vector = Eigen::Matrix<realtype, Eigen::Dynamic, 1>;

int main ()
{
  N_Vector v = NanoSim::N_VNew_Eigen<realtype>(2);

  N_Vector x = v->ops->nvcloneempty(v);
  // Check to make sure ops pointers are correct and that the content field of x is nullptr
  std::cout << std::boolalpha
  << (x->ops->nvgetlength == v->ops->nvgetlength) << std::endl
  << (x->ops->nvclone == v->ops->nvclone) << std::endl
  << (x->ops->nvcloneempty == v->ops->nvcloneempty) << std::endl
  << (x->ops->nvdestroy == v->ops->nvdestroy) << std::endl
  << (x->ops->nvspace == v->ops->nvspace) << std::endl
  << (x->ops->nvgetarraypointer == v->ops->nvgetarraypointer) << std::endl
  << (x->ops->nvsetarraypointer == v->ops->nvsetarraypointer) << std::endl
  << (x->ops->nvlinearsum == v->ops->nvlinearsum) << std::endl
  << (x->ops->nvconst == v->ops->nvconst) << std::endl
  << (x->ops->nvprod == v->ops->nvprod) << std::endl
  << (x->ops->nvdiv == v->ops->nvdiv) << std::endl
  << (x->ops->nvscale == v->ops->nvscale) << std::endl
  << (x->ops->nvabs == v->ops->nvabs) << std::endl
  << (x->ops->nvinv == v->ops->nvinv) << std::endl
  << (x->ops->nvaddconst == v->ops->nvaddconst) << std::endl
  << (x->ops->nvmaxnorm == v->ops->nvmaxnorm) << std::endl
  << (x->ops->nvwrmsnorm == v->ops->nvwrmsnorm) << std::endl
  << (x->ops->nvmin == v->ops->nvmin) << std::endl
  << (x->ops->nvminquotient == v->ops->nvminquotient) << std::endl
  << (x->ops->nvconstrmask == v->ops->nvconstrmask) << std::endl
  << (x->ops->nvcompare == v->ops->nvcompare) << std::endl
  << (x->ops->nvinvtest == v->ops->nvinvtest) << std::endl
  << (x->ops->nvlinearcombination == v->ops->nvlinearcombination) << std::endl
  << (x->ops->nvscaleaddmulti == v->ops->nvscaleaddmulti) << std::endl
  << (x->ops->nvdotprodmulti == v->ops->nvdotprodmulti) << std::endl
  << (x->ops->nvscalevectorarray == v->ops->nvscalevectorarray) << std::endl
  << (x->content == nullptr)
  << std::endl;

  v->ops->nvdestroy(v);
  x->ops->nvdestroy(x);
}