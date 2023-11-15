#ifndef NANOSIM_NVECTOREIGENDENSE_H
#define NANOSIM_NVECTOREIGENDENSE_H

#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>
#include <Eigen/Dense>

#include <cassert>



namespace NanoSim{
  template<typename Real>
  sunindextype
  nv_eigen_getlength(N_Vector v){
    return ( static_cast<Eigen::Vector<Real,Eigen::Dynamic>*>(v->content) )->size();
  }



  void  
  nv_eigen_copy_ops_pointers(N_Vector x, N_Vector y)
  {
    // Pointers are either NULL or point to functions, so there is not a concern about ownership.
    x->ops->nvgetlength         = y->ops->nvgetlength;
    x->ops->nvclone             = y->ops->nvclone;
    x->ops->nvcloneempty        = y->ops->nvcloneempty;
    x->ops->nvdestroy           = y->ops->nvdestroy;
    x->ops->nvspace             = y->ops->nvspace;
    x->ops->nvgetarraypointer   = y->ops->nvgetarraypointer;
    x->ops->nvsetarraypointer   = y->ops->nvsetarraypointer;
    x->ops->nvlinearsum         = y->ops->nvlinearsum;
    x->ops->nvconst             = y->ops->nvconst;
    x->ops->nvprod              = y->ops->nvprod;
    x->ops->nvdiv               = y->ops->nvdiv;
    x->ops->nvscale             = y->ops->nvscale;
    x->ops->nvabs               = y->ops->nvabs;
    x->ops->nvinv               = y->ops->nvinv;
    x->ops->nvaddconst          = y->ops->nvaddconst;
    x->ops->nvmaxnorm           = y->ops->nvmaxnorm;
    x->ops->nvwrmsnorm          = y->ops->nvwrmsnorm;
    x->ops->nvmin               = y->ops->nvmin;
    x->ops->nvminquotient       = y->ops->nvminquotient;
    x->ops->nvconstrmask        = y->ops->nvconstrmask;
    x->ops->nvcompare           = y->ops->nvcompare;
    x->ops->nvinvtest           = y->ops->nvinvtest;
    x->ops->nvlinearcombination = y->ops->nvlinearcombination;
    x->ops->nvscaleaddmulti     = y->ops->nvscaleaddmulti;
    x->ops->nvdotprodmulti      = y->ops->nvdotprodmulti;
    x->ops->nvscalevectorarray  = y->ops->nvscalevectorarray;
  }



  template <typename Real>
  N_Vector
  nv_eigen_clone( N_Vector w)
  {
    N_Vector v = N_VCloneEmpty(w);
    auto cloned = new Eigen::Vector<Real, Eigen::Dynamic>(w->ops->nvgetlength(w));
    cloned->setZero();
    v->content = cloned;
    return v;
  }



  N_Vector
  nv_eigen_clone_empty( N_Vector w)
  {
    sundials::Context cntxt;
    N_Vector v = N_VNewEmpty(cntxt);

    N_VCopyOps(w, v);
    v->content = nullptr;

    return v;
  }



  template<typename Real>
  void
  nv_eigen_destroy( N_Vector v)
  {

    if (v->content != nullptr)
    {
      Eigen::Vector<Real, Eigen::Dynamic> *content = static_cast<Eigen::Vector<Real, Eigen::Dynamic> *>(v->content);
      delete content;
      v->content = nullptr;
    }

    N_VFreeEmpty(v);
    v = nullptr;
  }



  void
  nv_eigen_space(N_Vector v, sunindextype* lrw, sunindextype* liw)
  {
    *lrw = v->ops->nvgetlength(v);
    *liw = 0;
  }



  template<typename Real>
  Real*
  nv_eigen_getarraypointer(N_Vector v)
  {
    Eigen::Vector<Real, Eigen::Dynamic>* vector = static_cast<Eigen::Vector<Real, Eigen::Dynamic>*>(v->content);
    return vector->data();
  }



  template<typename Real>
  void
  nv_eigen_setarraypointer(Real *v_data, N_Vector v)
  {
    assert(false);
  }



  template<typename Real>
  void
  nv_eigen_linearsum(Real a, N_Vector x, Real b, N_Vector y, N_Vector z)
  {
    auto x_vec = static_cast<Eigen::Vector<Real, Eigen::Dynamic>*>(x->content);
    auto y_vec = static_cast<Eigen::Vector<Real, Eigen::Dynamic>*>(y->content);
    auto z_vec = static_cast<Eigen::Vector<Real, Eigen::Dynamic>*>(z->content);

    *z_vec = a * (*x_vec) + b * (*y_vec);
  }



  template<typename Real>
  void
  nv_eigen_const(Real c, N_Vector z)
  {
    auto z_vec = static_cast<Eigen::Vector<Real, Eigen::Dynamic>*>(z->content);
    z_vec->setConstant(z_vec->size(), c);
  }



  template<typename Real>
  void
  nv_eigen_prod(N_Vector x, N_Vector y, N_Vector z)
  {
    auto x_vec = static_cast<Eigen::Vector<Real, Eigen::Dynamic>*>(x->content);
    auto y_vec = static_cast<Eigen::Vector<Real, Eigen::Dynamic>*>(y->content);
    auto z_vec = static_cast<Eigen::Vector<Real, Eigen::Dynamic>*>(z->content);

    *z_vec = ( x_vec->array() ) * ( y_vec->array() ) ;
  }



  template<typename Real>
  void
  nv_eigen_div(N_Vector x, N_Vector y, N_Vector z)
  {
    auto x_vec = static_cast<Eigen::Vector<Real, Eigen::Dynamic>*>(x->content);
    auto y_vec = static_cast<Eigen::Vector<Real, Eigen::Dynamic>*>(y->content);
    auto z_vec = static_cast<Eigen::Vector<Real, Eigen::Dynamic>*>(z->content);

    *z_vec = ( x_vec->array() ) / ( y_vec->array() ) ;
  }



  template<typename Real>
  void
  nv_eigen_scale(Real c, N_Vector x, N_Vector z)
  {
    auto x_vec = static_cast<Eigen::Vector<Real, Eigen::Dynamic>*>(x->content);
    auto z_vec = static_cast<Eigen::Vector<Real, Eigen::Dynamic>*>(z->content);

    *z_vec = c * (*x_vec);
  }



  template<typename Real>
  void
  nv_eigen_abs(N_Vector x, N_Vector z)
  {
    auto x_vec = static_cast<Eigen::Vector<Real, Eigen::Dynamic>*>(x->content);
    auto z_vec = static_cast<Eigen::Vector<Real, Eigen::Dynamic>*>(z->content);

    *z_vec = (*x_vec).array().abs();
  }



  template<typename Real>
  void
  nv_eigen_inv(N_Vector x, N_Vector z)
  {
    auto x_vec = static_cast<Eigen::Vector<Real, Eigen::Dynamic>*>(x->content);
    auto z_vec = static_cast<Eigen::Vector<Real, Eigen::Dynamic>*>(z->content);

    *z_vec = (*x_vec).array().inverse();
  }



  template<typename Real>
  void
  nv_eigen_addconst(N_Vector x, Real b, N_Vector z)
  {
    auto x_vec = static_cast<Eigen::Vector<Real, Eigen::Dynamic>*>(x->content);
    auto z_vec = static_cast<Eigen::Vector<Real, Eigen::Dynamic>*>(z->content);

    *z_vec = (*x_vec).array() + b;
  }



  template<typename Real>
  Real
  nv_eigen_maxnorm(N_Vector x)
  {
    auto x_vec = static_cast<Eigen::Vector<Real, Eigen::Dynamic>*>(x->content);
    return (*x_vec).array().abs().maxCoeff();
  }



  template<typename Real>
  Real
  nv_eigen_wrms_norm(N_Vector x, N_Vector w)
  {
    auto x_vec = static_cast<Eigen::Vector<Real, Eigen::Dynamic>*>(x->content);
    auto w_vec = static_cast<Eigen::Vector<Real, Eigen::Dynamic>*>(w->content);
    auto n = (*x_vec).size();
    auto sum = ( ( (*x_vec).array() * (*w_vec).array() ).pow(2) ).mean();
    return std::sqrt(sum);
  }



  template<typename Real>
  Real
  nv_eigen_min(N_Vector x)
  {
    auto x_vec = static_cast<Eigen::Vector<Real, Eigen::Dynamic>*>(x->content);
    return (*x_vec).minCoeff();
  }



  template<typename Real>
  Real
  nv_eigen_minquotient(N_Vector num, N_Vector denom)
  {
    auto num_vec = static_cast<Eigen::Vector<Real, Eigen::Dynamic>*>(num->content);
    auto denom_vec = static_cast<Eigen::Vector<Real, Eigen::Dynamic>*>(denom->content);

    Real result = BIG_REAL;
    for (unsigned int i = 0; i<num_vec->size(); ++i)
    {
      if ( (*denom_vec)(i) != 0 )
      {
        auto ratio = (*num_vec)(i) / (*denom_vec)(i);
        if (ratio < result)
        {
          result = ratio;
        }
      }
    }

    return result;
  }



  template<typename Real>
  booleantype
  nv_eigen_vconstrmask(N_Vector c, N_Vector x, N_Vector m)
  {
    auto c_vec = static_cast<Eigen::Vector<Real, Eigen::Dynamic>*>(c->content);
    auto x_vec = static_cast<Eigen::Vector<Real, Eigen::Dynamic>*>(x->content);
    auto m_vec = static_cast<Eigen::Vector<Real, Eigen::Dynamic>*>(m->content);

    booleantype result = SUNTRUE;

    for (unsigned int i=0; i<x->ops->nvgetlength(x); ++i)
    {
      auto constraint = (int)(*c_vec)(i);
      auto val = (*x_vec)(i);
      switch(constraint)
      {
        case 2:
          if (val > 0)
          {
            (*m_vec)(i) = 0;
          }
          else
          {
            (*m_vec)(i) = 1;
            result = SUNFALSE;
          }
          break;
        case 1:
          if (val >= 0)
          {
            (*m_vec)(i) = 0;
          }
          else
          {
            (*m_vec)(i) = 1;
            result = SUNFALSE;
          }
          break;
        case -1:
          if (val <= 0 )
          {
            (*m_vec)(i) = 0;
          }
          else
          {
            (*m_vec)(i) = 1;
            result = SUNFALSE;
          }
          break;
        case -2:
          if (val < 0)
          {
            (*m_vec)(i) = 0;
          }
          else
          {
            (*m_vec)(i) = 1;
            result = SUNFALSE;
          }
          break;
        default:
          (*m_vec)(i) = 0;
      }
    }

    return result;
  }



  template<typename Real>
  void
  nv_eigen_vcompare(Real c, N_Vector x, N_Vector z)
  {
    auto x_vec = static_cast<Eigen::Vector<Real, Eigen::Dynamic>*>(x->content);
    auto z_vec = static_cast<Eigen::Vector<Real, Eigen::Dynamic>*>(z->content);

    for (unsigned int i=0; i<x->ops->nvgetlength(x); ++i)
    {
      auto val = std::abs((*x_vec)(i));
      (val > c) ? (*z_vec)(i) = 1. : (*z_vec)(i) = 0.;
    }
  }



  template<typename Real>
  booleantype
  nv_eigen_vinvtest(N_Vector x, N_Vector z)
  {
    auto x_vec = static_cast<Eigen::Vector<Real, Eigen::Dynamic>*>(x->content);
    auto z_vec = static_cast<Eigen::Vector<Real, Eigen::Dynamic>*>(z->content);

    booleantype result = SUNTRUE;

    for (unsigned int i=0; i<x->ops->nvgetlength(x); ++i)
    {
      auto val = (*x_vec)(i);
      if (val != 0)
        (*z_vec)(i) = 1./val;
      else
        result = SUNFALSE;
    }
    return result;
  }



  template<typename Real>
  int
  nv_eigen_linearcombination(int nv, Real* c, N_Vector* X, N_Vector z)
  {
    // invalid number of vectors
    if (nv < 1) return -1;

    // should have called N_VScale in this case
    if (nv == 1)
    {
      nv_eigen_scale<Real>(c[0], X[0], z);
      return 0;
    }

    // should have called N_VLinearSum
    if (nv == 2)
    {
      nv_eigen_linearsum<Real>(c[0], X[0], c[1], X[1], z);
      return 0;
    }

    // when nv > 2, start with linear sum and then keep adding to z
    nv_eigen_linearsum<Real>(c[0], X[0], c[1], X[1], z);
    for (unsigned int i=2; i<nv; ++i)
      nv_eigen_linearsum<Real>(1, z, c[i], X[i], z);

    return 0;
  }



  template<typename Real>
  int
  nv_eigen_vscaleaddmulti(int nv, Real* c, N_Vector x, N_Vector* Y, N_Vector* Z)
  {
    // invalid number of vectors
    if (nv < 1) return -1;

    for (unsigned int i=0; i<nv; ++i)
    {
      nv_eigen_linearsum<Real>(c[i], x, 1, Y[i], Z[i]);
    }

    return 0;
  }



  template<typename Real>
  int
  nv_eigen_vdotprodmulti(int nv, N_Vector x, N_Vector* Y, Real* d)
  {
    // invalid number of vectors
    if (nv < 1) return -1;

    auto x_vec = static_cast<Eigen::Vector<Real, Eigen::Dynamic>*>(x->content);

    for (unsigned int i=0; i<nv; ++i)
    {
      auto y_vec = static_cast<Eigen::Vector<Real, Eigen::Dynamic>*>(Y[i]->content);
      d[i] = x_vec->dot(*y_vec);
    }
    return 0;
  }



  template<typename Real>
  int
  nv_eigen_vscalevectorarray(int nv, Real* c, N_Vector* X, N_Vector* Z)
  {
    // invalid number of vectors
    if (nv < 1) return -1;

    for (unsigned int i=0; i<nv; ++i)
    {
      nv_eigen_scale<Real>(c[i], X[i], Z[i]);
    }
    return 0;
  }









  template<typename Real>
  N_Vector
  N_VNew_Eigen(sunindextype N){
    // allocate eigen vector on heap
    auto vec = new Eigen::Vector<Real, Eigen::Dynamic>(N);
    vec->setZero();
    
    sundials::Context cntxt;
    N_Vector nv = N_VNewEmpty(cntxt);

    nv->content = (void*) vec;

    nv->ops->nvgetlength         = nv_eigen_getlength<Real>;
    nv->ops->nvclone             = nv_eigen_clone<Real>;
    nv->ops->nvcloneempty        = nv_eigen_clone_empty;
    nv->ops->nvdestroy           = nv_eigen_destroy<Real>;
    nv->ops->nvspace             = nv_eigen_space;
    nv->ops->nvgetarraypointer   = nv_eigen_getarraypointer<Real>;
    nv->ops->nvsetarraypointer   = nv_eigen_setarraypointer<Real>;
    nv->ops->nvlinearsum         = nv_eigen_linearsum<Real>;
    nv->ops->nvconst             = nv_eigen_const<Real>;
    nv->ops->nvprod              = nv_eigen_prod<Real>;
    nv->ops->nvdiv               = nv_eigen_div<Real>;
    nv->ops->nvscale             = nv_eigen_scale<Real>;
    nv->ops->nvabs               = nv_eigen_abs<Real>;
    nv->ops->nvinv               = nv_eigen_inv<Real>;
    nv->ops->nvaddconst          = nv_eigen_addconst<Real>;
    nv->ops->nvmaxnorm           = nv_eigen_maxnorm<Real>;
    nv->ops->nvwrmsnorm          = nv_eigen_wrms_norm<Real>;
    nv->ops->nvmin               = nv_eigen_min<Real>;
    nv->ops->nvminquotient       = nv_eigen_minquotient<Real>;
    nv->ops->nvconstrmask        = nv_eigen_vconstrmask<Real>;
    nv->ops->nvcompare           = nv_eigen_vcompare<Real>;
    nv->ops->nvinvtest           = nv_eigen_vinvtest<Real>;
    nv->ops->nvlinearcombination = nv_eigen_linearcombination<Real>;
    nv->ops->nvscaleaddmulti     = nv_eigen_vscaleaddmulti<Real>;
    nv->ops->nvdotprodmulti      = nv_eigen_vdotprodmulti<Real>;
    nv->ops->nvscalevectorarray  = nv_eigen_vscalevectorarray<Real>;

 
    return nv;
  }
}
#endif