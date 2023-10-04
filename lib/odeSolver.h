#ifndef NANOSIM_ODESOLVER_H
#define NANOSIM_ODESOLVER_H


#include<cvode/cvode.h>
#include "particleSystem.h"
#include "linearAlgebra.h"

#include <utility>

namespace NanoSim{
// TODO documentation
  template<typename Real>
  int rhs_function(Real time, N_Vector x, N_Vector x_dot, void* user_data){
    std::pair< NanoSim::particleSystem<Real>*, abstractLinearAlgebraOperations<Real>* >*
    data_pair = static_cast< std::pair< NanoSim::particleSystem<Real>*, abstractLinearAlgebraOperations<Real>* >* >(user_data);
    auto rxn = data_pair->first;
    auto lin_alg = data_pair->second;

    void * rhs_data = static_cast<void *>(lin_alg);
    auto rhs = rxn->composeRHSfunction();
    
    return rhs(time, x, x_dot, rhs_data);
  }



  template<typename Real>
  int jacobian_function(Real time, N_Vector x, N_Vector x_dot, SUNMatrix Jacobian, void * user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
    std::pair< NanoSim::particleSystem<Real>*, abstractLinearAlgebraOperations<Real>* >*
    data_pair = static_cast< std::pair< NanoSim::particleSystem<Real>*, abstractLinearAlgebraOperations<Real>* >* >(user_data);
    auto rxn = data_pair->first;
    auto lin_alg = data_pair->second;

    void * jac_data = static_cast<void *>(lin_alg);
    auto jac = rxn->composeJacobianfunction();
    
    return jac(time, x, x_dot, Jacobian, jac_data, tmp1, tmp2, tmp3);
  }



  template<typename Real>
  int sparse_jacobian_function(Real time, N_Vector x, N_Vector x_dot, SUNMatrix Jacobian, void * user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
    std::pair< NanoSim::particleSystem<Real>*, abstractLinearAlgebraOperations<Real>* >*
    data_pair = static_cast< std::pair< NanoSim::particleSystem<Real>*, abstractLinearAlgebraOperations<Real>* >* >(user_data);
    auto rxn = data_pair->first;
    auto lin_alg = data_pair->second;

    void * jac_data = static_cast<void *>(lin_alg);
    auto jac = rxn->composeSparseJacobianfunction();
    
    return jac(time, x, x_dot, Jacobian, jac_data, tmp1, tmp2, tmp3);
  }



  template<typename Real>
  struct cvodeOptions
  {
    cvodeOptions(std::string error_filename = "SUNDIALS_errors.txt",
      Real rtol = 1e-8,
      Real atol = 1e-14,
      long int max_steps = 1000,
      booleantype detect_bdf_stab = SUNTRUE)
      : rtol(rtol), 
        atol(atol), 
        max_steps(max_steps), 
        detect_bdf_stab(detect_bdf_stab) {
        err_file = fopen(error_filename.c_str(),"w");
    }

    ~cvodeOptions(){
      fclose(err_file);
    }

    FILE* err_file;
    Real rtol;
    Real atol;
    long int max_steps;
    booleantype detect_bdf_stab;
  };
  


// TODO check flag
  template<typename Real>
  void *
  prepareODESolver(N_Vector initial_condition,
    SUNMatrix template_matrix,
    SUNLinearSolver linear_solver,
    cvodeOptions<Real> & opts){
      sundials::Context cntxt;
      void * cvode_memory = CVodeCreate(CV_BDF, cntxt);

      auto flag = CVodeInit(cvode_memory, rhs_function<Real>, 0.0, initial_condition);

      flag = CVodeSetErrFile(cvode_memory, opts.err_file);

      flag = CVodeSStolerances(cvode_memory, opts.rtol, opts.atol);

      flag = CVodeSetLinearSolver(cvode_memory, linear_solver, template_matrix);

      flag = CVodeSetMaxNumSteps(cvode_memory, opts.max_steps);

      flag = CVodeSetStabLimDet(cvode_memory, opts.detect_bdf_stab);

      flag = CVodeSetJacFn(cvode_memory, jacobian_function<Real>);

      return cvode_memory;
    }
}


#endif 