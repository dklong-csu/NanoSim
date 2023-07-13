#ifndef NANOSIM_LINEARSOLVEREIGENDENSE_H
#define NANOSIM_LINEARSOLVEREIGENDENSE_H

#include "sundials/sundials_linearsolver.h"

#include "linearAlgebraEigen.h"
#include <Eigen/Sparse>

#include <iostream>

namespace NanoSim{
  SUNLinearSolver_Type getTypeEigenDenseSolver(SUNLinearSolver LS){
    return SUNLINEARSOLVER_DIRECT; // direct solver and this uses a direct solver always
  }


  SUNLinearSolver_Type getTypeEigenSparseSolver(SUNLinearSolver LS){
    return SUNLINEARSOLVER_MATRIX_ITERATIVE;
  }



  template<typename Real, typename SolverType>
  int setupEigenDenseSolver(SUNLinearSolver LS, SUNMatrix A){
    auto solver = static_cast< SolverType* >(LS->content);
    auto M = static_cast< Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>* >(A->content);
    solver->compute(*M);
    return 0; // indicates success to SUNDIALS
  }



  template<typename Real, typename SolverType>
  int setupEigenSparseSolver(SUNLinearSolver LS, SUNMatrix A){
    auto solver = static_cast< SolverType* >(LS->content);
    auto M = static_cast< Eigen::SparseMatrix<Real>* >(A->content);
    M->makeCompressed();
    solver->compute(*M);
    return 0;
  }



  template<typename Real, typename SolverType>
  int solveEigen(SUNLinearSolver LS, SUNMatrix A, N_Vector x, N_Vector b, realtype tol){
    auto xvec = static_cast< Eigen::Vector<Real, Eigen::Dynamic>* >(x->content);
    auto bvec = static_cast< Eigen::Vector<Real, Eigen::Dynamic>* >(b->content);

    auto solver = static_cast< SolverType* >(LS->content);
    *xvec = solver->solve(*bvec);
    return 0; // indicates success to SUNDIALS
  }



  template<typename Real, typename SolverType>
  int freeEigenSolver(SUNLinearSolver LS){
    auto solver = static_cast< SolverType* >(LS->content);
    delete solver;
    LS->content = nullptr;

    SUNLinSolFreeEmpty(LS);
    return 0;
  }



  template <typename Real, typename SolverType>
  SUNLinearSolver
  createLinearSolverEigenDense(){

    SUNLinearSolver LS = SUNLinSolNewEmpty();

    LS->ops->gettype = NanoSim::getTypeEigenDenseSolver;
    LS->ops->setup = NanoSim::setupEigenDenseSolver<Real, SolverType>;
    LS->ops->solve = NanoSim::solveEigen<Real, SolverType>;
    LS->ops->free  = NanoSim::freeEigenSolver<Real, SolverType>;

    auto solver = new SolverType;

    LS->content = (void*)solver;
    
    return LS;
  }



  template <typename Real, typename SolverType>
  SUNLinearSolver
  createLinearSolverEigenSparse(){

    SUNLinearSolver LS = SUNLinSolNewEmpty();

    LS->ops->gettype = NanoSim::getTypeEigenSparseSolver;
    LS->ops->setup = NanoSim::setupEigenSparseSolver<Real, SolverType>;
    LS->ops->solve = NanoSim::solveEigen<Real, SolverType>;
    LS->ops->free  = NanoSim::freeEigenSolver<Real, SolverType>;

    auto solver = new SolverType;

    LS->content = (void*)solver;
    
    return LS;
  }
}
#endif