#ifndef NANOSIM_LINEARSOLVERSCHURCOMPLEMENT_H
#define NANOSIM_LINEARSOLVERSCHURCOMPLEMENT_H

#include "sundials/sundials_linearsolver.h"
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_types.h>

#include <linearAlgebra.h>

#include <iostream>

namespace NanoSim{
  template<typename Real>
  struct schurSolverContent{
    abstractLinearAlgebraOperations<Real>* lin_alg;
    SUNLinearSolver solver;
    SUNMatrix AminusBDinvC;
    SUNMatrix A;
    SUNMatrix B;
    SUNMatrix C;
    N_Vector Diag;
    N_Vector SubDiag;
    N_Vector u;
    N_Vector v;
    N_Vector tmpU;
    N_Vector tmpV;
    N_Vector x;
    N_Vector y;
    sunindextype Adim;
    sunindextype Ddim;
  };



  



  template <typename Real>
  int
  schur_free(SUNLinearSolver LS){
    auto schur_content = static_cast< schurSolverContent<Real>* >(LS->content);

    // linear algebra object is handled outside of the schur solver, so just remove the pointer
    schur_content->lin_alg = nullptr;

    // SUNLinear solver has its own free routine
    int err = SUNLinSolFree(schur_content->solver);
    schur_content->solver = nullptr;

    // Free matrices
    SUNMatDestroy(schur_content->AminusBDinvC);
    schur_content->AminusBDinvC = nullptr;

    SUNMatDestroy(schur_content->A);
    schur_content->A = nullptr;

    SUNMatDestroy(schur_content->B);
    schur_content->B = nullptr;

    SUNMatDestroy(schur_content->C);
    schur_content->C = nullptr;

    // Free vectors
    N_VDestroy(schur_content->Diag);
    schur_content->Diag = nullptr;

    N_VDestroy(schur_content->SubDiag);
    schur_content->SubDiag = nullptr;

    N_VDestroy(schur_content->u);
    schur_content->u = nullptr;

    N_VDestroy(schur_content->v);
    schur_content->v = nullptr;

    N_VDestroy(schur_content->tmpU);
    schur_content->tmpU = nullptr;

    N_VDestroy(schur_content->tmpV);
    schur_content->tmpV = nullptr;

    N_VDestroy(schur_content->x);
    schur_content->x = nullptr;

    N_VDestroy(schur_content->y);
    schur_content->y = nullptr;

    // Now that the content is freed, we can delete it from the heap
    delete schur_content;
    LS->content = nullptr;

    // SUNDIALS can clean up the rest
    SUNLinSolFreeEmpty(LS);

    return err; 
  }



  template <typename Real>
  SUNLinearSolver_Type
  get_type_schur(SUNLinearSolver LS){
    auto schur_content = static_cast< schurSolverContent<Real>* >(LS->content);
    return SUNLinSolGetType(schur_content->solver);
  }



  template <typename Real>
  int
  schur_initialize(SUNLinearSolver LS){
    auto schur_content = static_cast< schurSolverContent<Real>* >(LS->content);
    if (schur_content->solver == NULL){
      return 0;
    } else {
      return schur_content->solver->ops->initialize(schur_content->solver);
    }
  }



  template <typename Real>
  int
  setup_schur_solver(SUNLinearSolver LS, SUNMatrix M){
    auto schur_content = static_cast< schurSolverContent<Real>* >(LS->content);
    // Decompose matrix
    get_diagonalD<Real>(M, 
      schur_content->Diag,
      schur_content->SubDiag,
      schur_content->lin_alg,
      schur_content->Adim,
      schur_content->Ddim);
    
    schur_content->lin_alg->getSubmatrix(M,
      schur_content->B,
      0, schur_content->Adim, 
      schur_content->Adim -1, schur_content->Adim + schur_content->Ddim -1);

    schur_content->lin_alg->getSubmatrix(M,
      schur_content->C,
      schur_content->Adim, 0,
      schur_content->Adim + schur_content->Ddim -1, schur_content->Adim -1);

    schur_content->lin_alg->getSubmatrix(M,
      schur_content->A,
      0, 0,
      schur_content->Adim -1, schur_content->Adim -1);

    // Need the matrix A - BD^{-1}C
    // Can form each column of E = D^{-1}C by with E_i = D^{-1}C_i

    // TODO
    //int err = SUNMatMatvecSetup(schur_content->B);
    int err = 0;
    for (sunindextype c=0;c<schur_content->Adim;++c){
      // Fill tmpV with column of C
      schur_content->lin_alg->matrixGetColumn(schur_content->C, schur_content->tmpV, c);

      // solver fills tmpV with column of E
      err = diag_subdiag_solve<Real>(schur_content->Diag,
        schur_content->SubDiag,
        schur_content->tmpV,
        schur_content->tmpV,
        schur_content->lin_alg);

      // tmpV now equals E(:, c)
      // For F = BD^{-1}C, we know that
      // F(:,c) = B*E(:,c)
      err = SUNMatMatvec(schur_content->B, schur_content->tmpV, schur_content->tmpU);


      // tmpU is now column c of BD^{-1}C and we can store it in AminusBDinvC
      schur_content->lin_alg->matrixInsertColumn(schur_content->AminusBDinvC, schur_content->tmpU, c);
    }

    // Form A - BD^{-1}C
    // AminusBDinvC = BDinvC right now
    // So multiply by -1.0 and add A gives correct result
    err = SUNMatScaleAdd(-1.0, schur_content->AminusBDinvC, schur_content->A);
  
    // Perform setup of provided linear solver
    err = SUNLinSolSetup(schur_content->solver, schur_content->AminusBDinvC);

    return err;
  }



  template<typename Real>
  int
  schur_solve(SUNLinearSolver solver, SUNMatrix M, N_Vector sol_vec, N_Vector b, Real tol){
    auto schur_content = static_cast< schurSolverContent<Real>* >(solver->content);

    /*
      Prototype system looks like
      | A  B | | x | | u | 
      |      |*|   |=|   |
      | C  D | | y | | v |

      Hence:
        Ax + By = u (1)
        Cx + Dy = v (2)
      Here we assume D has an easy inverse to compute. Thus, from (2) we get:
        y = Dinv(v - Cx) (3)
      Substitute (3) into (1) to get:
        Ax + B*Dinv*v - B*Dinv*C*x = u (4)
      and hence
        (A - BDinvC)x = u - BDinv*v (5)
      Therefore solution strategy is:
        solve (5) for x
        solve for y by substituting now known x into (3)
    */

    // Allocate working vectors
    schur_content->lin_alg->getSubvector(b, schur_content->u, 0);
    schur_content->lin_alg->getSubvector(b, schur_content->v, schur_content->Adim);

    // Form right hand side of (5)
    // Diag and SubDiag are formed in setup
    int err = diag_subdiag_solve<Real>(schur_content->Diag, 
      schur_content->SubDiag, 
      schur_content->v,
      schur_content->tmpV,
      schur_content->lin_alg);

    err = SUNMatMatvec(schur_content->B,
      schur_content->tmpV, 
      schur_content->tmpU);

    N_VLinearSum(1.0,schur_content->u,
      -1.0,schur_content->tmpU,
      schur_content->tmpU);


    // Solve for x
    err = SUNLinSolSolve(schur_content->solver,
      schur_content->AminusBDinvC,
      schur_content->x,
      schur_content->tmpU,
      tol);

    // Form vector in (3)
    err = SUNMatMatvec(schur_content->C,
      schur_content->x,
      schur_content->tmpV);

    N_VLinearSum(1.0,schur_content->v,
      -1.0,schur_content->tmpV,
      schur_content->tmpV);


    // Solve for y
    err = diag_subdiag_solve<Real>(schur_content->Diag,
      schur_content->SubDiag,
      schur_content->tmpV,
      schur_content->y,
      schur_content->lin_alg);

    
    // Stack x and y to form the solution
    schur_content->lin_alg->stackVectors(schur_content->x,
      schur_content->y,
      sol_vec);

    return err;
  }




  template<typename Real>
  void
  get_diagonalD(SUNMatrix M, N_Vector D, N_Vector Dlo, abstractLinearAlgebraOperations<Real> *lin_alg, sunindextype N, sunindextype nd){
    lin_alg->vectorInsert(D, lin_alg->matrixGetValue(M,N,N), 0);
    for (sunindextype i=N+1; i<N+nd; ++i){
      lin_alg->vectorInsert(D, lin_alg->matrixGetValue(M,i,i), i-N);
      lin_alg->vectorInsert(Dlo, lin_alg->matrixGetValue(M,i,i-1), i-N-1);
    }
  }



  template<typename Real>
  int
  diag_subdiag_solve(N_Vector Diag, N_Vector SubDiag, N_Vector b, N_Vector x, abstractLinearAlgebraOperations<Real> *lin_alg){
    Real val = lin_alg->vectorGetValue(b,0) / lin_alg->vectorGetValue(Diag,0);
    lin_alg->vectorInsert(x, val, 0);

    for (sunindextype i=1; i<N_VGetLength(x);++i){
      Real val = (lin_alg->vectorGetValue(b,i) - lin_alg->vectorGetValue(SubDiag,i-1)*lin_alg->vectorGetValue(x,i-1)) / lin_alg->vectorGetValue(Diag,i);
      lin_alg->vectorInsert(x, val, i);
    }

    return 0;
  }



  template <typename Real>
  SUNLinearSolver
  create_schur_complement_solver(abstractLinearAlgebraOperations<Real> *lin_alg, 
    SUNLinearSolver lin_solver, 
    sunindextype N, 
    sunindextype split){

    SUNLinearSolver solver = SUNLinSolNewEmpty();

    solver->ops->gettype = NanoSim::get_type_schur<Real>;
    solver->ops->setup = NanoSim::setup_schur_solver<Real>; 
    solver->ops->solve = NanoSim::schur_solve<Real>;
    solver->ops->free = NanoSim::schur_free<Real>; 
    solver->ops->initialize = NanoSim::schur_initialize<Real>;

    schurSolverContent<Real>* schur_content = new schurSolverContent<Real>;
    schur_content->lin_alg = lin_alg;
    schur_content->solver = lin_solver;

    // Create structures for the matrix decomposition
    /*
      Matrix M is NxN and decomposed as
          | A  B |
      M = | C  D |
      Split is the dimension where the schur decomposition occurs,
      so A is Split x Split
         B is Split x (N-Split)
         C is (N-Split) x Split
         D is (N-Split) x (N-Split)
    */
    sunindextype M = N - split;
    schur_content->AminusBDinvC = lin_alg->createNewMatrix(split, split);
    schur_content->A = lin_alg->createNewMatrix(split, split);
    schur_content->B = lin_alg->createNewMatrix(split, M);
    schur_content->C = lin_alg->createNewMatrix(M, split);
    schur_content->Diag = lin_alg->createNewVector(M);
    schur_content->SubDiag = lin_alg->createNewVector(M-1);
    schur_content->u = lin_alg->createNewVector(split);
    schur_content->v = lin_alg->createNewVector(M);
    schur_content->tmpU = lin_alg->createNewVector(split);
    schur_content->tmpV = lin_alg->createNewVector(M);
    schur_content->x = lin_alg->createNewVector(split);
    schur_content->y = lin_alg->createNewVector(M);

    schur_content->Adim = split;
    schur_content->Ddim = M;

    solver->content = (void*)schur_content;
    return solver;
  }
}
#endif