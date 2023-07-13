#include "particleSystem.h"
#include "odeSolver.h"
#include "linearAlgebraEigen.h"
#include "linearSolverEigen.h"

#include <vector>
#include <iostream>
#include <fstream>

/*
  This tutorial program shows how to construct a set of chemical reactions and then solve the ODEs for them
  TODO talk about more
*/

using Matrix = Eigen::Matrix<realtype, Eigen::Dynamic, Eigen::Dynamic>;
using Solver = Eigen::PartialPivLU< Matrix >;


int main(){
  NanoSim::particleSystem<double> my_rxns;
  /*
    Model a simple chemical reaction system
    A + B ->[k1] C
    A + C ->[k2] D
  */
  const realtype k1 = 2.0;
  const realtype k2 = 3.0;
  my_rxns.addReaction({{1,"A"},{1,"B"}},
    {{1,"C"}},
    k1);
  my_rxns.addReaction({{1,"A"}, {1,"C"}},
    {{1,"D"}},
    k2);
  /*
    Chemical reactions need to be finalized when all of them have been written
    This function will create an ordering of each species in the vectors used
    based on the order new species are added to the reaction system.

    Here, the logic is:
    First read the reaction from 
        my_rxns.addReaction({{1,"A"},{1,"B"}},
          {{1,"C"}},
          k1);
      Discover (in this order) chemical species
        A --> map to index 0 of vectors
        B --> map to index 1 of vectors
        C --> map to index 2 of vectors
    Then read the next reaction from
        my_rxns.addReaction({{1,"A"}, {1,"C"}},
          {{1,"D"}},
          k2);
      Discover (in this order) chemical species
        A --> already know about this
        C --> already know about this
        D --> map to index 3 of vectors
  */
  my_rxns.finalizeReactions();

  /*
    Set initial conditions for reaction:
    A(t=0) = 10.0
    B(t=0) = 5.0
    C(t=0) = 0.0
    D(t=0) = 0.0
  */
  NanoSim::eigenLinearAlgebraOperations<realtype, Matrix> lin_alg;
  const sunindextype length = my_rxns.getNumberOfSpecies();
  auto ic = lin_alg.createNewVector(length);
  N_VConst(0.0,ic);
  lin_alg.vectorInsert(ic, 10.0, 0); // Concentration of A
  lin_alg.vectorInsert(ic, 5.0, 1); // Concentration of B

  /*
    Prepare the ODE solver
  */
  auto template_matrix = lin_alg.createNewMatrix(length, length);

  auto lin_solve = NanoSim::createLinearSolverEigenDense<realtype, Solver>();

  NanoSim::cvodeOptions<realtype> opts("SUNDIALS_errors.txt", /* if SUNDIALS throws errors this file will display them */
    1e-8, /* relative tolerance of solver */
    1e-14, /* absolute tolerance of solver */
    1000, /* max number of iterations */
    SUNFALSE /* use algorithm to detect BDF stability? */);
  auto cvode_mem = NanoSim::prepareODESolver<realtype>(ic, 
    template_matrix, 
    lin_solve, 
    opts);

  std::pair< NanoSim::particleSystem<realtype>*, NanoSim::abstractLinearAlgebraOperations<realtype>* >
  data_pair = {&my_rxns, &lin_alg};
  void * user_data = static_cast<void *>(&data_pair);
  auto flag = CVodeSetUserData(cvode_mem, user_data);

  /*
    Solve the ODEs for a few different times
  */
  std::vector<realtype> times;
  const realtype dt = 0.01;
  const realtype endtime = 2.0;
  realtype t = 0.0;
  while (t < endtime){
    t += dt;
    times.push_back(t);
  }

  std::ofstream outputfile; // This will be used to store results to visualize in e.g. Matlab
  // first output all times since we know them already
  outputfile.open("step01-output-times.txt");
  outputfile << 0.0 << "\n"; // for initial condition
  for (const auto t : times){
    outputfile << t << "\n";
  }
  outputfile.close();

  // Now open the file for the concentrations we solve for
  outputfile.open("step01-output-concentrations.txt"); 
  outputfile << lin_alg.vectorGetValue(ic, 0)
    << "    " << lin_alg.vectorGetValue(ic, 1)
    << "    " << lin_alg.vectorGetValue(ic, 2)
    << "    " << lin_alg.vectorGetValue(ic, 3)
    << "\n";

  for (const auto solve_time : times){
    // Need to create a vector to hold the solution
    N_Vector sol = ic->ops->nvclone(ic);
    // You just need this for the ODE solver -- reasons are unimportant for now
    realtype t;
    // This solves the ODE at time solve_time
    // err is an error code the CVode will provide if something went wrong
    // meanings of the error code can be found at: TODO

    auto err = CVode(cvode_mem, solve_time, sol, &t, CV_NORMAL);

    // Write the solution to a file to plot
    outputfile << lin_alg.vectorGetValue(sol, 0)
      << "    " << lin_alg.vectorGetValue(sol, 1)
      << "    " << lin_alg.vectorGetValue(sol, 2)
      << "    " << lin_alg.vectorGetValue(sol, 3)
      << "\n";

    sol->ops->nvdestroy(sol); // We have to manage this memory manually, unfortunately
  }
  outputfile.close(); // Need to close the file at the end

  // We have to manage this memory manually, unfortunately
  CVodeFree(&cvode_mem);
  N_VDestroy(ic);
  SUNLinSolFree(lin_solve);
  SUNMatDestroy(template_matrix);
}