#include "particleSystem.h"
#include "odeSolver.h"
#include "linearAlgebraEigen.h"
#include "linearSolverEigen.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <random>

using Matrix = Eigen::Matrix<realtype, Eigen::Dynamic, Eigen::Dynamic>;
using Solver = Eigen::PartialPivLU< Matrix >;

/*
  This tutorial shows how to construct a cost function that is dependent on:
    Input: data, parameters
    Performs:
      (1) simulates Rxns based on parameters
      (2) compares simulation results to data
    Output: quantitative value corresponding to comparisons
*/


// TODO -- explain stuff
std::vector< N_Vector > simulateRxns(const Eigen::Vector<realtype, Eigen::Dynamic> & parameters){
  // Similar to the function in step02 but stores the simulation results
  // in a vector instead of writing them to a file
  NanoSim::particleSystem<double> my_rxns;

  const realtype k1 = parameters(0);
  const realtype k2 = parameters(1); 
  my_rxns.addReaction({{1,"A"},{1,"B"}},
    {{1,"C"}},
    k1);
  my_rxns.addReaction({{1,"A"}, {1,"C"}},
    {{1,"D"}},
    k2);

  my_rxns.finalizeReactions();


  NanoSim::eigenLinearAlgebraOperations<realtype, Matrix> lin_alg;
  const sunindextype length = my_rxns.getNumberOfSpecies();
  auto ic = lin_alg.createNewVector(length);
  ic->ops->nvconst(0.0,ic);
  lin_alg.vectorInsert(ic, 10.0, 0); 
  lin_alg.vectorInsert(ic, 5.0, 1); 


  auto template_matrix = lin_alg.createNewMatrix(length, length);

  auto lin_solve = NanoSim::createLinearSolverEigenDense<realtype, Solver>();

  NanoSim::cvodeOptions<realtype> opts("SUNDIALS_errors.txt", 
    1e-8, 
    1e-14, 
    1000, 
    SUNFALSE );
  auto cvode_mem = NanoSim::prepareODESolver<realtype>(ic, 
    template_matrix, 
    lin_solve, 
    opts);

  std::pair< NanoSim::particleSystem<realtype>*, NanoSim::abstractLinearAlgebraOperations<realtype>* >
  data_pair = {&my_rxns, &lin_alg};
  void * user_data = static_cast<void *>(&data_pair);
  auto flag = CVodeSetUserData(cvode_mem, user_data);


  std::vector<realtype> times;
  const realtype dt = 0.05;
  const realtype endtime = 2.0;
  realtype t = 0.0;
  while (t < endtime){
    t += dt;
    times.push_back(t);
  }

  std::vector< N_Vector > solutions;
  std::ofstream outputfile;
  outputfile.open("step03-times.txt");
  for (const auto solve_time : times){
    outputfile << solve_time << "\n";
    N_Vector sol = ic->ops->nvclone(ic);
    realtype t;

    auto err = CVode(cvode_mem, solve_time, sol, &t, CV_NORMAL);

    solutions.push_back(sol);
    // We need to use sol later so we no longer delete sol now
    // we MUST delete sol later though!!! 
  }
  outputfile.close();

  CVodeFree(&cvode_mem);
  N_VDestroy(ic);
  SUNLinSolFree(lin_solve);
  SUNMatDestroy(template_matrix);

  return solutions;
}


// TODO -- explain stuff
std::vector < N_Vector > makeFakeData(){
  // Simulate reactions with "known" parameters
  Eigen::Vector<realtype, Eigen::Dynamic> true_parameters(2);
  true_parameters << 2.0, 3.0;
  auto measurements = simulateRxns(true_parameters);

  // Add some noise to the solutions to simulate an actual data measurement
  std::mt19937 gen(1); // seed for reproducibility
  std::normal_distribution<realtype> distr(0, 0.1); // noise follows normal distribution
  NanoSim::eigenLinearAlgebraOperations<realtype, Matrix> lin_alg;
  auto noise = lin_alg.createNewVector(N_VGetLength(measurements[0]));
  for (auto vec : measurements){
    for (sunindextype idx=0; idx<N_VGetLength(noise); ++idx){
      lin_alg.vectorInsert(noise, distr(gen), idx);
    }
    N_VLinearSum(1.0, vec, 1.0, noise, vec);
  }
  return measurements;
}


// TODO -- explain stuff
realtype costFunction(const Eigen::Vector<realtype, Eigen::Dynamic> & parameters,
  std::vector< N_Vector > & data){
  // First perform a simulation with the parameters
  auto sim = simulateRxns(parameters);

  // Check against data
  realtype cost = 0.0;
  NanoSim::eigenLinearAlgebraOperations<realtype, Matrix> lin_alg;
  const sunindextype length = N_VGetLength(sim[0]);
  N_Vector compare = lin_alg.createNewVector(length);
  N_Vector weights = lin_alg.createNewVector(length);
  N_VConst(1.0, weights); // all weights = 1 so we just do a sum of least squares
  for (unsigned int i=0; i< data.size(); ++i){
    N_VLinearSum(1.0, data[i], -1.0, sim[i], compare);
    cost += N_VWrmsNorm(compare, weights);
  }

  
  // Now we are done with the N_Vector's in "sim" so we
  // delete their memory!!
  for (auto v : sim){
    N_VDestroy(v);
  }
  N_VDestroy(compare);
  N_VDestroy(weights);

  return cost;
}

int main(){
  auto data = makeFakeData();
  // Output data for visualization purposes
  NanoSim::eigenLinearAlgebraOperations<realtype, Matrix> lin_alg;
  std::ofstream outputfile;
  outputfile.open("step03-data-conc.txt");
  for (auto v : data){
    outputfile << lin_alg.vectorGetValue(v, 0)
      << "    " << lin_alg.vectorGetValue(v, 1)
      << "    " << lin_alg.vectorGetValue(v, 2)
      << "    " << lin_alg.vectorGetValue(v, 3)
      << "\n";
  }
  outputfile.close();

  // Check cost when the parameters are far from "true"
  Eigen::Vector<realtype, Eigen::Dynamic> bad_param(2);
  bad_param << 10.0, 0.2;
  auto bad_cost = costFunction(bad_param, data);
  std::cout << "Results of analysis\n"
    << "_____________________________________________________\n"
    << "Bad parameters:    "
    << "k1 = " << bad_param(0)
    << "  k2 = " << bad_param(1)
    << "  Cost = " << bad_cost << "\n";
  // Output simulation for visualization purposes
  auto bad_sim = simulateRxns(bad_param);
  outputfile.open("step03-badsim-conc.txt");
  for (auto v : bad_sim){
    outputfile << lin_alg.vectorGetValue(v, 0)
      << "    " << lin_alg.vectorGetValue(v, 1)
      << "    " << lin_alg.vectorGetValue(v, 2)
      << "    " << lin_alg.vectorGetValue(v, 3)
      << "\n";
  }
  outputfile.close();


  // Check cost when the parameters are exact
  Eigen::Vector<realtype, Eigen::Dynamic> good_param(2);
  good_param << 2.0, 3.0;
  auto good_cost = costFunction(good_param, data);
  std::cout << "Good parameters:    "
    << "k1 = " << good_param(0)
    << "  k2 = " << good_param(1)
    << "  Cost = " << good_cost << "\n";
  // Output simulation for visualization purposes
  auto good_sim = simulateRxns(good_param);
  outputfile.open("step03-goodsim-conc.txt");
  for (auto v : good_sim){
    outputfile << lin_alg.vectorGetValue(v, 0)
      << "    " << lin_alg.vectorGetValue(v, 1)
      << "    " << lin_alg.vectorGetValue(v, 2)
      << "    " << lin_alg.vectorGetValue(v, 3)
      << "\n";
  }
  outputfile.close();


  // Delete data at the end manually
  for (auto v : data){
    N_VDestroy(v);
  }
  for (auto v : bad_sim){
    N_VDestroy(v);
  }
  for (auto v : good_sim){
    N_VDestroy(v);
  }
}