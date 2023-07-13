#include "particleSystem.h"
#include "odeSolver.h"
#include "linearAlgebraEigen.h"
#include "linearSolverEigen.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <random>

#include <omp.h>

#define OPTIM_ENABLE_EIGEN_WRAPPERS
#include "optim.hpp"


using Matrix = Eigen::Matrix<realtype, Eigen::Dynamic, Eigen::Dynamic>;
using Solver = Eigen::PartialPivLU< Matrix >;


// TODO -- explain stuff
std::vector< N_Vector > simulateRxns(const Eigen::Vector<realtype, Eigen::Dynamic> & parameters){
  // Similar to the function in step02 but stores the simulation results
  // in a vector instead of writing them to a file
  NanoSim::particleSystem<double> my_rxns;

  const realtype k1 = parameters[0];
  const realtype k2 = parameters[1]; 
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

// TODO -- explain stuff
double costFunctionOptimLib(const Eigen::Vector<realtype, Eigen::Dynamic> &parameters,
  Eigen::Vector<realtype, Eigen::Dynamic>* grad,
  void* f_data){

  auto data = makeFakeData();

  auto cost = costFunction(parameters, data);
  for (auto v : data){
    N_VDestroy(v);
  }

  return cost;
}



int main(){
  // Choose starting parameters at random
  Eigen::Vector<realtype, Eigen::Dynamic> param(2);
  std::random_device rng;
  std::mt19937 gen(rng());
  std::uniform_real_distribution<> unif(0.0, 10.0);
  param(0) = unif(gen);
  param(1) = unif(gen);

  std::cout << "Starting guess: "
    << "k1 = " << param[0]
    << "   k2 = " << param[1]
    << "\n";

  Eigen::Vector<realtype, Eigen::Dynamic> lb(2);
  lb << 0.0, 0.0;
  Eigen::Vector<realtype, Eigen::Dynamic> ub(2);
  ub << 10.0, 10.0;

  optim::algo_settings_t options;
  // options.vals_bound = true;
  options.lower_bounds = lb;
  options.lower_bounds = ub;
  options.print_level = 2;
  // options.iter_max = 10;
  options.de_settings.n_gen = 100;
  options.de_settings.initial_lb = lb;
  options.de_settings.initial_ub = ub;

  bool success = optim::de(param, costFunctionOptimLib, nullptr, options);
  if (success){
    std::cout << "Optimum is: f(" << param(0) << ", " << param(1) << ") = "
      << options.opt_fn_value << "\n";
  }
  else {
    std::cout << "Failed :(\n";
  }

}