#include "particleSystem.h"
#include "odeSolver.h"
#include "linearAlgebraEigen.h"
#include "linearSolverEigen.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <chrono>
#include <cmath>
#include <string>

#include <Eigen/Sparse>

/*
  This tutorial program shows how to construct a set of chemical reactions and then solve the ODEs for them
  TODO talk about more
*/

using Matrix = Eigen::Matrix<realtype, Eigen::Dynamic, Eigen::Dynamic>;
using SparseMatrix = Eigen::SparseMatrix<realtype>;
using Solver = Eigen::PartialPivLU< Matrix >;
using SolverSparse = Eigen::BiCGSTAB< SparseMatrix, Eigen::IncompleteLUT<realtype> >;


std::vector<double> simulateRxnsDense(const Eigen::Vector<realtype, Eigen::Dynamic> & parameters){
  const realtype factor = parameters(3);
  const realtype kb = 1.37e5 * factor;
  const realtype kf = 5e-7 * kb * factor;
  const realtype k1 = 7.69e4 * factor;
  const realtype k2 = 1.40e4 * factor;
  const realtype k3 = 7.15e3 * factor;
  const realtype k4 = 1.74e3 * factor;
  const unsigned int M = 111;
  const realtype S = 11.3;
  const realtype T = 4.838;

  const realtype icA = parameters(0);
  const unsigned int max_size = parameters(1);
  const unsigned int n_solves = parameters(2);

  
  NanoSim::particleSystem<double> my_rxns;

  my_rxns.defineParticle(3, max_size, M);

  my_rxns.addReaction({{1,"A"}},
    {{1,"As"}, {1,"L"}},
    kf * S * S);
  my_rxns.addReaction({{1,"As"}, {1,"L"}},
    {{1,"A"}},
    kb);

  my_rxns.addNucleation({{1, "A"}, {2,"As"}},
    {{1,"L"}},
    k1);

  const std::function<realtype(const unsigned int)> growth_kernel = 
    [&](const unsigned int size){
      return 2.677 * std::pow(size * 1.0, 0.72) * (size <= M ? k2 : k3);
    };
  my_rxns.addGrowth({{1, "A"}},
    {{1,"L"}},
    growth_kernel);

  const std::function<realtype(const unsigned int, const unsigned int)> agglom_kernel =
    [&](const unsigned int size1, const unsigned int size2){
      return 2.677 * std::pow(size1 * 1.0, 0.72) * 2.677 * std::pow(size2 * 1.0, 0.72) * k4;
    };
  my_rxns.addAgglomeration({},{},agglom_kernel);

  my_rxns.finalizeReactions();

  const auto length = my_rxns.getNumberOfSpecies();
  NanoSim::eigenLinearAlgebraOperations<realtype, Matrix> lin_alg;
  auto ic = lin_alg.createNewVector(length);
  N_VConst(0.0, ic);
  lin_alg.vectorInsert(ic, icA, 0);

  auto template_matrix = lin_alg.createNewMatrix(length, length);
  auto lin_solve = NanoSim::createLinearSolverEigenDense<realtype, Solver>();

  NanoSim::cvodeOptions<realtype> opts("SUNDIALS_errors.txt", 
    1e-8, 
    1e-14, 
    1000, 
    SUNTRUE );


  // solve it a few times to benchmark the time it takes
  std::vector<double> solve_times;
  std::string prev_out = "";
  std::string out_text = "";

  out_text = "Solved 0/" + std::to_string(n_solves);
  std::cout << out_text << std::flush;

  for (unsigned int i=0; i<n_solves;++i){
    auto cvode_mem = NanoSim::prepareODESolver<realtype>(ic, 
      template_matrix, 
      lin_solve, 
      opts);

    std::pair< NanoSim::particleSystem<realtype>*, NanoSim::abstractLinearAlgebraOperations<realtype>* >
    data_pair = {&my_rxns, &lin_alg};
    void * user_data = static_cast<void *>(&data_pair);
    auto flag = CVodeSetUserData(cvode_mem, user_data);

    auto sol = N_VClone(ic);
    realtype t;
    // start timing
    auto start = std::chrono::high_resolution_clock::now();
    auto err = CVode(cvode_mem, T, sol, &t, CV_NORMAL);
    // end timing
    auto end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
    solve_times.push_back(elapsed);
    

    if (i == 0){
      // output solution
      std::ofstream outputfile;
      std::string filename = "step05-densesolve-" + std::to_string(max_size) + "-particles-PSD.txt";
      outputfile.open(filename.c_str());
      // just PSD -- starts at index 3
      for (sunindextype idx=3;idx<=max_size;++idx){
        outputfile << lin_alg.vectorGetValue(sol, idx) << "\n";
      }
      outputfile.close();
    }
    N_VDestroy(sol);
    CVodeFree(&cvode_mem);
    prev_out = out_text;
    out_text = "Solved " + std::to_string(i+1) + "/" + std::to_string(n_solves);
    std::cout << "\r" << out_text << std::flush;
  }
  std::cout << "\n";

  N_VDestroy(ic);
  SUNLinSolFree(lin_solve);
  SUNMatDestroy(template_matrix);

  std::ofstream outputfile;
  std::string filename = "step05-densesolve-" + std::to_string(max_size) + "-particles-solvetimes.txt";
  outputfile.open(filename.c_str());
  for (auto t : solve_times){
    outputfile << t << "\n";
  }
  outputfile.close();

  return solve_times;
}

std::vector<double> simulateRxnsSparse(const Eigen::Vector<realtype, Eigen::Dynamic> & parameters){
  const realtype factor = parameters(3);
  const realtype kb = 1.37e5 ;
  const realtype kf = 5e-7 * kb;
  const realtype k1 = 7.69e4 ;
  const realtype k2 = 1.40e4 * factor;
  const realtype k3 = 7.15e3 * factor;
  const realtype k4 = 1.74e3 * factor;
  const unsigned int M = 111;
  const realtype S = 11.3;
  const realtype T = 10.0;

  const realtype icA = parameters(0) * factor;
  const unsigned int max_size = parameters(1);
  const unsigned int n_solves = parameters(2);

  
  NanoSim::particleSystem<double> my_rxns;

  my_rxns.defineParticle(3, max_size, M);

  my_rxns.addReaction({{1,"A"}},
    {{1,"As"}, {1,"L"}},
    kf * S * S);
  my_rxns.addReaction({{1,"As"}, {1,"L"}},
    {{1,"A"}},
    kb);

  my_rxns.addNucleation({{1, "A"}, {2,"As"}},
    {{1,"L"}},
    k1);

  const std::function<realtype(const unsigned int)> growth_kernel = 
    [&](const unsigned int size){
      return 2.677 * std::pow(size * 1.0, 0.72) * (size <= M * factor ? k2 : k3);
    };
  my_rxns.addGrowth({{1, "A"}},
    {{1,"L"}},
    growth_kernel);

  const std::function<realtype(const unsigned int, const unsigned int)> agglom_kernel =
    [&](const unsigned int size1, const unsigned int size2){
      return 2.677 * std::pow(size1 * 1.0, 0.72) * 2.677 * std::pow(size2 * 1.0, 0.72) * k4;
    };
  my_rxns.addAgglomeration({},{},agglom_kernel);

  my_rxns.finalizeReactions();

  const auto length = my_rxns.getNumberOfSpecies();
  NanoSim::eigenLinearAlgebraOperations<realtype, SparseMatrix> lin_alg;
  auto ic = lin_alg.createNewVector(length);
  N_VConst(0.0, ic);
  lin_alg.vectorInsert(ic, icA, 0);

  auto template_matrix = lin_alg.createNewMatrix(length, length);
  auto lin_solve = NanoSim::createLinearSolverEigenSparse<realtype, SolverSparse>();

  NanoSim::cvodeOptions<realtype> opts("SUNDIALS_errors.txt", 
    1e-8, 
    1e-14, 
    1000, 
    SUNTRUE );


  // solve it a few times to benchmark the time it takes
  std::vector<double> solve_times;
  std::string prev_out = "";
  std::string out_text = "";

  out_text = "Solved 0/" + std::to_string(n_solves);
  std::cout << out_text << std::flush;

  for (unsigned int i=0; i<n_solves;++i){
    auto cvode_mem = NanoSim::prepareODESolver<realtype>(ic, 
      template_matrix, 
      lin_solve, 
      opts);

    std::pair< NanoSim::particleSystem<realtype>*, NanoSim::abstractLinearAlgebraOperations<realtype>* >
    data_pair = {&my_rxns, &lin_alg};
    void * user_data = static_cast<void *>(&data_pair);
    auto flag = CVodeSetUserData(cvode_mem, user_data);

    auto sol = N_VClone(ic);
    realtype t;
    // start timing
    auto start = std::chrono::high_resolution_clock::now();
    auto err = CVode(cvode_mem, T, sol, &t, CV_NORMAL);
    // end timing
    auto end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
    solve_times.push_back(elapsed);
    

    if (i == 0){
      // output solution
      std::ofstream outputfile;
      std::string filename = "step05-sparsesolve-" + std::to_string(max_size) + "-particles-PSD.txt";
      outputfile.open(filename.c_str());
      // just PSD -- starts at index 3
      for (sunindextype idx=3;idx<=max_size;++idx){
        outputfile << lin_alg.vectorGetValue(sol, idx) << "\n";
      }
      outputfile.close();
    }
    N_VDestroy(sol);
    CVodeFree(&cvode_mem);
    prev_out = out_text;
    out_text = "Solved " + std::to_string(i+1) + "/" + std::to_string(n_solves);
    std::cout << "\r" << out_text << std::flush;
  }
  std::cout << "\n";

  N_VDestroy(ic);
  SUNLinSolFree(lin_solve);
  SUNMatDestroy(template_matrix);

  std::ofstream outputfile;
  std::string filename = "step05-sparsesolve-" + std::to_string(max_size) + "-particles-solvetimes.txt";
  outputfile.open(filename.c_str());
  for (auto t : solve_times){
    outputfile << t << "\n";
  }
  outputfile.close();

  return solve_times;
}


int main(int argc, char **argv){
  // default settings
  realtype ic = 0.0012;
  realtype max_size = 2500;
  realtype n_solves = 1;
  realtype factor = 1.0;
  // check if command line arguments were provided
  if (argc > 1){
    ic = std::atof(argv[1]);
  }
  if (argc > 2){
    max_size = std::atof(argv[2]);
  }
  if (argc > 3){
    n_solves = std::atof(argv[3]);
  }
  if (argc > 4){
    factor = std::atof(argv[4]);
  }
  Eigen::Vector<realtype, Eigen::Dynamic> settings(4);
  settings(0) = ic; // initial concentration of precursor
  settings(1) = max_size;   // largest particle size to track
  settings(2) = n_solves;      // number of times to solve ODE to benchmark time it takes
  settings(3) = factor;

  // std::cout << "Dense matrix solve:\n";
  // auto dense_times = simulateRxnsDense(settings);

  std::cout << "Sparse matrix solve:\n";
  auto sparse_times = simulateRxnsSparse(settings);

}