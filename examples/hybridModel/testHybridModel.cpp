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
#include <numeric>

#include <Eigen/Dense>
#include <Eigen/Sparse>



using Matrix = Eigen::Matrix<realtype, Eigen::Dynamic, Eigen::Dynamic>;
using SparseMatrix = Eigen::SparseMatrix<realtype>;
using Solver = Eigen::PartialPivLU< Matrix >;
using SolverSparse = Eigen::BiCGSTAB< SparseMatrix, Eigen::IncompleteLUT<realtype> >;

//-----------------------------------------------------------------------------------
//  Global parameters
//-----------------------------------------------------------------------------------
const realtype _G_kb = 1.37e5 ;
const realtype _G_kf = 5e-7 * _G_kb;
const realtype _G_k1 = 7.69e4 ;
const realtype _G_k2 = 1.40e4 ;
const realtype _G_k3 = 7.15e3 ;
const realtype _G_k4 = 1.74e3 ;
const unsigned int _G_M = 111;
const realtype _G_S = 11.3;
const realtype _G_T = 10.0;


template<typename T>
void
output2File(const std::string filename,
  const T & data
)
{
  std::ofstream file;
  file.open(filename);
  for (auto v : data){
    file << v << "\n";
  }
  file.close();
}


realtype
rFcn(const realtype size)
{
  return 8./3. * std::pow(size, 2./3.);
}


NanoSim::particleSystem<double>
defineUnfinalizedParticleSystem(const realtype ic, 
  const unsigned int maxSize, 
  const realtype factor,
  const realtype reductionFactor
)
{
  //--------------------------------
  //  Various constants used
  //--------------------------------
  const realtype kb = _G_kb ;
  const realtype kf = _G_kf;
  const realtype k1 = _G_k1 ;
  const realtype k2 = _G_k2 * factor;
  const realtype k3 = _G_k3 * factor;
  const realtype k4 = _G_k4 * factor;
  const unsigned int M = _G_M;
  const realtype S = _G_S;

  const realtype icA = ic * factor;

  //--------------------------------
  //  Setup for the system
  //--------------------------------
  NanoSim::particleSystem<double> my_rxns;

  //--------------------------------
  //  Function converting between number of atoms and diameter (assuming spheres of Iridium with bulk density)
  //--------------------------------
  std::function<realtype(const int)> atoms2diameter 
    = [](const int atoms){ return 0.3 * std::cbrt(1.0*atoms);};

  //--------------------------------
  //  Define a particle by specifying min/max size, size conversion, and if the model is reduced
  //--------------------------------
  my_rxns.defineParticle(3, maxSize, atoms2diameter, reductionFactor);

  //--------------------------------
  //  Add chemical reactions
  //--------------------------------
  my_rxns.addReaction({{1,"A"}},
    {{1,"As"}, {1,"L"}},
    kf * S * S);
  my_rxns.addReaction({{1,"As"}, {1,"L"}},
    {{1,"A"}},
    kb);

  //--------------------------------
  //  Nucleation is the reaction that creates the first particle
  //--------------------------------
  my_rxns.addNucleation({{1, "A"}, {2,"As"}},
    {{1,"L"}},
    k1);

  //--------------------------------
  //  Add particle growth
  //--------------------------------
  const std::function<realtype(const unsigned int)> growth_kernel = 
    [&](const unsigned int size){
      return rFcn(size) * (size <= M ? k2 : k3);
    };
  my_rxns.addGrowth({{1, "A"}},
    {{1,"L"}},
    growth_kernel);
  
  //--------------------------------
  //  Add particle agglomeration
  //--------------------------------
  const std::function<realtype(const unsigned int, const unsigned int)> agglom_kernel =
    [&](const unsigned int size1, const unsigned int size2){
      if (size1 > M || size2 > M){
        return 0.0;
      } else {
        return rFcn(size1) * rFcn(size2) * k4 * 0.0;
      }
    };
  my_rxns.addAgglomeration({},{},agglom_kernel);

  return my_rxns;
}


// Return
//  { simulation times,
//    simulation diameters,
//    simulation concentrations}
std::tuple<
  std::vector<double>,
  std::vector<realtype>,
  std::vector<realtype>
>
simulateRxnsSparse(const Eigen::Vector<realtype, Eigen::Dynamic> & parameters){
  const realtype factor = parameters(3);
  const realtype icA = parameters(0) * factor;
  const unsigned int max_size = parameters(1);

  auto my_rxns = defineUnfinalizedParticleSystem(icA, max_size, factor, 0.0);
  my_rxns.finalizeReactions();

  std::cout << "Full model has " << my_rxns.getNumberOfParticleBins() << " particles\n";
  std::cout << "Full model has " << my_rxns.getNumberOfReactions() << " reactions\n";
  
  
  const unsigned int n_solves = parameters(2);
  
  
  //  Solve ODEs
  const auto length = my_rxns.getNumberOfSpecies();
  NanoSim::eigenLinearAlgebraOperations<realtype, SparseMatrix> lin_alg;
  auto ic = lin_alg.createNewVector(length);
  N_VConst(0.0, ic);
  const auto idx_A = my_rxns.species_to_index_map["A"].vector_index;
  lin_alg.vectorInsert(ic, icA, idx_A);

  auto template_matrix = lin_alg.createNewMatrix(length, length);
  auto lin_solve = NanoSim::createLinearSolverEigenSparse<realtype, SolverSparse>();

  NanoSim::cvodeOptions<realtype> opts("SUNDIALS_errors.txt", 
    1e-8, 
    1e-14, 
    10000,
    CV_BDF, 
    SUNTRUE );


  // solve it a few times to benchmark the time it takes
  std::vector<double> solve_times;
  std::string prev_out = "";
  std::string out_text = "";

  out_text = "Solved 0/" + std::to_string(n_solves);
  std::cout << out_text << std::flush;

  // Vectors to hold the solution for comparison
  std::vector<realtype> diams;
  std::vector<realtype> concs;

  const realtype T = _G_T;

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
    // std::cout << err << "\n"
    //   << t << "\n";
    // end timing
    auto end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
    solve_times.push_back(elapsed);
    

    if (i == 0){
      // just PSD -- starts at index 3
      for (sunindextype bin=0; bin<=my_rxns.getNumberOfParticleBins(); ++bin){
        auto bin_info = my_rxns.species_to_index_map["__PARTICLE__" + std::to_string(bin)];
        auto idx = bin_info.vector_index;
        auto smallest_size = bin_info.smallest_size;
        auto n_binned = bin_info.n_binned_particles;

        const realtype conc = lin_alg.vectorGetValue(sol, idx)/n_binned;

        const int avg_size = smallest_size + (n_binned-1)/2;

        const realtype sizeL = 0.3 * std::pow(avg_size - 0.5, 1./3.);
        const realtype sizeC = 0.3 * std::pow(avg_size * 1.0, 1./3.);
        const realtype sizeR = 0.3 * std::pow(avg_size + 0.5, 1./3.); 
        const realtype dens = conc / (sizeR - sizeL);
        concs.push_back(dens);
        diams.push_back(sizeC); 
      }
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

  return {solve_times, diams, concs};
}



// Return
//  { simulation times,
//    simulation diameters,
//    simulation concentrations}
std::tuple<
  std::vector<double>,
  std::vector<realtype>,
  std::vector<realtype>
> 
simulateRxnsReduced(const Eigen::Vector<realtype, Eigen::Dynamic> & parameters){
  const realtype factor = parameters(3);
  const realtype kb = _G_kb ;
  const realtype kf = _G_kf;
  const realtype k1 = _G_k1 ;
  const realtype k2 = _G_k2 * factor;
  const realtype k3 = _G_k3 * factor;
  const realtype k4 = _G_k4 * factor;
  const unsigned int M = _G_M;
  const realtype S = _G_S;
  const realtype T = _G_T;

  const realtype icA = parameters(0) * factor;
  const unsigned int max_size = parameters(1);
  const unsigned int n_solves = parameters(2);

  
  NanoSim::particleSystem<double> my_rxns;

  std::function<realtype(const int)> atoms2diameter 
    = [](const int atoms){ return 0.3 * std::cbrt(1.0*atoms);};

  my_rxns.defineParticle(3, max_size, atoms2diameter, 0.01);

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
      return 8./3. * std::pow(size * 1.0, 2./3.) * (size <= M ? k2 : k3);
    };
  my_rxns.addGrowth({{1, "A"}},
    {{1,"L"}},
    growth_kernel);

  const std::function<realtype(const unsigned int, const unsigned int)> agglom_kernel =
    [&](const unsigned int size1, const unsigned int size2){
      if (size1 > M || size2 > M){
        return 0.0;
      } else {
        return 8./3. * std::pow(size1 * 1.0, 2./3.) * 8./3. * std::pow(size2 * 1.0, 2./3.) * k4;
      }
    };
  my_rxns.addAgglomeration({},{},agglom_kernel);

  my_rxns.finalizeReactions();

  std::cout << "Reduced model has " << my_rxns.getNumberOfParticleBins() << " particles\n";
  std::cout << "Reduced model has " << my_rxns.getNumberOfReactions() << " reactions\n";

  const auto length = my_rxns.getNumberOfSpecies();
  NanoSim::eigenLinearAlgebraOperations<realtype, Matrix> lin_alg;
  auto ic = lin_alg.createNewVector(length);
  N_VConst(0.0, ic);
  lin_alg.vectorInsert(ic, icA, 0);

  auto template_matrix = lin_alg.createNewMatrix(length, length);
  // std::cout << length << std::endl;
  auto lin_solve = NanoSim::createLinearSolverEigenDense<realtype, Solver>();

  NanoSim::cvodeOptions<realtype> opts("SUNDIALS_errors.txt", 
    1e-8, 
    1e-14, 
    10000,
    CV_BDF, 
    SUNTRUE );


  // solve it a few times to benchmark the time it takes
  std::vector<double> solve_times;
  std::string prev_out = "";
  std::string out_text = "";

  // Solution to ODE
  std::vector<realtype> diams;
  std::vector<realtype> concs;

  out_text = "Solved 0/" + std::to_string(n_solves);
  // std::cout << out_text << std::flush;

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
      // just PSD -- starts at index 3
      for (sunindextype bin=0; bin<=my_rxns.getNumberOfParticleBins(); ++bin){
        auto bin_info = my_rxns.species_to_index_map["__PARTICLE__" + std::to_string(bin)];
        auto idx = bin_info.vector_index;
        auto smallest_size = bin_info.smallest_size;
        auto n_binned = bin_info.n_binned_particles;

        const realtype conc = lin_alg.vectorGetValue(sol, idx)/n_binned;

        const int avg_size = smallest_size + (n_binned-1)/2;

        const realtype sizeL = 0.3 * std::pow(avg_size - 0.5, 1./3.);
        const realtype sizeC = 0.3 * std::pow(avg_size * 1.0, 1./3.);
        const realtype sizeR = 0.3 * std::pow(avg_size + 0.5, 1./3.); 
        const realtype dens = conc / (sizeR - sizeL);
        concs.push_back(dens);
        diams.push_back(sizeC); 
      }
    }
    N_VDestroy(sol);
    CVodeFree(&cvode_mem);
    prev_out = out_text;
    out_text = "Solved " + std::to_string(i+1) + "/" + std::to_string(n_solves);
    // std::cout << "\r" << out_text << std::flush;
  }
  // std::cout << "\n";

  N_VDestroy(ic);
  SUNLinSolFree(lin_solve);
  SUNMatDestroy(template_matrix);

  // std::ofstream outputfile;
  // std::string filename = "step05-densesolve-" + std::to_string(max_size) + "-particles-solvetimes.txt";
  // outputfile.open(filename.c_str());
  // for (auto t : solve_times){
  //   outputfile << t << "\n";
  // }
  // outputfile.close();

  return {solve_times, diams, concs};
}

/*fixme add return*/void 
simulateRxnseMoM(const Eigen::Vector<realtype, Eigen::Dynamic> & parameters){
  const realtype factor = parameters(3);
  const realtype kb = _G_kb ;
  const realtype kf = _G_kf;
  const realtype k1 = _G_k1 ;
  const realtype k2 = _G_k2 * factor;
  const realtype k3 = _G_k3 * factor;
  const realtype k4 = _G_k4 * factor;
  const unsigned int M = _G_M;
  const realtype S = _G_S;
  const realtype T = _G_T;
  const realtype deltaT = 0.00001; // interval for which to approximate density from moments

  /*
    eMoM variables for clarity
  */
  const realtype k_inflow = k3;
  const realtype k_eMoM = k3;
  const realtype delx = 0.3;
  const realtype xM = delx * std::pow( 2.0 * M, 1./3.); // FIXME

  const realtype icA = parameters(0) * factor;
//   const unsigned int max_size_filename = parameters(1);
  const unsigned int max_size = 2*M; // FIXME 2*M;
  const unsigned int n_solves = parameters(2);

  NanoSim::particleSystem<double> my_rxns;

  std::function<realtype(const int)> atoms2diameter 
    = [](const int atoms){ return 0.3 * std::cbrt(1.0*atoms);};

  my_rxns.defineParticle(3, max_size, atoms2diameter);

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
      return 8./3. * std::pow(size * 1.0, 2./3.) * (size <= M ? k2 : k3);
    };
  my_rxns.addGrowth({{1, "A"}},
    {{1,"L"}},
    growth_kernel);

  const std::function<realtype(const unsigned int, const unsigned int)> agglom_kernel =
    [&](const unsigned int size1, const unsigned int size2){
      if (size1 > M || size2 > M){
        return 0.0;
      } else {
        return 8./3. * std::pow(size1 * 1.0, 2./3.) * 8./3. * std::pow(size2 * 1.0, 2./3.) * k4;
      }
    };
  my_rxns.addAgglomeration({},{},agglom_kernel);
  // FIXME
  my_rxns.addeMoMGrowth("A",
    {{1, "L"}},
    k_inflow,
    k_eMoM);

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
    10000, 
    CV_BDF,
    SUNTRUE );


  // solve it a few times to benchmark the time it takes
  std::vector<double> solve_times;
  std::string prev_out = "";
  std::string out_text = "";

  out_text = "Solved 0/" + std::to_string(n_solves);
  // std::cout << out_text << std::flush;

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
    // for eMoM it is necessary to save output along the way to reconstruct the PSD
    // so for a fair comparison, that should be done here
    realtype t_iter = deltaT;
    std::vector<realtype> ui;
    std::vector<realtype> ai;

    // know the initial conditions
    ui.push_back(0.0);
    ai.push_back(8. * delx / 9. * k_eMoM * icA * deltaT);
    // TODO
    while (t_iter <= T){
      auto err = CVode(cvode_mem, t_iter, sol, &t, CV_NORMAL);
      // std::cout << "t = " << t ;
      // std::cout << "    err = " << err << "\n";

      const auto BM = lin_alg.vectorGetValue(sol, N_VGetLength(sol)-1); // last index
      const auto A  = lin_alg.vectorGetValue(sol, 0); // index 0 for this problem
      ui.push_back(
        3. * k_inflow * std::pow(max_size, 2./3.) * BM / delx / k_eMoM
      );

      ai.push_back(
        8. * delx / 9. * k_eMoM * A * deltaT
      );
      t_iter += deltaT;
    }
    // end timing
    auto end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
    solve_times.push_back(elapsed);
    

    if (i == 0){
      /*
        Output density from the full ODE part:
          Have # atoms = 3, 4, ..., max_size
          diameter = 0.3 * (# atoms) ^ (1/3)
        Approximate continuous density from discrete concentrations:
          Conci = concentration of cluster with i atoms
          Assume ODE system is discretization around (i - 1/2, i + 1/2)
          densityi = Conci / ( 0.3 * (i + 1/2)^(1/3) - 0.3 * (i - 1/2)^(1/3) )
      */
      for (sunindextype size=3;size<=max_size;++size){
        // outputfileSize << 0.3 * std::pow(1.0 * size, 1./3.) << "\n";
        
        const auto sizeL = 0.3 * std::pow(size - 0.5, 1./3.);
        const auto sizeR = 0.3 * std::pow(size + 0.5, 1./3.);
        auto idx = my_rxns.species_to_index_map.find("__PARTICLE__" + std::to_string(size))->second;
        // const auto conc = lin_alg.vectorGetValue(sol, idx.first); 
        // outputfilePSD << conc / (sizeR - sizeL) << "\n";
      }
      /*
        Output density from the eMoM part
      */
      
      // Doing a reverse cumulative sum of ai gives the particle diameters
      std::vector<realtype> xx(ai.size());
      std::partial_sum(ai.rbegin(), ai.rend(), xx.rbegin());
      // Now xM + xx = particles sizes corresponding to densities ui
      // Both are in descending order of particle size so write to file in reverse
      for (auto it = xx.rbegin(); it != xx.rend(); ++it){
        // std::cout << "size = xM + xx = " << xM << " + " << *it << "\n";s
        // outputfileSize << xM + *it << "\n";
      }
      for (auto it = ui.rbegin(); it != ui.rend(); ++it){
        // outputfilePSD << *it << "\n";
      }

    //   outputfilePSD.close();
    //   outputfileSize.close();
    }
    N_VDestroy(sol);
    CVodeFree(&cvode_mem);
    prev_out = out_text;
    out_text = "Solved " + std::to_string(i+1) + "/" + std::to_string(n_solves);
    // std::cout << "\r" << out_text << std::flush;
  }
  // std::cout << "\n";

  N_VDestroy(ic);
  SUNLinSolFree(lin_solve);
  SUNMatDestroy(template_matrix);

  // std::ofstream outputfile;
  // std::string filename = "step05-densesolve-eMoM-" + std::to_string(max_size_filename) + "-particles-solvetimes.txt";
  // outputfile.open(filename.c_str());
  // for (auto t : solve_times){
  //   outputfile << t << "\n";
  // }
  // outputfile.close();

//   return solve_times;
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


  auto sparse_results = simulateRxnsSparse(settings);

  auto sparse_times = std::get<0>(sparse_results);
  auto sparse_diams = std::get<1>(sparse_results);
  auto sparse_concs = std::get<2>(sparse_results);

  output2File("sparse_times.txt", sparse_times);
  output2File("sparse_diams.txt", sparse_diams);
  output2File("sparse_concs.txt", sparse_concs);

  for (auto t : sparse_times){
    std::cout << t << "  ";
  }
  std::cout << std::endl;

  // auto reduced_results = simulateRxnsReduced(settings);

  // auto reduced_times = std::get<0>(reduced_results);
  // auto reduced_diams = std::get<1>(reduced_results);
  // auto reduced_concs = std::get<2>(reduced_results);


  // std::cout << "Reduced solve times: ";
  // for (auto t : reduced_times){
  //   std::cout << t << "  ";
  // }
  // std::cout << std::endl;


  // simulateRxnseMoM(settings);

}