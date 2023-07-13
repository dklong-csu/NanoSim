#include "particleSystem.h"
#include "odeSolver.h"
#include "linearAlgebraEigen.h"
#include "linearSolverEigen.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <string>

using Matrix = Eigen::Matrix<realtype, Eigen::Dynamic, Eigen::Dynamic>;
using Solver = Eigen::PartialPivLU< Matrix >;

/*
  This tutorial shows how to wrap simulating chemical reactions into a 
  function that is parametrized in some manner. In this case, the function
  will take reaction rates as parameters.
  TODO add more stuff
*/


// TODO -- explain stuff
void simulateRxns(const Eigen::Vector<realtype, Eigen::Dynamic> & parameters,
  std::string filename /* just for the purposes of this tutorial
                          you probably wouldn't want to print to a file in this function in a real application*/
                          ){
  // Most of this function is copy-pasted from step01-solveRxnODE.cpp
  // The documentation in that file is removed here and only changes are documented
  // so that the differences are more obvious.
  NanoSim::particleSystem<double> my_rxns;

  // Parameters passed are the two reaction rates!
  const realtype k1 = parameters(0); // This is problem dependent!
  const realtype k2 = parameters(1); // This is problem dependent!
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
  const realtype dt = 0.01;
  const realtype endtime = 2.0;
  realtype t = 0.0;
  while (t < endtime){
    t += dt;
    times.push_back(t);
  }

  std::ofstream outputfile;
  std::string times_filename = filename + "-output-times.txt"; 
  outputfile.open(times_filename.c_str());
  outputfile << 0.0 << "\n"; 
  for (const auto t : times){
    outputfile << t << "\n";
  }
  outputfile.close();

  std::string conc_filename = filename + "-output-concentrations.txt";
  outputfile.open(conc_filename.c_str()); 
  outputfile << lin_alg.vectorGetValue(ic, 0)
    << "    " << lin_alg.vectorGetValue(ic, 1)
    << "    " << lin_alg.vectorGetValue(ic, 2)
    << "    " << lin_alg.vectorGetValue(ic, 3)
    << "\n";

  for (const auto solve_time : times){
    
    N_Vector sol = ic->ops->nvclone(ic);
    realtype t;

    auto err = CVode(cvode_mem, solve_time, sol, &t, CV_NORMAL);

    outputfile << lin_alg.vectorGetValue(sol, 0)
      << "    " << lin_alg.vectorGetValue(sol, 1)
      << "    " << lin_alg.vectorGetValue(sol, 2)
      << "    " << lin_alg.vectorGetValue(sol, 3)
      << "\n";

    sol->ops->nvdestroy(sol); 
  }
  outputfile.close(); 


  CVodeFree(&cvode_mem);
  N_VDestroy(ic);
  SUNLinSolFree(lin_solve);
  SUNMatDestroy(template_matrix);
}


int main(){
  Eigen::Vector<realtype, Eigen::Dynamic> prm1(2);
  prm1 << 2.0, 3.0; // same as in step01
  std::string filename1 = "step02-01";
  simulateRxns(prm1, filename1);

  Eigen::Vector<realtype, Eigen::Dynamic> prm2(2);
  prm2 << 0.2, 0.3; // Slower, but at same ratio
  std::string filename2 = "step02-02";
  simulateRxns(prm2, filename2);
}