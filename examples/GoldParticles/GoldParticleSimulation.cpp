// -------------------------------------------------------------------------------
//  Outside libraries
// -------------------------------------------------------------------------------

#include "particleSystem.h"
#include "odeSolver.h"
#include "linearAlgebraEigen.h"
#include "linearSolverEigen.h"

#include <Eigen/Sparse>
#include <omp.h>

using Matrix = Eigen::Matrix<realtype, Eigen::Dynamic, Eigen::Dynamic>;
using Solver = Eigen::PartialPivLU< Matrix >;


// -------------------------------------------------------------------------------
//  C++ standard library
// -------------------------------------------------------------------------------

#include <vector>
#include <iostream>
#include <fstream>
#include <chrono>
#include <cmath>
#include <string>
#include <numeric>


// -------------------------------------------------------------------------------
//  Global settings to test the effects of different options
//  TODO -- explain
// -------------------------------------------------------------------------------

constexpr realtype reduction_threshold = 0.05;
constexpr int min_particle_atoms = 1;
constexpr int max_particle_atoms = 60000;


// -------------------------------------------------------------------------------
//  Custom class that gives meaning to the parameters vector
//  Also makes it easier to adjust the parameter vector without having to change
//  other functions
// -------------------------------------------------------------------------------
template <typename VectorType>
class myParameters
{
    public:
    // Require that parameters are provided
    myParameters() = delete;
    myParameters(const VectorType & parameters)
    : reduction_rate(parameters[0]),
      prob_a(parameters[1]),
      prob_b(parameters[2])
    {}

    // ------------------------------------
    //  How fast does the reduction reaction occur
    //      This cause particle nucleation
    // ------------------------------------
    const realtype reduction_rate;

    // ------------------------------------
    //  Constants in adhesion prob
    // ------------------------------------
    const realtype prob_a;
    const realtype prob_b;
};



// -------------------------------------------------------------------------------
//  Custom function that converts number of atoms in a particle to a diameter
// -------------------------------------------------------------------------------
realtype atoms_to_diameter(
    const realtype n_atoms
)
{
    // From Negishi paper on Gold particles
    // Basically just use the bulk density of gold and the mass of a single gold atom
    // to convert
    const double factor = 0.3186854328;
    return factor * std::pow( 1.0*n_atoms, 1./3.);
}

// -------------------------------------------------------------------------------
//  Custom function for the growth kernel
//      Just return a value of 0 if no growth is desired
//      Can use parameters or not
//          Just make sure parameter indexing is consistent with optimization algorithm
// -------------------------------------------------------------------------------
template <typename VectorType>
realtype my_growth_kernel(
    const unsigned int n_atoms,
    const VectorType & parameters
)
{
    return 0.0;
}


// -------------------------------------------------------------------------------
//  Custom function for the agglomeration kernel
//      Just return a value of 0 if no agglomeration is desired
//      Can use parameters or not
//          Just make sure parameter indexing is consistent with optimization algorithm
// -------------------------------------------------------------------------------
template <typename VectorType>
realtype my_agglomeration_kernel(
    const unsigned int n_atoms1,
    const unsigned int n_atoms2,
    const VectorType & parameters
)
{
    const auto diam1 = atoms_to_diameter(n_atoms1);
    // const auto volume1 = 4.0 / 3.0 * M_PI * std::pow(diam1/2.0, 3.0);

    const auto diam2 = atoms_to_diameter(n_atoms2);
    // const auto volume2 = 4.0 / 3.0 * M_PI * std::pow(diam2/2.0, 3.0);

    myParameters prm(parameters);

    // Probability-based model
    //  prob = a * exp(-b*min(d1^2,d2^2))
    // or min(1,prob) ????
    const realtype prob_a = prm.prob_a;
    const realtype prob_b = prm.prob_b;
    const realtype adhesion_prob = prob_a * std::exp( -prob_b * std::min(diam1*diam1, diam2*diam2) );

    const realtype kernel_pow = 1.0;
    const realtype brownian = (std::pow(diam1, 1./kernel_pow) + std::pow(diam2, 1./kernel_pow))
        * (std::pow(diam1, -1./kernel_pow) + std::pow(diam2, -1./kernel_pow));

    // ---------------------------------------------
    //  Compute the prefactor on the brownian kernel
    //      PF = kB*Temp/viscosity * AvoNum
    // ---------------------------------------------
    const double kB = 1.380649e-23;
    const double Temp = 273.15 + 20; // 20 degrees Celcius
    const double viscosity = 1.0016e-3; // Viscosity of water
    const double AvoNumber = 6.02214076e23;

    const realtype physical_constants = kB*Temp*AvoNumber/viscosity;

    return physical_constants * brownian * std::min(1.0, adhesion_prob);
}


// -------------------------------------------------------------------------------
//  Creates the chemical reaction system
// -------------------------------------------------------------------------------
template <typename VectorType>
NanoSim::particleSystem<realtype>
create_particle_system(
    const VectorType & parameters, 
    const realtype model_reduction_tolerance
)
{
    const std::function<realtype(const int)> atoms2diameter = [](const int atoms){return atoms_to_diameter(atoms);};
    const std::function<realtype(const unsigned int)> growth_kernel
        = [&](const unsigned int size){
            return my_growth_kernel(size, parameters);
        };
    const std::function<realtype(const unsigned int, const unsigned int)> agglom_kernel
        = [&](const unsigned int size1, const unsigned int size2){
            return my_agglomeration_kernel(size1, size2, parameters);
        };

    const myParameters model_prm(parameters);

    NanoSim::particleSystem<realtype> my_rxns;
    // For gold particles, 60000 atoms is ~12nm which is a find place to stop tracking
    // based on the data we have
    my_rxns.defineParticle(min_particle_atoms, max_particle_atoms, atoms2diameter, model_reduction_tolerance);

    my_rxns.addNucleation(
        {{1,"NaAuCl4"}, {1,"NaBH4"}},
        {{1,"BoricAcid"}},
        model_prm.reduction_rate
    );

    my_rxns.addGrowth(
        {{1,"NaAuCl4"}},
        {},
        growth_kernel
    );

    my_rxns.addAgglomeration(
        {},
        {},
        agglom_kernel
    );

    my_rxns.finalizeReactions();
    return my_rxns;
}


// -------------------------------------------------------------------------------
//  Splits a string based on a delimiter
//  Used to import data for optimization
// -------------------------------------------------------------------------------
std::vector<std::string> split(std::string s, std::string delimiter) {
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    std::string token;
    std::vector<std::string> res;

    while ((pos_end = s.find(delimiter, pos_start)) != std::string::npos) {
        token = s.substr (pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back (token);
    }

    res.push_back (s.substr (pos_start));
    return res;
}



// -------------------------------------------------------------------------------
//  FIXME
// -------------------------------------------------------------------------------
std::vector< std::vector<realtype> > solve_odes(
    const std::vector<double> &X,
    const std::vector<realtype> & solve_times
)
{
    // GLOBAL SETTING
    auto my_rxns = create_particle_system(X, reduction_threshold);

    // Set up ODE solver
    const auto dim = my_rxns.getNumberOfSpecies();
    NanoSim::eigenLinearAlgebraOperations<realtype, Matrix> lin_alg;
    auto ic = lin_alg.createNewVector(dim);
    const auto idx_NaAuCl4 = my_rxns.species_to_index_map["NaAuCl4"].vector_index;
    const auto idx_NaBH4   = my_rxns.species_to_index_map["NaBH4"].vector_index;
    lin_alg.vectorInsert(ic, 0.0001, idx_NaAuCl4);
    lin_alg.vectorInsert(ic, 0.0003, idx_NaBH4);

    auto matrix = lin_alg.createNewMatrix(dim,dim);
    auto lin_solver = NanoSim::createLinearSolverEigenDense<realtype, Solver>();

    NanoSim::cvodeOptions<realtype> opts(
        "SUNDIALS_errors.txt",
        1e-8,
        1e-14,
        10000,
        CV_BDF,
        SUNTRUE
    );

    auto cvode_mem = NanoSim::prepareODESolver<realtype>(
        ic,
        matrix,
        lin_solver,
        opts
    );

    std::pair<
        NanoSim::particleSystem<realtype>*,
        NanoSim::abstractLinearAlgebraOperations<realtype>*
    > data_pair = {&my_rxns, &lin_alg};
    void* user_data = static_cast<void*>(&data_pair);
    auto flag = CVodeSetUserData(cvode_mem, user_data);

    auto sol = N_VClone(ic);
    realtype t;

    // Create vector with all of the particle diameters
    std::vector< realtype > simulation_diameters;
    for (int idx=0; idx<my_rxns.getNumberOfParticleBins(); ++idx){
        // Gather metadata about the particle bin
        auto bin_info = my_rxns.species_to_index_map["__PARTICLE__"+std::to_string(idx)];
        auto vec_idx = bin_info.vector_index;
        auto smallest_size = bin_info.smallest_size;
        auto n_binned = bin_info.n_binned_particles;

        // Compute the average particle size (# atoms) in the bin
        const realtype bin_conc = lin_alg.vectorGetValue(sol, vec_idx);
        const realtype conc = std::max(0.0, bin_conc / n_binned);
        const realtype avg_size = smallest_size + (n_binned-1.)/2.;

        // Convert the average particle size to a diameter and add to the diameter vector
        const realtype diam = atoms_to_diameter(avg_size);
        simulation_diameters.push_back(diam);
    }


    // Solve ODE
    std::vector< std::vector<realtype> > fcn_result;
    fcn_result.push_back(simulation_diameters);
    for (auto T : solve_times){
        auto err = CVode(cvode_mem, T, sol, &t, CV_NORMAL);

        // Retrieve solution
        // FIXME -- make this a function in NanoSim
        std::vector<realtype> sim_diam;
        std::vector<realtype> sim_dens;
        for (int idx=0; idx<my_rxns.getNumberOfParticleBins();++idx){
            auto bin_info = my_rxns.species_to_index_map["__PARTICLE__"+std::to_string(idx)];
            auto vec_idx = bin_info.vector_index;
            auto smallest_size = bin_info.smallest_size;
            auto n_binned = bin_info.n_binned_particles;

            const realtype bin_conc = lin_alg.vectorGetValue(sol, vec_idx);
            const realtype conc = std::max(0.0, bin_conc / n_binned);
            const realtype avg_size = smallest_size + (n_binned-1.)/2.;

            const realtype diamL = atoms_to_diameter(avg_size - 0.5);
            const realtype diamR = atoms_to_diameter(avg_size + 0.5);
            const realtype diamC = atoms_to_diameter(avg_size);
            // Also volume weight!
            const realtype dens = std::pow(diamC,3.0) * conc / (diamR - diamL);
            sim_dens.push_back(dens);
        }
        // Add simulated density to the final result
        fcn_result.push_back(sim_dens);
    }

    return fcn_result;
}


// -------------------------------------------------------------------------------
//  Main function to be executed
// -------------------------------------------------------------------------------
int main(int argc, char **argv){
    // Solve ODEs for an example set of parameters
    auto start = std::chrono::high_resolution_clock::now();

    // Command line arguments
    std::vector<std::string> argList(argv + 0, argv + argc);

    std::string parameter_file_name = "parameters.inp";
    std::string diameter_file_name = "diameters.out";
    std::string density_file_name = "density.out";
    if (argc > 1){
        parameter_file_name = argList[1];
    }
    if (argc > 2){
        diameter_file_name = argList[2];
    }
    if (argc > 3){
        density_file_name = argList[3];
    }

    std::cout << "Model parameters read from: " << parameter_file_name << "\n";
    std::cout << "Simulated diameters written to: " << diameter_file_name << "\n";
    std::cout << "Simulated densities written to: " << density_file_name << "\n";

    // Read parameters from specified file
    std::vector< realtype > prm;
    std::ifstream infile(parameter_file_name);
    std::string line;
    while (std::getline(infile, line)){
        prm.push_back(
            std::stod(line)
        );
    }

    for (auto p : prm){
        std::cout << p << "\n";
    }

    auto sim = solve_odes(
        prm,
        {2,5,7,10,20,30,60,120,180,240,300} /* FIXME */
    );

    std::ofstream diam_file;
    diam_file.open(diameter_file_name);
    for (auto d : sim[0]){
        diam_file << d << "\n";
    }
    diam_file.close();

    std::ofstream q_file;
    q_file.open(density_file_name);

    for (int i=1; i<sim.size();++i){
        for (auto q : sim[i]){
            q_file << q << "\t";
        }
        q_file << "\n";
    }
    q_file.close();

    auto end = std::chrono::high_resolution_clock::now();
    const double runtime = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
    std::cout << "ODE solve time: "
        << runtime 
        << " ms\n";
    
    return 0;
}

