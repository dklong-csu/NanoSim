// -------------------------------------------------------------------------------
//  Outside libraries
// -------------------------------------------------------------------------------

#include "particleSystem.h"
#include "odeSolver.h"
#include "linearAlgebraEigen.h"
#include "linearSolverEigen.h"

#include <Eigen/Sparse>
#include <omp.h>
#include <matplot/matplot.h>
#include <boost/math/interpolators/makima.hpp>

using Matrix = Eigen::Matrix<realtype, Eigen::Dynamic, Eigen::Dynamic>;
using Solver = Eigen::PartialPivLU< Matrix >;

#include <nlopt.hpp>

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
      agglomeration_rate(parameters[1]),
      prob_a(parameters[2]),
      prob_b(parameters[3])
    {}

    // ------------------------------------
    //  How fast does the reduction reaction occur
    //      This cause particle nucleation
    // ------------------------------------
    const realtype reduction_rate;

    // ------------------------------------
    //  Additional factor multiplied to agglomeration rate
    // ------------------------------------
    const realtype agglomeration_rate;

    // ------------------------------------
    //  Constants in adhesion prob
    // ------------------------------------
    const realtype prob_a;
    const realtype prob_b;
};


// -------------------------------------------------------------------------------
//  Struct to hold measurement data
// -------------------------------------------------------------------------------
struct MeasurementData
{
    std::vector<realtype> measurement_times;
    std::vector< std::vector<realtype> > measured_psds;
    std::vector< realtype > measured_diameters;
    bool compare_cdf = false;
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
    // if ( (n_atoms1 > 1 && n_atoms1 < 10) || (n_atoms2 > 1 && n_atoms2 < 10) ){
    //     return 0.0;
    // } else{
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


        const auto RR = prm.agglomeration_rate;

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

        // Correction to convert from density to concentration??
        // const realtype conv_factor = std::pow(M_PI / 186., -2.0);


        return physical_constants * brownian * std::min(1.0, adhesion_prob);
    // }
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
    my_rxns.defineParticle(1, 60000, atoms2diameter, model_reduction_tolerance);

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
//  Loads data for optimization
// -------------------------------------------------------------------------------

MeasurementData importData(
    const std::string filename, 
    const int skip_first_N_lines,
    const realtype min_diam = std::numeric_limits<realtype>::lowest(),
    const realtype max_diam = std::numeric_limits<realtype>::max(),
    const bool negative_density_to_zero = true
)
{
    MeasurementData my_data;

    std::ifstream infile(filename);
    std::string line;

    // Skip the first few lines
    for (int l=0;l<skip_first_N_lines;++l){
        std::getline(infile, line);
    }

    // Assume the first line after skipped give times
    std::getline(infile, line);
    const auto times_string = split(line, "s\t");
    for (int i=0; i<times_string.size();++i){
        auto t = times_string[i];
        my_data.measurement_times.push_back( std::stod(t) );
    }

    // Data in the format
    // diam density density density ...
    const auto n_times = my_data.measurement_times.size();
    std::vector< std::vector<realtype> > data_psds(n_times);
    std::vector< realtype > data_diams;

    while (std::getline(infile,line)){
        const auto split_line = split(line, "\t");
        // First is the diameter
        // Only import data between specified diameters
        const auto diam = std::stod(split_line[0]);
        
        if (diam >= min_diam && diam <= max_diam){
            data_diams.push_back( diam );

            // The rest are densities corresponding to measurement times
            for (int i=1; i<split_line.size(); ++i){
                auto dens = std::stod( split_line[i] );
                if (negative_density_to_zero && dens < 0){
                    dens = 0.0;
                }
                data_psds[i-1].push_back( dens );
            }
        }
    }

    my_data.measured_diameters = data_diams;
    my_data.measured_psds = data_psds;

    return my_data;
}


// -------------------------------------------------------------------------------
//  Trapezoidal rule for integration
// -------------------------------------------------------------------------------
realtype trapz(
    const std::vector<realtype> & x,
    const std::vector<realtype> & y
)
{
    realtype integral = 0.0;
    for (int i=0;i<x.size()-1;++i){
        const realtype deltaX = x[i+1] - x[i];
        const realtype avg_y = (y[i+1] + y[i])/2.0;
        integral += avg_y*deltaX;
    }
    return integral;
}



// -------------------------------------------------------------------------------
//  Optimization function
// -------------------------------------------------------------------------------
int __OPTIMIZATION_COUNTER__ = 0;
double obj_fcn(
    const std::vector<double> &X,
    std::vector<double> &grad,
    void* f_data
)
{
    auto start = std::chrono::high_resolution_clock::now();

    const MeasurementData* measured_data = static_cast<MeasurementData*>(f_data);

    auto my_rxns = create_particle_system(X, 0.05);

    // Set up ODE solver
    const auto dim = my_rxns.getNumberOfSpecies();
    NanoSim::eigenLinearAlgebraOperations<realtype, Matrix> lin_alg;
    auto ic = lin_alg.createNewVector(dim);
    const auto idx_NaAuCl4 = my_rxns.species_to_index_map["NaAuCl4"].vector_index;
    const auto idx_NaBH4   = my_rxns.species_to_index_map["NaBH4"].vector_index;
    const auto idx_Au1    = my_rxns.species_to_index_map["__PARTICLE__0"].vector_index;
    const auto idx_Au10    = my_rxns.species_to_index_map["__PARTICLE__9"].vector_index;
    // FIXME -- nucleation event creates size 10 (happens immediately)
    lin_alg.vectorInsert(ic, 0.0001, idx_NaAuCl4);
    lin_alg.vectorInsert(ic, 0.0003, idx_NaBH4);
    // lin_alg.vectorInsert(ic, 0.005, idx_Au1);
    // lin_alg.vectorInsert(ic, 0.001, idx_Au10);

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

    // Solve at each time FIXME #1
    realtype obj_value = 0.0;
    for (int time = 1; time < measured_data->measurement_times.size();++time){
        const auto T = measured_data->measurement_times[time];
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
            sim_diam.push_back(diamC);
        }

        std::vector<realtype> sim_cdf(sim_dens.size());
        sim_cdf[0] = 0.0;
        for (int idx=0;idx<sim_cdf.size()-1;++idx){
            // integrate one trapezoid
            const auto integral = trapz(
                {sim_diam[idx], sim_diam[idx+1]},
                {sim_dens[idx], sim_dens[idx+1]}
            );
            sim_cdf[idx+1] = integral;
        }
        std::partial_sum(sim_cdf.begin(), sim_cdf.end(), sim_cdf.begin());
        const auto max_denssum = sim_cdf.back();
        for (int idx=0;idx<sim_cdf.size();++idx){
            sim_cdf[idx] /= max_denssum;
            sim_dens[idx] /= max_denssum;
        }

        std::vector<realtype> data_cdf(measured_data->measured_diameters.size());
        data_cdf[0] = 0.0;
        for (int idx=0;idx<data_cdf.size()-1;++idx){
            // Integrate one trapezoid
            const auto integral = trapz(
                {measured_data->measured_diameters[idx], measured_data->measured_diameters[idx+1]},
                {measured_data->measured_psds[time][idx], measured_data->measured_psds[time][idx+1]}
            );
            data_cdf[idx+1] = integral;
        }
        std::partial_sum(data_cdf.begin(), data_cdf.end(), data_cdf.begin());
        const auto max_datacdf = data_cdf.back();
        for (int idx=0; idx<data_cdf.size(); ++idx){
            data_cdf[idx] /= max_datacdf;
        }

        // ---------------------------------------------
        //  Wasserstein is for CDF
        //  L2 is for PDF
        // ---------------------------------------------
        if (measured_data->compare_cdf){

            auto sim_interp = boost::math::interpolators::makima(
                std::move(sim_diam),
                std::move(sim_cdf)
            );

            std::vector<realtype> data_cdf(measured_data->measured_diameters.size());
            data_cdf[0] = 0.0;
            for (int idx=0;idx<data_cdf.size()-1;++idx){
                // Integrate one trapezoid
                const auto integral = trapz(
                    {measured_data->measured_diameters[idx], measured_data->measured_diameters[idx+1]},
                    {measured_data->measured_psds[time][idx], measured_data->measured_psds[time][idx+1]}
                );
                data_cdf[idx+1] = integral;
            }
            std::partial_sum(data_cdf.begin(), data_cdf.end(), data_cdf.begin());
            const auto max_datacdf = data_cdf.back();
            for (int idx=0; idx<data_cdf.size(); ++idx){
                data_cdf[idx] /= max_datacdf;
            }

            std::vector<realtype> cdf_diff;
            for (int idx=0; idx<measured_data->measured_diameters.size();++idx){
                const realtype diam = measured_data->measured_diameters[idx];
                const realtype sim_cdf_value = sim_interp(diam);
                const realtype data_cdf_value = data_cdf[idx];
                const realtype diff = std::abs(sim_cdf_value - data_cdf_value);
                cdf_diff.push_back(diff);
            }

            const auto wasserstein = trapz(
                measured_data->measured_diameters,
                cdf_diff
            );
            obj_value += wasserstein;
        } else {
            auto sim_interp = boost::math::interpolators::makima(
                std::move(sim_diam),
                std::move(sim_dens)
            );

            realtype L2 = 0.0;
            for (int i=0; i<measured_data->measured_diameters.size(); ++i){
                const realtype diam = measured_data->measured_diameters[i];
                const realtype interp_sim_dense = sim_interp(diam);
                const realtype data_psd_value = measured_data->measured_psds[time][i] / max_datacdf;
                L2 += std::pow(data_psd_value - interp_sim_dense, 2.0);
            }
            obj_value += L2;
        }
    }
    // Delete SUNDIALS stuff
    N_VDestroy(ic);
    N_VDestroy(sol);
    SUNLinSolFree(lin_solver);
    SUNMatDestroy(matrix);


    auto end = std::chrono::high_resolution_clock::now();
    const double runtime = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
        
    
    ++__OPTIMIZATION_COUNTER__;
    std::cout << "Fcn eval #" << __OPTIMIZATION_COUNTER__ << ": " 
        << "f(";
    for (auto prmval : X){
        std::cout << prmval << ",";
    }
    std::cout << ") = " << obj_value << " (" << runtime/1000 << "s)" << std::endl;
    return obj_value;
}


// -------------------------------------------------------------------------------
//  Plot solutions
// -------------------------------------------------------------------------------

void plot_data_vs_sim(
    const std::vector<double> &X,
    const MeasurementData* measured_data,
    const std::string filename 
)
{
    auto my_rxns = create_particle_system(X, 0.01);


    // Set up ODE solver
    const auto dim = my_rxns.getNumberOfSpecies();
    NanoSim::eigenLinearAlgebraOperations<realtype, Matrix> lin_alg;
    auto ic = lin_alg.createNewVector(dim);
    const auto idx_NaAuCl4 = my_rxns.species_to_index_map["NaAuCl4"].vector_index;
    const auto idx_NaBH4   = my_rxns.species_to_index_map["NaBH4"].vector_index;
    const auto idx_Au1    = my_rxns.species_to_index_map["__PARTICLE__0"].vector_index;
    const auto idx_Au10    = my_rxns.species_to_index_map["__PARTICLE__9"].vector_index;
    // FIXME -- nucleation event creates size 10 (happens immediately)
    lin_alg.vectorInsert(ic, 0.0001, idx_NaAuCl4);
    lin_alg.vectorInsert(ic, 0.0003, idx_NaBH4);
    // lin_alg.vectorInsert(ic, 0.005, idx_Au1);
    // lin_alg.vectorInsert(ic, 0.001, idx_Au10);

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

    // Solve at each time
    for (int time = 1; time < measured_data->measurement_times.size();++time){
        const auto T = measured_data->measurement_times[time];
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
            // std::cout << "__Particle__" << idx
            //     << ": " << smallest_size << " atoms - "
            //     << smallest_size+n_binned-1 << " atoms\n";

            const realtype bin_conc = lin_alg.vectorGetValue(sol, vec_idx);
            const realtype conc = std::max(bin_conc / n_binned,0.0);
            const realtype avg_size = smallest_size + (n_binned-1.)/2.;

            const realtype diamL = atoms_to_diameter(avg_size - 0.5);
            const realtype diamR = atoms_to_diameter(avg_size + 0.5);
            const realtype diamC = atoms_to_diameter(avg_size);
            // Also volume weight!
            const realtype dens = std::pow(diamC,3.0) * conc / (diamR - diamL);
            sim_dens.push_back(dens);
            sim_diam.push_back(diamC);
        }      

        std::vector<realtype> sim_cdf(sim_dens.size());
        sim_cdf[0] = 0.0;
        for (int idx=0;idx<sim_cdf.size()-1;++idx){
            // integrate one trapezoid
            const auto integral = trapz(
                {sim_diam[idx], sim_diam[idx+1]},
                {sim_dens[idx], sim_dens[idx+1]}
            );
            sim_cdf[idx+1] = integral;
        }
        std::partial_sum(sim_cdf.begin(), sim_cdf.end(), sim_cdf.begin());
        const auto max_denssum = sim_cdf.back();
        for (int idx=0;idx<sim_cdf.size();++idx){
            sim_cdf[idx] /= max_denssum;
            sim_dens[idx] /= max_denssum;
        }

        std::vector<realtype> data_cdf(measured_data->measured_diameters.size());
        data_cdf[0] = 0.0;
        for (int idx=0;idx<data_cdf.size()-1;++idx){
            // Integrate one trapezoid
            const auto integral = trapz(
                {measured_data->measured_diameters[idx], measured_data->measured_diameters[idx+1]},
                {measured_data->measured_psds[time][idx], measured_data->measured_psds[time][idx+1]}
            );
            data_cdf[idx+1] = integral;
        }
        std::partial_sum(data_cdf.begin(), data_cdf.end(), data_cdf.begin());
        const auto max_datacdf = data_cdf.back();
        std::vector<realtype> data_dens = measured_data->measured_psds[time];
        for (int idx=0; idx<data_cdf.size(); ++idx){
            data_cdf[idx] /= max_datacdf;
            data_dens[idx] /= max_datacdf;
        }


        const auto max_dens_in_sim = std::max_element(sim_dens.begin(), sim_dens.end());
        // std::cout << "----------------------------------------\n"
        //     << "Simulation information at time "
        //     << (int)T << " seconds\n"
        //     << "----------------------------------------\n";
        // std::cout << "[NaAuCl4]   = " << lin_alg.vectorGetValue(sol, idx_NaAuCl4) << "\n";
        // std::cout << "[NaBH4]     = " << lin_alg.vectorGetValue(sol, idx_NaBH4) << "\n";
        // std::cout << "max density = " << *max_dens_in_sim << "\n";
        // std::cout << "----------------------------------------\n\n";



        // ------------------------------------------------------------------
        //  CDF plots
        // ------------------------------------------------------------------

        auto plot_sim = matplot::plot(sim_diam, sim_cdf, "o");
        plot_sim->display_name("Simulation");

        matplot::hold(true);

        auto plot_data = matplot::plot(
            measured_data->measured_diameters,
            data_cdf,
            "-");
        plot_data->display_name("Data");

        auto lgd = matplot::legend(true);
        lgd->location(matplot::legend::general_alignment::bottomright);
        matplot::title("PSD at " + std::to_string((int)T) + "s");
        matplot::xlabel("Diameter / nm");
        matplot::ylabel("Cumulative Number Density / ");

        matplot::hold(false);

        matplot::save(filename + ".cdf." + std::to_string(time) + ".png");

        // ------------------------------------------------------------------
        //  Density plots
        // ------------------------------------------------------------------

        auto plot_psd_sim = matplot::plot(sim_diam, sim_dens,"o");
        plot_psd_sim->display_name("Simulation");
        matplot::hold(true);

        auto plot_psd_data = matplot::plot(
            measured_data->measured_diameters,
            data_dens,
            "-"
        );
        plot_psd_data->display_name("Data");

        auto lgd2 = matplot::legend(true);
        lgd2->location(matplot::legend::general_alignment::topright);
        matplot::title("PSD at " + std::to_string((int)T) + "s");
        matplot::xlabel("Diameter / nm");
        matplot::ylabel("Volume-Weighted Number Density / nm^{-1}");

        matplot::ylim({0,1});

        matplot::save(filename + ".psd." + std::to_string(time) + ".png");
        matplot::hold(false);
    }
    // Delete SUNDIALS stuff
    N_VDestroy(ic);
    N_VDestroy(sol);
    SUNLinSolFree(lin_solver);
    SUNMatDestroy(matrix);
}


void export_sim(
    const std::vector<double> &X,
    const realtype ratio,
    const std::string filename 
)
{
    auto my_rxns = create_particle_system(X, 0.01);

    // Set up ODE solver
    const auto dim = my_rxns.getNumberOfSpecies();
    NanoSim::eigenLinearAlgebraOperations<realtype, Matrix> lin_alg;
    auto ic = lin_alg.createNewVector(dim);
    const auto idx_NaAuCl4 = my_rxns.species_to_index_map["NaAuCl4"].vector_index;
    const auto idx_NaBH4   = my_rxns.species_to_index_map["NaBH4"].vector_index;
    lin_alg.vectorInsert(ic, 0.01, idx_NaAuCl4);
    lin_alg.vectorInsert(ic, ratio*0.01, idx_NaBH4);

    // std::cout << "NaBH4 ic: " << lin_alg.vectorGetValue(ic, idx_NaBH4) << std::endl;

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

    // Solve to 5 minutes (300 seconds)
    const realtype T = 300;
    auto err = CVode(cvode_mem, T, sol, &t, CV_NORMAL);

    // Retrieve simulation diameters and concentrations
    std::vector<realtype> sim_diam;
    std::vector<realtype> sim_dens;

    for (int idx=0; idx<my_rxns.getNumberOfParticleBins(); ++idx){
        auto bin_info = my_rxns.species_to_index_map["__PARTICLE__"+std::to_string(idx)];
        auto vec_idx = bin_info.vector_index;
        auto smallest_size = bin_info.smallest_size;
        auto n_binned = bin_info.n_binned_particles;

        const realtype bin_conc = lin_alg.vectorGetValue(sol, vec_idx);
        const realtype conc = std::max(bin_conc / n_binned,0.0);
        const realtype avg_size = smallest_size + (n_binned-1.)/2.;

        const realtype diamL = atoms_to_diameter(avg_size - 0.5);
        const realtype diamR = atoms_to_diameter(avg_size + 0.5);
        const realtype diamC = atoms_to_diameter(avg_size);
        // Also volume weight!
        const realtype dens = std::pow(diamC,3.0) * conc / (diamR - diamL);
        sim_dens.push_back(dens);
        sim_diam.push_back(diamC);
    }

    // Integrate density so we can normalize the area to 1
    const auto area = trapz(sim_diam, sim_dens);
    for (auto & d : sim_dens){
        d /= area;
    }

    // Output to txt files
    const std::string fname_diam = filename + ".diameters.txt";
    const std::string fname_dens = filename + ".density.txt";

    std::ofstream diam_file;
    diam_file.open(fname_diam);
    for (auto d : sim_diam){
        diam_file << d << "\n";
    }
    diam_file.close();

    std::ofstream dens_file;
    dens_file.open(fname_dens);
    for (auto d : sim_dens){
        dens_file << d << "\n";
    }
    dens_file.close();

    
    // Delete SUNDIALS stuff
    N_VDestroy(ic);
    N_VDestroy(sol);
    SUNLinSolFree(lin_solver);
    SUNMatDestroy(matrix);
}

// -------------------------------------------------------------------------------
//  Main function to be executed
// -------------------------------------------------------------------------------
int main(){
    const realtype min_size = atoms_to_diameter(1);
    std::cout << "min diam=" << min_size;
    const realtype max_size = atoms_to_diameter(60000);
    std::cout << "\tmax diam=" << max_size << std::endl;
    // auto psd_data = importData("data.txt",2,min_size,max_size);
    auto psd_data = importData("Au_quench_corediameter_hplc.txt",0,min_size,max_size);

    // ---------------------------------------------------------------------
    //  false --> Computes L2 norm between particle size distribution
    //  true  --> Computes Wasserstein distance between cumulative size distribution
    // ---------------------------------------------------------------------
    psd_data.compare_cdf = false;

    void* opt_data = static_cast<void*>(&psd_data);


    const int n_dims = 4;
    nlopt::opt global_optimizer(nlopt::GN_ESCH, n_dims);
    global_optimizer.set_min_objective(obj_fcn, opt_data);
    
    // ---------------------------------------------------------------------
    //  Parameters: {reduct, R, a, b}
    // ---------------------------------------------------------------------
    const std::vector<double> lb = {0, 1, 0, 0};
    const std::vector<double> ub = {1e10, 1, 1, 1.5};

    // FIXME #3
    const int max_evals = 100;
    const double max_time_seconds = 36000;

    global_optimizer.set_lower_bounds(lb);
    global_optimizer.set_upper_bounds(ub);

    global_optimizer.set_maxeval( max_evals );
    global_optimizer.set_maxtime( max_time_seconds );

    global_optimizer.set_ftol_abs(1e-8);

    // ---------------------------------------------------------------
    //  BEST OBSERVED SO FAR
    //      For fitting to data times 8-11
    //          Minimum found: f(5.83245e+06,470.49,1,0.0623937,) = 6.0226
    //      For fitting to all data times
    //          Minimum found: f(8.54958e+06,11981.2,1,0.393163,) = 121.624
    //      For forcing agglom rate less than diffusion (all times)
    // ---------------------------------------------------------------

    // FIXME #2
    std::vector<double> x = {1.33529e+09,1,0.484272,0.394759};
    double fcn_val;

    std::cout << "\n============================================================================\n"
        <<         "                             GLOBAL OPTIMIZATION\n"
        <<         "============================================================================\n\n";
    global_optimizer.optimize(x, fcn_val);

    // ---------------------------------------------------------------
    //  Now polish global optimization witih local optimization
    // ---------------------------------------------------------------
    nlopt::opt local_optimizer(nlopt::LN_BOBYQA, n_dims);
    local_optimizer.set_min_objective(obj_fcn, opt_data);
    
    // ---------------------------------------------------------------------
    //  Parameters: {reduction rate, agglom rate, logistic_k, logistic_x0, min_prob}
    // ---------------------------------------------------------------------

    local_optimizer.set_lower_bounds(lb);
    local_optimizer.set_upper_bounds(ub);

    local_optimizer.set_maxeval( max_evals );
    local_optimizer.set_maxtime( max_time_seconds );

    local_optimizer.set_ftol_abs(1e-8);

    std::cout << "\n============================================================================\n"
        <<         "                             LOCAL OPTIMIZATION\n"
        <<         "============================================================================\n\n";
    local_optimizer.optimize(x, fcn_val);


    std::cout << "\n============================================================================\n"
        <<         "                                   RESULT\n"
        <<         "============================================================================\n\n";


    std::cout << "Minimum found: f(";
    for (auto v : x){
        std::cout << v << ",";
    }
    std::cout << ") = " << fcn_val << "\n";

    std::cout << "\nFor Matlab plot\n";
    std::vector<std::string> var_list = {"R", "a", "b"};
    for (int i=1;i<x.size();++i){
        std::cout << var_list[i-1] << " = " << x[i] << ";" << std::endl;
    }

    plot_data_vs_sim(x, &psd_data, "psd_fit");

    // export_sim(x, 3.0, "SIMULATION_PSD_ratio3");
    // export_sim(x, 5.0, "SIMULATION_PSD_ratio5");
    // export_sim(x, 10.0, "SIMULATION_PSD_ratio10");
    // export_sim(x, 15.0, "SIMULATION_PSD_ratio15");
    // export_sim(x, 20.0, "SIMULATION_PSD_ratio20");

}

