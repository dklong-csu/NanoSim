#ifndef MEPBM_PARTICLESYSTEM_H
#define MEPBM_PARTICLESYSTEM_H

#include<vector> // needed to store the list of chemical reactions
#include<functional> // needed to pass a function as an argument
#include<unordered_map> // needed to automatically index chemical species
#include<utility> // needed to have coefficient-chemical pair for reactions
#include<string> // needed to have text descriptions of chemical species
#include<iostream> // needed to output
#include<cmath> // needed to raise a number to a power

// includes for SUNDIALS to solve ODEs
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_matrix.h>

namespace NanoSim{
  template<typename Real>
  class particleSystem {
    public:
    /// Constructor



    /// Adds an individual chemical reaction
    void addReaction(const std::vector< std::pair<int,std::string> > & reactants,
      const std::vector< std::pair<int,std::string> > & products,
      const Real reaction_rate);



    /// Define characteristics of a particle
    void defineParticle(const int fewest_atoms, const int most_atoms, const int largest_agglomeration_size);



    /// Adds in particle nucleation
    void addNucleation(const std::vector< std::pair<int, std::string> > & reactants,
      const std::vector< std::pair<int, std::string> > & products,
      const Real reaction_rate,
      const int n_particles_created=1);



    /// Adds in particle growth
    void addGrowth(const std::vector< std::pair<int, std::string> > & reactants,
      const std::vector< std::pair<int, std::string> > & products,
      const std::function<Real(const unsigned int)> & growth_rate_fcn,
      const int particle_size_increase=1);



    /// Adds in particle agglomeration
    void addAgglomeration(const std::vector< std::pair<int, std::string> > & reactants,
      const std::vector< std::pair<int, std::string> > & products,
      const std::function<Real(const unsigned int, const unsigned int)> & agglomeration_rate_fcn);



    /// Finalizes chemical system so it is prepared for simulations
    void finalizeReactions();
    
    
    
    /// Prints the chemical reaction system
    void printChemicalReactions();



    /// Print ODE right hand side equations
    void printODErhs();



    /// Print ODE Jacobian equations
    void printODEJacobian();



    /// Computes the rhs vector for the ODEs
    std::function<int(Real, N_Vector, N_Vector, void*)>
    composeRHSfunction() const;

    /// Computes the Jacobian for the ODEs
    std::function<int(Real, N_Vector, N_Vector, SUNMatrix, void *, N_Vector, N_Vector, N_Vector)>
    composeJacobianfunction() const;
    
    /// Computes the PSD by solving the ODEs

    /// Extracts particles from the PSD vector based on internal indexing

    /// Provides particles sizes corresponding to the same vector indices in extract particles

    /// Options for the ODE solver in non-default settings are desired



    private:
    /*==============================
      Variables
      =============================*/
    /// Flag to prevent reactions being added after finalization has occurred 
    bool finalized = false;



    /// Flag indicating nucleation is present
    bool has_nucleation = false;



    /// Flag indicating growth is present
    bool has_growth = false;



    /// Flag indicating agglomeration is present
    bool has_agglomeration = false;



    /// Flag indicating whether a particle is defined
    bool has_particle = false;



    /// Stores a map of all chemicals species to their indices in vectors
    std::unordered_map<std::string, int> species_to_index_map;



    /// Tracks how many chemical species are in the reaction network
    int n_species = 0;



    /// Stores a map of all chemical reactions for simple conversion to ODEs
    std::vector< std::tuple< std::vector< std::pair<int, int> >, std::vector< std::pair<int, int> >, Real> > all_rxns;



    /// Stores information about nucleation in preparation for the final indexing operation
    std::tuple< std::vector< std::pair<int, std::string> >, std::vector< std::pair<int, std::string> >, Real, int> nucleation_info;



    /// Stores information about particle growth in preparation for the final indexing operation
    std::tuple< std::vector< std::pair<int, std::string> >, std::vector< std::pair<int, std::string> >, std::function<Real(const unsigned int)>, int> growth_info;



    /// Stores information about particle agglomeration in preparation for the final indexing operation
    std::tuple< std::vector< std::pair<int, std::string> >, std::vector< std::pair<int, std::string> >, std::function<Real(const unsigned int, const unsigned int)> > agglomeration_info;



    /// Particle size range and last size to track agglomeration for
    std::tuple< int, int, int> particle_info;



    /// Internal object to perform ODE integration



  };




  template<typename Real>
  void
  particleSystem<Real>::addReaction(const std::vector< std::pair<int,std::string> > & reactants,
    const std::vector< std::pair<int,std::string> > & products,
    const Real reaction_rate){
    // TODO error handling for particle keyword
    // TODO error handle if reactions are finalized


    // Vectors to put reaction into format for creating ODEs internally
    std::vector< std::pair<int, int> > R;
    std::vector< std::pair<int, int> > P;

    // Check hash table to see if a chemical has been included yet
    for (const auto & r : reactants){
      if (!species_to_index_map.contains(r.second)){
        species_to_index_map[r.second] = n_species;
        n_species += 1; // Vector indexes from 0 so increment after
      }
      R.push_back({r.first, species_to_index_map[r.second]});
    }

    for (const auto & p : products){
      if (!species_to_index_map.contains(p.second)){
        species_to_index_map[p.second] = n_species;
        n_species += 1; // Vector indexes from 0 so increment after
      }
      P.push_back({p.first, species_to_index_map[p.second]});
    }

    // 
    all_rxns.push_back({R, P, reaction_rate});
  }



  template<typename Real>
  void
  particleSystem<Real>::defineParticle(const int fewest_atoms, const int most_atoms, const int largest_agglomeration_size){
    // TODO do some error handling on sizes
    particle_info= {fewest_atoms, most_atoms, largest_agglomeration_size};
    has_particle = true;
  }



  template<typename Real>
  void
  particleSystem<Real>::addNucleation(const std::vector< std::pair<int, std::string> > & reactants,
    const std::vector< std::pair<int, std::string> > & products,
    const Real reaction_rate,
    const int n_particles_created){
    // TODO error handle if reactions are finalized
    nucleation_info = {reactants, products, reaction_rate, n_particles_created};
    has_nucleation = true;
  }



  template<typename Real>
  void
  particleSystem<Real>::addGrowth(const std::vector< std::pair<int, std::string> > & reactants,
    const std::vector< std::pair<int, std::string> > & products,
    const std::function<Real(const unsigned int)> & growth_rate_fcn,
    const int particle_size_increase){
      // TODO error handle if reactions are finalized
      growth_info = {reactants, products, growth_rate_fcn, particle_size_increase};
      has_growth = true;
  }



  template<typename Real>
  void
  particleSystem<Real>::addAgglomeration(const std::vector< std::pair<int, std::string> > & reactants,
    const std::vector< std::pair<int, std::string> > & products,
    const std::function<Real(const unsigned int, const unsigned int)> & agglomeration_rate_fcn){
      // TODO error handle if reactions are finalized
      agglomeration_info = {reactants, products, agglomeration_rate_fcn};
      has_agglomeration = true;
  }



  template<typename Real>
  void
  particleSystem<Real>::finalizeReactions(){
    // TODO error handle if already called


    // n_species now indicates the number of non-particle species and no more of these will be added
    // Hence we are free to create the vector indices for the particle species now
    if (has_particle){
      const unsigned int fewest_atoms = std::get<0>(particle_info);
      const unsigned int most_atoms = std::get<1>(particle_info);
      const unsigned int max_agglomeration_atoms = std::get<2>(particle_info);
      for (unsigned int n_atoms = fewest_atoms; n_atoms <= most_atoms; ++n_atoms){
        // Define a key name for the particle that's easy to construct to lookup in the hash table
        std::string particle_key = "__PARTICLE__" + std::to_string(n_atoms);
        
        // Add the new species to the hash table
        species_to_index_map[particle_key] = n_species;
        n_species += 1; // Vector indexes from 0 so increment after
      }

      // Construct individual reaction for nucleation
      // It is assumed that the provided smallest size is the size that is created
      if (has_nucleation)
      {
        std::vector< std::pair<int, std::string> > reactants = std::get<0>(nucleation_info);
        std::vector< std::pair<int, std::string> > products = std::get<1>(nucleation_info);

        std::string nucleated_particle = "__PARTICLE__" + std::to_string(fewest_atoms);

        const unsigned int n_particles_created = std::get<3>(nucleation_info);
        products.push_back({n_particles_created, nucleated_particle});

        const Real reaction_rate = std::get<2>(nucleation_info);

        addReaction(reactants, products, reaction_rate);
      }

      // Construct individual reactions for particle growth
      // Largest particle does not grow so that we conserve mass (this is an assumption that is fine if the largest atom is sufficiently large)
      if (has_growth){
        for (unsigned int n_atoms = fewest_atoms; n_atoms < most_atoms; ++n_atoms){
          const unsigned int created_atoms = n_atoms + std::get<3>(growth_info);
          if (created_atoms <= most_atoms){
            std::vector< std::pair<int, std::string> > reactants = std::get<0>(growth_info);
            std::vector< std::pair<int, std::string> > products  = std::get<1>(growth_info);

            const std::string consumed_particle = "__PARTICLE__" + std::to_string(n_atoms);
            const std::string created_particle  = "__PARTICLE__" + std::to_string(created_atoms);

            reactants.push_back({1, consumed_particle});
            products.push_back({1, created_particle});

            const auto reaction_rate_function = std::get<2>(growth_info);
            const Real reaction_rate = reaction_rate_function(n_atoms);

            addReaction(reactants, products, reaction_rate);
          }
        }
      }

      // Construct individual reactions for particle agglomeration
      if (has_agglomeration){
        for (unsigned int n_atoms1 = fewest_atoms; n_atoms1 <= max_agglomeration_atoms; ++n_atoms1){
          for (unsigned int n_atoms2 = n_atoms1; n_atoms2 <= max_agglomeration_atoms; ++n_atoms2){
            const unsigned int created_atoms = n_atoms1 + n_atoms2;
            if (created_atoms <= most_atoms){
              std::vector< std::pair<int, std::string> > reactants = std::get<0>(agglomeration_info);
              std::vector< std::pair<int, std::string> > products  = std::get<1>(agglomeration_info);

              const std::string consumed_particle1 = "__PARTICLE__" + std::to_string(n_atoms1);
              const std::string consumed_particle2 = "__PARTICLE__" + std::to_string(n_atoms2);
              const std::string created_particle  = "__PARTICLE__" + std::to_string(created_atoms);

              if (n_atoms1 == n_atoms2){
                reactants.push_back({2, consumed_particle1});
              } else {
                reactants.push_back({1, consumed_particle1});
                reactants.push_back({1, consumed_particle2});
              }
              products.push_back({1, created_particle});

              const auto reaction_rate_function = std::get<2>(agglomeration_info);
              const Real reaction_rate = reaction_rate_function(n_atoms1, n_atoms2);

              addReaction(reactants, products, reaction_rate);
            }
          }
        }
      }
    }
    finalized = true;
  }
    
    
    
  template<typename Real>
  void
  particleSystem<Real>::printChemicalReactions(){
    // TODO error handling for finalization
    for (const auto & rxn : all_rxns){
      const auto reactants = std::get<0>(rxn);
      const auto products = std::get<1>(rxn);
      const auto rate = std::get<2>(rxn);

      for (unsigned int i=0; i<reactants.size();++i){
        const auto r = reactants[i];
        std::cout << r.first << "X(" << r.second << ")";
        if (i<reactants.size()-1){
          std::cout << " + ";
        }
      }
      std::cout << " ->[" << rate << "] ";
      for (unsigned int i=0; i<products.size();++i){
        const auto p = products[i];
        std::cout << p.first << "X(" << p.second << ")";
        if (i<products.size()-1){
          std::cout << " + ";
        }
      }
      std::cout << "\n";
    }
  }



  template<typename Real>
  void
  particleSystem<Real>::printODErhs(){
    std::cout << "TODO\n";
  }



  template<typename Real>
  void
  particleSystem<Real>::printODEJacobian(){
    std::cout << "TODO\n";
  }



  template<typename Real>
  std::function<int(Real, N_Vector, N_Vector, void*)>
  particleSystem<Real>::composeRHSfunction() const
  {
    auto fcn = [&](Real time, N_Vector x, N_Vector x_dot, void* user_data)
    {
      // Reset the x_dot array to all zeros
      x_dot->ops->nvconst(0.0,x_dot);

      // Retrieve underlying data
      const Real* x_data = x->ops->nvgetarraypointer(x);
      Real* x_dot_data = x_dot->ops->nvgetarraypointer(x_dot);

      for (const auto & rxn : all_rxns){
        const std::vector< std::pair<int,int> > reactants = std::get<0>(rxn);
        const std::vector< std::pair<int,int> > products = std::get<1>(rxn);
        const Real reaction_rate = std::get<2>(rxn);

        // All effected species have derivative contributions proportional to
        // reaction_rate * [conc1]^coeff1 * [conc2]^coeff2 * ... 
        Real dx = reaction_rate;
        for (const auto & r : reactants){
          const unsigned int coeff = r.first;
          const unsigned int idx = r.second;

          dx *= std::pow(x_data[idx], coeff);
        }

        // Contributions to reactants (all negative!)
        for (const auto & r : reactants){
          const unsigned int coeff = r.first;
          const unsigned int idx = r.second;

          x_dot_data[idx] -= coeff * dx;
        }

        // Contributions to products (all positive!)
        for (const auto & p : products){
          const unsigned int coeff = p.first;
          const unsigned int idx = p.second;

          x_dot_data[idx] += coeff * dx;
        }
      }

      return 0;
    };
    return fcn;
  }



  template<typename Real>
  std::function<int(Real, N_Vector, N_Vector, SUNMatrix, void *, N_Vector, N_Vector, N_Vector)>
  particleSystem<Real>::composeJacobianfunction() const
  {
    auto fcn = [&](Real time, N_Vector x, N_Vector x_dot, SUNMatrix Jacobian, void * user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
    {
      // Zero out the Jacobian
      Jacobian->ops->zero(Jacobian);
 
      // Efficient insertion into the Jacobian is implementation specific so a function to do this this provided in the user_data
      const auto jacobianInsertAdd = *static_cast<std::function<int(SUNMatrix,sunindextype,sunindextype,Real)>*>(user_data);
 
      // Retrieve underlying data
      const Real* x_data = x->ops->nvgetarraypointer(x);

      for (const auto & rxn : all_rxns){
        const std::vector< std::pair<int,int> > reactants = std::get<0>(rxn);
        const std::vector< std::pair<int,int> > products = std::get<1>(rxn);
        const Real reaction_rate = std::get<2>(rxn);

        // Take derivative with respect to each reactant
        for (unsigned int r=0; r<reactants.size(); ++r){
          Real dx = reaction_rate;
          const unsigned int col = reactants[r].second;
 
          for (unsigned int r2=0; r2<reactants.size();++r2){
            const unsigned int coeff = reactants[r2].first;
            const unsigned int idx = reactants[r2].second;
    
            if (r2==r){
              // Derivative of x^n = n*x^(n-1)
              dx *= coeff * std::pow(x_data[idx], coeff-1); 
            } else {
              dx *= std::pow(x_data[idx], coeff);
            }
          }
          // Derivative wrt r computed up to +/- and coeff
          for (unsigned int r2=0; r2<reactants.size(); ++r2){
            // r = column, r2 = row
            const Real coeff = reactants[r2].first;
            const unsigned int idx = reactants[r2].second;
            // Scale by negative since reactant
            jacobianInsertAdd(Jacobian, idx, col, -coeff*dx);
          }

          for (unsigned int p=0; p<products.size(); ++p){
            // r = column, p = row
            const Real coeff = products[p].first;
            const unsigned int idx = products[p].second;
            // Scale by positive since reactant
            jacobianInsertAdd(Jacobian, idx, col, coeff*dx);
          }
        }
      }

      return 0;
    };
  
    return fcn;
  }
}
#endif