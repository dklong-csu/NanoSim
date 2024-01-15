#ifndef NANOSIM_PARTICLESYSTEM_H
#define NANOSIM_PARTICLESYSTEM_H

#include <vector> // needed to store the list of chemical reactions
#include <functional> // needed to pass a function as an argument
#include <unordered_map> // needed to automatically index chemical species
#include <utility> // needed to have coefficient-chemical pair for reactions
#include <string> // needed to have text descriptions of chemical species
#include <iostream> // needed to output
#include <cmath> // needed to raise a number to a power
#include <cassert> // needed to assert certain conditions for error checking
#include <stdexcept> // needed for error handling
#include <fstream> // needed to write sparsity pattern to file
#include <chrono>

// includes for SUNDIALS to solve ODEs
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_types.h>


#include "linearAlgebra.h"

#include <omp.h>

namespace NanoSim{

  /**
   * TODO documentatioin
  */
  struct speciesInfo {
    int vector_index = -1;

    bool is_particle = false;

    int smallest_size = 0;
    int n_binned_particles = 0;
  };



  /**
   * TODO documentation
   */ 
  template<typename Real>
  class particleSystem {
    public:

      // ------------------------------------------------------------------------------
      //  Adds a reaction to the list of all reactions
      //    Input:
      //      reactants     --> vector where each element is {stoichiometry coefficient, name of chemical species}
      //                        Note: the stoichiometry coefficient must be an integer for reactants for proper law of mass action differential equations
      //      products      --> vector where each element is {stoichiometry coefficient, name of chemical species}
      //                        Note: the stoichiometry coefficient is a real number here
      //      reaction_rate --> the rate at which this individual reaction takes place
      //      
      // ------------------------------------------------------------------------------
      void addReaction(const std::vector< std::pair<int,std::string> > & reactants,
        const std::vector< std::pair<Real,std::string> > & products,
        const Real reaction_rate);


      // ------------------------------------------------------------------------------
      //  Defines the necessary information for representing a particle in chemical reactions
      //    Input:
      //      fewest_atoms        --> the fewest number of atoms that a particle contains
      //                              Note: this is the size particle created through nucleation
      //      most_atoms          --> the largest number of atoms in a tracked particle
      //      atoms2diameter      --> a function that computes a particle diameter based on the number of atoms
      //      reduction_tolerance --> Particles are binned based on how close their diameters are together
      //                              This provides a percentage value that decides whether particles are binned
      //                              If reduction_tolerance = 0 then all individual particles are tracked
      //                              If, for example, reduction_tolerance = 0.01, then particles are binned together
      //                              such that a particle bin spans a 1% change in particle diameter.
      // ------------------------------------------------------------------------------
      void defineParticle(const int fewest_atoms, 
        const int most_atoms, 
        const std::function<Real(const int)> atoms2diameter,
        const Real reduction_tolerance = 0.0);


      // ------------------------------------------------------------------------------
      //  Adds a nucleation event following the chemical reaction
      //    Reactants --> Products + smallest particle
      //    Input:
      //      reactants     --> vector where each element is {stoichiometry coefficient, name of chemical species}
      //                        Note: the stoichiometry coefficient must be an integer for reactants for proper law of mass action differential equations
      //      products      --> vector where each element is {stoichiometry coefficient, name of chemical species}
      //                        Note: the stoichiometry coefficient is a real number here and the smallest particle will be added to this list
      //      reaction_rate --> the rate at which this individual reaction takes place
      // ------------------------------------------------------------------------------
      void addNucleation(const std::vector< std::pair<int, std::string> > & reactants,
        const std::vector< std::pair<Real, std::string> > & products,
        const Real reaction_rate);


      // ------------------------------------------------------------------------------
      //  Adds particle growth
      //    Input:
      //      reactants       --> vector where each element is {stoichiometry coefficient, name of chemical species}
      //                          Note: the stoichiometry coefficient must be an integer for reactants for proper law of mass action differential equations
      //                          Note: This list should not include any particles because those are added later
      //      products        --> vector where each element is {stoichiometry coefficient, name of chemical species}
      //                          Note: the stoichiometry coefficient is a real number here
      //                          Note: This list should not include any particles because those are added later
      //      growth_rate_fcn --> A function that takes in the number of atoms in a particle and returns the associated reaction rate for the growth reaction
      // ------------------------------------------------------------------------------
      void addGrowth(const std::vector< std::pair<int, std::string> > & reactants,
        const std::vector< std::pair<Real, std::string> > & products,
        const std::function<Real(const unsigned int)> & growth_rate_fcn);


      // ------------------------------------------------------------------------------
      //  Adds particle agglomeration
      //    Input:
      //      reactants               --> vector where each element is {stoichiometry coefficient, name of chemical species}
      //                                  Note: the stoichiometry coefficient must be an integer for reactants for proper law of mass action differential equations
      //                                  Note: This list should not include any particles because those are added later
      //      products                --> vector where each element is {stoichiometry coefficient, name of chemical species}
      //                                  Note: the stoichiometry coefficient is a real number here
      //                                  Note: This list should not include any particles because those are added later
      //      agglomeration_rate_fcn  --> A function that takes in the number of atoms in each of the two particles agglomerating and returns the associated reaction rate for that agglomeration reaction
      // ------------------------------------------------------------------------------
      void addAgglomeration(const std::vector< std::pair<int, std::string> > & reactants,
        const std::vector< std::pair<Real, std::string> > & products,
        const std::function<Real(const unsigned int, const unsigned int)> & agglomeration_rate_fcn);


      // ------------------------------------------------------------------------------
      //  Adds growth with exact Method of Moments which can simulate arbitrarily large particles with little computational effort
      //    Input:
      //      growth_precursor
      //      products
      //      rate_inflow
      //      rate_eMoM
      //      conversion_factor
      //  FIXME -- FINISH THIS DOCUMENTATION
      // ------------------------------------------------------------------------------
      void addeMoMGrowth(const std::string & growth_precursor,
        const std::vector< std::pair<Real, std::string> > & products,
        const Real rate_inflow,
        const Real rate_eMoM,
        const Real conversion_factor = 0.3 /* FIXME */);


      // ------------------------------------------------------------------------------
      //  Finalizes all of the reactions that occur in this particle system
      //    - Constructs a list of every chemical species and particle
      //    - Bins particles together to reduce the number of equations
      //    - Adds the nucleation reaction (if present)
      //    - Adds every reaction corresponding to particle growth (if present)
      //    - Adds every reaction corresponding to particle agglomeration (if present)
      //    - Adds equations to track the Moments needs for exact Method of Moments (if present)
      // ------------------------------------------------------------------------------
      void finalizeReactions();


      // ------------------------------------------------------------------------------
      //  Debugging tool used to print a list of every chemical species in the system
      // ------------------------------------------------------------------------------
      void printChemicalSpecies();


      // ------------------------------------------------------------------------------
      //  Debugging tool used to print a list of every reaction included in the system
      // ------------------------------------------------------------------------------
      void printChemicalReactions();


      // ------------------------------------------------------------------------------
      //  Creates a function corresponding to the right hand side of the ODEs that 
      //  describe this chemical reaction system.
      //  In a format that is compatible with the SUNDIALS software that solves the ODEs
      // ------------------------------------------------------------------------------
      std::function<int(Real, N_Vector, N_Vector, void*)>
      composeRHSfunction() const;


      // ------------------------------------------------------------------------------
      //  Creates a function corresponding to the Jacobian matrix of the ODEs that
      //  describe this chemical reaction system.
      //  In a format that is compatible with the SUNDIALS software that solves the ODEs
      // ------------------------------------------------------------------------------
      std::function<int(Real, N_Vector, N_Vector, SUNMatrix, void *, N_Vector, N_Vector, N_Vector)>
      composeJacobianfunction() const;


      // ------------------------------------------------------------------------------
      //  Returns the number of chemical species in this chemical reaction system.
      //  A set of binned particles counts as 1 species.
      //  In other words, this provides the dimension of the vector that tracks concentrations in the ODEs
      // ------------------------------------------------------------------------------
      unsigned int getNumberOfSpecies() const;


      // ------------------------------------------------------------------------------
      //  Returns the number of particle bins in this chemical reaction system.
      // ------------------------------------------------------------------------------
      unsigned int getNumberOfParticleBins() const;


      // ------------------------------------------------------------------------------
      //  Returns the number of reactions present in the system.
      // ------------------------------------------------------------------------------
      unsigned int getNumberOfReactions() const;


      // ------------------------------------------------------------------------------
      //  Debugging tool used to print the sparsity pattern of the Jacobian matrix
      // ------------------------------------------------------------------------------
      void printJacobianSparsityPattern(std::string out_filename, abstractLinearAlgebraOperations<Real> *lin_alg);

      // ------------------------------------------------------------------------------
      // FIXME -- probably rename this so it reflects all the information contained in the hash
      // ------------------------------------------------------------------------------
      std::unordered_map< std::string, speciesInfo > species_to_index_map;

    private:

      // ------------------------------------------------------------------------------
      //  Whether the set of chemical reactions has been finalized or not
      //  Once finalized, one cannot add more reactions to this system
      // ------------------------------------------------------------------------------
      bool finalized = false;


      // ------------------------------------------------------------------------------
      //  Whether particle nucleation is included in this reaction system
      // ------------------------------------------------------------------------------
      bool has_nucleation = false;


      // ------------------------------------------------------------------------------
      //  Whether particle growth is included in this reaction system
      // ------------------------------------------------------------------------------
      bool has_growth = false;


      // ------------------------------------------------------------------------------
      //  Whether particle agglomeration is included in this reaction system
      // ------------------------------------------------------------------------------
      bool has_agglomeration = false;


      // ------------------------------------------------------------------------------
      //  Whether this reaction system contains particles
      //    If not then one can still solve the dynamics of a set of chemical reactions
      // ------------------------------------------------------------------------------
      bool has_particle = false;


      // ------------------------------------------------------------------------------
      //  Whether this reaction system contains growth modeled by exact Method of Moments
      // ------------------------------------------------------------------------------
      bool has_eMoM = false;


      // ------------------------------------------------------------------------------
      //  The tolerance for binning particles together
      //    If zero, then all particles are tracked individually
      //    If nonzero, then Particles are bins containing 1 or more particles
      //      Each bin has a smallest size N_S and contain particles with N atoms satisfying
      //        diameter(N_S) <= diameter(N) <= 100%*(diameter_tolerance)*diameter(N_S)
      //    The bin boundaries N_S are determined by starting with the smallest particle size
      //    and iteratively constructing the bins
      // ------------------------------------------------------------------------------
      Real diameter_tolerance = 0.0;


      // ------------------------------------------------------------------------------
      //  The number of chemical species in this chemical reaction system
      //  Binned particles are considered a single species even if a bin contains
      //  particles of multiple sizes
      // ------------------------------------------------------------------------------
      int n_species = 0;


      // ------------------------------------------------------------------------------
      //  The number of binned particles in this chemical reaction system
      // ------------------------------------------------------------------------------
      int n_particles = 0;


      // ------------------------------------------------------------------------------
      //  A list of all reactions
      //  Each element of the vector is a tuple containing
      //    1) List of all reactants in the form {stoichiometry coefficient, vector index of chemical species}
      //    2) List of all products in the form {stoichiometry coefficient, vector index of chemical species}
      //    3) Reaction rate of this particular reaction
      // ------------------------------------------------------------------------------
      std::vector< 
        std::tuple< 
          std::vector< std::pair<int, int > >, 
          std::vector< std::pair<Real, int > >, 
          Real
        > 
      > all_rxns;


      // ------------------------------------------------------------------------------
      //  Necessary information to describe nucleation
      //    1) List of all reactants in the form {stoichiometry coefficient, vector index of chemical species}
      //    2) List of all products in the form {stoichiometry coefficient, vector index of chemical species}
      //    3) Reaction rate of nucleation
      // ------------------------------------------------------------------------------
      std::tuple< 
        std::vector< std::pair<int, std::string> >, 
        std::vector< std::pair<Real, std::string> >, 
        Real
      > nucleation_info;


      // ------------------------------------------------------------------------------
      //  Necessary information to describe particle growth
      //    1) List of all reactants in the form {stoichiometry coefficient, vector index of chemical species}
      //    2) List of all products in the form {stoichiometry coefficient, vector index of chemical species}
      //    3) Function that provides the growth reaction rate based on the number of atoms in the particle
      // ------------------------------------------------------------------------------
      std::tuple< 
        std::vector< std::pair<int, std::string> >, 
        std::vector< std::pair<Real, std::string> >, 
        std::function<Real(const unsigned int)>
      > growth_info;


      // ------------------------------------------------------------------------------
      //  Necessary information to describe particle agglomeration
      //    1) List of all reactants in the form {stoichiometry coefficient, vector index of chemical species}
      //    2) List of all products in the form {stoichiometry coefficient, vector index of chemical species}
      //    3) Function that provides the agglomeration reaction rate based on the number of atoms in each particle
      // ------------------------------------------------------------------------------
      std::tuple< 
        std::vector< std::pair<int, std::string> >, 
        std::vector< std::pair<Real, std::string> >, 
        std::function<Real(const unsigned int, const unsigned int)> 
      > agglomeration_info;


      // ------------------------------------------------------------------------------
      //   FIXME -- documentation
      // ------------------------------------------------------------------------------
      std::tuple< 
        std::string, 
        std::vector< std::pair<Real, std::string> >, 
        Real, 
        Real, 
        Real
      > eMoM_info;


      // ------------------------------------------------------------------------------
      //  FIXME -- documentation
      // ------------------------------------------------------------------------------
      std::tuple< int, int> particle_info;


      // ------------------------------------------------------------------------------
      //  FIXME -- documentation
      // ------------------------------------------------------------------------------
      std::function<Real(const int)> fcn_atoms2diameter;


      // ------------------------------------------------------------------------------
      //  Checks if a chemical species has been added to the list yet and adds it if not
      //    Input:
      //      key --> Name of the chemical species
      //    If key is already in the map, nothing is done
      //    Otherwise, key is added to the map and the vector index is automatically calculated
      //    In addition, the number of chemical species in the system is also increased accordingly
      // ------------------------------------------------------------------------------
      void addSpeciesToMap(const std::string & key);
  };



  template<typename Real>
  void
  particleSystem<Real>::addSpeciesToMap(const std::string & key)
  {
    const bool in_map = species_to_index_map.contains(key);
    if (!in_map){
      speciesInfo species_info;
      species_info.vector_index = n_species;
      species_to_index_map[key] = species_info;
      n_species += 1;
    }
  }




  template<typename Real>
  void
  particleSystem<Real>::addReaction(
    const std::vector< std::pair<int,std::string> > & reactants,
    const std::vector< std::pair<Real,std::string> > & products,
    const Real reaction_rate
  )
  {

    if (finalized){
      throw std::logic_error(std::string("ERROR: Attempt to add a reaction after the chemical system has been finalized.\nThe method 'addReaction()' MUST be called prior to 'finalizeReactions()'.\n"));
    }

    if (reaction_rate > 0){
      // Vectors to put reaction into format for creating ODEs internally
      std::vector< std::pair<int, int> > R;
      std::vector< std::pair<Real, int>  > P;

      // Check hash table to see if a chemical has been included yet
      for (const auto & r : reactants){
        const auto key = r.second;
        addSpeciesToMap(key);
        // if (!species_to_index_map.contains(r.second)){
        //   // If the species is unknown at this point, assume it is not a particle
        //   speciesInfo my_chemical;
        //   my_chemical.vector_index = n_species;
        //   species_to_index_map[r.second] = my_chemical;
        //   n_species += 1; // Vector indexes from 0 so increment after
        // }
        R.push_back({r.first, species_to_index_map[key].vector_index});
      }

      for (const auto & p : products){
        const auto key = p.second;
        addSpeciesToMap(key);
        // if (!species_to_index_map.contains(p.second)){
        //   // If the species is unknown at this point, assume it is not a particle
        //   speciesInfo my_chemical;
        //   my_chemical.vector_index = n_species;
        //   species_to_index_map[p.second] = my_chemical;
        //   n_species += 1; // Vector indexes from 0 so increment after
        // }
        P.push_back({p.first, species_to_index_map[key].vector_index});
      }

      // 
      all_rxns.push_back({R, P, reaction_rate});
    }
  }



  template<typename Real>
  void
  particleSystem<Real>::defineParticle(const int fewest_atoms, 
    const int most_atoms, 
    const std::function<Real(const int)> atoms2diameter,
    const Real reduction_tolerance){

    assert(fewest_atoms <= most_atoms);
    particle_info= {fewest_atoms, most_atoms};
    has_particle = true;
    fcn_atoms2diameter = atoms2diameter;
    diameter_tolerance = reduction_tolerance;
  }



  template<typename Real>
  void
  particleSystem<Real>::addNucleation(const std::vector< std::pair<int, std::string> > & reactants,
    const std::vector< std::pair<Real, std::string> > & products,
    const Real reaction_rate){
     if (finalized){
      throw std::logic_error(std::string("ERROR: Attempt to add particle nucleation after the chemical system has been finalized.\nThe method 'addNucleation()' MUST be called prior to 'finalizeReactions()'.\n"));
    }
    nucleation_info = {reactants, products, reaction_rate};
    has_nucleation = true;
  }



  template<typename Real>
  void
  particleSystem<Real>::addGrowth(const std::vector< std::pair<int, std::string> > & reactants,
    const std::vector< std::pair<Real, std::string> > & products,
    const std::function<Real(const unsigned int)> & growth_rate_fcn){
      if (finalized){
        throw std::logic_error(std::string("ERROR: Attempt to add particle growth after the chemical system has been finalized.\nThe method 'addGrowth()' MUST be called prior to 'finalizeReactions()'.\n"));
      }
      growth_info = {reactants, products, growth_rate_fcn};
      has_growth = true;
  }



  template<typename Real>
  void
  particleSystem<Real>::addAgglomeration(const std::vector< std::pair<int, std::string> > & reactants,
    const std::vector< std::pair<Real, std::string> > & products,
    const std::function<Real(const unsigned int, const unsigned int)> & agglomeration_rate_fcn){
      if (finalized){
        throw std::logic_error(std::string("ERROR: Attempt to add particle agglomeration after the chemical system has been finalized.\nThe method 'addAgglomeration()' MUST be called prior to 'finalizeReactions()'.\n"));
      }
      agglomeration_info = {reactants, products, agglomeration_rate_fcn};
      has_agglomeration = true;
  }



  template<typename Real>
  void
  particleSystem<Real>::addeMoMGrowth(const std::string & growth_precursor,
    const std::vector< std::pair<Real, std::string> > & products,
    const Real rate_inflow,
    const Real rate_eMoM,
    const Real conversion_factor){
    if (finalized){
        throw std::logic_error(std::string("ERROR: Attempt to add particle growth after the chemical system has been finalized.\nThe method 'addeMoMGrowth()' MUST be called prior to 'finalizeReactions()'.\n"));
      }
    eMoM_info = {growth_precursor, products, rate_inflow, rate_eMoM, conversion_factor};
    has_eMoM = true;
  }


  template<typename Real>
  void
  particleSystem<Real>::finalizeReactions(){
    if (finalized){
      throw std::logic_error(std::string("ERROR: Attempt to finalize the chemical system a second time.\nThe method 'finalizeReactions()' can only be called one time. Otherwise chemical reactions will be duplicated and the resulting simulation will be incorrect.\n"));
    }

    // Need to see if nucleation has any new chemical species to keep vector ordering correct
    if (has_nucleation){
      std::vector< std::pair<int, std::string> > reactants = std::get<0>(nucleation_info);
      std::vector< std::pair<Real, std::string> > products = std::get<1>(nucleation_info);

      for (const auto & r : reactants){
        const auto key = r.second;
        addSpeciesToMap(key);
        // if (!species_to_index_map.contains(r.second)){
        //   // A normal chemical species does not have a particle diameter, so set it to zero
        //   species_to_index_map[r.second] = {n_species, 0.0};
        //   n_species += 1; // Vector indexes from 0 so increment after
        // }
      }

      for (const auto & p : products){
        const auto key = p.second;
        addSpeciesToMap(key);
        // if (!species_to_index_map.contains(p.second)){
        //   // A normal chemical species does not have a particle diameter, so set it to zero
        //   species_to_index_map[p.second] = {n_species, 0.0};
        //   n_species += 1; // Vector indexes from 0 so increment after
        // }
      }
    }

    // Need to see if growth has any new chemical species to keep vector ordering correct
    if (has_growth){
      std::vector< std::pair<int, std::string> > reactants = std::get<0>(growth_info);
      std::vector< std::pair<Real, std::string> > products = std::get<1>(growth_info);

      for (const auto & r : reactants){
        const auto key = r.second;
        addSpeciesToMap(key);
      }

      for (const auto & p : products){
        const auto key = p.second;
        addSpeciesToMap(key);
      }
    }

    // Need to see if agglomeration has any new chemical species to keep vector ordering correct
    if (has_agglomeration){
      std::vector< std::pair<int, std::string> > reactants = std::get<0>(agglomeration_info);
      std::vector< std::pair<Real, std::string> > products = std::get<1>(agglomeration_info);

      for (const auto & r : reactants){
        const auto key = r.second;
        addSpeciesToMap(key);
      }

      for (const auto & p : products){
        const auto key = p.second;
        addSpeciesToMap(key);
      }
    }

    // TODO
    // Need to see if eMoM growth has any new chemical species to keep vector ordering correct
    if (has_eMoM){
      auto precursor = std::get<0>(eMoM_info);
      auto products = std::get<1>(eMoM_info);
      addSpeciesToMap(precursor);

      for (const auto & p : products){
        const auto key = p.second;
        addSpeciesToMap(key);
      }

      // Add in "species" for the moments that eMoM models
      // TODO throw error if moments already defined
      addSpeciesToMap("__MOMENT3__");
      addSpeciesToMap("__MOMENT2__");
      addSpeciesToMap("__MOMENT1__");
      addSpeciesToMap("__MOMENT0__");
    }


    // n_species now indicates the number of non-particle species and no more of these will be added
    // Hence we are free to create the vector indices for the particle species now
    if (has_particle){
      const unsigned int fewest_atoms = std::get<0>(particle_info);
      const unsigned int most_atoms = std::get<1>(particle_info);

      // Loop through all particle sizes and bin them
      
      // Always track the first particle as an individual
      n_particles = 0; // No particles at this point
      std::string key = "__PARTICLE__" + std::to_string(n_particles);
      addSpeciesToMap(key);
      // ---------------------
      //  Map defaults speciesInfo to non-particle
      //  Update this information now
      //    Need to modify:
      //      is_particle = true
      //      smallest_size = fewest_atoms
      //      n_binned_particles = 1
      // ---------------------
      species_to_index_map[key].is_particle = true;
      species_to_index_map[key].smallest_size = fewest_atoms;
      species_to_index_map[key].n_binned_particles = 1;
      // One particle is now included to update this counter
      n_particles += 1;

      // information for binning particles
      int bin_size1 = fewest_atoms+1;
      int bin_num_particles = 1;
      auto bin_diam1 = fcn_atoms2diameter(bin_size1);
      // Loop through all other particles except for the largest size
      for (unsigned int n_atoms = fewest_atoms+2; n_atoms < most_atoms; ++n_atoms){
        const auto diam = fcn_atoms2diameter(n_atoms);
        const auto diff = std::abs(diam - bin_diam1);
        const bool in_bin = diff < diameter_tolerance * bin_diam1;
        if (in_bin){
          // Add this particle to the bin
          bin_num_particles += 1;
        } else {
          // finalize bin
          key = "__PARTICLE__" + std::to_string(n_particles);
          addSpeciesToMap(key);
          species_to_index_map[key].is_particle = true;
          species_to_index_map[key].smallest_size = bin_size1;
          species_to_index_map[key].n_binned_particles = bin_num_particles;
          n_particles += 1;

          // Set data for next particle bin
          bin_size1 = n_atoms;
          bin_num_particles = 1;
          bin_diam1 = fcn_atoms2diameter(bin_size1);
        }
      }
      // Last bin is left unfinalized
      key = "__PARTICLE__" + std::to_string(n_particles);
      addSpeciesToMap(key);
      species_to_index_map[key].is_particle = true;
      species_to_index_map[key].smallest_size = bin_size1;
      species_to_index_map[key].n_binned_particles = bin_num_particles;
      n_particles += 1;

      // Very last particle size is its own bin
      key = "__PARTICLE__" + std::to_string(n_particles);
      addSpeciesToMap(key);
      species_to_index_map[key].is_particle = true;
      species_to_index_map[key].smallest_size = most_atoms;
      species_to_index_map[key].n_binned_particles = 1;
      n_particles += 1;



      // Construct individual reaction for nucleation
      // It is assumed that the provided smallest size is the size that is created
      if (has_nucleation)
      {
        std::vector< std::pair<int, std::string> > reactants = std::get<0>(nucleation_info);
        std::vector< std::pair<Real, std::string> > products = std::get<1>(nucleation_info);

        std::string nucleated_particle = "__PARTICLE__0";

        products.push_back({1, nucleated_particle});

        const Real reaction_rate = std::get<2>(nucleation_info);

        addReaction(reactants, products, reaction_rate);
      }

      // Construct individual reactions for particle growth
      // Largest particle does not grow so that we conserve mass (this is an assumption that is fine if the largest atom is sufficiently large)
      if (has_growth){
        /*
            Loop through all binned particles
            If bin I contains N particles, then the growth of the first N-1 particles
            creates a particle in bin I. The Nth particle grows into a particle in bin J.
            The reaction of the binned particle is shorthand for the N reactions occurring within the bin
            where we assume the kinetics behave all in the same way. 
        */
        // std::cout << "Adding growth reactions...";
        for (unsigned int bin=0; bin<n_particles-1; ++bin){
          const std::string bin_key = "__PARTICLE__" + std::to_string(bin);
          const std::string next_bin_key = "__PARTICLE__" + std::to_string(bin+1);

          std::vector< std::pair<int, std::string> > reactants = std::get<0>(growth_info);
          std::vector< std::pair<Real, std::string> > products  = std::get<1>(growth_info);

          // This bin is part of reactants
          reactants.push_back( {1, bin_key} );

          // FIXME
          const Real N = species_to_index_map[bin_key].n_binned_particles;
          const Real weight1 = (N-1.0)/N; // percentage of reactions that remain in same bin
          const Real weight2 = 1.0/N; // percentage of reactions that go to next bin

          if (weight1 > 0){
            products.push_back({weight1, bin_key}); 
          } 

          products.push_back( {weight2, next_bin_key} ); 

          // Binned particle acts like average size
          const auto min_size = species_to_index_map[bin_key].smallest_size;
          const auto n_binned = species_to_index_map[bin_key].n_binned_particles;

          const int avg_size = min_size + (n_binned-1)/2;

          const auto reaction_rate_fcn = std::get<2>(growth_info);
          const auto reaction_rate = reaction_rate_fcn(avg_size);
          addReaction(reactants, products, reaction_rate);
        }
      }

      // Construct individual reactions for particle agglomeration
      if (has_agglomeration){
        /*
            Loop through all combinations of bin I agglomerating with bin J
            Compute the proportion that goes into bin K
            If the agglomerated size is bigger than the largest tracked particle
            then default to putting it in the last bin
            Ignore agglomeration with the last bin and just treat it as a "catch all" for larger particles
        */
      //  std::cout << "Adding agglomeration reactions..." << std::flush;
        auto start = std::chrono::high_resolution_clock::now();
        for (unsigned int bin1=0; bin1 < n_particles-1; ++bin1){
          for (unsigned int bin2=bin1; bin2 < n_particles-1; ++bin2){
            std::vector< std::pair<int, std::string> > reactants = std::get<0>(agglomeration_info);
            std::vector< std::pair<Real, std::string> > products  = std::get<1>(agglomeration_info);

            const std::string key_bin1 = "__PARTICLE__" + std::to_string(bin1);
            const std::string key_bin2 = "__PARTICLE__" + std::to_string(bin2);
            // Add particles to reactant list
            if (bin1 == bin2){
              reactants.push_back( {2, key_bin1} );
            } else {
              reactants.push_back( {1, key_bin1} );
              reactants.push_back( {1, key_bin2} );
            }

            int first_bin_check = bin2; // Don't need to check if created particle is in smaller bins
            std::vector<int> weights(n_particles);

            const auto smallest_bin1 = species_to_index_map[key_bin1].smallest_size;
            const auto N_bin1 = species_to_index_map[key_bin1].n_binned_particles;
            const auto smallest_bin2 = species_to_index_map[key_bin2].smallest_size;
            const auto N_bin2 = species_to_index_map[key_bin2].n_binned_particles;
            const Real n_rxns = N_bin1*N_bin2;

            // Bin I acts like its average size -- fixme
            const int avg_size_bin1 = smallest_bin1 + (N_bin1 - 1)/2;
            // Bin J acts like its average size
            const int avg_size_bin2 = smallest_bin2 + (N_bin2 - 1)/2;

            const auto reaction_rate_function = std::get<2>(agglomeration_info);
            const Real reaction_rate = reaction_rate_function(avg_size_bin1, avg_size_bin2);

            if (reaction_rate > 0){
              // -------------------------------------------------------------------------
              //  Binning algorithm
              //    - Start with the bigger bin between Bin I and Bin J (because can't create a particle smaller than the larger of the two bins)
              //    - Bin K = max(Bin I, Bin J)
              //    - Created particles = [min(Bin I) + min(Bin J), max(Bin I) + max(Bin J)] = [lo, hi]
              //    - N = size(Bin I), M = size(Bin J)
              //    - Slope up = [lo, lo + min(N,M) - 1] = [sup1, sup2]
              //    - Slope down = [hi - min(N,M) + 1, hi] = [sdown1, sdown2]
              //    - 
              //    - Total created = |Bin I| * |Bin J|
              //    - check bin = Bin K
              //    - working bin = [lo, hi]
              //    - while checking bins
              //      - in check bin = [lo, min( hi, max(check bin) )] = [a, b]
              //      - in slope up = []
              // -------------------------------------------------------------------------
              for (unsigned int size1=smallest_bin1; size1<smallest_bin1+N_bin1;++size1){
                int size2 = smallest_bin2;
                const int max_size2 = smallest_bin2+N_bin2-1;
                while (size2 <= max_size2){
                  // Check which bin size1+size2 is in
                  const int created_size = size1+size2;
                  int search_bin = first_bin_check;
                  bool bin_found = false;
                  while (!bin_found){
                    if (search_bin == n_particles-1){
                      // If searching reaches the last bin, the just put it in the last bin
                      bin_found = true;
                    } else {
                      const std::string search_key = "__PARTICLE__" + std::to_string(search_bin);
                      const int search_bin_smallest = species_to_index_map[search_key].smallest_size;
                      const int search_bin_largest = species_to_index_map[search_key].n_binned_particles - 1 + search_bin_smallest;
                      const bool in_bin = (created_size >= search_bin_smallest) && (created_size <= search_bin_largest);
                      if (in_bin){
                        bin_found = true;
                      } else {
                        search_bin += 1;
                      }
                    }
                  } // end search for bin
                  // Now check how many more particles can be placed in this same bin
                  // If the found bin is the last bin, then all remaining size2 particles will also go here
                  if (search_bin == n_particles-1){
                    const int particles_left = max_size2 - size2 + 1;
                    weights[search_bin] += particles_left;
                    size2 = max_size2+1; // Can end the while loop over size2 now
                  } else {
                    int num_in_bin = 1;
                    // Determine how many more consecutive particles can go into search_bin
                    const std::string search_key = "__PARTICLE__" + std::to_string(search_bin);
                    const int search_bin_smallest = species_to_index_map[search_key].smallest_size;
                    const int search_bin_largest = species_to_index_map[search_key].n_binned_particles - 1 + search_bin_smallest;
                    const int remaining_in_bin = search_bin_largest - created_size;
                    const int remaining_size2 = max_size2 - size2;
                    const int additional_in_bin = std::min(remaining_in_bin, remaining_size2);
                    num_in_bin += additional_in_bin;

                    // Add created particles to the bin
                    weights[search_bin] += num_in_bin;

                    // Increment size2
                    size2 += num_in_bin;
                  }
                } // Finished loop over size2
              } // Finished loop over size1

              // Add bins to product list
              for (int product_bin=first_bin_check; product_bin<n_particles; ++product_bin){
                if (weights[product_bin]>0){
                  const std::string product_key = "__PARTICLE__" + std::to_string(product_bin);
                  // FIXME -- add weight to product
                  const Real wij = 1.0*weights[product_bin] / n_rxns;
                  products.push_back( {wij, product_key} );
                }
              } // Finished adding particles to the product list
              addReaction(reactants, products, reaction_rate);
            }
            
          }
        }
        auto end = std::chrono::high_resolution_clock::now();
        const double runtime = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
        // std::cout << "Agglomeration loop took " << runtime << " ms\n";
        // std::cout << "done!\n";
      }

      // Add the last growth reaction for eMoM cases!
      if (has_eMoM){
        std::vector< std::pair<int, std::string> > reactants = std::get<0>(growth_info);
        std::vector< std::pair<Real, std::string> > products  = std::get<1>(growth_info);

        const std::string consumed_particle = "__PARTICLE__" + std::to_string(n_particles-1);
        reactants.push_back({1,consumed_particle});

        const auto reaction_rate_function = std::get<2>(growth_info);
        const Real reaction_rate = reaction_rate_function(most_atoms);

        addReaction(reactants, products, reaction_rate);
      }
    }
    finalized = true;
  }
    
    
 

  template<typename Real>
  void
  particleSystem<Real>::printChemicalSpecies(){

    for (const auto & kv : species_to_index_map){
      std::cout << kv.first << "  index " << kv.second.vector_index << "\n";
    }
  }



  template<typename Real>
  void
  particleSystem<Real>::printChemicalReactions(){
    if (!finalized){
      std::cout << "WARNING: Chemical reactions not finalized. Final results may be different from the printed chemical reactions.\n";
    }
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
  std::function<int(Real, N_Vector, N_Vector, void*)>
  particleSystem<Real>::composeRHSfunction() const
  {
    auto fcn = [&](Real time, N_Vector x, N_Vector x_dot, void* user_data)
    {
      // user_data is a pointer to an object containing linear algebra operations
      abstractLinearAlgebraOperations<Real> *lin_algebra = static_cast< abstractLinearAlgebraOperations<Real>* >(user_data);

      // Reset the x_dot array to all zeros
      x_dot->ops->nvconst(0.0,x_dot);

      // std::cout << "All rxns loop\n";
      for (const auto & rxn : all_rxns){
        const std::vector< std::pair< int, int > > reactants = std::get<0>(rxn);
        const std::vector< std::pair< Real, int > > products = std::get<1>(rxn);
        const Real reaction_rate = std::get<2>(rxn);

        // All affected species have derivative contributions proportional to
        // reaction_rate * [conc1]^coeff1 * [conc2]^coeff2 * ... 
        Real dx = reaction_rate;
        for (const auto & r : reactants){
          const unsigned int coeff = r.first;
          const int idx = r.second;
          Real conc = lin_algebra->vectorGetValue(x,idx);
          dx *= std::pow(conc, coeff);
        }

        // Contributions to reactants (all negative!)
        for (const auto & r : reactants){
          const unsigned int coeff = r.first;
          const int idx = r.second;
          lin_algebra->vectorInsertAdd(x_dot, -1.0*coeff*dx, idx);
        }

        // Contributions to products (all positive!)
        for (const auto & p : products){
          const Real coeff = p.first;
          const int idx = p.second;
          lin_algebra->vectorInsertAdd(x_dot, 1.0*coeff*dx, idx);
        }
      }
      // std::cout << "All rxns loop done\n";

      // add in eMoM if present!
      // FIXME
      if (has_eMoM){
        /*
          Values used
        */
        const Real delx = std::get<4>(eMoM_info); 
        const unsigned int M = std::get<1>(particle_info); // largest tracked cluster which flows into eMoM
        const Real xM = delx * std::pow(1.0 * M, 1.0 / 3.0); // continuous diameter size of largest cluster
        const Real k_inflow = std::get<2>(eMoM_info); // Growth rate of cluster of size M growing into eMoM
        const Real k_eMoM   = std::get<3>(eMoM_info); // Growth rate of eMoM

        const auto M_idx = species_to_index_map.find("__PARTICLE__" + std::to_string(M))->second.vector_index;

        /*
          Concentrations
        */
        const auto A_idx = species_to_index_map.find(std::get<0>(eMoM_info))->second.vector_index;
        const auto A = lin_algebra->vectorGetValue(x, A_idx);

        const auto M3_idx = species_to_index_map.find("__MOMENT3__")->second.vector_index;
        const auto M2_idx = species_to_index_map.find("__MOMENT2__")->second.vector_index;
        const auto M1_idx = species_to_index_map.find("__MOMENT1__")->second.vector_index;
        const auto M0_idx = species_to_index_map.find("__MOMENT0__")->second.vector_index;

        const auto M3 = lin_algebra->vectorGetValue(x, M3_idx);
        const auto M2 = lin_algebra->vectorGetValue(x, M2_idx);
        const auto M1 = lin_algebra->vectorGetValue(x, M1_idx);
        const auto M0 = lin_algebra->vectorGetValue(x, M0_idx);

        /*
          Flux:
          q = 3 * k_inflow * M^(2/3) * [BM] / delx / k_eMoM
        */
        const Real q = 3 * k_inflow * std::pow( M * 1.0, 2.0 / 3.0) * lin_algebra->vectorGetValue(x, M_idx) / delx / k_eMoM;

        /*
          Precursor contribution
          dA/dt += -8*delx/9 * [A] * k_eMoM * xn^3 * q - k_eMoM * 8 / (delx^2 * 3) * [A] * [M2]
        */
        lin_algebra->vectorInsertAdd(x_dot,
          -8. * delx / 9. * A * k_eMoM * std::pow(xM, 3.0) * q - k_eMoM * 8. / std::pow(delx, 2.0) / 3. * A * M2,
          A_idx);

        /*
          Products contribution
          dP/dt += 8*delx/9 * [A] * k_eMoM * xn^3 * q + k_eMoM * 8 / (delx^2 * 3) * [A] * [M2]
        */
        const auto products = std::get<1>(eMoM_info);
        for (auto p : products){
          const auto coeff = p.first;
          const auto p_ID = p.second;
          const auto p_idx = species_to_index_map.find(p_ID)->second.vector_index;
          lin_algebra->vectorInsertAdd(x_dot, 
            coeff * (8. * delx / 9. * A * k_eMoM * std::pow(xM, 3.0) * q + k_eMoM * 8. / std::pow(delx, 2.0) / 3. * A * M2),
            p_idx); 
        }


       /*
        Moment 3
        dM3/dt = 8 * delx / 9 * [A] * k_eMoM * xM^3 * [q] + 3 * 8 * k_eMoM * delx / 9 * [A] * [M2]
       */
      lin_algebra->vectorInsert(x_dot,
        8. * delx / 9. * A * k_eMoM * std::pow(xM, 3.0) * q + 3. * 8. * k_eMoM * delx / 9. * A * M2,
        M3_idx);

      /*
        Moment 2
        dM2/dt = 8 * delx / 9 * [A] * k_eMoM * xM^2 * [q] + 2 * 8 * k_eMoM * delx / 9 * [A] * [M1]
       */
      lin_algebra->vectorInsert(x_dot,
        8. * delx / 9. * A * k_eMoM * std::pow(xM, 2.0) * q + 2. * 8. * k_eMoM * delx / 9. * A * M1,
        M2_idx);

      /*
        Moment 1
        dM1/dt = 8 * delx / 9 * [A] * k_eMoM * xM^1 * [q] + 1 * 8 * k_eMoM * delx / 9 * [A] * [M0]
       */
      lin_algebra->vectorInsert(x_dot,
        8. * delx / 9. * A * k_eMoM * std::pow(xM, 1.0) * q + 1. * 8. * k_eMoM * delx / 9. * A * M0,
        M1_idx);

      /*
        Moment 0
        dM0/dt = 8 * delx / 9 * [A] * k_eMoM * xM^0 * [q] 
       */
      lin_algebra->vectorInsert(x_dot,
        8. * delx / 9. * A * k_eMoM * q,
        M0_idx);
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
      // user_data is a pointer to an object containing linear algebra operations
      abstractLinearAlgebraOperations<Real> *lin_algebra = static_cast< abstractLinearAlgebraOperations<Real>* >(user_data);

      // Zero out the Jacobian
      Jacobian->ops->zero(Jacobian);

      for (const auto & rxn : all_rxns){
        const std::vector< std::pair<int, int > > reactants = std::get<0>(rxn);
        const std::vector< std::pair<Real, int > > products = std::get<1>(rxn);
        const Real reaction_rate = std::get<2>(rxn);

        // Take derivative with respect to each reactant
        for (unsigned int r=0; r<reactants.size(); ++r){
          Real dx = reaction_rate;
          const int col = reactants[r].second;
          for (unsigned int r2=0; r2<reactants.size();++r2){
            const unsigned int coeff = reactants[r2].first;
            const int row = reactants[r2].second;
            
            Real conc = lin_algebra->vectorGetValue(x,row);
            if (r2==r){
              dx *= coeff * std::pow(conc, coeff-1);
            } else {
              dx *= std::pow(conc, coeff);
            }
          }

          for (unsigned int r2=0; r2<reactants.size(); ++r2){
            const Real coeff = reactants[r2].first;
            const int row = reactants[r2].second;
            lin_algebra->matrixInsertAdd(Jacobian, -coeff*dx, row, col);
          }

          for (unsigned int p=0; p<products.size(); ++p){
            const Real coeff = products[p].first;
            const int row = products[p].second;
            lin_algebra->matrixInsertAdd(Jacobian, coeff*dx, row, col);
          }
        }
      }


      // add in eMoM if present!
      // FIXME (has_eMoM)
      if (has_eMoM){
        /*
          Values used
        */
        const Real delx = std::get<4>(eMoM_info); 
        const unsigned int M = std::get<1>(particle_info); // largest tracked cluster which flows into eMoM
        const Real xM = delx * std::pow(1.0 * M, 1.0 / 3.0); // continuous diameter size of largest cluster
        const Real k_inflow = std::get<2>(eMoM_info); // Growth rate of cluster of size M growing into eMoM
        const Real k_eMoM   = std::get<3>(eMoM_info); // Growth rate of eMoM

        const auto M_idx = species_to_index_map.find("__PARTICLE__" + std::to_string(M))->second.vector_index;

        /*
          Concentrations
        */
        const auto A_idx = species_to_index_map.find(std::get<0>(eMoM_info))->second.vector_index;
        const auto A = lin_algebra->vectorGetValue(x, A_idx);

        const auto M3_idx = species_to_index_map.find("__MOMENT3__")->second.vector_index;
        const auto M2_idx = species_to_index_map.find("__MOMENT2__")->second.vector_index;
        const auto M1_idx = species_to_index_map.find("__MOMENT1__")->second.vector_index;
        const auto M0_idx = species_to_index_map.find("__MOMENT0__")->second.vector_index;

        const auto M3 = lin_algebra->vectorGetValue(x, M3_idx);
        const auto M2 = lin_algebra->vectorGetValue(x, M2_idx);
        const auto M1 = lin_algebra->vectorGetValue(x, M1_idx);
        const auto M0 = lin_algebra->vectorGetValue(x, M0_idx);

        /*
          Flux:
          q = 3 * k_inflow * M^(2/3) * [BM] / delx / k_eMoM
          Derivative
          dq/dBM = 3 * k_inflow * M^(2/3) / delx / k_eMoM
        */
        const Real q = 3 * k_inflow * std::pow( M * 1.0, 2.0 / 3.0) * lin_algebra->vectorGetValue(x, M_idx) / delx / k_eMoM;
        const Real dqdB = 3 * k_inflow * std::pow( M * 1.0, 2.0 / 3.0)  / delx / k_eMoM;

        /*
          Precursor contribution
          dA/dt += -8*delx/9 * [A] * k_eMoM * xn^3 * q - k_eMoM * 8 / (delx^2 * 3) * [A] * [M2]
          Derivatives
          dA'/dA += -8*delx/9 * k_eMoM * xn^3 * q - k_eMoM * 8 / (delx^2 * 3) * [M2]
          dA'/dBM += -8*delx/9 * k_eMoM * xn^3 * [A] * dq/dBM
          dA'/dM2 += -k_eMoM / delx^2 * 8 / 3 * [A]
        */
        lin_algebra->matrixInsertAdd(Jacobian, 
          -8. * delx / 9. * k_eMoM * std::pow(xM, 3.0) * q - k_eMoM * 8. / std::pow(delx, 2.0) / 3. * M2,
          A_idx,
          A_idx); // dA'/dA

        lin_algebra->matrixInsertAdd(Jacobian,
          -8. * delx / 9. * k_eMoM * std::pow(xM, 3.0) * A * dqdB,
          A_idx,
          M_idx); // dA'/dBM

        lin_algebra->matrixInsertAdd(Jacobian,
          -k_eMoM / std::pow(delx, 2.0) * 8. / 3. * A,
          A_idx,
          M2_idx); // dA'/dM2

        /*
          Products contribution
          dP/dt += 8*delx/9 * [A] * k_eMoM * xn^3 * q + k_eMoM * 8 / (delx^2 * 3) * [A] * [M2]
          Derivatives
          dP'/dA += 8*delx/9 * k_eMoM * xn^3 * q + k_eMoM * 8 / (delx^2 * 3) * [M2]
          dP'/dBM += 8*delx/9 * k_eMoM * xn^3 * [A] * dq/dBM
          dP'/dM2 += k_eMoM / delx^2 * 8 / 3 * [A]
        */
        const auto products = std::get<1>(eMoM_info);
        for (auto p : products){
          const auto coeff = p.first;
          const auto p_ID = p.second;
          const auto p_idx = species_to_index_map.find(p_ID)->second.vector_index;
          lin_algebra->matrixInsertAdd(Jacobian, 
            coeff * (8. * delx / 9. * k_eMoM * std::pow(xM, 3.0) * q + k_eMoM * 8. / std::pow(delx, 2.0) / 3. * M2),
            p_idx,
            A_idx); // dP'/dA

          lin_algebra->matrixInsertAdd(Jacobian,
            coeff * 8. * delx / 9. * k_eMoM * std::pow(xM, 3.0) * A * dqdB,
            p_idx,
            M_idx); // dP'/dBM

          lin_algebra->matrixInsertAdd(Jacobian,
            coeff * k_eMoM / std::pow(delx, 2.0) * 8. / 3. * A,
            p_idx,
            M2_idx); // dP'/dM2
        }


       /*
        Moment 3
        dM3/dt = 8 * delx / 9 * [A] * k_eMoM * xM^3 * [q] + 3 * 8 * k_eMoM * delx / 9 * [A] * [M2]
        Derivatives
        dM3'/dA = 8 * delx / 9 * k_eMoM * xM^3 * [q] + 3 * 8 * k_eMoM * delx / 9 * [M2]
        dM3'/dBM = 8 * delx / 9 * [A] * k_eMoM * xM^3 * dqdB
        dM3'/dM2 = 3 * 8 * k_eMoM * delx / 9 * [A]
       */
      lin_algebra->matrixInsert(Jacobian,
        8. * delx / 9. * k_eMoM * std::pow(xM, 3.0) * q + 3. * 8. * k_eMoM * delx / 9. * M2,
        M3_idx,
        A_idx); // dM3'/dA

      lin_algebra->matrixInsert(Jacobian,
        8. * delx / 9. * A *  k_eMoM * std::pow(xM, 3.0) * dqdB,
        M3_idx,
        M_idx); // dM3'/dBM

      lin_algebra->matrixInsert(Jacobian,
        3. * 8. * k_eMoM * delx / 9. * A,
        M3_idx,
        M2_idx); // dM3'/dM2

      /*
        Moment 2
        dM2/dt = 8 * delx / 9 * [A] * k_eMoM * xM^2 * [q] + 2 * 8 * k_eMoM * delx / 9 * [A] * [M1]
        Derivatives
        dM2'/dA = 8 * delx / 9 * k_eMoM * xM^2 * [q] + 2 * 8 * k_eMoM * delx / 9 * [M1]
        dM2'/dBM = 8 * delx / 9 * [A] * k_eMoM * xM^2 * dqdB
        dM2'/dM2 = 2 * 8 * k_eMoM * delx / 9 * [A]
       */
      lin_algebra->matrixInsert(Jacobian,
        8. * delx / 9. * k_eMoM * std::pow(xM, 2.0) * q + 2. * 8. * k_eMoM * delx / 9. * M1,
        M2_idx,
        A_idx); // dM2'/dA

      lin_algebra->matrixInsert(Jacobian,
        8. * delx / 9. * A * k_eMoM * std::pow(xM, 2.0) * dqdB,
        M2_idx,
        M_idx); // dM2'/dBM

      lin_algebra->matrixInsert(Jacobian,
        2. * 8. * k_eMoM * delx / 9. * A,
        M2_idx,
        M1_idx); // dM2'/dM1

      /*
        Moment 1
        dM1/dt = 8 * delx / 9 * [A] * k_eMoM * xM^1 * [q] + 1 * 8 * k_eMoM * delx / 9 * [A] * [M0]
        Derivatives
        dM1'/dA = 8 * delx / 9 * k_eMoM * xM^1 * [q] + 1 * 8 * k_eMoM * delx / 9 * [M0]
        dM1'/dBM = 8 * delx / 9 * [A] * k_eMoM * xM^1 * dqdB
        dM1'/dM0 = 1 * 8 * k_eMoM * delx / 9 * [A]
       */
      lin_algebra->matrixInsert(Jacobian,
        8. * delx / 9. * k_eMoM * xM * q + 8. * k_eMoM * delx / 9. * M0,
        M1_idx,
        A_idx); // dM1'/dA

      lin_algebra->matrixInsert(Jacobian,
        8. * delx / 9. * A * k_eMoM * xM * dqdB,
        M1_idx,
        M_idx); // dM1'/dBM

      lin_algebra->matrixInsert(Jacobian,
        8. * k_eMoM * delx / 9. * A,
        M1_idx,
        M0_idx); // dM1'/dM0

      /*
        Moment 0
        dM0/dt = 8 * delx / 9 * [A] * k_eMoM * xM^0 * [q] 
        Derivatives
        dM0'/dA = 8 * delx / 9 * k_eMoM * xM^0 * [q] 
        dM0'/dBM = 8 * delx / 9 * [A] * k_eMoM * xM^0 * dqdB
       */
      lin_algebra->matrixInsert(Jacobian,
        8. * delx / 9. * k_eMoM * q,
        M0_idx,
        A_idx); // dM0'/dA

      lin_algebra->matrixInsert(Jacobian,
        8. * delx / 9. * A * k_eMoM * dqdB,
        M0_idx,
        M_idx); // dM0'/dBM
        

      }

      return 0;
    };
  
    return fcn;
  }

  

  template<typename Real>
  unsigned int
  particleSystem<Real>::getNumberOfSpecies() const {
    return n_species;
  }



  template<typename Real>
  unsigned int
  particleSystem<Real>::getNumberOfParticleBins() const {
    return n_particles;
  }



  template<typename Real>
  unsigned int
  particleSystem<Real>::getNumberOfReactions() const {
    return all_rxns.size();
  }

  template<typename Real>
  void
  particleSystem<Real>::printJacobianSparsityPattern(const std::string out_filename, abstractLinearAlgebraOperations<Real> *lin_alg){
    std::ofstream myfile;
    myfile.open(out_filename);
    auto jFcn = composeJacobianfunction();
    
    N_Vector x = lin_alg->createNewVector(n_species);
    N_VConst(1.0, x);
    N_Vector x_dot = N_VClone(x);
    N_Vector tmp1 = N_VClone(x);
    N_Vector tmp2 = N_VClone(x);
    N_Vector tmp3 = N_VClone(x);

    SUNMatrix J = lin_alg->createNewMatrix(n_species, n_species);
    void * user_data = static_cast<void *>(lin_alg);
    jFcn(0.0, x, x_dot, J, user_data, tmp1, tmp2, tmp3);

    for(unsigned int row=0; row<n_species;++row){
      for (unsigned int col=0; col<n_species;++col){
        auto jval = lin_alg->matrixGetValue(J, row, col);
        if (jval != 0.0){
          myfile << "#";
        } else {
          myfile << " ";
        }
        if (col == n_species-1){
          myfile << "\n";
        }
      }
    }

    myfile.close();
    SUNMatDestroy(J);
    N_VDestroy(x);
    N_VDestroy(x_dot);
    N_VDestroy(tmp1);
    N_VDestroy(tmp2);
    N_VDestroy(tmp3);
  }


 

}
#endif