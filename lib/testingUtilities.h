#ifndef NANOSIM_TESTINGUTILITIES_H
#define NANOSIM_TESTINGUTILITIES_H

#include "particleSystem.h"
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_types.h>

namespace NanoSim{
  namespace Testing{
    template<typename Real, typename LinearAlgebraType>
    void
    printParticleSystemOutput(NanoSim::particleSystem<Real> & my_rxns,
      N_Vector x,
      N_Vector x_dot,
      SUNMatrix J,
      LinearAlgebraType lin_alg){
        my_rxns.printChemicalReactions();

        const sunindextype length = x->ops->nvgetlength(x);

        std::cout << "x = \n";
        for (sunindextype i=0; i<length; ++i){
          std::cout << lin_alg.vectorGetValue(x,i) << "\n";
        }

        std::cout << "dx = \n";
        for (sunindextype i=0; i<length; ++i){
          std::cout << lin_alg.vectorGetValue(x_dot,i) << "\n";
        }

        std::cout << "J =\n";
        for (sunindextype row=0; row<length; ++row){
          for (sunindextype col=0; col<length; ++col){
            std::cout << lin_alg.matrixGetValue(J, row, col) << " ";
            if (col == length-1){
              std::cout << "\n";
            }
          }
        }
      }
  }
}

#endif