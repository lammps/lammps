/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator

   Original Version:
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   See the README file in the top-level LAMMPS directory.

   -----------------------------------------------------------------------

   USER-CUDA Package and associated modifications:
   https://sourceforge.net/projects/lammpscuda/

   Christian Trott, christian.trott@tu-ilmenau.de
   Lars Winterfeld, lars.winterfeld@tu-ilmenau.de
   Theoretical Physics II, University of Technology Ilmenau, Germany

   See the README file in the USER-CUDA directory.

   This software is distributed under the GNU General Public License.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(tersoff/zbl/cuda,PairTersoffZBLCuda)

#else

#ifndef PAIR_TERSOFF_ZBL_CUDA_H
#define PAIR_TERSOFF_ZBL_CUDA_H

#include "pair_tersoff_cuda.h"

namespace LAMMPS_NS {

class PairTersoffZBLCuda : public PairTersoffCuda
{
        public:
                PairTersoffZBLCuda(class LAMMPS *);
         private:
          double global_a_0;    // Bohr radius for Coulomb repulsion
          double global_epsilon_0;  // permittivity of vacuum for Coulomb repulsion
          double global_e;    // proton charge (negative of electron charge)

          void read_file(char *);
          void coeff(int narg, char **arg);
};

}

#endif
#endif
