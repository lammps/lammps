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

#ifdef INTEGRATE_CLASS

IntegrateStyle(verlet/cuda,VerletCuda)

#else


#ifndef LMP_VERLET_CUDA_H
#define LMP_VERLET_CUDA_H
#include "verlet.h"
#include "modify_cuda.h"

namespace LAMMPS_NS {

class VerletCuda : public Verlet
{
        public:
                VerletCuda(class LAMMPS *, int, char **);
                void setup();
                 void setup_minimal(int);
                  void run(int);

                void test_atom(int atom,char* astring); //debugging purpose
                int dotestatom;        //debugging purpose

        protected:
                class Cuda *cuda;
                void force_clear();
            double time_pair;
            double time_kspace;
            double time_comm;
            double time_modify;
            double time_fulliterate;
            ModifyCuda* modify_cuda;
};

}

#endif
#endif
