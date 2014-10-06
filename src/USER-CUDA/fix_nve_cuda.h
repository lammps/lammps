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

#ifdef FIX_CLASS

FixStyle(nve/cuda,FixNVECuda)

#else

#ifndef LMP_FIX_NVE_CUDA_H
#define LMP_FIX_NVE_CUDA_H

#include "fix.h"
#include "cuda_precision.h"

namespace LAMMPS_NS {

class FixNVECuda : public Fix
{
        public:
                FixNVECuda(class LAMMPS *, int, char **);
                int setmask();
                virtual void init();
                virtual void initial_integrate(int);
                virtual void final_integrate();
                void initial_integrate_respa(int, int, int);
                void final_integrate_respa(int, int);
                void reset_dt();

                X_CFLOAT triggerneighsq;

        protected:
                class Cuda *cuda;
                double dtv, dtf;
                double *step_respa;
                int mass_require;

};

}

#endif
#endif
