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

PairStyle(lj/gromacs/coul/gromacs/cuda,PairLJGromacsCoulGromacsCuda)

#else

#ifndef LMP_PAIR_LJ_GROMACS_COUL_GROMACS_CUDA_H
#define LMP_PAIR_LJ_GROMACS_COUL_GROMACS_CUDA_H

#include "pair_lj_gromacs_coul_gromacs.h"
#include "cuda_data.h"

namespace LAMMPS_NS {

class PairLJGromacsCoulGromacsCuda : public PairLJGromacsCoulGromacs
{
        public:
                PairLJGromacsCoulGromacsCuda(class LAMMPS *);
                void compute(int, int);
                void settings(int, char **);
                void coeff(int, char **);
                void init_list(int, class NeighList *);
                void init_style();
                void ev_setup(int eflag, int vflag);
        protected:
                class Cuda *cuda;
                void allocate();
                bool allocated2;
                class CudaNeighList* cuda_neigh_list;
                cCudaData<double  , F_CFLOAT , x >* cu_lj1_gm;
                cCudaData<double  , F_CFLOAT , x >* cu_lj2_gm;
                cCudaData<double  , F_CFLOAT , x >* cu_lj3_gm;
                cCudaData<double  , F_CFLOAT , x >* cu_lj4_gm;
                cCudaData<double  , F_CFLOAT , x >* cu_ljsw1_gm;
                cCudaData<double  , F_CFLOAT , x >* cu_ljsw2_gm;
                cCudaData<double  , F_CFLOAT , x >* cu_ljsw3_gm;
                cCudaData<double  , F_CFLOAT , x >* cu_ljsw4_gm;
                cCudaData<double  , F_CFLOAT , x >* cu_ljsw5_gm;

};

}

#endif
#endif
