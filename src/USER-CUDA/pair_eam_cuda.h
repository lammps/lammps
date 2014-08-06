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

/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
#ifdef PAIR_CLASS

PairStyle(eam/cuda,PairEAMCuda)

#else

#ifndef PAIR_EAM_CUDA_H
#define PAIR_EAM_CUDA_H

#include "cuda_data.h"
#include "pair_eam.h"

namespace LAMMPS_NS {

class PairEAMCuda : public PairEAM
{
        public:
                PairEAMCuda(class LAMMPS *);
                void compute(int, int);
                void settings(int, char **);
                void coeff(int, char **);
                void init_list(int, class NeighList *);
                void init_style();
                void array2spline();
                int pack_forward_comm(int n, int *iswap, double *buf, 
                                      int pbc_flag, int *pbc);
                void unpack_forward_comm(int n, int first, double *buf);
        protected:
                class Cuda *cuda;
                void allocate();
                bool allocated2;
                virtual void ev_setup(int eflag, int vflag);
                class CudaNeighList* cuda_neigh_list;
                cCudaData<double, F_FLOAT, x>* cu_rho;
                cCudaData<double, F_FLOAT, x>* cu_fp;
            cCudaData<double, F_FLOAT, xyz>* cu_rhor_spline;
            cCudaData<double, F_FLOAT, xyz>* cu_z2r_spline;
            cCudaData<double, F_FLOAT, xyz>* cu_frho_spline;

};

}

#endif
#endif
