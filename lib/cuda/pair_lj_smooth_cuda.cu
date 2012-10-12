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

#include <stdio.h>

#define _lj1 MY_AP(coeff1_gm)
#define _lj2 MY_AP(coeff2_gm)
#define _lj3 MY_AP(coeff3_gm)
#define _lj4 MY_AP(coeff4_gm)
#define _ljsw1 MY_AP(coeff5_gm)
#define _ljsw2 MY_AP(coeff6_gm)
#define _ljsw3 MY_AP(coeff7_gm)
#define _ljsw4 MY_AP(coeff8_gm)
#define _ljsw0 MY_AP(coeff9_gm)

#include "pair_lj_smooth_cuda_cu.h"
#include "pair_lj_smooth_cuda_kernel_nc.cu"

#include <time.h>

void Cuda_PairLJSmoothCuda_Init(cuda_shared_data* sdata)
{
  Cuda_Pair_Init_AllStyles(sdata, 9, false, true, true);
}



void Cuda_PairLJSmoothCuda(cuda_shared_data* sdata, cuda_shared_neighlist* sneighlist, int eflag, int vflag,
                           int eflag_atom, int vflag_atom)
{
  // initialize only on first call
  static  short init = 0;

  if(! init) {
    init = 1;
    Cuda_PairLJSmoothCuda_Init(sdata);
  }

  dim3 grid, threads;
  int sharedperproc;

  Cuda_Pair_PreKernel_AllStyles(sdata, sneighlist, eflag, vflag, grid, threads, sharedperproc, false, 192);

  cudaStream_t* streams = (cudaStream_t*) CudaWrapper_returnStreams();

  if(sdata->pair.use_block_per_atom)
    Pair_Kernel_BpA<PAIR_LJ_SMOOTH, COUL_NONE, DATA_NONE>
    <<< grid, threads, sharedperproc* sizeof(ENERGY_FLOAT)*threads.x, streams[1]>>> (eflag, vflag, eflag_atom, vflag_atom);
  else
    Pair_Kernel_TpA<PAIR_LJ_SMOOTH, COUL_NONE, DATA_NONE>
    <<< grid, threads, sharedperproc* sizeof(ENERGY_FLOAT)*threads.x, streams[1]>>> (eflag, vflag, eflag_atom, vflag_atom);

  Cuda_Pair_PostKernel_AllStyles(sdata, grid, sharedperproc, eflag, vflag);
}

#undef _lj1
#undef _lj2
#undef _lj3
#undef _lj4
#undef _ljsw1
#undef _ljsw2
#undef _ljsw3
#undef _ljsw4
#undef _ljsw0
