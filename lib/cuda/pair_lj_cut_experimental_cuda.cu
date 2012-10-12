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

#define _lj1 MY_AP(coeff1)
#define _lj2 MY_AP(coeff2)
#define _lj3 MY_AP(coeff3)
#define _lj4 MY_AP(coeff4)

#include "pair_lj_cut_experimental_cuda_cu.h"

#include <time.h>

void Cuda_PairLJCutExperimentalCuda_Init(cuda_shared_data* sdata)
{
  Cuda_Pair_Init_AllStyles(sdata, 4);
}

void Cuda_PairLJCutExperimentalCuda(cuda_shared_data* sdata, cuda_shared_neighlist* sneighlist, int eflag, int vflag, int eflag_atom, int vflag_atom)
{


  static  short init = 0;

  if(! init) {
    init = 1;
    Cuda_PairLJCutExperimentalCuda_Init(sdata);
  }

  dim3 grid, threads;
  int sharedperproc;

  //int maxthreads=192*sizeof(double)/sizeof(F_FLOAT);
  //if(CUDA_ARCH==20) maxthreads*=2;
  //cudaFuncSetCacheConfig(Pair_Kernel_TpA_opt<PAIR_LJ_CUT,COUL_NONE,DATA_NONE>,cudaFuncCachePreferL1);
  Cuda_Pair_PreKernel_AllStyles(sdata, sneighlist, eflag, vflag, grid, threads, sharedperproc, false, 192);

  if(sharedperproc == 0) sharedperproc++;

  //printf("comm_phase: %i\n",sdata->comm.comm_phase);

  cudaStream_t* streams = (cudaStream_t*) CudaWrapper_returnStreams();

  if(sdata->pair.use_block_per_atom)
    Pair_Kernel_BpA<PAIR_LJ_CUT, COUL_NONE, DATA_NONE>
    <<< grid, threads, sharedperproc* sizeof(ENERGY_FLOAT)*threads.x, streams[1]>>> (eflag, vflag, eflag_atom, vflag_atom);
  else
    Pair_Kernel_TpA_opt<PAIR_LJ_CUT, COUL_NONE, DATA_NONE>
    <<< grid, threads, sharedperproc* sizeof(ENERGY_FLOAT)*threads.x, streams[1]>>> (eflag, vflag, eflag_atom, vflag_atom, sdata->comm.comm_phase);
  Cuda_Pair_PostKernel_AllStyles(sdata, grid, sharedperproc, eflag, vflag);
}


#undef _lj1
#undef _lj2
#undef _lj3
#undef _lj4
