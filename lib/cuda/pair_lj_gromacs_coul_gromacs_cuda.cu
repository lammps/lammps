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
#define _ljsw5 MY_AP(coeff9_gm)

#define _cut_coul_inner_global MY_AP(cut_coul_inner_global)
#define _coulsw1 MY_AP(coulsw1)
#define _coulsw2 MY_AP(coulsw2)
#define _coulsw5 MY_AP(coulsw5)
__device__ __constant__ F_CFLOAT _cut_coul_inner_global;
__device__ __constant__ F_CFLOAT _coulsw1;
__device__ __constant__ F_CFLOAT _coulsw2;
__device__ __constant__ F_CFLOAT _coulsw5;


#include "pair_lj_gromacs_coul_gromacs_cuda_cu.h"
#include "pair_lj_gromacs_coul_gromacs_cuda_kernel_nc.cu"

#include <time.h>

void Cuda_PairLJGromacsCoulGromacsCuda_Init(cuda_shared_data* sdata, F_CFLOAT cut_coul_inner, F_CFLOAT coulsw1, F_CFLOAT coulsw2, F_CFLOAT coulsw5)
{
  Cuda_Pair_Init_AllStyles(sdata, 9, true, true, true);
  cudaMemcpyToSymbol(MY_AP(cut_coul_inner_global) , &cut_coul_inner  , sizeof(F_CFLOAT));
  cudaMemcpyToSymbol(MY_AP(coulsw1) , &coulsw1  , sizeof(F_CFLOAT));
  cudaMemcpyToSymbol(MY_AP(coulsw2) , &coulsw2  , sizeof(F_CFLOAT));
  cudaMemcpyToSymbol(MY_AP(coulsw5) , &coulsw5  , sizeof(F_CFLOAT));

  return;
}



void Cuda_PairLJGromacsCoulGromacsCuda(cuda_shared_data* sdata, cuda_shared_neighlist* sneighlist, int eflag, int vflag,
                                       int eflag_atom, int vflag_atom, F_CFLOAT cut_coul_inner, F_CFLOAT coulsw1, F_CFLOAT coulsw2, F_CFLOAT coulsw5)
{
  static  short init = 0;

  if(! init) {
    init = 1;
    Cuda_PairLJGromacsCoulGromacsCuda_Init(sdata, cut_coul_inner, coulsw1, coulsw2, coulsw5);
  }

  dim3 grid, threads;
  int sharedperproc;

  Cuda_Pair_PreKernel_AllStyles(sdata, sneighlist, eflag, vflag, grid, threads, sharedperproc, true, 192);

  cudaStream_t* streams = (cudaStream_t*) CudaWrapper_returnStreams();

  if(sdata->pair.use_block_per_atom)
    Pair_Kernel_BpA<PAIR_LJ_GROMACS, COUL_GROMACS, DATA_NONE>
    <<< grid, threads, sharedperproc* sizeof(ENERGY_CFLOAT)*threads.x, streams[1]>>> (eflag, vflag, eflag_atom, vflag_atom);
  else
    Pair_Kernel_TpA<PAIR_LJ_GROMACS, COUL_GROMACS, DATA_NONE>
    <<< grid, threads, sharedperproc* sizeof(ENERGY_CFLOAT)*threads.x, streams[1]>>> (eflag, vflag, eflag_atom, vflag_atom);

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
#undef _ljsw5
#undef _cut_coul_inner_global
#undef _coulsw1
#undef _coulsw2
#undef _coulsw5
