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

#define _r0 MY_AP(coeff1)
#define _alpha MY_AP(coeff2)
#define _morse1 MY_AP(coeff3)
#define _d0 MY_AP(coeff4)
#define _c0 MY_AP(coeff5)

#include "pair_morse_coul_long_cuda_cu.h"
#include "pair_morse_coul_long_cuda_kernel_nc.cu"

#include <time.h>

void Cuda_PairMorseCoulLongCuda_Init(cuda_shared_data* sdata)
{
  Cuda_Pair_Init_AllStyles(sdata, 5, true);
}

void Cuda_PairMorseCoulLongCuda(cuda_shared_data* sdata, cuda_shared_neighlist* sneighlist, int eflag, int vflag, int eflag_atom, int vflag_atom)
{


  static  short init = 0;

  if(! init) {
    init = 1;
    Cuda_PairMorseCoulLongCuda_Init(sdata);
  }

  dim3 grid, threads;
  int sharedperproc;


  Cuda_Pair_PreKernel_AllStyles(sdata, sneighlist, eflag, vflag, grid, threads, sharedperproc, true, 192);

  cudaStream_t* streams = (cudaStream_t*) CudaWrapper_returnStreams();

  if(sdata->pair.use_block_per_atom)
    Pair_Kernel_BpA<PAIR_MORSE_R6, COUL_LONG, DATA_NONE>
    <<< grid, threads, sharedperproc* sizeof(ENERGY_CFLOAT)*threads.x, streams[1]>>> (eflag, vflag, eflag_atom, vflag_atom);
  else
    Pair_Kernel_TpA<PAIR_MORSE_R6, COUL_LONG, DATA_NONE>
    <<< grid, threads, sharedperproc* sizeof(ENERGY_CFLOAT)*threads.x, streams[1]>>> (eflag, vflag, eflag_atom, vflag_atom);

  Cuda_Pair_PostKernel_AllStyles(sdata, grid, sharedperproc, eflag, vflag);
}



#undef _rhoinv
#undef _sigma
#undef _a
#undef _c
#undef _d
#undef _c0

