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

#include "pair_morse_cuda_cu.h"
#include "pair_morse_cuda_kernel_nc.cu"
#include <time.h>



void Cuda_PairMorseCuda_Init(cuda_shared_data* sdata)
{
  Cuda_Pair_Init_AllStyles(sdata, 4);
}




void Cuda_PairMorseCuda(cuda_shared_data* sdata, cuda_shared_neighlist* sneighlist, int eflag, int vflag, int eflag_atom, int vflag_atom)
{

  // initialize only on first call
  static  short init = 0;

  if(! init) {
    init = 1;
    Cuda_PairMorseCuda_Init(sdata);
  }

  dim3 grid, threads;
  int sharedperproc;

  Cuda_Pair_PreKernel_AllStyles(sdata, sneighlist, eflag, vflag, grid, threads, sharedperproc, false, 256);

  cudaStream_t* streams = (cudaStream_t*) CudaWrapper_returnStreams();

  if(sdata->pair.use_block_per_atom)
    Pair_Kernel_BpA<PAIR_MORSE, COUL_NONE, DATA_NONE>
    <<< grid, threads, sharedperproc* sizeof(ENERGY_CFLOAT)*threads.x, streams[1]>>> (eflag, vflag, eflag_atom, vflag_atom);
  else
    Pair_Kernel_TpA<PAIR_MORSE, COUL_NONE, DATA_NONE>
    <<< grid, threads, sharedperproc* sizeof(ENERGY_CFLOAT)*threads.x, streams[1]>>> (eflag, vflag, eflag_atom, vflag_atom);

  Cuda_Pair_PostKernel_AllStyles(sdata, grid, sharedperproc, eflag, vflag);
}

#undef _r0
#undef _alpha
#undef _morse1
#undef _d0


