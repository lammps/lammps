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

#include "pair_lj_cut_cuda_cu.h"
#include "pair_lj_cut_cuda_kernel_nc.cu"

#include <time.h>
#include <sys/time.h>

static unsigned long long int *d_nb_blocks_done=NULL;
static bool g_science_init=false;

void Cuda_PairLJCutCuda_Init(cuda_shared_data* sdata)
{
  Cuda_Pair_Init_AllStyles(sdata, 4);
}

void Cuda_PairLJCutCuda(cuda_shared_data* sdata, cuda_shared_neighlist* sneighlist, int eflag, int vflag, int eflag_atom, int vflag_atom)
{


  static  short init = 0;

  if(! init) {
    init = 1;
    Cuda_PairLJCutCuda_Init(sdata);
  }

  dim3 grid, threads;
  int sharedperproc;

  //int maxthreads=192*sizeof(double)/sizeof(F_CFLOAT);
  //if(CUDA_ARCH==20) maxthreads*=2;
  //cudaFuncSetCacheConfig(Pair_Kernel_TpA_opt<PAIR_LJ_CUT,COUL_NONE,DATA_NONE>,cudaFuncCachePreferL1);
  Cuda_Pair_PreKernel_AllStyles(sdata, sneighlist, eflag, vflag, grid, threads, sharedperproc, false, 192);
  cudaStream_t* streams = (cudaStream_t*) CudaWrapper_returnStreams();

  if(!g_science_init){
	  cudaMalloc(&d_nb_blocks_done, sizeof(unsigned long long int));
	  g_science_init=true;
  }

  if(sdata->pair.use_block_per_atom){
    Pair_Kernel_BpA<PAIR_LJ_CUT, COUL_NONE, DATA_NONE>
    <<< grid, threads, sharedperproc* sizeof(ENERGY_CFLOAT)*threads.x, streams[1]>>> (eflag, vflag, eflag_atom, vflag_atom);
  }else{
    //Pair_Kernel_TpA<PAIR_LJ_CUT, COUL_NONE, DATA_NONE>
    //<<< grid, threads, sharedperproc* sizeof(ENERGY_CFLOAT)*threads.x, streams[1]>>> (eflag, vflag, eflag_atom, vflag_atom);
    //struct timeval start, end;
    unsigned long long int nb_blocks_todo = (grid.x * grid.y * grid.z), tmp=0;
    cudaMemcpyAsync(d_nb_blocks_done, &tmp, sizeof(unsigned long long int), cudaMemcpyHostToDevice, streams[1]);
    //cudaDeviceSynchronize();
    //gettimeofday(&start, NULL);
    Pair_Kernel_TpA_sw<PAIR_LJ_CUT, COUL_NONE, DATA_NONE>
    <<< 240, threads, sharedperproc* sizeof(ENERGY_CFLOAT)*threads.x, streams[1]>>> (d_nb_blocks_done, nb_blocks_todo, grid, eflag, vflag, eflag_atom, vflag_atom); //240 => K40 # of SMs x Max # of blocks/SM; TODO: compute from cudaGetDeviceProperties(...)
    //cudaDeviceSynchronize();
    //gettimeofday(&end, NULL);
    //cudaMemcpyAsync(&tmp, d_nb_blocks_done, sizeof(unsigned long long int), cudaMemcpyDeviceToHost, streams[1]);
    //cudaDeviceSynchronize();
    //printf("Pair_Kernel_TpA_sw nb_blocks_todo = %llu, nb_blocks_done = %llu, duration(us) = %llu\n", nb_blocks_todo, tmp-240, (end.tv_sec-start.tv_sec)*1000000+(end.tv_usec-start.tv_usec)); //Each block increments the global counter once extra
  }

  Cuda_Pair_PostKernel_AllStyles(sdata, grid, sharedperproc, eflag, vflag);
}


#undef _lj1
#undef _lj2
#undef _lj3
#undef _lj4
