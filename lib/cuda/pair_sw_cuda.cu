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

#include "pair_sw_cuda_cu.h"
__device__ __constant__ ParamSW_Float params_sw[MANYBODY_NPAIR* MANYBODY_NPAIR* MANYBODY_NPAIR];

#include "pair_sw_cuda_kernel_nc.cu"

#include <time.h>


void Cuda_PairSWCuda_Init(cuda_shared_data* sdata, ParamSW_Float* params_host, void* map_host, void* elem2param_host, int nelements_h)
{
  unsigned cuda_ntypes = sdata->atom.ntypes + 1;
  X_FLOAT box_size[3] = {
    sdata->domain.subhi[0] - sdata->domain.sublo[0],
    sdata->domain.subhi[1] - sdata->domain.sublo[1],
    sdata->domain.subhi[2] - sdata->domain.sublo[2]
  };

  cudaMemcpyToSymbol(MY_AP(box_size)     , box_size                      , sizeof(X_FLOAT) * 3);
  cudaMemcpyToSymbol(MY_AP(cuda_ntypes)  , &cuda_ntypes                   , sizeof(unsigned));
  cudaMemcpyToSymbol(MY_AP(virial)       , &sdata->pair.virial.dev_data   , sizeof(ENERGY_FLOAT*));
  cudaMemcpyToSymbol(MY_AP(eng_vdwl)     , &sdata->pair.eng_vdwl.dev_data , sizeof(ENERGY_FLOAT*));
  cudaMemcpyToSymbol(MY_AP(periodicity)  , sdata->domain.periodicity     , sizeof(int) * 3);
  cudaMemcpyToSymbol(MY_AP(collect_forces_later), &sdata->pair.collect_forces_later  , sizeof(int));
  cudaMemcpyToSymbol(params_sw, params_host  , sizeof(ParamSW_Float)*nelements_h * nelements_h * nelements_h);
  cudaMemcpyToSymbol(elem2param, elem2param_host  , sizeof(int)*nelements_h * nelements_h * nelements_h);
  cudaMemcpyToSymbol(map, map_host  , sizeof(int)*cuda_ntypes);
  cudaMemcpyToSymbol(nelements, &nelements_h, sizeof(int));
}

void Cuda_PairSWCuda(cuda_shared_data* sdata, cuda_shared_neighlist* sneighlist, int eflag, int vflag, int eflag_atom, int vflag_atom)
{
  static int glob_ij_size = 0;
  static F_FLOAT4* glob_r_ij = NULL;
  static int* glob_numneigh_red = NULL;
  static int* glob_neighbors_red = NULL;
  static int* glob_neightype_red = NULL;

  if(glob_ij_size < sdata->atom.nall * sneighlist->maxneighbors * sizeof(F_FLOAT)) {
    glob_ij_size = sdata->atom.nall * sneighlist->maxneighbors * sizeof(F_FLOAT);
    cudaFree(glob_r_ij);
    cudaFree(glob_numneigh_red);
    cudaFree(glob_neighbors_red);
    cudaFree(glob_neightype_red);
    cudaMalloc(&glob_r_ij, glob_ij_size * 4);
    cudaMalloc(&glob_numneigh_red, sdata->atom.nall * sizeof(int));
    cudaMalloc(&glob_neighbors_red, sdata->atom.nall * sneighlist->maxneighbors * sizeof(int));
    cudaMalloc(&glob_neightype_red, sdata->atom.nall * sneighlist->maxneighbors * sizeof(int));
    cudaMemcpyToSymbol(_glob_numneigh_red, &glob_numneigh_red  , sizeof(int*));
    cudaMemcpyToSymbol(_glob_neighbors_red, &glob_neighbors_red  , sizeof(int*));
    cudaMemcpyToSymbol(_glob_neightype_red, &glob_neightype_red  , sizeof(int*));
    cudaMemcpyToSymbol(_glob_r_ij, &glob_r_ij  , sizeof(F_FLOAT4*));
  }

  dim3 grid, threads;
  int sharedperproc;

  Cuda_Pair_PreKernel_AllStyles(sdata, sneighlist, eflag, vflag, grid, threads, sharedperproc, false, 64);
  cudaStream_t* streams = (cudaStream_t*) CudaWrapper_returnStreams();



  dim3 grid2;

  if(sdata->atom.nall <= 256 * 64000) {
    grid2.x = (sdata->atom.nall + 255) / 256;
    grid2.y = 1;
  } else {
    grid2.x = (sdata->atom.nall + 256 * 128 - 1) / (256 * 128);
    grid2.y = 128;
  }

  grid2.z = 1;
  dim3 threads2;
  threads2.x = 256;
  threads2.y = 1;
  threads2.z = 1;

  timespec time1, time2;

  //pre-calculate all neighbordistances and zeta_ij
  clock_gettime(CLOCK_REALTIME, &time1);
  Pair_SW_Kernel_TpA_RIJ <<< grid2, threads2, 0, streams[1]>>>();
  cudaThreadSynchronize();
  clock_gettime(CLOCK_REALTIME, &time2);
  sdata->cuda_timings.test1 +=
    time2.tv_sec - time1.tv_sec + 1.0 * (time2.tv_nsec - time1.tv_nsec) / 1000000000;
  clock_gettime(CLOCK_REALTIME, &time1);

  //actual force calculation
  unsigned int sharedsize = (sharedperproc * sizeof(ENERGY_FLOAT) + 4 * sizeof(F_FLOAT)) * threads.x; //extra 4 floats per thread used to reduce register pressure

  if(eflag) {
    if(vflag)
      Pair_SW_Kernel_TpA<1, 1> <<< grid, threads, sharedsize, streams[1]>>>
      (eflag_atom, vflag_atom);
    else
      Pair_SW_Kernel_TpA<1, 0> <<< grid, threads, sharedsize, streams[1]>>>
      (eflag_atom, vflag_atom);
  } else {
    if(vflag)
      Pair_SW_Kernel_TpA<0, 1> <<< grid, threads, sharedsize, streams[1]>>>
      (eflag_atom, vflag_atom);
    else
      Pair_SW_Kernel_TpA<0, 0> <<< grid, threads, sharedsize, streams[1]>>>
      (eflag_atom, vflag_atom);
  }
  cudaThreadSynchronize();
  clock_gettime(CLOCK_REALTIME, &time2);
  sdata->cuda_timings.test2 +=
    time2.tv_sec - time1.tv_sec + 1.0 * (time2.tv_nsec - time1.tv_nsec) / 1000000000;

  Cuda_Pair_PostKernel_AllStyles(sdata, grid, sharedperproc, eflag, vflag);
}

