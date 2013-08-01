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
#include <time.h>
#define MY_PREFIX neighbor
#define IncludeCommonNeigh
#include "cuda_shared.h"
#include "cuda_common.h"
#include "crm_cuda_utils.cu"
#include "cuda_wrapper_cu.h"

#define _cutneighsq     MY_AP(cutneighsq)
#define _ex_type     	MY_AP(ex_type)
#define _nex_type     	MY_AP(nex_type)
#define _ex1_bit     	MY_AP(ex1_bit)
#define _ex2_bit     	MY_AP(ex2_bit)
#define _nex_group     	MY_AP(nex_group)
#define _ex_mol_bit     MY_AP(ex_mol_bit)
#define _nex_mol     	MY_AP(nex_mol)
__device__ __constant__ CUDA_FLOAT* _cutneighsq;
__device__ __constant__ int* _ex_type;
__device__ __constant__ int _nex_type;
__device__ __constant__ int* _ex1_bit;
__device__ __constant__ int* _ex2_bit;
__device__ __constant__ int _nex_group;
__device__ __constant__ int* _ex_mol_bit;
__device__ __constant__ int _nex_mol;

#include "neighbor_cu.h"
#include "neighbor_kernel.cu"

void Cuda_Neighbor_UpdateBuffer(cuda_shared_data* sdata, cuda_shared_neighlist* sneighlist)
{
  CUT_CHECK_ERROR("Cuda_PairLJCutCuda: before updateBuffer failed");

  int size = (unsigned)(sizeof(int) * 20 + sneighlist->bin_dim[0] * sneighlist->bin_dim[1] * sneighlist->bin_dim[2] * (sizeof(int) + sneighlist->bin_nmax * 3 * sizeof(CUDA_FLOAT)));

  if(sdata->buffersize < size) {
    MYDBG(printf("Cuda_Neighbor Resizing Buffer at %p with %i kB to\n", sdata->buffer, sdata->buffersize);)

    if(sdata->buffer != NULL) CudaWrapper_FreeCudaData(sdata->buffer, sdata->buffersize);

    sdata->buffer = CudaWrapper_AllocCudaData(size);
    sdata->buffersize = size;
    sdata->buffer_new++;
    MYDBG(printf("New buffer at %p with %i kB\n", sdata->buffer, sdata->buffersize);)
  }

  cudaMemcpyToSymbol(MY_AP(buffer), & sdata->buffer, sizeof(int*));
  CUT_CHECK_ERROR("Cuda_PairLJCutCuda: updateBuffer failed");
}

int Cuda_BinAtoms(cuda_shared_data* sdata, cuda_shared_neighlist* sneighlist)
{
  if(sdata->buffer_new)
    Cuda_Neighbor_UpdateBuffer(sdata, sneighlist);

  // initialize only on first call
  CUDA_FLOAT rez_bin_size[3] = {
    (1.0 * sneighlist->bin_dim[0] - 4.0) / (sdata->domain.subhi[0] - sdata->domain.sublo[0]),
    (1.0 * sneighlist->bin_dim[1] - 4.0) / (sdata->domain.subhi[1] - sdata->domain.sublo[1]),
    (1.0 * sneighlist->bin_dim[2] - 4.0) / (sdata->domain.subhi[2] - sdata->domain.sublo[2])
  };

  short init = 0;

  if(! init) {
    init = 0;
    cudaMemcpyToSymbol(MY_AP(x)              , & sdata->atom.x         .dev_data, sizeof(X_FLOAT*));
    cudaMemcpyToSymbol(MY_AP(nall)         , & sdata->atom.nall                    , sizeof(unsigned));
    cudaMemcpyToSymbol(MY_AP(nmax)           , & sdata->atom.nmax                    , sizeof(unsigned));
    cudaMemcpyToSymbol(MY_AP(sublo)          ,   sdata->domain.sublo                 , sizeof(X_FLOAT) * 3);
  }


  int3 layout = getgrid(sdata->atom.nall); // sneighlist->inum
  dim3 threads(layout.z, 1, 1);
  dim3 grid(layout.x, layout.y, 1);

  my_times starttime, endtime;
  my_gettime(CLOCK_REALTIME, &starttime);

  cudaMemset((int*)(sdata->buffer), 0, sizeof(int) * (20 + (sneighlist->bin_dim[0]) * (sneighlist->bin_dim[1]) * (sneighlist->bin_dim[2])) + 3 * sizeof(CUDA_FLOAT) * (sneighlist->bin_dim[0]) * (sneighlist->bin_dim[1]) * (sneighlist->bin_dim[2]) * (sneighlist->bin_nmax));

  Binning_Kernel <<< grid, threads>>> (sneighlist->binned_id, sneighlist->bin_nmax, sneighlist->bin_dim[0], sneighlist->bin_dim[1], sneighlist->bin_dim[2], rez_bin_size[0], rez_bin_size[1], rez_bin_size[2]);
  cudaThreadSynchronize();

  my_gettime(CLOCK_REALTIME, &endtime);
  sdata->cuda_timings.neigh_bin +=
    endtime.tv_sec - starttime.tv_sec + 1.0 * (endtime.tv_nsec - starttime.tv_nsec) / 1000000000;


  int binning_error;
  cudaMemcpy((void*) &binning_error, (void*) sdata->buffer, 1 * sizeof(int), cudaMemcpyDeviceToHost);

  if(binning_error) {
    sneighlist->bin_extraspace += 0.05;
  } else {
    MYDBG(printf("CUDA: binning successful\n");)
  }
  CUT_CHECK_ERROR("Cuda_Binning: binning Kernel execution failed");
  return binning_error;
}

int Cuda_NeighborBuildFullBin(cuda_shared_data* sdata, cuda_shared_neighlist* sneighlist)
{
  //Cuda_Neighbor_UpdateBuffer(sdata,sneighlist);
  CUDA_FLOAT globcutoff = -1.0;

  short init = 0;

  if(! init) {
    init = 1;

    // !! LAMMPS indexes atom types starting with 1 !!

    unsigned cuda_ntypes = sdata->atom.ntypes + 1;

    unsigned nx = sizeof(CUDA_FLOAT) * cuda_ntypes * cuda_ntypes;

    CUDA_FLOAT* acutneighsq = (CUDA_FLOAT*) malloc(nx);
    //printf("Allocate: %i\n",nx);
    sneighlist->cu_cutneighsq = (CUDA_FLOAT*) CudaWrapper_AllocCudaData(nx);

    if(sneighlist->cutneighsq) {
      int cutoffsdiffer = 0;
      double cutoff0 = sneighlist->cutneighsq[1][1];

      for(int i = 1; i <= sdata->atom.ntypes; ++i) {
        for(int j = 1; j <= sdata->atom.ntypes; ++j) {
          acutneighsq[i * cuda_ntypes + j] = (CUDA_FLOAT)(sneighlist->cutneighsq[i][j]);

          if((sneighlist->cutneighsq[i][j] - cutoff0) * (sneighlist->cutneighsq[i][j] - cutoff0) > 1e-6) cutoffsdiffer++;
        }
      }

      if(not cutoffsdiffer) globcutoff = (CUDA_FLOAT) cutoff0;
    } else {
      MYEMUDBG(printf("# CUDA: Cuda_NeighborBuild: cutneighsq == NULL\n");)
      return 0;
    }

    int size = 100;

    if(sdata->buffersize < size) {
      MYDBG(printf("Cuda_NeighborBuild Resizing Buffer at %p with %i kB to\n", sdata->buffer, sdata->buffersize);)
      CudaWrapper_FreeCudaData(sdata->buffer, sdata->buffersize);
      sdata->buffer = CudaWrapper_AllocCudaData(size);
      sdata->buffersize = size;
      sdata->buffer_new++;
      MYDBG(printf("New buffer at %p with %i kB\n", sdata->buffer, sdata->buffersize);)
    }

    CudaWrapper_UploadCudaData(acutneighsq, sneighlist->cu_cutneighsq, nx);
    cudaMemcpyToSymbol(MY_AP(cutneighsq)       , &sneighlist->cu_cutneighsq       , sizeof(CUDA_FLOAT*));

    cudaMemcpyToSymbol(MY_AP(cuda_ntypes)      , & cuda_ntypes                    , sizeof(unsigned));
    cudaMemcpyToSymbol(MY_AP(special_flag)     , sdata->atom.special_flag         , 4 * sizeof(int));
    cudaMemcpyToSymbol(MY_AP(molecular)        , & sdata->atom.molecular          , sizeof(int));
  }

  cudaMemcpyToSymbol(MY_AP(neighbor_maxlocal), & sneighlist->firstneigh.dim[0]  , sizeof(unsigned));
  //cudaMemcpyToSymbol(MY_AP(firstneigh)       , & sneighlist->firstneigh.dev_data, sizeof(int*)     );
  cudaMemcpyToSymbol(MY_AP(ilist)            , & sneighlist->ilist     .dev_data, sizeof(int*));
  cudaMemcpyToSymbol(MY_AP(inum)             , & sneighlist->inum               , sizeof(int));
  cudaMemcpyToSymbol(MY_AP(nlocal)           , & sdata->atom.nlocal             , sizeof(int));
  cudaMemcpyToSymbol(MY_AP(nall)             , & sdata->atom.nall            , sizeof(int));
  cudaMemcpyToSymbol(MY_AP(numneigh)         , & sneighlist->numneigh  .dev_data, sizeof(int*));
  cudaMemcpyToSymbol(MY_AP(type)             , & sdata->atom.type      .dev_data, sizeof(int*));
  cudaMemcpyToSymbol(MY_AP(mask)             , & sdata->atom.mask      .dev_data, sizeof(int*));
  cudaMemcpyToSymbol(MY_AP(tag)              , & sdata->atom.tag       .dev_data, sizeof(int*));
  cudaMemcpyToSymbol(MY_AP(special)          , & sdata->atom.special   .dev_data, sizeof(int*));
  cudaMemcpyToSymbol(MY_AP(maxspecial)       , & sdata->atom.maxspecial         , sizeof(int));
  cudaMemcpyToSymbol(MY_AP(nspecial)         , & sdata->atom.nspecial  .dev_data, sizeof(int*));
  cudaMemcpyToSymbol(MY_AP(maxneighbors)     , & sneighlist->maxneighbors	 , sizeof(int));
  cudaMemcpyToSymbol(MY_AP(debugdata)        , & sdata->debugdata	 , sizeof(int*));
  cudaMemcpyToSymbol(MY_AP(overlap_comm)     , & sdata->overlap_comm, sizeof(int));
  cudaMemcpyToSymbol(MY_AP(neighbors) 		  , & sneighlist->neighbors.dev_data, sizeof(int*));
  cudaMemcpyToSymbol(MY_AP(ex_type) 		  , & sneighlist->ex_type.dev_data, sizeof(int*));
  cudaMemcpyToSymbol(MY_AP(ex1_bit) 		  , & sneighlist->ex1_bit.dev_data, sizeof(int*));
  cudaMemcpyToSymbol(MY_AP(ex2_bit) 		  , & sneighlist->ex2_bit.dev_data, sizeof(int*));
  cudaMemcpyToSymbol(MY_AP(ex_mol_bit) 	  , & sneighlist->ex_mol_bit.dev_data, sizeof(int*));
  cudaMemcpyToSymbol(MY_AP(nex_type)     	  , & sneighlist->nex_type, sizeof(int));
  cudaMemcpyToSymbol(MY_AP(nex_group)     	  , & sneighlist->nex_group, sizeof(int));
  cudaMemcpyToSymbol(MY_AP(nex_mol)     	  , & sneighlist->nex_mol, sizeof(int));

  if(sdata->overlap_comm) {
    cudaMemcpyToSymbol(MY_AP(numneigh_border)  , & sneighlist->numneigh_border .dev_data, sizeof(int*));
    cudaMemcpyToSymbol(MY_AP(numneigh_inner)   , & sneighlist->numneigh_inner  .dev_data, sizeof(int*));
    cudaMemcpyToSymbol(MY_AP(neighbors_border) , & sneighlist->neighbors_border.dev_data, sizeof(int*));
    cudaMemcpyToSymbol(MY_AP(neighbors_inner)  , & sneighlist->neighbors_inner .dev_data, sizeof(int*));
    cudaMemcpyToSymbol(MY_AP(ilist_border)     , & sneighlist->ilist_border    .dev_data, sizeof(int*));
    cudaMemcpyToSymbol(MY_AP(inum_border)      , & sneighlist->inum_border     .dev_data, sizeof(int*));
  }

  //dim3 threads(sneighlist->bin_nmax,1,1);
  dim3 threads(MIN(128, sneighlist->bin_nmax), 1, 1);
  dim3 grid(sneighlist->bin_dim[0]*sneighlist->bin_dim[1], sneighlist->bin_dim[2], 1);

  //printf("Configuration: %i %i %i %i %i\n",grid.x,grid.y,threads.x,(sizeof(int)+3*sizeof(X_FLOAT))*threads.x,sneighlist->bin_nmax);
  int buffer[20];
  buffer[0] = 1;
  buffer[1] = 0;
  CudaWrapper_UploadCudaData(buffer, sdata->buffer, 2 * sizeof(int));
  CUT_CHECK_ERROR("Cuda_NeighborBuild: pre neighbor build kernel error");
  //cudaMemset(sdata->debugdata,0,100*sizeof(int));
  unsigned int shared_size = (sizeof(int) + 3 * sizeof(CUDA_FLOAT)) * threads.x;
  MYDBG(printf("Configuration: %i %i %i %u %i\n", grid.x, grid.y, threads.x, shared_size, sneighlist->bin_nmax);)
  //shared_size=2056;
  my_times starttime, endtime;
  my_gettime(CLOCK_REALTIME, &starttime);
  //for(int i=0;i<100;i++)
  {
    if(sdata->overlap_comm)
      NeighborBuildFullBin_OverlapComm_Kernel <<< grid, threads, shared_size>>>
      (sneighlist->binned_id, sneighlist->bin_nmax, sneighlist->bin_dim[0], sneighlist->bin_dim[1], globcutoff, sdata->pair.use_block_per_atom);
    else {
      int exclude = sneighlist->nex_mol | sneighlist->nex_group | sneighlist->nex_type;

      if(exclude)
        NeighborBuildFullBin_Kernel<1> <<< grid, threads, shared_size>>>
        (sneighlist->binned_id, sneighlist->bin_nmax, sneighlist->bin_dim[0], sneighlist->bin_dim[1], globcutoff, sdata->pair.use_block_per_atom, sdata->pair.neighall);
      else
        NeighborBuildFullBin_Kernel<0> <<< grid, threads, shared_size>>>
        (sneighlist->binned_id, sneighlist->bin_nmax, sneighlist->bin_dim[0], sneighlist->bin_dim[1], globcutoff, sdata->pair.use_block_per_atom, sdata->pair.neighall);
    }
    //NeighborBuildFullBin_Kernel_Restrict<<<grid,threads,(2*sizeof(int)+3*sizeof(X_FLOAT))*threads.x+sizeof(int)>>>
    //	(sneighlist->binned_id,sneighlist->bin_nmax,sneighlist->bin_dim[0],sneighlist->bin_dim[1],globcutoff);

    cudaThreadSynchronize();
    CUT_CHECK_ERROR("Cuda_NeighborBuild: neighbor build kernel execution failed");
    my_gettime(CLOCK_REALTIME, &endtime);
    sdata->cuda_timings.neigh_build +=
      endtime.tv_sec - starttime.tv_sec + 1.0 * (endtime.tv_nsec - starttime.tv_nsec) / 1000000000;
    //dim3 threads,grid;
    CudaWrapper_DownloadCudaData(buffer, sdata->buffer, sizeof(int));

    if(buffer[0] >= 0 && true && sdata->atom.molecular) {
      //printf("Find Special: %i %i\n",sneighlist->inum,sdata->atom.nall);
      my_gettime(CLOCK_REALTIME, &starttime);
      int3 layout = getgrid(sdata->atom.nlocal, 0, 512);
      threads.x = layout.z;
      threads.y = 1;
      threads.z = 1;
      grid.x = layout.x;
      grid.y = layout.y;
      grid.z = 1;
      FindSpecial <<< grid, threads>>>(sdata->pair.use_block_per_atom);
      cudaThreadSynchronize();
      CUT_CHECK_ERROR("Cuda_NeighborBuild: FindSpecial kernel execution failed");
      my_gettime(CLOCK_REALTIME, &endtime);
      sdata->cuda_timings.neigh_special +=
        endtime.tv_sec - starttime.tv_sec + 1.0 * (endtime.tv_nsec - starttime.tv_nsec) / 1000000000;
    }
  }
  //printf("Neightime: %lf\n",sdata->cuda_timings.test1);
  CUT_CHECK_ERROR("Cuda_NeighborBuild: neighbor build kernel execution failed");

  //CudaWrapper_DownloadCudaData(buffer, sneighlist->numneigh_border .dev_data, sizeof(int));

  MYDBG(printf("Cuda_NeighborBuildFullBin build neighbor list ... end\n");)
  return buffer[0];
}

int Cuda_NeighborBuildFullNsq(cuda_shared_data* sdata, cuda_shared_neighlist* sneighlist)
{
  MYDBG(printf("Cuda_NeighborBuildFullNsq build neighbor list ... start\n");)
  // initialize only on first call
  /*static*/ short init = 0;

  if(! init) {
    init = 1;

    // !! LAMMPS indexes atom types starting with 1 !!

    unsigned cuda_ntypes = sdata->atom.ntypes + 1;

    if(cuda_ntypes * cuda_ntypes > CUDA_MAX_TYPES2)
      printf("# CUDA: Cuda_PairLJCutCuda_Init: you need %u types. this is more than %u "
             "(assumed at compile time). re-compile with -DCUDA_MAX_TYPES_PLUS_ONE=32 "
             "or ajust this in cuda_common.h\n", cuda_ntypes, CUDA_MAX_TYPES2);

    unsigned nx = sizeof(CUDA_FLOAT) * cuda_ntypes * cuda_ntypes;
    CUDA_FLOAT* acutneighsq = (CUDA_FLOAT*) malloc(nx);

    if(sneighlist->cutneighsq) {
      for(int i = 1; i <= sdata->atom.ntypes; ++i) {
        for(int j = 1; j <= sdata->atom.ntypes; ++j) {
          acutneighsq[i * cuda_ntypes + j] = (CUDA_FLOAT)(sneighlist->cutneighsq[i][j]);
          //printf("CUTOFFS: %i %i %i %e\n",i,j,cuda_ntypes,acutneighsq[i * cuda_ntypes + j]);
        }
      }
    } else {
      MYEMUDBG(printf("# CUDA: Cuda_NeighborBuild: cutneighsq == NULL\n");)
      return 0;
    }

    int size = 100;

    if(sdata->buffersize < size) {
      MYDBG(printf("Cuda_NeighborBuild Resizing Buffer at %p with %i kB to\n", sdata->buffer, sdata->buffersize);)
      CudaWrapper_FreeCudaData(sdata->buffer, sdata->buffersize);
      sdata->buffer = CudaWrapper_AllocCudaData(size);
      sdata->buffersize = size;
      sdata->buffer_new++;
      MYDBG(printf("New buffer at %p with %i kB\n", sdata->buffer, sdata->buffersize);)
    }

    cudaMemcpyToSymbol(MY_AP(buffer)           , & sdata->buffer                  , sizeof(int*));
    cudaMemcpyToSymbol(MY_AP(cuda_ntypes)      , & cuda_ntypes                    , sizeof(unsigned));
    cudaMemcpyToSymbol(MY_AP(cutneighsq)       , acutneighsq                    , nx);
    cudaMemcpyToSymbol(MY_AP(neighbor_maxlocal), & sneighlist->firstneigh.dim[0]  , sizeof(unsigned));
    cudaMemcpyToSymbol(MY_AP(firstneigh)       , & sneighlist->firstneigh.dev_data, sizeof(int*));
    cudaMemcpyToSymbol(MY_AP(ilist)            , & sneighlist->ilist     .dev_data, sizeof(int*));
    cudaMemcpyToSymbol(MY_AP(inum)             , & sneighlist->inum               , sizeof(int));
    cudaMemcpyToSymbol(MY_AP(nlocal)           , & sdata->atom.nlocal             , sizeof(int));
    cudaMemcpyToSymbol(MY_AP(nall)             , & sdata->atom.nall               , sizeof(int));
    cudaMemcpyToSymbol(MY_AP(nmax)             , & sdata->atom.nmax               , sizeof(int));
    cudaMemcpyToSymbol(MY_AP(numneigh)         , & sneighlist->numneigh  .dev_data, sizeof(int*));
    cudaMemcpyToSymbol(MY_AP(type)             , & sdata->atom.type      .dev_data, sizeof(int*));
    cudaMemcpyToSymbol(MY_AP(x)                , & sdata->atom.x         .dev_data, sizeof(X_FLOAT*));
    cudaMemcpyToSymbol(MY_AP(maxneighbors)     , & sneighlist->maxneighbors	 , sizeof(int));

    free(acutneighsq);
  }

  int3 layout = getgrid(sdata->atom.nlocal); // sneighlist->inum
  dim3 threads(layout.z, 1, 1);
  dim3 grid(layout.x, layout.y, 1);

  int return_value = 1;
  CudaWrapper_UploadCudaData(& return_value, sdata->buffer, sizeof(int));

  CUT_CHECK_ERROR("Cuda_NeighborBuild: pre neighbor build kernel execution failed");
  NeighborBuildFullNsq_Kernel <<< grid, threads>>> ();
  cudaThreadSynchronize();
  CUT_CHECK_ERROR("Cuda_NeighborBuild: neighbor build kernel execution failed");

  int buffer[20];
  CudaWrapper_DownloadCudaData(buffer, sdata->buffer, sizeof(int) * 20);
  MYDBG(printf("Cuda_NeighborBuildFullNSQ build neighbor list ... end\n");)
  return return_value = buffer[0];
}
