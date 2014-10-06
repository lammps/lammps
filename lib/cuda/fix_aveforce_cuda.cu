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
#define MY_PREFIX fix_ave_force_cuda
#include "cuda_shared.h"
#include "cuda_common.h"

#include "crm_cuda_utils.cu"

#include "fix_aveforce_cuda_cu.h"
#include "fix_aveforce_cuda_kernel.cu"

void Cuda_FixAveForceCuda_UpdateBuffer(cuda_shared_data* sdata)
{
  int3 layout = getgrid(sdata->atom.nlocal, 4 * sizeof(F_CFLOAT));
  dim3 threads(layout.z, 1, 1);
  dim3 grid(layout.x, layout.y, 1);
  int size = (unsigned)(layout.z * layout.y * layout.x) * 4 * sizeof(F_CFLOAT);

  if(sdata->buffersize < size) {
    MYDBG(printf("Cuda_FixAveForceCuda Resizing Buffer at %p with %i kB to\n", sdata->buffer, sdata->buffersize);)
    CudaWrapper_FreeCudaData(sdata->buffer, sdata->buffersize);
    sdata->buffer = CudaWrapper_AllocCudaData(size);
    sdata->buffersize = size;
    sdata->buffer_new++;
    MYDBG(printf("New buffer at %p with %i kB\n", sdata->buffer, sdata->buffersize);)
  }

  cudaMemcpyToSymbol(MY_AP(buffer), & sdata->buffer, sizeof(int*));
}

void Cuda_FixAveForceCuda_UpdateNmax(cuda_shared_data* sdata)
{
  cudaMemcpyToSymbol(MY_AP(mask)    , & sdata->atom.mask .dev_data, sizeof(int*));
  cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal        , sizeof(int));
  cudaMemcpyToSymbol(MY_AP(nmax)    , & sdata->atom.nmax          , sizeof(int));
  cudaMemcpyToSymbol(MY_AP(x)       , & sdata->atom.x    .dev_data, sizeof(X_CFLOAT*));
  cudaMemcpyToSymbol(MY_AP(f)       , & sdata->atom.f    .dev_data, sizeof(F_CFLOAT*));
}

void Cuda_FixAveForceCuda_Init(cuda_shared_data* sdata)
{
  Cuda_FixAveForceCuda_UpdateNmax(sdata);
}

void Cuda_FixAveForceCuda_PostForce_FOrg(cuda_shared_data* sdata, int groupbit, F_CFLOAT* aforiginal)
{
  if(sdata->atom.update_nmax)
    Cuda_FixAveForceCuda_UpdateNmax(sdata);

  if(sdata->atom.update_nlocal)
    cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal        , sizeof(int));

  if(sdata->buffer_new)
    Cuda_FixAveForceCuda_UpdateBuffer(sdata);

  int3 layout = getgrid(sdata->atom.nlocal, 4 * sizeof(F_CFLOAT));
  dim3 threads(layout.z, 1, 1);
  dim3 grid(layout.x, layout.y, 1);


  Cuda_FixAveForceCuda_PostForce_FOrg_Kernel <<< grid, threads, threads.x* 4* sizeof(F_CFLOAT)>>> (groupbit);
  cudaThreadSynchronize();
  CUT_CHECK_ERROR("Cuda_FixAveForceCuda_PostForce: fix ave_force post_force Compute Kernel execution failed");

  int oldgrid = grid.x;
  grid.x = 4;
  threads.x = 512;
  Cuda_FixAveForceCuda_reduce_foriginal <<< grid, threads, threads.x* sizeof(F_CFLOAT)>>> (oldgrid, aforiginal);
  cudaThreadSynchronize();
  CUT_CHECK_ERROR("Cuda_FixAveForceCuda_PostForce: fix ave_force post_force Reduce Kernel execution failed");

}

void Cuda_FixAveForceCuda_PostForce_Set(cuda_shared_data* sdata, int groupbit, int xflag, int yflag, int zflag, F_CFLOAT axvalue, F_CFLOAT ayvalue, F_CFLOAT azvalue)
{
  int3 layout = getgrid(sdata->atom.nlocal);
  dim3 threads(layout.z, 1, 1);
  dim3 grid(layout.x, layout.y, 1);


  Cuda_FixAveForceCuda_PostForce_Set_Kernel <<< grid, threads, 0>>> (groupbit, xflag, yflag, zflag, axvalue, ayvalue, azvalue);
  cudaThreadSynchronize();
  CUT_CHECK_ERROR("Cuda_FixAveForceCuda_PostForce_Set: fix ave_force post_force Compute Kernel execution failed");

}
