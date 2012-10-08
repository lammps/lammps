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
#define MY_PREFIX compute_temp_cuda
#include "cuda_shared.h"
#include "cuda_common.h"

#include "crm_cuda_utils.cu"

#include "compute_temp_cuda_cu.h"
#include "compute_temp_cuda_kernel.cu"

void Cuda_ComputeTempCuda_UpdateBuffer(cuda_shared_data* sdata)
{
  int size = (unsigned)((sdata->atom.nlocal + 63) / 64.0) * 6 * sizeof(ENERGY_FLOAT);

  if(sdata->buffersize < size) {
    MYDBG(printf("Cuda_ComputeTempCuda Resizing Buffer at %p with %i kB to\n", sdata->buffer, sdata->buffersize);)
    CudaWrapper_FreeCudaData(sdata->buffer, sdata->buffersize);
    sdata->buffer = CudaWrapper_AllocCudaData(size);
    sdata->buffersize = size;
    sdata->buffer_new++;
    MYDBG(printf("New buffer at %p with %i kB\n", sdata->buffer, sdata->buffersize);)
  }

  cudaMemcpyToSymbol(MY_AP(buffer), & sdata->buffer, sizeof(int*));
}

void Cuda_ComputeTempCuda_UpdateNmax(cuda_shared_data* sdata)
{
  cudaMemcpyToSymbol(MY_AP(mask)    , & sdata->atom.mask .dev_data, sizeof(int*));
  cudaMemcpyToSymbol(MY_AP(mass)    , & sdata->atom.mass .dev_data, sizeof(V_FLOAT*));

  if(sdata->atom.rmass_flag)
    cudaMemcpyToSymbol(MY_AP(rmass)   , & sdata->atom.rmass.dev_data, sizeof(V_FLOAT*));

  cudaMemcpyToSymbol(MY_AP(rmass_flag)   , & sdata->atom.rmass_flag, sizeof(int));
  cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal        , sizeof(int));
  cudaMemcpyToSymbol(MY_AP(nmax)    , & sdata->atom.nmax          , sizeof(int));
  cudaMemcpyToSymbol(MY_AP(v)       , & sdata->atom.v    .dev_data, sizeof(V_FLOAT*));
  cudaMemcpyToSymbol(MY_AP(type)       , & sdata->atom.type    .dev_data, sizeof(int*));
}

void Cuda_ComputeTempCuda_Init(cuda_shared_data* sdata)
{
  Cuda_ComputeTempCuda_UpdateNmax(sdata);
}


void Cuda_ComputeTempCuda_Vector(cuda_shared_data* sdata, int groupbit, ENERGY_FLOAT* t)
{
  //if(sdata->atom.update_nmax) //is most likely not called every timestep, therefore update of constants is necessary
  Cuda_ComputeTempCuda_UpdateNmax(sdata);
  //if(sdata->atom.update_nlocal)
  cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal        , sizeof(int));
  //if(sdata->buffer_new)
  Cuda_ComputeTempCuda_UpdateBuffer(sdata);

  int3 layout = getgrid(sdata->atom.nlocal);
  dim3 threads(layout.z, 1, 1);
  dim3 grid(layout.x, layout.y, 1);

  if(sdata->atom.nlocal > 0) {
    Cuda_ComputeTempCuda_Vector_Kernel <<< grid, threads, threads.x* 6* sizeof(ENERGY_FLOAT)>>> (groupbit);
    cudaThreadSynchronize();
    CUT_CHECK_ERROR("Cuda_ComputeTempCuda_Vector: compute_vector Kernel execution failed");

    int oldgrid = grid.x * grid.y;
    grid.x = 6;
    grid.y = 1;
    threads.x = 512;
    Cuda_ComputeTempCuda_Reduce_Kernel <<< grid, threads, threads.x* sizeof(ENERGY_FLOAT)>>> (oldgrid, t);
    cudaThreadSynchronize();
    CUT_CHECK_ERROR("Cuda_ComputeTempCuda_Vector: reduce_vector Kernel execution failed");
  }
}

void Cuda_ComputeTempCuda_Scalar(cuda_shared_data* sdata, int groupbit, ENERGY_FLOAT* t)
{
  //if(sdata->atom.update_nmax) //is most likely not called every timestep, therefore update of constants is necessary
  Cuda_ComputeTempCuda_UpdateNmax(sdata);
  //if(sdata->atom.update_nlocal)
  cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal        , sizeof(int));
  //if(sdata->buffer_new)
  Cuda_ComputeTempCuda_UpdateBuffer(sdata);
  MYDBG(printf("#CUDA ComputeTempCuda_Scalar: %i\n", sdata->atom.nlocal);)
  int3 layout = getgrid(sdata->atom.nlocal);
  dim3 threads(layout.z, 1, 1);
  dim3 grid(layout.x, layout.y, 1);

  if(sdata->atom.nlocal > 0) {
    CUT_CHECK_ERROR("Cuda_ComputeTempCuda_Scalar: pre compute_scalar Kernel");
    Cuda_ComputeTempCuda_Scalar_Kernel <<< grid, threads, threads.x* sizeof(ENERGY_FLOAT)>>> (groupbit);
    cudaThreadSynchronize();
    CUT_CHECK_ERROR("Cuda_ComputeTempCuda_Scalar: compute_scalar Kernel execution failed");

    int oldgrid = grid.x * grid.y;
    grid.x = 1;
    grid.y = 1;
    threads.x = 512;
    Cuda_ComputeTempCuda_Reduce_Kernel <<< grid, threads, threads.x* sizeof(ENERGY_FLOAT)>>> (oldgrid, t);
    cudaThreadSynchronize();
    CUT_CHECK_ERROR("Cuda_ComputeTempCuda_Scalar: reduce_scalar Kernel execution failed");
  }
}
