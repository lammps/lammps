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
#define MY_PREFIX fix_gravity_cuda
#include "cuda_shared.h"
#include "cuda_common.h"
#include "crm_cuda_utils.cu"

#include "fix_gravity_cuda_cu.h"
#include "fix_gravity_cuda_kernel.cu"

void Cuda_FixGravityCuda_UpdateBuffer(cuda_shared_data* sdata)
{
  int3 layout = getgrid(sdata->atom.nlocal, 3 * sizeof(F_CFLOAT));
  dim3 threads(layout.z, 1, 1);
  dim3 grid(layout.x, layout.y, 1);
  int size = (unsigned)(layout.z * layout.y * layout.x) * 3 * sizeof(F_CFLOAT);

  if(sdata->buffersize < size) {
    MYDBG(printf("Cuda_FixGravityCuda Resizing Buffer at %p with %i kB to\n", sdata->buffer, sdata->buffersize);)
    CudaWrapper_FreeCudaData(sdata->buffer, sdata->buffersize);
    sdata->buffer = CudaWrapper_AllocCudaData(size);
    sdata->buffersize = size;
    sdata->buffer_new++;
    MYDBG(printf("New buffer at %p with %i kB\n", sdata->buffer, sdata->buffersize);)

  }

  cudaMemcpyToSymbol(MY_AP(buffer) , & sdata->buffer, sizeof(int*));
}

void Cuda_FixGravityCuda_UpdateNmax(cuda_shared_data* sdata)
{
  cudaMemcpyToSymbol(MY_AP(mask)    , & sdata->atom.mask .dev_data, sizeof(int*));
  cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal        , sizeof(int));
  cudaMemcpyToSymbol(MY_AP(nmax)    , & sdata->atom.nmax          , sizeof(int));
  cudaMemcpyToSymbol(MY_AP(x)       , & sdata->atom.x    .dev_data, sizeof(X_CFLOAT*));
  cudaMemcpyToSymbol(MY_AP(f)       , & sdata->atom.f    .dev_data, sizeof(F_CFLOAT*));
  cudaMemcpyToSymbol(MY_AP(type)       , & sdata->atom.type    .dev_data, sizeof(int*));
  cudaMemcpyToSymbol(MY_AP(rmass_flag)       , & sdata->atom.rmass_flag, sizeof(int));
  cudaMemcpyToSymbol(MY_AP(rmass)       , & sdata->atom.rmass    .dev_data, sizeof(V_CFLOAT*));
  cudaMemcpyToSymbol(MY_AP(mass)       , & sdata->atom.mass    .dev_data, sizeof(V_CFLOAT*));
}

void Cuda_FixGravityCuda_Init(cuda_shared_data* sdata)
{
  Cuda_FixGravityCuda_UpdateNmax(sdata);

}


void Cuda_FixGravityCuda_PostForce(cuda_shared_data* sdata, int groupbit, F_CFLOAT xacc, F_CFLOAT yacc, F_CFLOAT zacc)
{
  if(sdata->atom.update_nmax)
    Cuda_FixGravityCuda_UpdateNmax(sdata);

  if(sdata->atom.update_nlocal)
    cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal        , sizeof(int));

  if(sdata->buffer_new)
    Cuda_FixGravityCuda_UpdateBuffer(sdata);


  int3 layout = getgrid(sdata->atom.nlocal);
  dim3 threads(layout.z, 1, 1);
  dim3 grid(layout.x, layout.y, 1);

  Cuda_FixGravityCuda_PostForce_Kernel <<< grid, threads>>> (groupbit, xacc, yacc, zacc);
  cudaThreadSynchronize();
  CUT_CHECK_ERROR("Cuda_Cuda_FixGravityCuda_PostForce: fix add_force post_force compute Kernel execution failed");
}
