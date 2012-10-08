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
#define MY_PREFIX fix_viscous_cuda
#include "cuda_shared.h"
#include "cuda_common.h"
#include "crm_cuda_utils.cu"

#include "fix_viscous_cuda_cu.h"
#include "fix_viscous_cuda_kernel.cu"

void Cuda_FixViscousCuda_UpdateNmax(cuda_shared_data* sdata)
{
  cudaMemcpyToSymbol(MY_AP(mask)    , & sdata->atom.mask .dev_data, sizeof(int*));
  cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal        , sizeof(int));
  cudaMemcpyToSymbol(MY_AP(nmax)    , & sdata->atom.nmax          , sizeof(int));
  cudaMemcpyToSymbol(MY_AP(v)       , & sdata->atom.x    .dev_data, sizeof(X_FLOAT*));
  cudaMemcpyToSymbol(MY_AP(f)       , & sdata->atom.f    .dev_data, sizeof(F_FLOAT*));
  cudaMemcpyToSymbol(MY_AP(type)    , & sdata->atom.type .dev_data, sizeof(int*));
}

void Cuda_FixViscousCuda_Init(cuda_shared_data* sdata)
{
  Cuda_FixViscousCuda_UpdateNmax(sdata);

}


void Cuda_FixViscousCuda_PostForce(cuda_shared_data* sdata, int groupbit, void* gamma)
{
  if(sdata->atom.update_nmax)
    Cuda_FixViscousCuda_UpdateNmax(sdata);

  if(sdata->atom.update_nlocal)
    cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal        , sizeof(int));


  int3 layout = getgrid(sdata->atom.nlocal, 0);
  dim3 threads(layout.z, 1, 1);
  dim3 grid(layout.x, layout.y, 1);

  Cuda_FixViscousCuda_PostForce_Kernel <<< grid, threads, 0>>> (groupbit, (F_FLOAT*) gamma);
  cudaThreadSynchronize();
  CUT_CHECK_ERROR("Cuda_Cuda_FixViscousCuda_PostForce: Kernel execution failed");

}
