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
#define MY_PREFIX fix_temp_berendsen_cuda
#include "cuda_shared.h"
#include "cuda_common.h"
#include "crm_cuda_utils.cu"

#include "fix_temp_berendsen_cuda_cu.h"
#include "fix_temp_berendsen_cuda_kernel.cu"


void Cuda_FixTempBerendsenCuda_UpdateNmax(cuda_shared_data* sdata)
{
  cudaMemcpyToSymbol(MY_AP(mask)    , & sdata->atom.mask .dev_data, sizeof(int*));
  cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal        , sizeof(int));
  cudaMemcpyToSymbol(MY_AP(nmax)    , & sdata->atom.nmax          , sizeof(int));
  cudaMemcpyToSymbol(MY_AP(v)       , & sdata->atom.v    .dev_data, sizeof(X_CFLOAT*));
}

void Cuda_FixTempBerendsenCuda_Init(cuda_shared_data* sdata)
{
  Cuda_FixTempBerendsenCuda_UpdateNmax(sdata);

}


void Cuda_FixTempBerendsenCuda_EndOfStep(cuda_shared_data* sdata, int groupbit, double afactor)
{
  V_CFLOAT factor = afactor;

  if(sdata->atom.update_nmax)
    Cuda_FixTempBerendsenCuda_UpdateNmax(sdata);

  if(sdata->atom.update_nlocal)
    cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal        , sizeof(int));

  int3 layout = getgrid(sdata->atom.nlocal);
  dim3 threads(layout.z, 1, 1);
  dim3 grid(layout.x, layout.y, 1);

  Cuda_FixTempBerendsenCuda_EndOfStep_Kernel <<< grid, threads, 0>>> (groupbit, factor);
  cudaThreadSynchronize();
  CUT_CHECK_ERROR("Cuda_Cuda_FixTempBerendsenCuda_PostForce: fix add_force post_force compute Kernel execution failed");
}
