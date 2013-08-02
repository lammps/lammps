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
#define MY_PREFIX fix_nh_cuda
#define IncludeCommonNeigh
#include "cuda_shared.h"
#include "cuda_common.h"
#include "crm_cuda_utils.cu"
#include "fix_nh_cuda_cu.h"
#include "fix_nh_cuda_kernel.cu"

void Cuda_FixNHCuda_UpdateNmax(cuda_shared_data* sdata)
{
  cudaMemcpyToSymbol(MY_AP(f)       , & sdata->atom.f    .dev_data, sizeof(F_FLOAT*));
  cudaMemcpyToSymbol(MY_AP(mask)    , & sdata->atom.mask .dev_data, sizeof(int*));
  cudaMemcpyToSymbol(MY_AP(tag)     , & sdata->atom.tag  .dev_data, sizeof(int*));
  cudaMemcpyToSymbol(MY_AP(debugdata)     , & sdata->debugdata, sizeof(int*));
  cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal        , sizeof(int));
  cudaMemcpyToSymbol(MY_AP(nmax)    , & sdata->atom.nmax          , sizeof(int));
  cudaMemcpyToSymbol(MY_AP(rmass)   , & sdata->atom.rmass.dev_data, sizeof(V_FLOAT*));
  cudaMemcpyToSymbol(MY_AP(mass)    , & sdata->atom.mass.dev_data, sizeof(V_FLOAT*));
  cudaMemcpyToSymbol(MY_AP(type)    , & sdata->atom.type .dev_data, sizeof(int*));
  cudaMemcpyToSymbol(MY_AP(v)       , & sdata->atom.v    .dev_data, sizeof(V_FLOAT*));
  cudaMemcpyToSymbol(MY_AP(x)       , & sdata->atom.x    .dev_data, sizeof(X_FLOAT*));
  cudaMemcpyToSymbol(MY_AP(xhold)   , & sdata->atom.xhold.dev_data, sizeof(X_FLOAT*));  //might be moved to a neighbor record in sdata
  cudaMemcpyToSymbol(MY_AP(maxhold)   , & sdata->atom.maxhold, sizeof(int));  //might be moved to a neighbor record in sdata
  cudaMemcpyToSymbol(MY_AP(reneigh_flag), & sdata->buffer, sizeof(int*));  //might be moved to a neighbor record in sdata
  cudaMemcpyToSymbol(MY_AP(triggerneighsq), & sdata->atom.triggerneighsq, sizeof(X_FLOAT)); //might be moved to a neighbor record in sdata
}

void Cuda_FixNHCuda_UpdateBuffer(cuda_shared_data* sdata)
{
  int size = (unsigned)10 * sizeof(int);

  if(sdata->buffersize < size) {
    MYDBG(printf("Cuda_FixNHCuda Resizing Buffer at %p with %i kB to\n", sdata->buffer, sdata->buffersize);)
    CudaWrapper_FreeCudaData(sdata->buffer, sdata->buffersize);
    sdata->buffer = CudaWrapper_AllocCudaData(size);
    sdata->buffersize = size;
    sdata->buffer_new++;
    MYDBG(printf("New buffer at %p with %i kB\n", sdata->buffer, sdata->buffersize);)

  }

  cudaMemcpyToSymbol(MY_AP(buffer) , & sdata->buffer, sizeof(int*));
  cudaMemcpyToSymbol(MY_AP(reneigh_flag), & sdata->buffer, sizeof(int*));  //might be moved to a neighbor record in sdata
}

void Cuda_FixNHCuda_Init(cuda_shared_data* sdata, X_FLOAT dtv, V_FLOAT dtf)
{
  cudaMemcpyToSymbol(MY_AP(mass)    , & sdata->atom.mass.dev_data , sizeof(V_FLOAT*));
  cudaMemcpyToSymbol(MY_AP(dtf)     , & dtf                       		, sizeof(V_FLOAT));
  cudaMemcpyToSymbol(MY_AP(dtv)     , & dtv                            , sizeof(X_FLOAT));
  cudaMemcpyToSymbol(MY_AP(triggerneighsq), &sdata->atom.triggerneighsq, sizeof(X_FLOAT));
  cudaMemcpyToSymbol(MY_AP(dist_check), & sdata->atom.dist_check       , sizeof(int));
  cudaMemcpyToSymbol(MY_AP(rmass_flag), & sdata->atom.rmass_flag       , sizeof(int));       //
  Cuda_FixNHCuda_UpdateNmax(sdata);
}


void Cuda_FixNHCuda_nh_v_press(cuda_shared_data* sdata, int groupbit, double* factor_h, int mynlocal, int p_triclinic) //mynlocal can be nfirst if firstgroup==igroup  see cpp
{
  my_times atime1, atime2;
  my_gettime(CLOCK_REALTIME, &atime1);

  if(sdata->atom.update_nmax)
    Cuda_FixNHCuda_UpdateNmax(sdata);

  if(sdata->atom.update_nlocal)
    cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal , sizeof(int));

  my_gettime(CLOCK_REALTIME, &atime2);
  sdata->cuda_timings.test1 +=
    atime2.tv_sec - atime1.tv_sec + 1.0 * (atime2.tv_nsec - atime1.tv_nsec) / 1000000000;

  if(sdata->buffer_new)
    Cuda_FixNHCuda_UpdateBuffer(sdata);

  F_FLOAT3 factor = {factor_h[0], factor_h[1], factor_h[2]};
  F_FLOAT3 factor2;

  if(p_triclinic) {
    factor2.x = factor_h[3], factor2.y = factor_h[4];
    factor2.z = factor_h[5];
  }

  int3 layout = getgrid(mynlocal);
  dim3 threads(layout.z, 1, 1);
  dim3 grid(layout.x, layout.y, 1);
  FixNHCuda_nh_v_press_Kernel <<< grid, threads>>> (groupbit, factor, p_triclinic, factor2);
  cudaThreadSynchronize();
  CUT_CHECK_ERROR("FixNHCuda: fix nh v_press Kernel execution failed");

}

void Cuda_FixNHCuda_nh_v_press_and_nve_v_NoBias(cuda_shared_data* sdata, int groupbit, double* factor_h, int mynlocal, int p_triclinic) //mynlocal can be nfirst if firstgroup==igroup  see cpp
{
  if(sdata->atom.update_nmax)
    Cuda_FixNHCuda_UpdateNmax(sdata);

  if(sdata->atom.update_nlocal)
    cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal , sizeof(int));

  if(sdata->buffer_new)
    Cuda_FixNHCuda_UpdateBuffer(sdata);

  F_FLOAT3 factor = {factor_h[0], factor_h[1], factor_h[2]};
  F_FLOAT3 factor2;

  if(p_triclinic) {
    factor2.x = factor_h[3], factor2.y = factor_h[4];
    factor2.z = factor_h[5];
  }

  int3 layout = getgrid(mynlocal);
  dim3 threads(layout.z, 1, 1);
  dim3 grid(layout.x, layout.y, 1);
  CUT_CHECK_ERROR("FixNHCuda: fix nh v_press pre Kernel execution failed");
  FixNHCuda_nh_v_press_and_nve_v_NoBias_Kernel <<< grid, threads>>> (groupbit, factor, p_triclinic, factor2);
  cudaThreadSynchronize();
  CUT_CHECK_ERROR("FixNHCuda: fix nh v_press Kernel execution failed");

}

void Cuda_FixNHCuda_nh_v_temp(cuda_shared_data* sdata, int groupbit, F_FLOAT factor_eta, int mynlocal) //mynlocal can be nfirst if firstgroup==igroup  see cpp
{
  my_times atime1, atime2;
  my_gettime(CLOCK_REALTIME, &atime1);

  if(sdata->atom.update_nmax)
    Cuda_FixNHCuda_UpdateNmax(sdata);

  if(sdata->atom.update_nlocal)
    cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal , sizeof(int));

  my_gettime(CLOCK_REALTIME, &atime2);
  sdata->cuda_timings.test1 +=
    atime2.tv_sec - atime1.tv_sec + 1.0 * (atime2.tv_nsec - atime1.tv_nsec) / 1000000000;

  if(sdata->buffer_new)
    Cuda_FixNHCuda_UpdateBuffer(sdata);

  int3 layout = getgrid(mynlocal);
  dim3 threads(layout.z, 1, 1);
  dim3 grid(layout.x, layout.y, 1);
  FixNHCuda_nh_v_temp_Kernel <<< grid, threads>>> (groupbit, factor_eta);
  cudaThreadSynchronize();
  CUT_CHECK_ERROR("FixNHCuda: fix nh v_temp Kernel execution failed");

}
void Cuda_FixNHCuda_nve_v(cuda_shared_data* sdata, int groupbit, int mynlocal) //mynlocal can be nfirst if firstgroup==igroup  see cpp
{
  my_times atime1, atime2;
  my_gettime(CLOCK_REALTIME, &atime1);

  if(sdata->atom.update_nmax)
    Cuda_FixNHCuda_UpdateNmax(sdata);

  if(sdata->atom.update_nlocal)
    cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal , sizeof(int));

  my_gettime(CLOCK_REALTIME, &atime2);
  sdata->cuda_timings.test1 +=
    atime2.tv_sec - atime1.tv_sec + 1.0 * (atime2.tv_nsec - atime1.tv_nsec) / 1000000000;

  if(sdata->buffer_new)
    Cuda_FixNHCuda_UpdateBuffer(sdata);

  int3 layout = getgrid(mynlocal);
  dim3 threads(layout.z, 1, 1);
  dim3 grid(layout.x, layout.y, 1);
  FixNHCuda_nve_v_Kernel <<< grid, threads>>> (groupbit);
  cudaThreadSynchronize();
  CUT_CHECK_ERROR("FixNHCuda: nve_v Kernel execution failed");
}


void Cuda_FixNHCuda_nve_x(cuda_shared_data* sdata, int groupbit, int mynlocal) //mynlocal can be nfirst if firstgroup==igroup  see cpp
{
  my_times atime1, atime2;
  my_gettime(CLOCK_REALTIME, &atime1);

  if(sdata->atom.update_nmax)
    Cuda_FixNHCuda_UpdateNmax(sdata);

  if(sdata->atom.update_nlocal)
    cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal , sizeof(int));

  my_gettime(CLOCK_REALTIME, &atime2);
  sdata->cuda_timings.test1 +=
    atime2.tv_sec - atime1.tv_sec + 1.0 * (atime2.tv_nsec - atime1.tv_nsec) / 1000000000;

  if(sdata->buffer_new)
    Cuda_FixNHCuda_UpdateBuffer(sdata);

  int3 layout = getgrid(mynlocal);
  dim3 threads(layout.z, 1, 1);
  dim3 grid(layout.x, layout.y, 1);
  cudaMemset(sdata->buffer, 0, sizeof(int));
  FixNHCuda_nve_x_Kernel <<< grid, threads>>> (groupbit);
  cudaThreadSynchronize();
  int reneigh_flag;
  cudaMemcpy((void*)(&reneigh_flag), sdata->buffer, sizeof(int), cudaMemcpyDeviceToHost);
  sdata->atom.reneigh_flag += reneigh_flag;
  CUT_CHECK_ERROR("FixNHCuda: nve_x Kernel execution failed");
}

void Cuda_FixNHCuda_nve_v_and_nh_v_press_NoBias(cuda_shared_data* sdata, int groupbit, double* factor_h, int mynlocal, int p_triclinic) //mynlocal can be nfirst if firstgroup==igroup  see cpp
{
  if(sdata->atom.update_nmax)
    Cuda_FixNHCuda_UpdateNmax(sdata);

  if(sdata->atom.update_nlocal)
    cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal , sizeof(int));

  if(sdata->buffer_new)
    Cuda_FixNHCuda_UpdateBuffer(sdata);

  F_FLOAT3 factor = {factor_h[0], factor_h[1], factor_h[2]};
  F_FLOAT3 factor2;

  if(p_triclinic) {
    factor2.x = factor_h[3], factor2.y = factor_h[4];
    factor2.z = factor_h[5];
  }

  int3 layout = getgrid(mynlocal);
  dim3 threads(layout.z, 1, 1);
  dim3 grid(layout.x, layout.y, 1);
  FixNHCuda_nve_v_and_nh_v_press_NoBias_Kernel <<< grid, threads>>> (groupbit, factor, p_triclinic, factor2);
  cudaThreadSynchronize();
  CUT_CHECK_ERROR("FixNHCuda__nve_v_and_nh_v_press_NoBias:   Kernel execution failed");
}

