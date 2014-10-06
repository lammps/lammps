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
#define MY_PREFIX domain
#include "cuda_shared.h"
#include "cuda_common.h"

#include "crm_cuda_utils.cu"

#include "domain_cu.h"
#include "domain_kernel.cu"

void Cuda_Domain_UpdateBuffer(cuda_shared_data* sdata, int size)
{
  if(sdata->buffersize < size) {
    MYDBG(printf("Cuda_Domain Resizing Buffer at %p with %i kB to\n", sdata->buffer, sdata->buffersize);)
    CudaWrapper_FreeCudaData(sdata->buffer, sdata->buffersize);
    sdata->buffer = CudaWrapper_AllocCudaData(size);
    sdata->buffersize = size;
    sdata->buffer_new++;
    MYDBG(printf("New buffer at %p with %i kB\n", sdata->buffer, sdata->buffersize);)
  }

  cudaMemcpyToSymbol(MY_AP(buffer), & sdata->buffer, sizeof(int*));
}

void Cuda_Domain_UpdateNmax(cuda_shared_data* sdata)
{
  cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal        , sizeof(int));
  cudaMemcpyToSymbol(MY_AP(nmax)    , & sdata->atom.nmax          , sizeof(int));
  cudaMemcpyToSymbol(MY_AP(x)       , & sdata->atom.x    .dev_data, sizeof(X_CFLOAT*));
  cudaMemcpyToSymbol(MY_AP(v)       , & sdata->atom.v    .dev_data, sizeof(V_CFLOAT*));
  cudaMemcpyToSymbol(MY_AP(mask)    , & sdata->atom.mask .dev_data, sizeof(int*));
  cudaMemcpyToSymbol(MY_AP(tag)    , & sdata->atom.tag .dev_data, sizeof(int*));
  cudaMemcpyToSymbol(MY_AP(image)   , & sdata->atom.image.dev_data, sizeof(int*));
}

void Cuda_Domain_UpdateDomain(cuda_shared_data* sdata)
{
  cudaMemcpyToSymbol(MY_AP(boxlo)   ,  sdata->domain.boxlo       , 3 * sizeof(X_CFLOAT));
  cudaMemcpyToSymbol(MY_AP(boxhi)   ,  sdata->domain.boxhi       , 3 * sizeof(X_CFLOAT));
  cudaMemcpyToSymbol(MY_AP(sublo)   ,  sdata->domain.sublo       , 3 * sizeof(X_CFLOAT));
  cudaMemcpyToSymbol(MY_AP(subhi)   ,  sdata->domain.subhi       , 3 * sizeof(X_CFLOAT));
  cudaMemcpyToSymbol(MY_AP(prd)     ,  sdata->domain.prd         , 3 * sizeof(X_CFLOAT));
  cudaMemcpyToSymbol(MY_AP(periodicity)   ,   sdata->domain.periodicity , 3 * sizeof(int));
  cudaMemcpyToSymbol(MY_AP(triclinic)     , & sdata->domain.triclinic   , sizeof(int));
  cudaMemcpyToSymbol(MY_AP(boxlo_lamda)   ,   sdata->domain.boxlo_lamda , 3 * sizeof(X_CFLOAT));
  cudaMemcpyToSymbol(MY_AP(boxhi_lamda)   ,   sdata->domain.boxhi_lamda , 3 * sizeof(X_CFLOAT));
  cudaMemcpyToSymbol(MY_AP(prd_lamda)	   ,   sdata->domain.prd_lamda   , 3 * sizeof(X_CFLOAT));
  cudaMemcpyToSymbol(MY_AP(h)	   	 ,   sdata->domain.h   		  , 6 * sizeof(X_CFLOAT));
  cudaMemcpyToSymbol(MY_AP(h_inv)	 ,   sdata->domain.h_inv   	  , 6 * sizeof(X_CFLOAT));
  cudaMemcpyToSymbol(MY_AP(h_rate)	 ,   sdata->domain.h_rate     , 6 * sizeof(V_CFLOAT));
  cudaMemcpyToSymbol(MY_AP(flag)	 ,   &sdata->flag     , sizeof(int*));
  cudaMemcpyToSymbol(MY_AP(debugdata)	 ,   &sdata->debugdata     , sizeof(int*));
}

void Cuda_Domain_Init(cuda_shared_data* sdata)
{
  Cuda_Domain_UpdateNmax(sdata);
  Cuda_Domain_UpdateDomain(sdata);
}

void Cuda_Domain_PBC(cuda_shared_data* sdata, int deform_remap, int deform_groupbit, double* extent)
{
  Cuda_Domain_UpdateNmax(sdata);
  //if(sdata->domain.update)
  Cuda_Domain_UpdateDomain(sdata);
  cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal        , sizeof(int));

  int box_change = 0;

  if(extent) box_change = 1;

  int sharedmem = 0;

  if(box_change) sharedmem = 6 * sizeof(X_CFLOAT);

  int3 layout = getgrid(sdata->atom.nlocal, sharedmem);
  dim3 threads(layout.z, 1, 1);
  dim3 grid(layout.x, layout.y, 1);
  sharedmem *= threads.x;

  if((box_change) && (sdata->buffer_new or (6 * sizeof(X_CFLOAT)*grid.x * grid.y > sdata->buffersize)))
    Cuda_Domain_UpdateBuffer(sdata, layout.x * layout.y * 6 * sizeof(X_CFLOAT));


  Domain_PBC_Kernel <<< grid, threads, sharedmem>>>(deform_remap, deform_groupbit, box_change);
  cudaThreadSynchronize();

  CUT_CHECK_ERROR("Cuda_Domain_PBC: Kernel execution failed");

  if(box_change) {
    X_CFLOAT buf2[6 * layout.x * layout.y];
    X_CFLOAT* buf = buf2;
    int flag;
    cudaMemcpy(buf, sdata->buffer, 6 * layout.x * layout.y * sizeof(X_CFLOAT), cudaMemcpyDeviceToHost);
    cudaMemcpy(&flag, sdata->flag, sizeof(int), cudaMemcpyDeviceToHost);
    //printf("Flag: %i\n",flag);
    X_CFLOAT min, max;
    min = 1.0 * BIG;
    max = -1.0 * BIG;

    for(int i = 0; i < layout.x * layout.y; i++) {
      if(buf[i] < min) min = buf[i];

      if(buf[i + layout.x * layout.y] > max) max = buf[i + layout.x * layout.y];
    }

    extent[0] = min;
    extent[1] = max;

    buf += 2 * layout.x * layout.y;
    min = 1.0 * BIG;
    max = -1.0 * BIG;

    for(int i = 0; i < layout.x * layout.y; i++) {
      if(buf[i] < min) min = buf[i];

      if(buf[i + layout.x * layout.y] > max) max = buf[i + layout.x * layout.y];
    }

    extent[2] = min;
    extent[3] = max;

    buf += 2 * layout.x * layout.y;
    min = 1.0 * BIG;
    max = -1.0 * BIG;

    for(int i = 0; i < layout.x * layout.y; i++) {
      if(buf[i] < min) min = buf[i];

      if(buf[i + layout.x * layout.y] > max) max = buf[i + layout.x * layout.y];
    }

    extent[4] = min;
    extent[5] = max;
    //printf("Extent: %lf %lf %lf %lf %lf %lf\n",extent[0],extent[1],extent[2],extent[3],extent[4],extent[5]);
    /*	   int n=grid.x*grid.y;
    	   if(n<128) threads.x=32;
    	   else if(n<256) threads.x=64;
    	   else threads.x=128;
    	   sharedmem=n*sizeof(X_CFLOAT);
    	   grid.x=6;
    	   grid.y=1;
    	   Domain_reduceBoxExtent<<<grid, threads,sharedmem>>>(extent,n);
    	   cudaThreadSynchronize();
    	   CUT_CHECK_ERROR("Cuda_Domain_reduceBoxExtent: Kernel execution failed");*/
  }
}

void Cuda_Domain_lamda2x(cuda_shared_data* sdata, int n)
{
  Cuda_Domain_UpdateNmax(sdata);
  //if(sdata->domain.update)
  Cuda_Domain_UpdateDomain(sdata);
  cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal        , sizeof(int));

  int3 layout = getgrid(n);
  dim3 threads(layout.z, 1, 1);
  dim3 grid(layout.x, layout.y, 1);

  Domain_lamda2x_Kernel <<< grid, threads, 0>>>(n);
  cudaThreadSynchronize();
  CUT_CHECK_ERROR("Cuda_Domain_lamda2x: Kernel execution failed");
}

void Cuda_Domain_x2lamda(cuda_shared_data* sdata, int n)
{
  Cuda_Domain_UpdateNmax(sdata);
  //if(sdata->domain.update)
  Cuda_Domain_UpdateDomain(sdata);
  cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal        , sizeof(int));

  int3 layout = getgrid(n);
  dim3 threads(layout.z, 1, 1);
  dim3 grid(layout.x, layout.y, 1);

  Domain_x2lamda_Kernel <<< grid, threads, 0>>>(n);
  cudaThreadSynchronize();
  CUT_CHECK_ERROR("Cuda_Domain_x2lamda: Kernel execution failed");
}
