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
#define MY_PREFIX comm_cuda
#include "cuda_shared.h"
#include "cuda_common.h"

#include "crm_cuda_utils.cu"

#include "comm_cuda_cu.h"
#include "comm_cuda_kernel.cu"
#include <ctime>

void Cuda_CommCuda_UpdateBuffer(cuda_shared_data* sdata, int n)
{
  int size = n * 3 * sizeof(X_CFLOAT);

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


void Cuda_CommCuda_UpdateNmax(cuda_shared_data* sdata)
{
  cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal        , sizeof(int));
  cudaMemcpyToSymbol(MY_AP(nmax)    , & sdata->atom.nmax          , sizeof(int));
  cudaMemcpyToSymbol(MY_AP(x)       , & sdata->atom.x    .dev_data, sizeof(X_CFLOAT*));
  cudaMemcpyToSymbol(MY_AP(v)       , & sdata->atom.v    .dev_data, sizeof(X_CFLOAT*));
  cudaMemcpyToSymbol(MY_AP(f)       , & sdata->atom.f    .dev_data, sizeof(F_CFLOAT*));
  cudaMemcpyToSymbol(MY_AP(type)    , & sdata->atom.type .dev_data, sizeof(int*));
}


void Cuda_CommCuda_Init(cuda_shared_data* sdata)
{
  Cuda_CommCuda_UpdateNmax(sdata);
  int ntypesp = sdata->atom.ntypes + 1;
  cudaMemcpyToSymbol(MY_AP(cuda_ntypes)   , &ntypesp, sizeof(int));
  cudaMemcpyToSymbol(MY_AP(prd)   , sdata->domain.prd, 3 * sizeof(X_CFLOAT));
  cudaMemcpyToSymbol(MY_AP(flag)  , &sdata->flag, sizeof(int*));
  cudaMemcpyToSymbol(MY_AP(debugdata)  , &sdata->debugdata, sizeof(int*));
}

int Cuda_CommCuda_PackComm(cuda_shared_data* sdata, int n, int iswap, void* buf_send, int* pbc, int pbc_flag)
{

  my_times time1, time2;

  if(sdata->atom.update_nmax)
    Cuda_CommCuda_UpdateNmax(sdata);

  if(sdata->atom.update_nlocal)
    cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal        , sizeof(int));

  int size = n * 3 * sizeof(X_CFLOAT);

  if(sdata->buffer_new or (size > sdata->buffersize))
    Cuda_CommCuda_UpdateBuffer(sdata, n);

  X_CFLOAT dx = 0.0;
  X_CFLOAT dy = 0.0;
  X_CFLOAT dz = 0.0;

  if(pbc_flag != 0) {
    if(sdata->domain.triclinic == 0) {
      dx = pbc[0] * sdata->domain.prd[0];
      dy = pbc[1] * sdata->domain.prd[1];
      dz = pbc[2] * sdata->domain.prd[2];
    } else {
      dx = pbc[0] * sdata->domain.prd[0] + pbc[5] * sdata->domain.xy + pbc[4] * sdata->domain.xz;
      dy = pbc[1] * sdata->domain.prd[1] + pbc[3] * sdata->domain.yz;
      dz = pbc[2] * sdata->domain.prd[2];
    }
  }

  int3 layout = getgrid(n);
  dim3 threads(layout.z, 1, 1);
  dim3 grid(layout.x, layout.y, 1);

  if(sdata->atom.nlocal > 0) {
    cudaMemset(sdata->flag, 0, sizeof(int));

    my_gettime(CLOCK_REALTIME, &time1);

    void* buf = sdata->overlap_comm ? sdata->comm.buf_send_dev[iswap] : sdata->buffer;
    Cuda_CommCuda_PackComm_Kernel <<< grid, threads, 0>>>((int*) sdata->comm.sendlist.dev_data, n
        , sdata->comm.maxlistlength, iswap, dx, dy, dz, buf);
    cudaThreadSynchronize();

    my_gettime(CLOCK_REALTIME, &time2);
    sdata->cuda_timings.comm_forward_kernel_pack +=
      time2.tv_sec - time1.tv_sec + 1.0 * (time2.tv_nsec - time1.tv_nsec) / 1000000000;

    CUT_CHECK_ERROR("Cuda_CommCuda_PackComm: Kernel execution failed");

    if(not sdata->overlap_comm)
      cudaMemcpy(buf_send, sdata->buffer, n * 3 * sizeof(X_CFLOAT), cudaMemcpyDeviceToHost);
    //cudaMemcpy(buf_send, sdata->comm.buf_send_dev[iswap], n*3*sizeof(X_CFLOAT), cudaMemcpyDeviceToHost);

    my_gettime(CLOCK_REALTIME, &time1);
    sdata->cuda_timings.comm_forward_download +=
      time1.tv_sec - time2.tv_sec + 1.0 * (time1.tv_nsec - time2.tv_nsec) / 1000000000;

    int aflag;
    cudaMemcpy(&aflag, sdata->flag, sizeof(int), cudaMemcpyDeviceToHost);
    if(aflag != 0) printf("aflag PackComm: %i\n", aflag);
    CUT_CHECK_ERROR("Cuda_CommCuda_PackComm: Kernel execution failed");

  }

  return 3 * n;
}

int Cuda_CommCuda_PackCommVel(cuda_shared_data* sdata, int n, int iswap, void* buf_send, int* pbc, int pbc_flag)
{

  my_times time1, time2;

  if(sdata->atom.update_nmax)
    Cuda_CommCuda_UpdateNmax(sdata);

  if(sdata->atom.update_nlocal)
    cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal        , sizeof(int));

  int size = n * 6 * sizeof(X_CFLOAT);

  if(sdata->buffer_new or (size > sdata->buffersize))
    Cuda_CommCuda_UpdateBuffer(sdata, n);

  X_CFLOAT dx = 0.0;
  X_CFLOAT dy = 0.0;
  X_CFLOAT dz = 0.0;

  if(pbc_flag != 0) {
    if(sdata->domain.triclinic == 0) {
      dx = pbc[0] * sdata->domain.prd[0];
      dy = pbc[1] * sdata->domain.prd[1];
      dz = pbc[2] * sdata->domain.prd[2];
    } else {
      dx = pbc[0] * sdata->domain.prd[0] + pbc[5] * sdata->domain.xy + pbc[4] * sdata->domain.xz;
      dy = pbc[1] * sdata->domain.prd[1] + pbc[3] * sdata->domain.yz;
      dz = pbc[2] * sdata->domain.prd[2];
    }
  }

  int3 layout = getgrid(n);
  dim3 threads(layout.z, 1, 1);
  dim3 grid(layout.x, layout.y, 1);

  if(sdata->atom.nlocal > 0) {
    cudaMemset(sdata->flag, 0, sizeof(int));

    my_gettime(CLOCK_REALTIME, &time1);

    void* buf = sdata->overlap_comm ? sdata->comm.buf_send_dev[iswap] : sdata->buffer;
    Cuda_CommCuda_PackComm_Kernel <<< grid, threads, 0>>>((int*) sdata->comm.sendlist.dev_data, n
        , sdata->comm.maxlistlength, iswap, dx, dy, dz, buf);
    cudaThreadSynchronize();

    my_gettime(CLOCK_REALTIME, &time2);
    sdata->cuda_timings.comm_forward_kernel_pack +=
      time2.tv_sec - time1.tv_sec + 1.0 * (time2.tv_nsec - time1.tv_nsec) / 1000000000;

    CUT_CHECK_ERROR("Cuda_CommCuda_PackComm: Kernel execution failed");

    if(not sdata->overlap_comm)
      cudaMemcpy(buf_send, sdata->buffer, n * 6 * sizeof(X_CFLOAT), cudaMemcpyDeviceToHost);
    //cudaMemcpy(buf_send, sdata->comm.buf_send_dev[iswap], n*3*sizeof(X_CFLOAT), cudaMemcpyDeviceToHost);

    my_gettime(CLOCK_REALTIME, &time1);
    sdata->cuda_timings.comm_forward_download +=
      time1.tv_sec - time2.tv_sec + 1.0 * (time1.tv_nsec - time2.tv_nsec) / 1000000000;

    int aflag;
    cudaMemcpy(&aflag, sdata->flag, sizeof(int), cudaMemcpyDeviceToHost);
    if(aflag != 0) printf("aflag PackComm: %i\n", aflag);
    CUT_CHECK_ERROR("Cuda_CommCuda_PackComm: Kernel execution failed");

  }

  return 6 * n;
}

int Cuda_CommCuda_PackComm_Self(cuda_shared_data* sdata, int n, int iswap, int first, int* pbc, int pbc_flag)
{
  MYDBG(printf(" # CUDA: CommCuda_PackComm_Self\n");)
  my_times time1, time2;

  if(sdata->atom.update_nmax)
    Cuda_CommCuda_UpdateNmax(sdata);

  if(sdata->atom.update_nlocal)
    cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal        , sizeof(int));

  int size = n * 3 * sizeof(X_CFLOAT);

  if(sdata->buffer_new or (size > sdata->buffersize))
    Cuda_CommCuda_UpdateBuffer(sdata, n);

  static int count = -1;
  count++;
  X_CFLOAT dx = 0.0;
  X_CFLOAT dy = 0.0;
  X_CFLOAT dz = 0.0;

  if(pbc_flag != 0) {
    if(sdata->domain.triclinic == 0) {
      dx = pbc[0] * sdata->domain.prd[0];
      dy = pbc[1] * sdata->domain.prd[1];
      dz = pbc[2] * sdata->domain.prd[2];
    } else {
      dx = pbc[0] * sdata->domain.prd[0] + pbc[5] * sdata->domain.xy + pbc[4] * sdata->domain.xz;
      dy = pbc[1] * sdata->domain.prd[1] + pbc[3] * sdata->domain.yz;
      dz = pbc[2] * sdata->domain.prd[2];
    }
  }



  int3 layout = getgrid(n);
  dim3 threads(layout.z, 1, 1);
  dim3 grid(layout.x, layout.y, 1);

  if(sdata->atom.nlocal > 0) {

    my_gettime(CLOCK_REALTIME, &time1);

    Cuda_CommCuda_PackComm_Self_Kernel <<< grid, threads, 0>>>((int*) sdata->comm.sendlist.dev_data, n, sdata->comm.maxlistlength, iswap, dx, dy, dz, first);
    cudaThreadSynchronize();

    my_gettime(CLOCK_REALTIME, &time2);
    sdata->cuda_timings.comm_forward_kernel_self +=
      time2.tv_sec - time1.tv_sec + 1.0 * (time2.tv_nsec - time1.tv_nsec) / 1000000000;

    CUT_CHECK_ERROR("Cuda_CommCuda_PackComm_Self: Kernel execution failed");
  }

  return 3 * n;
}

int Cuda_CommCuda_PackCommVel_Self(cuda_shared_data* sdata, int n, int iswap, int first, int* pbc, int pbc_flag)
{
  MYDBG(printf(" # CUDA: CommCuda_PackComm_Self\n");)
  my_times time1, time2;

  if(sdata->atom.update_nmax)
    Cuda_CommCuda_UpdateNmax(sdata);

  if(sdata->atom.update_nlocal)
    cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal        , sizeof(int));

  int size = n * 6 * sizeof(X_CFLOAT);

  if(sdata->buffer_new or (size > sdata->buffersize))
    Cuda_CommCuda_UpdateBuffer(sdata, n);

  static int count = -1;
  count++;
  X_CFLOAT dx = 0.0;
  X_CFLOAT dy = 0.0;
  X_CFLOAT dz = 0.0;

  if(pbc_flag != 0) {
    if(sdata->domain.triclinic == 0) {
      dx = pbc[0] * sdata->domain.prd[0];
      dy = pbc[1] * sdata->domain.prd[1];
      dz = pbc[2] * sdata->domain.prd[2];
    } else {
      dx = pbc[0] * sdata->domain.prd[0] + pbc[5] * sdata->domain.xy + pbc[4] * sdata->domain.xz;
      dy = pbc[1] * sdata->domain.prd[1] + pbc[3] * sdata->domain.yz;
      dz = pbc[2] * sdata->domain.prd[2];
    }
  }



  int3 layout = getgrid(n);
  dim3 threads(layout.z, 1, 1);
  dim3 grid(layout.x, layout.y, 1);

  if(sdata->atom.nlocal > 0) {

    my_gettime(CLOCK_REALTIME, &time1);

    Cuda_CommCuda_PackComm_Self_Kernel <<< grid, threads, 0>>>((int*) sdata->comm.sendlist.dev_data, n, sdata->comm.maxlistlength, iswap, dx, dy, dz, first);
    cudaThreadSynchronize();

    my_gettime(CLOCK_REALTIME, &time2);
    sdata->cuda_timings.comm_forward_kernel_self +=
      time2.tv_sec - time1.tv_sec + 1.0 * (time2.tv_nsec - time1.tv_nsec) / 1000000000;

    CUT_CHECK_ERROR("Cuda_CommCuda_PackComm_Self: Kernel execution failed");
  }

  return 6 * n;
}

void Cuda_CommCuda_UnpackComm(cuda_shared_data* sdata, int n, int first, void* buf_recv, int iswap)
{
  my_times time1, time2;

  if(sdata->atom.update_nmax)
    Cuda_CommCuda_UpdateNmax(sdata);

  if(sdata->atom.update_nlocal)
    cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal        , sizeof(int));

  int size = n * 3 * sizeof(X_CFLOAT);

  if(sdata->buffer_new or (size > sdata->buffersize))
    Cuda_CommCuda_UpdateBuffer(sdata, n);

  int3 layout = getgrid(n);
  dim3 threads(layout.z, 1, 1);
  dim3 grid(layout.x, layout.y, 1);

  if(sdata->atom.nlocal > 0) {
    my_gettime(CLOCK_REALTIME, &time1);

    if(not sdata->overlap_comm || iswap < 0)
      cudaMemcpy(sdata->buffer, (void*)buf_recv, n * 3 * sizeof(X_CFLOAT), cudaMemcpyHostToDevice);

    my_gettime(CLOCK_REALTIME, &time2);
    sdata->cuda_timings.comm_forward_upload +=
      time2.tv_sec - time1.tv_sec + 1.0 * (time2.tv_nsec - time1.tv_nsec) / 1000000000;
    void* buf = (sdata->overlap_comm && iswap >= 0) ? sdata->comm.buf_recv_dev[iswap] : sdata->buffer;
    Cuda_CommCuda_UnpackComm_Kernel <<< grid, threads, 0>>>(n, first, buf);
    cudaThreadSynchronize();

    my_gettime(CLOCK_REALTIME, &time1);
    sdata->cuda_timings.comm_forward_kernel_unpack +=
      time1.tv_sec - time2.tv_sec + 1.0 * (time1.tv_nsec - time2.tv_nsec) / 1000000000;

    CUT_CHECK_ERROR("Cuda_CommCuda_UnpackComm: Kernel execution failed");

  }
}

void Cuda_CommCuda_UnpackCommVel(cuda_shared_data* sdata, int n, int first, void* buf_recv, int iswap)
{
  my_times time1, time2;

  if(sdata->atom.update_nmax)
    Cuda_CommCuda_UpdateNmax(sdata);

  if(sdata->atom.update_nlocal)
    cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal        , sizeof(int));

  int size = n * 6 * sizeof(X_CFLOAT);

  if(sdata->buffer_new or (size > sdata->buffersize))
    Cuda_CommCuda_UpdateBuffer(sdata, n);

  int3 layout = getgrid(n);
  dim3 threads(layout.z, 1, 1);
  dim3 grid(layout.x, layout.y, 1);

  if(sdata->atom.nlocal > 0) {
    my_gettime(CLOCK_REALTIME, &time1);

    if(not sdata->overlap_comm || iswap < 0)
      cudaMemcpy(sdata->buffer, (void*)buf_recv, n * 6 * sizeof(X_CFLOAT), cudaMemcpyHostToDevice);

    my_gettime(CLOCK_REALTIME, &time2);
    sdata->cuda_timings.comm_forward_upload +=
      time2.tv_sec - time1.tv_sec + 1.0 * (time2.tv_nsec - time1.tv_nsec) / 1000000000;
    void* buf = (sdata->overlap_comm && iswap >= 0) ? sdata->comm.buf_recv_dev[iswap] : sdata->buffer;
    Cuda_CommCuda_UnpackComm_Kernel <<< grid, threads, 0>>>(n, first, buf);
    cudaThreadSynchronize();

    my_gettime(CLOCK_REALTIME, &time1);
    sdata->cuda_timings.comm_forward_kernel_unpack +=
      time1.tv_sec - time2.tv_sec + 1.0 * (time1.tv_nsec - time2.tv_nsec) / 1000000000;

    CUT_CHECK_ERROR("Cuda_CommCuda_UnpackComm: Kernel execution failed");

  }
}

int Cuda_CommCuda_PackReverse(cuda_shared_data* sdata, int n, int first, void* buf_send)
{
  if(sdata->atom.update_nmax)
    Cuda_CommCuda_UpdateNmax(sdata);

  if(sdata->atom.update_nlocal)
    cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal        , sizeof(int));

  int size = n * 3 * sizeof(F_CFLOAT);

  if(sdata->buffer_new or (size > sdata->buffersize))
    Cuda_CommCuda_UpdateBuffer(sdata, n);


  F_CFLOAT* buf = (F_CFLOAT*)buf_send;
  F_CFLOAT* f_dev = (F_CFLOAT*)sdata->atom.f.dev_data;
  f_dev += first;
  cudaMemcpy(buf, f_dev, n * sizeof(F_CFLOAT), cudaMemcpyDeviceToHost);
  buf += n;
  f_dev += sdata->atom.nmax;
  cudaMemcpy(buf, f_dev, n * sizeof(F_CFLOAT), cudaMemcpyDeviceToHost);
  buf += n;
  f_dev += sdata->atom.nmax;
  cudaMemcpy(buf, f_dev, n * sizeof(F_CFLOAT), cudaMemcpyDeviceToHost);
  return 	n * 3;
}


void Cuda_CommCuda_UnpackReverse(cuda_shared_data* sdata, int n, int iswap, void* buf_recv)
{
  if(sdata->atom.update_nmax)
    Cuda_CommCuda_UpdateNmax(sdata);

  if(sdata->atom.update_nlocal)
    cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal        , sizeof(int));

  int size = n * 3 * sizeof(F_CFLOAT);

  if(sdata->buffer_new or (size > sdata->buffersize))
    Cuda_CommCuda_UpdateBuffer(sdata, n);


  int3 layout = getgrid(n);
  dim3 threads(layout.z, 1, 1);
  dim3 grid(layout.x, layout.y, 1);

  if(sdata->atom.nlocal > 0) {
    cudaMemcpy(sdata->buffer, buf_recv, size, cudaMemcpyHostToDevice);
    Cuda_CommCuda_UnpackReverse_Kernel <<< grid, threads, 0>>>((int*) sdata->comm.sendlist.dev_data, n, sdata->comm.maxlistlength, iswap);
    cudaThreadSynchronize();
    CUT_CHECK_ERROR("Cuda_CommCuda_UnpackReverse: Kernel execution failed");
  }
}

void Cuda_CommCuda_UnpackReverse_Self(cuda_shared_data* sdata, int n, int iswap, int first)
{
  if(sdata->atom.update_nmax)
    Cuda_CommCuda_UpdateNmax(sdata);

  if(sdata->atom.update_nlocal)
    cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal        , sizeof(int));

  int size = n * 3 * sizeof(X_CFLOAT);

  if(sdata->buffer_new or (size > sdata->buffersize))
    Cuda_CommCuda_UpdateBuffer(sdata, n);

  int3 layout = getgrid(n);
  dim3 threads(layout.z, 1, 1);
  dim3 grid(layout.x, layout.y, 1);

  if(sdata->atom.nlocal > 0) {
    Cuda_CommCuda_UnpackReverse_Self_Kernel <<< grid, threads, 0>>>((int*) sdata->comm.sendlist.dev_data, n, sdata->comm.maxlistlength, iswap, first);
    cudaThreadSynchronize();
    CUT_CHECK_ERROR("Cuda_CommCuda_PackReverse_Self: Kernel execution failed");

  }
}


int Cuda_CommCuda_BuildSendlist(cuda_shared_data* sdata, int bordergroup, int ineed, int style, int atom_nfirst, int nfirst, int nlast, int dim, int iswap)
{
  MYDBG(printf(" # CUDA: CommCuda_BuildSendlist\n");)
  my_times time1, time2;

  if(sdata->atom.update_nmax)
    Cuda_CommCuda_UpdateNmax(sdata);

  if(sdata->atom.update_nlocal)
    cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal        , sizeof(int));

  if(sdata->buffer_new or (80 > sdata->buffersize))
    Cuda_CommCuda_UpdateBuffer(sdata, 10);

  int n;

  if(!bordergroup || ineed >= 2)
    n = nlast - nfirst + 1;
  else {
    n = atom_nfirst;

    if(nlast - sdata->atom.nlocal + 1 > n) n = nlast - sdata->atom.nlocal + 1;
  }

  int3 layout = getgrid(n, 0, 512, true);
  dim3 threads(layout.z, 1, 1);
  dim3 grid(layout.x + 1, layout.y, 1);


  cudaMemset((int*)(sdata->buffer), 0, sizeof(int));

  my_gettime(CLOCK_REALTIME, &time1);

  if(style == 1)
    Cuda_CommCuda_BuildSendlist_Single <<< grid, threads, (threads.x + 1)*sizeof(int) >>> (bordergroup, ineed, atom_nfirst, nfirst, nlast, dim, iswap, (X_CFLOAT*) sdata->comm.slablo.dev_data, (X_CFLOAT*) sdata->comm.slabhi.dev_data, (int*) sdata->comm.sendlist.dev_data, sdata->comm.maxlistlength);
  else
    Cuda_CommCuda_BuildSendlist_Multi <<< grid, threads, (threads.x + 1)*sizeof(int) >>> (bordergroup, ineed, atom_nfirst, nfirst, nlast, dim, iswap, (X_CFLOAT*) sdata->comm.multilo.dev_data, (X_CFLOAT*) sdata->comm.multihi.dev_data, (int*) sdata->comm.sendlist.dev_data, sdata->comm.maxlistlength);

  cudaThreadSynchronize();
  my_gettime(CLOCK_REALTIME, &time2);
  sdata->cuda_timings.comm_border_kernel_buildlist +=
    time2.tv_sec - time1.tv_sec + 1.0 * (time2.tv_nsec - time1.tv_nsec) / 1000000000;

  CUT_CHECK_ERROR("Cuda_CommCuda_BuildSendlist: Kernel execution failed");
  int nsend;
  cudaMemcpy(&nsend, sdata->buffer, sizeof(int), cudaMemcpyDeviceToHost);
  return nsend;


}

