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
#define MY_PREFIX atom_vec_cuda
#include "cuda_shared.h"
#include "cuda_common.h"
#include "cuda_wrapper_cu.h"
#include "crm_cuda_utils.cu"

#include "atom_vec_cuda_kernel.cu"

int AtomVecCuda_CountDataItems(unsigned int data_mask)
{
  int n = 0;

  if(data_mask & X_MASK) n += 3;

  if(data_mask & V_MASK) n += 3;

  if(data_mask & F_MASK) n += 3;

  if(data_mask & TAG_MASK) n++;

  if(data_mask & TYPE_MASK) n++;

  if(data_mask & MASK_MASK) n++;

  if(data_mask & IMAGE_MASK) n++;

  if(data_mask & Q_MASK) n++;

  if(data_mask & MOLECULE_MASK) n++;

  if(data_mask & RMASS_MASK) n++;

  if(data_mask & RADIUS_MASK) n++;

  if(data_mask & DENSITY_MASK) n++;

  if(data_mask & OMEGA_MASK) n += 3;

  if(data_mask & TORQUE_MASK) n++;

  //if(data_mask & NSPECIAL_MASK) n+=3;
  return n;
}

void Cuda_AtomVecCuda_UpdateBuffer(cuda_shared_data* sdata, int size)
{
  if(sdata->buffersize < size) {
    MYDBG(printf("Cuda_AtomVecCuda Resizing Buffer at %p with %i kB to\n", sdata->buffer, sdata->buffersize);)
    CudaWrapper_FreeCudaData(sdata->buffer, sdata->buffersize);
    sdata->buffer = CudaWrapper_AllocCudaData(size);
    sdata->buffersize = size;
    sdata->buffer_new++;
    MYDBG(printf("New buffer at %p with %i kB\n", sdata->buffer, sdata->buffersize);)
  }

  cudaMemcpyToSymbol(MY_AP(buffer), & sdata->buffer, sizeof(int*));
}

template <const unsigned int data_mask>
void Cuda_AtomVecCuda_UpdateNmax(cuda_shared_data* sdata)
{
  cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal        , sizeof(int));
  cudaMemcpyToSymbol(MY_AP(nmax)    , & sdata->atom.nmax          , sizeof(int));
  cudaMemcpyToSymbol(MY_AP(x)       , & sdata->atom.x    .dev_data, sizeof(X_CFLOAT*));
  cudaMemcpyToSymbol(MY_AP(v)       , & sdata->atom.v    .dev_data, sizeof(V_CFLOAT*));
  cudaMemcpyToSymbol(MY_AP(f)       , & sdata->atom.f    .dev_data, sizeof(F_CFLOAT*));
  cudaMemcpyToSymbol(MY_AP(tag)     , & sdata->atom.tag  .dev_data, sizeof(int*));
  cudaMemcpyToSymbol(MY_AP(type)    , & sdata->atom.type .dev_data, sizeof(int*));
  cudaMemcpyToSymbol(MY_AP(mask)    , & sdata->atom.mask .dev_data, sizeof(int*));
  cudaMemcpyToSymbol(MY_AP(image)   , & sdata->atom.image.dev_data, sizeof(int*));

  if(data_mask & Q_MASK) cudaMemcpyToSymbol(MY_AP(q)       , & sdata->atom.q    .dev_data, sizeof(F_CFLOAT*));

  if(data_mask & MOLECULE_MASK) cudaMemcpyToSymbol(MY_AP(molecule)   , & sdata->atom.molecule.dev_data, sizeof(int*));

  if(data_mask & RADIUS_MASK) cudaMemcpyToSymbol(MY_AP(radius)   , & sdata->atom.radius.dev_data, sizeof(int*));

  if(data_mask & DENSITY_MASK) cudaMemcpyToSymbol(MY_AP(density)   , & sdata->atom.density.dev_data, sizeof(int*));

  if(data_mask & RMASS_MASK) cudaMemcpyToSymbol(MY_AP(rmass)   , & sdata->atom.rmass.dev_data, sizeof(int*));

  if(data_mask & OMEGA_MASK) cudaMemcpyToSymbol(MY_AP(omega)   , & sdata->atom.omega.dev_data, sizeof(int*));

  //if(data_mask & NSPECIAL_MASK) cudaMemcpyToSymbol(MY_AP(nspecial)   , & sdata->atom.nspecial.dev_data, sizeof(int*) );
  cudaMemcpyToSymbol(MY_AP(flag)    , & sdata->flag, sizeof(int*));
}

template <const unsigned int data_mask>
void Cuda_AtomVecCuda_Init(cuda_shared_data* sdata)
{
  MYDBG(printf("# CUDA: Cuda_AtomVecCuda_Init ... start\n");)

  if(sdata->atom.update_nmax)
    Cuda_AtomVecCuda_UpdateNmax<data_mask>(sdata);

  if(sdata->atom.update_nlocal)
    cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal        , sizeof(int));

  MYDBG(printf("# CUDA: Cuda_AtomVecCuda_Init ... post Nmax\n");)
  cudaMemcpyToSymbol(MY_AP(prd)   , sdata->domain.prd, 3 * sizeof(X_CFLOAT));
  cudaMemcpyToSymbol(MY_AP(sublo)   , & sdata->domain.sublo, 3 * sizeof(X_CFLOAT));
  cudaMemcpyToSymbol(MY_AP(subhi)   , & sdata->domain.subhi, 3 * sizeof(X_CFLOAT));
  cudaMemcpyToSymbol(MY_AP(flag)   , & sdata->flag, sizeof(int*));
  cudaThreadSynchronize();
  MYDBG(printf("# CUDA: Cuda_AtomVecCuda_Init ... end\n");)
}


template <const unsigned int data_mask>
int Cuda_AtomVecCuda_PackComm(cuda_shared_data* sdata, int n, int iswap, void* buf_send, int* pbc, int pbc_flag)
{

  my_times time1, time2;

  if(sdata->atom.update_nmax)
    Cuda_AtomVecCuda_UpdateNmax<data_mask>(sdata);

  if(sdata->atom.update_nlocal)
    cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal        , sizeof(int));

  int n_data_items = AtomVecCuda_CountDataItems(data_mask);
  int size = (n * n_data_items) * sizeof(X_CFLOAT);

  if(sdata->buffer_new or (size > sdata->buffersize))
    Cuda_AtomVecCuda_UpdateBuffer(sdata, size);

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
    Cuda_AtomVecCuda_PackComm_Kernel<data_mask> <<< grid, threads, 0>>>((int*) sdata->comm.sendlist.dev_data, n
        , sdata->comm.maxlistlength, iswap, dx, dy, dz, buf);
    cudaThreadSynchronize();

    my_gettime(CLOCK_REALTIME, &time2);
    sdata->cuda_timings.comm_forward_kernel_pack +=
      time2.tv_sec - time1.tv_sec + 1.0 * (time2.tv_nsec - time1.tv_nsec) / 1000000000;

    CUT_CHECK_ERROR("Cuda_AtomVecCuda_PackComm: Kernel execution failed");

    if(not sdata->overlap_comm)
      cudaMemcpy(buf_send, sdata->buffer, n* n_data_items* sizeof(X_CFLOAT), cudaMemcpyDeviceToHost);
    //cudaMemcpy(buf_send, sdata->comm.buf_send_dev[iswap], n*3*sizeof(X_CFLOAT), cudaMemcpyDeviceToHost);

    my_gettime(CLOCK_REALTIME, &time1);
    sdata->cuda_timings.comm_forward_download +=
      time1.tv_sec - time2.tv_sec + 1.0 * (time1.tv_nsec - time2.tv_nsec) / 1000000000;

    int aflag;
    cudaMemcpy(&aflag, sdata->flag, sizeof(int), cudaMemcpyDeviceToHost);
    if(aflag != 0) printf("aflag PackComm: %i\n", aflag);
    CUT_CHECK_ERROR("Cuda_AtomVecCuda_PackComm: Kernel execution failed");

  }

  return n_data_items * n;
}


template <const unsigned int data_mask>
int Cuda_AtomVecCuda_PackComm_Self(cuda_shared_data* sdata, int n, int iswap, int first, int* pbc, int pbc_flag)
{
  MYDBG(printf(" # CUDA: AtomVecCuda_PackComm_Self\n");)
  my_times time1, time2;

  if(sdata->atom.update_nmax)
    Cuda_AtomVecCuda_UpdateNmax<data_mask>(sdata);

  if(sdata->atom.update_nlocal)
    cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal        , sizeof(int));

  int n_data_items = AtomVecCuda_CountDataItems(data_mask);
  int size = (n * n_data_items) * sizeof(X_CFLOAT);

  if(sdata->buffer_new or (size > sdata->buffersize))
    Cuda_AtomVecCuda_UpdateBuffer(sdata, size);

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
    CUT_CHECK_ERROR("Cuda_AtomVecCuda_PackComm_Self:Pre Kernel execution failed");

    Cuda_AtomVecCuda_PackComm_Self_Kernel<data_mask> <<< grid, threads, 0>>>((int*) sdata->comm.sendlist.dev_data, n, sdata->comm.maxlistlength, iswap, dx, dy, dz, first);
    cudaThreadSynchronize();

    my_gettime(CLOCK_REALTIME, &time2);
    sdata->cuda_timings.comm_forward_kernel_self +=
      time2.tv_sec - time1.tv_sec + 1.0 * (time2.tv_nsec - time1.tv_nsec) / 1000000000;

    CUT_CHECK_ERROR("Cuda_AtomVecCuda_PackComm_Self: Kernel execution failed");
  }

  return n_data_items * n;
}


template <const unsigned int data_mask>
void Cuda_AtomVecCuda_UnpackComm(cuda_shared_data* sdata, int n, int first, void* buf_recv, int iswap)
{
  my_times time1, time2;

  if(sdata->atom.update_nmax)
    Cuda_AtomVecCuda_UpdateNmax<data_mask>(sdata);

  if(sdata->atom.update_nlocal)
    cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal        , sizeof(int));

  int n_data_items = AtomVecCuda_CountDataItems(data_mask);
  int size = (n * n_data_items) * sizeof(X_CFLOAT);

  if(sdata->buffer_new or (size > sdata->buffersize))
    Cuda_AtomVecCuda_UpdateBuffer(sdata, size);

  int3 layout = getgrid(n);
  dim3 threads(layout.z, 1, 1);
  dim3 grid(layout.x, layout.y, 1);

  if(sdata->atom.nlocal > 0) {
    my_gettime(CLOCK_REALTIME, &time1);

    if(not sdata->overlap_comm || iswap < 0)
      cudaMemcpy(sdata->buffer, (void*)buf_recv, n_data_items * n * sizeof(X_CFLOAT), cudaMemcpyHostToDevice);

    my_gettime(CLOCK_REALTIME, &time2);
    sdata->cuda_timings.comm_forward_upload +=
      time2.tv_sec - time1.tv_sec + 1.0 * (time2.tv_nsec - time1.tv_nsec) / 1000000000;
    void* buf = (sdata->overlap_comm && iswap >= 0) ? sdata->comm.buf_recv_dev[iswap] : sdata->buffer;
    Cuda_AtomVecCuda_UnpackComm_Kernel<data_mask> <<< grid, threads, 0>>>(n, first, buf);
    cudaThreadSynchronize();

    my_gettime(CLOCK_REALTIME, &time1);
    sdata->cuda_timings.comm_forward_kernel_unpack +=
      time1.tv_sec - time2.tv_sec + 1.0 * (time1.tv_nsec - time2.tv_nsec) / 1000000000;

    CUT_CHECK_ERROR("Cuda_AtomVecCuda_UnpackComm: Kernel execution failed");

  }
}

template <const unsigned int data_mask>
int Cuda_AtomVecCuda_PackExchangeList(cuda_shared_data* sdata, int n, int dim, void* buf_send)
{
  MYDBG(printf("# CUDA: Cuda_AtomVecCuda_PackExchangeList ... start dim %i \n", dim);)
  CUT_CHECK_ERROR("Cuda_AtomVecCuda_PackExchangeList: pre Kernel execution failed");
  cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal        , sizeof(int));
  Cuda_AtomVecCuda_Init<data_mask>(sdata);
  int size = n * sizeof(double);

  if(sdata->buffer_new or (size > sdata->buffersize))
    Cuda_AtomVecCuda_UpdateBuffer(sdata, size);

  cudaMemset((int*)(sdata->buffer), 0, sizeof(int));

  int3 layout = getgrid(sdata->atom.nlocal, sizeof(int), 256, true);
  dim3 threads(layout.z, 1, 1);
  dim3 grid(layout.x, layout.y, 1);

  my_times time1, time2;
  my_gettime(CLOCK_REALTIME, &time1);

  Cuda_AtomVecCuda_PackExchangeList_Kernel <<< grid, threads, (threads.x + 1)*sizeof(int) >>> (n - 1, dim);
  cudaThreadSynchronize();
  CUT_CHECK_ERROR("Cuda_AtomVecCuda_PackExchangeList: Kernel execution failed");

  my_gettime(CLOCK_REALTIME, &time2);
  sdata->cuda_timings.comm_exchange_kernel_pack +=
    time2.tv_sec - time1.tv_sec + 1.0 * (time2.tv_nsec - time1.tv_nsec) / 1000000000;

  cudaMemcpy(buf_send, sdata->buffer, sizeof(double), cudaMemcpyDeviceToHost);
  int return_value = ((int*) buf_send)[0];

  if(n > 1 + return_value)
    cudaMemcpy(buf_send, sdata->buffer, (1 + return_value)*sizeof(double), cudaMemcpyDeviceToHost);

  CUT_CHECK_ERROR("Cuda_AtomVecCuda_PackExchangeList: return copy failed");

  my_gettime(CLOCK_REALTIME, &time1);
  sdata->cuda_timings.comm_exchange_download +=
    time1.tv_sec - time2.tv_sec + 1.0 * (time1.tv_nsec - time2.tv_nsec) / 1000000000;

  MYDBG(printf("# CUDA: Cuda_AtomVecCuda_PackExchangeList ... done\n");)
  return return_value;
}

template <const unsigned int data_mask>
int Cuda_AtomVecCuda_PackExchange(cuda_shared_data* sdata, int nsend, void* buf_send, void* copylist)
{
  MYDBG(printf("# CUDA: Cuda_AtomVecCuda_PackExchange ... start \n");)

  if(sdata->atom.update_nmax)
    Cuda_AtomVecCuda_UpdateNmax<data_mask>(sdata);

  //if(sdata->atom.update_nlocal)
  cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal        , sizeof(int));

  int n_data_items = AtomVecCuda_CountDataItems(data_mask) + 1;
  int size = (nsend * n_data_items + 1) * sizeof(double);

  if(sdata->buffer_new or (size > sdata->buffersize))
    Cuda_AtomVecCuda_UpdateBuffer(sdata, size);

  cudaMemset((int*)(sdata->buffer), 0, sizeof(int));

  int3 layout = getgrid(nsend, 0);
  dim3 threads(layout.z, 1, 1);
  dim3 grid(layout.x, layout.y, 1);

  my_times time1, time2;
  my_gettime(CLOCK_REALTIME, &time1);

  Cuda_AtomVecCuda_PackExchange_Kernel<data_mask> <<< grid, threads, 0>>>(nsend, (int*) copylist);
  cudaThreadSynchronize();
  CUT_CHECK_ERROR("Cuda_AtomVecCuda_PackExchange: Kernel execution failed");

  my_gettime(CLOCK_REALTIME, &time2);
  sdata->cuda_timings.comm_exchange_kernel_pack +=
    time2.tv_sec - time1.tv_sec + 1.0 * (time2.tv_nsec - time1.tv_nsec) / 1000000000;

  cudaMemcpy(buf_send, sdata->buffer, size, cudaMemcpyDeviceToHost);

  my_gettime(CLOCK_REALTIME, &time1);
  sdata->cuda_timings.comm_exchange_download +=
    time1.tv_sec - time2.tv_sec + 1.0 * (time1.tv_nsec - time2.tv_nsec) / 1000000000;

  MYDBG(printf("# CUDA: Cuda_AtomVecCuda_PackExchange ... done\n");)
  return nsend * n_data_items + 1;
}


template <const unsigned int data_mask>
int Cuda_AtomVecCuda_UnpackExchange(cuda_shared_data* sdata, int nsend, void* buf_send, void* copylist)
{
  Cuda_AtomVecCuda_UpdateNmax<data_mask>(sdata);
  cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal        , sizeof(int));
  int n_data_items = AtomVecCuda_CountDataItems(data_mask) + 1;

  int size = (nsend * n_data_items + 1) * sizeof(double);

  if(sdata->buffer_new or (size > sdata->buffersize))
    Cuda_AtomVecCuda_UpdateBuffer(sdata, size);

  cudaMemcpyToSymbol(MY_AP(flag)   , & sdata->flag, sizeof(int*));

  cudaMemset((int*)(sdata->flag), 0, sizeof(int));

  if(nsend) {
    int3 layout = getgrid(nsend, 0);
    dim3 threads(layout.z, 1, 1);
    dim3 grid(layout.x, layout.y, 1);

    if(sdata->atom.nlocal > 0) {
      my_times time1, time2;
      my_gettime(CLOCK_REALTIME, &time1);

      cudaMemcpy(sdata->buffer, buf_send , size, cudaMemcpyHostToDevice);

      my_gettime(CLOCK_REALTIME, &time2);
      sdata->cuda_timings.comm_exchange_upload +=
        time2.tv_sec - time1.tv_sec + 1.0 * (time2.tv_nsec - time1.tv_nsec) / 1000000000;

      Cuda_AtomVecCuda_UnpackExchange_Kernel<data_mask> <<< grid, threads, 0>>>(sdata->exchange_dim, nsend, (int*) copylist);
      cudaThreadSynchronize();

      my_gettime(CLOCK_REALTIME, &time1);
      sdata->cuda_timings.comm_exchange_kernel_unpack +=
        time1.tv_sec - time2.tv_sec + 1.0 * (time1.tv_nsec - time2.tv_nsec) / 1000000000;

      CUT_CHECK_ERROR("Cuda_AtomVecCuda_UnpackExchange: Kernel execution failed");
    }
  }

  int naccept;
  cudaMemcpy((void*)&naccept, sdata->flag, sizeof(int), cudaMemcpyDeviceToHost);

  return naccept;
}

template <const unsigned int data_mask>
int Cuda_AtomVecCuda_PackBorder(cuda_shared_data* sdata, int nsend, int iswap, void* buf_send, int* pbc, int pbc_flag)
{
  my_times atime1, atime2;
  my_gettime(CLOCK_REALTIME, &atime1);

  if(sdata->atom.update_nmax)
    Cuda_AtomVecCuda_UpdateNmax<data_mask>(sdata);

  if(sdata->atom.update_nlocal)
    cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal        , sizeof(int));

  my_gettime(CLOCK_REALTIME, &atime2);
  sdata->cuda_timings.test1 +=
    atime2.tv_sec - atime1.tv_sec + 1.0 * (atime2.tv_nsec - atime1.tv_nsec) / 1000000000;

  int n_data_items = AtomVecCuda_CountDataItems(data_mask);

  int size = nsend * n_data_items * sizeof(X_CFLOAT);

  if(sdata->buffer_new or (size > sdata->buffersize))
    Cuda_AtomVecCuda_UpdateBuffer(sdata, size);

  X_CFLOAT dx = 0.0;
  X_CFLOAT dy = 0.0;
  X_CFLOAT dz = 0.0;

  if(pbc_flag != 0) {
    if(sdata->domain.triclinic == 0) {
      dx = pbc[0] * sdata->domain.prd[0];
      dy = pbc[1] * sdata->domain.prd[1];
      dz = pbc[2] * sdata->domain.prd[2];
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
  }

  int3 layout = getgrid(nsend);
  dim3 threads(layout.z, 1, 1);
  dim3 grid(layout.x, layout.y, 1);

  if(sdata->atom.nlocal > 0) {
    my_times time1, time2;
    my_gettime(CLOCK_REALTIME, &time1);

    Cuda_AtomVecCuda_PackBorder_Kernel<data_mask> <<< grid, threads, 0>>>((int*) sdata->comm.sendlist.dev_data, nsend, sdata->comm.maxlistlength, iswap, dx, dy, dz);
    cudaThreadSynchronize();

    my_gettime(CLOCK_REALTIME, &time2);
    sdata->cuda_timings.comm_border_kernel_pack +=
      time2.tv_sec - time1.tv_sec + 1.0 * (time2.tv_nsec - time1.tv_nsec) / 1000000000;

    cudaMemcpy(buf_send, sdata->buffer, size, cudaMemcpyDeviceToHost);
    CUT_CHECK_ERROR("Cuda_AtomVecCuda_PackBorder: Kernel execution failed");

    my_gettime(CLOCK_REALTIME, &time1);
    sdata->cuda_timings.comm_border_download +=
      time1.tv_sec - time2.tv_sec + 1.0 * (time1.tv_nsec - time2.tv_nsec) / 1000000000;

  }

  return nsend * n_data_items;
}

template <const unsigned int data_mask>
int Cuda_AtomVecCuda_PackBorder_Self(cuda_shared_data* sdata, int n, int iswap, int first, int* pbc, int pbc_flag)
{
  if(sdata->atom.update_nmax)
    Cuda_AtomVecCuda_UpdateNmax<data_mask>(sdata);

  if(sdata->atom.update_nlocal)
    cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal        , sizeof(int));

  int n_data_items = AtomVecCuda_CountDataItems(data_mask);

  int size = n * n_data_items * sizeof(X_CFLOAT);

  if(sdata->buffer_new or (size > sdata->buffersize))
    Cuda_AtomVecCuda_UpdateBuffer(sdata, size);

  X_CFLOAT dx = 0.0;
  X_CFLOAT dy = 0.0;
  X_CFLOAT dz = 0.0;

  if(pbc_flag != 0) {
    if(sdata->domain.triclinic == 0) {
      dx = pbc[0] * sdata->domain.prd[0];
      dy = pbc[1] * sdata->domain.prd[1];
      dz = pbc[2] * sdata->domain.prd[2];
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
  }

  int3 layout = getgrid(n);
  dim3 threads(layout.z, 1, 1);
  dim3 grid(layout.x, layout.y, 1);

  if(sdata->atom.nlocal > 0) {
    my_times time1, time2;
    my_gettime(CLOCK_REALTIME, &time1);

    Cuda_AtomVecCuda_PackBorder_Self_Kernel<data_mask> <<< grid, threads, 0>>>((int*) sdata->comm.sendlist.dev_data, n, sdata->comm.maxlistlength, iswap, dx, dy, dz, first);
    cudaThreadSynchronize();

    my_gettime(CLOCK_REALTIME, &time2);
    sdata->cuda_timings.comm_border_kernel_self +=
      time2.tv_sec - time1.tv_sec + 1.0 * (time2.tv_nsec - time1.tv_nsec) / 1000000000;

    CUT_CHECK_ERROR("Cuda_AtomVecCuda_PackBorder_Self: Kernel execution failed");

  }

  return n * n_data_items;
}


template <const unsigned int data_mask>
int Cuda_AtomVecCuda_UnpackBorder(cuda_shared_data* sdata, int n, int first, void* buf_recv)
{
  my_times atime1, atime2;
  my_gettime(CLOCK_REALTIME, &atime1);

  if(sdata->atom.update_nmax)
    Cuda_AtomVecCuda_UpdateNmax<data_mask>(sdata);

  if(sdata->atom.update_nlocal)
    cudaMemcpyToSymbol(MY_AP(nlocal)  , & sdata->atom.nlocal        , sizeof(int));

  my_gettime(CLOCK_REALTIME, &atime2);
  sdata->cuda_timings.test1 +=
    atime2.tv_sec - atime1.tv_sec + 1.0 * (atime2.tv_nsec - atime1.tv_nsec) / 1000000000;

  int n_data_items = AtomVecCuda_CountDataItems(data_mask);

  int size = n * n_data_items * sizeof(X_CFLOAT);

  if(sdata->buffer_new or (size > sdata->buffersize))
    Cuda_AtomVecCuda_UpdateBuffer(sdata, size);

  int3 layout = getgrid(n);
  dim3 threads(layout.z, 1, 1);
  dim3 grid(layout.x, layout.y, 1);

  if(sdata->atom.nlocal > 0) {
    my_times time1, time2;
    my_gettime(CLOCK_REALTIME, &time1);

    cudaMemset((int*)(sdata->flag), 0, sizeof(int));
    cudaMemcpy(sdata->buffer, (void*)buf_recv, size, cudaMemcpyHostToDevice);

    my_gettime(CLOCK_REALTIME, &time2);
    sdata->cuda_timings.comm_border_upload +=
      time2.tv_sec - time1.tv_sec + 1.0 * (time2.tv_nsec - time1.tv_nsec) / 1000000000;

    Cuda_AtomVecCuda_UnpackBorder_Kernel<data_mask> <<< grid, threads, 0>>>(n, first);
    cudaThreadSynchronize();

    my_gettime(CLOCK_REALTIME, &time1);
    sdata->cuda_timings.comm_border_kernel_unpack +=
      time1.tv_sec - time2.tv_sec + 1.0 * (time1.tv_nsec - time2.tv_nsec) / 1000000000;

    cudaMemcpy(&sdata->comm.grow_flag, sdata->flag, sizeof(int), cudaMemcpyDeviceToHost);

    CUT_CHECK_ERROR("Cuda_AtomVecCuda_UnpackBorder: Kernel execution failed");

  }

  return sdata->comm.grow_flag;
}


#include "atom_vec_angle_cuda.cu"
#include "atom_vec_atomic_cuda.cu"
#include "atom_vec_charge_cuda.cu"
#include "atom_vec_full_cuda.cu"
//#include "atom_vec_granular_cuda.cu"
