/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author (triclinic) : Pieter in 't Veld (SNL)
------------------------------------------------------------------------- */

#ifdef LAMMPS_BIGBIG
#error LAMMPS_BIGBIG not supported by this file
#endif

#include "mpi.h"
#include <cmath>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include "comm_cuda.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "pair.h"
#include "domain.h"
#include "neighbor.h"
#include "modify.h"
#include "fix.h"
#include "group.h"
#include "compute.h"
#include "user_cuda.h"
#include "error.h"
#include "memory.h"
#include "comm_cuda_cu.h"

using namespace LAMMPS_NS;

#define BUFFACTOR 1.5
#define BUFMIN 1000
#define BUFEXTRA 1000



#define BIG 1.0e20

enum{SINGLE,MULTI};

/* ----------------------------------------------------------------------
   setup MPI and allocate buffer space
------------------------------------------------------------------------- */

CommCuda::CommCuda(LAMMPS *lmp) : CommBrick(lmp)
{
  cuda = lmp->cuda;
   if(cuda == NULL)
        error->all(FLERR,"You cannot use a /cuda class, without activating 'cuda' acceleration. Provide '-c on' as command-line argument to LAMMPS..");

  cu_pbc=NULL;
  cu_slablo=NULL;
  cu_slabhi=NULL;
  cu_multilo=NULL;
  cu_multihi=NULL;
  cu_sendlist=NULL;


  memory->sfree(buf_send);
  memory->sfree(buf_recv);
  buf_send = NULL;
  buf_recv = NULL;

  CommBrick::free_swap();
  allocate_swap(maxswap);
}

/* ---------------------------------------------------------------------- */

CommCuda::~CommCuda()
{
  delete cu_sendlist;
  if(cuda->pinned)
  {
    CudaWrapper_FreePinnedHostData((void*)buf_send);
    CudaWrapper_FreePinnedHostData((void*)buf_recv);
  }
  else
  {
    memory->sfree(buf_send);
    memory->sfree(buf_recv);
  }
  buf_send=NULL;
  buf_recv=NULL;
}

/* ---------------------------------------------------------------------- */

void CommCuda::init()
{
  int factor = 1;
  if(cuda->shared_data.overlap_comm) factor=maxswap;
  if(not buf_send)
  grow_send(maxsend,0);
  if(not buf_recv)
  grow_recv(maxrecv);
  if(not cu_sendlist)
  {
    cu_sendlist=new cCudaData<int, int, xy> ((int*)sendlist,maxswap,BUFMIN);
    cuda->shared_data.comm.sendlist.dev_data=cu_sendlist->dev_data();
    cuda->shared_data.comm.maxswap=maxswap;
    cuda->shared_data.comm.maxlistlength=BUFMIN;
    cu_sendlist->upload();
  }
  delete cu_pbc;
  cu_pbc=new cCudaData<int, int, xy> ((int*)pbc,cuda->shared_data.comm.maxswap,6);
  cu_pbc->upload();

  delete cu_slablo;
  cu_slablo = new cCudaData<double, X_CFLOAT,x>(slablo,cuda->shared_data.comm.maxswap);
  cu_slablo->upload();

  delete cu_slabhi;
  cu_slabhi = new cCudaData<double, X_CFLOAT,x>(slabhi,cuda->shared_data.comm.maxswap);
  cu_slabhi->upload();

  cuda->shared_data.comm.pbc.dev_data=cu_pbc->dev_data();
  cuda->shared_data.comm.slablo.dev_data=cu_slablo->dev_data();
  cuda->shared_data.comm.slabhi.dev_data=cu_slabhi->dev_data();

  CommBrick::init();
}

/* ----------------------------------------------------------------------
   setup spatial-decomposition communication patterns
   function of neighbor cutoff(s) & cutghostuser & current box size
   single style sets slab boundaries (slablo,slabhi) based on max cutoff
   multi style sets type-dependent slab boundaries (multilo,multihi)
------------------------------------------------------------------------- */

void CommCuda::setup()
{
  if(cuda->shared_data.pair.neighall) cutghostuser = MAX(2.0*neighbor->cutneighmax,cutghostuser);
  CommBrick::setup();

  //upload changed geometry to device
    if(style == SINGLE)
    {
            if(cu_slablo) cu_slablo->upload();
            if(cu_slabhi) cu_slabhi->upload();
    }
        else
    {
            if(cu_multilo) cu_multilo->upload();
            if(cu_multihi) cu_multihi->upload();
    }
}

/* ----------------------------------------------------------------------
   forward communication of atom coords every timestep
   other per-atom attributes may also be sent via pack/unpack routines
------------------------------------------------------------------------- */

void CommCuda::forward_comm(int mode)
{
  if(mode==0) return forward_comm_cuda();
  if(mode==1) return forward_comm_pack_cuda();
  if(mode==2) return forward_comm_transfer_cuda();
  if(mode==3) return forward_comm_unpack_cuda();
}


void CommCuda::forward_comm_cuda()
{
  static int count=0;
  static double kerneltime=0.0;
  static double copytime=0.0;
  my_times time1,time2,time3;

  int n;
  MPI_Request request;
  MPI_Status status;
  AtomVec *avec = atom->avec;
  double **x = atom->x;

  cuda->shared_data.domain.xy=domain->xy;
  cuda->shared_data.domain.xz=domain->xz;
  cuda->shared_data.domain.yz=domain->yz;
  cuda->shared_data.domain.prd[0]=domain->prd[0];
  cuda->shared_data.domain.prd[1]=domain->prd[1];
  cuda->shared_data.domain.prd[2]=domain->prd[2];
  cuda->shared_data.domain.triclinic=domain->triclinic;
  if(not comm_x_only && not avec->cudable)
  {
          cuda->downloadAll();
    CommBrick::forward_comm();
    cuda->uploadAll();
    return;
  }

  // exchange data with another proc
  // if other proc is self, just copy
  // if comm_x_only set, exchange or copy directly to x, don't unpack

  for (int iswap = 0; iswap < nswap; iswap++) {
    if (sendproc[iswap] != me)
    {
      if (comm_x_only)
      {

        int size_forward_recv_now=0;

        if((sizeof(X_CFLOAT)!=sizeof(double)) && size_forward_recv[iswap]) //some complicated way to safe some transfer size if single precision is used
          size_forward_recv_now=(size_forward_recv[iswap]+1)*sizeof(X_CFLOAT)/sizeof(double);
        else
          size_forward_recv_now=size_forward_recv[iswap];
my_gettime(CLOCK_REALTIME,&time1);

        MPI_Irecv(buf_recv,size_forward_recv_now,MPI_DOUBLE,
                 recvproc[iswap],0,world,&request);
        n = Cuda_CommCuda_PackComm(&cuda->shared_data,sendnum[iswap],iswap,(void*) buf_send,pbc[iswap],pbc_flag[iswap]);

my_gettime(CLOCK_REALTIME,&time2);

        if((sizeof(X_CFLOAT)!=sizeof(double)) && n) //some complicated way to safe some transfer size if single precision is used
          n=(n+1)*sizeof(X_CFLOAT)/sizeof(double);

                //printf("RecvSize: %i SendSize: %i\n",size_forward_recv_now,n);
            MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap],0,world);
            MPI_Wait(&request,&status);

my_gettime(CLOCK_REALTIME,&time3);
cuda->shared_data.cuda_timings.comm_forward_mpi_upper+=
      time3.tv_sec-time1.tv_sec+1.0*(time3.tv_nsec-time1.tv_nsec)/1000000000;
cuda->shared_data.cuda_timings.comm_forward_mpi_lower+=
      time3.tv_sec-time2.tv_sec+1.0*(time3.tv_nsec-time2.tv_nsec)/1000000000;

        Cuda_CommCuda_UnpackComm(&cuda->shared_data,recvnum[iswap],firstrecv[iswap],(void*)buf_recv,iswap); //Unpack for cpu exchange happens implicitely since buf==x[firstrecv]

      }
      else if (ghost_velocity)
      {
            MPI_Irecv(buf_recv,size_forward_recv[iswap],MPI_DOUBLE,
                  recvproc[iswap],0,world,&request);

        if(avec->cudable)
          n = avec->pack_comm_vel(sendnum[iswap],&iswap,
                           buf_send,pbc_flag[iswap],pbc[iswap]);
        else
              n = avec->pack_comm_vel(sendnum[iswap],sendlist[iswap],
                                buf_send,pbc_flag[iswap],pbc[iswap]);

            MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap],0,world);
            MPI_Wait(&request,&status);
            avec->unpack_comm_vel(recvnum[iswap],firstrecv[iswap],buf_recv);
      }
      else
      {
            MPI_Irecv(buf_recv,size_forward_recv[iswap],MPI_DOUBLE,
                  recvproc[iswap],0,world,&request);

        if(avec->cudable)
          n = avec->pack_comm(sendnum[iswap],&iswap,
                           buf_send,pbc_flag[iswap],pbc[iswap]);
        else
              n = avec->pack_comm(sendnum[iswap],sendlist[iswap],
                            buf_send,pbc_flag[iswap],pbc[iswap]);

            MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap],0,world);
            MPI_Wait(&request,&status);
            avec->unpack_comm(recvnum[iswap],firstrecv[iswap],buf_recv);
      }

    }
    else  //sendproc == me
    {
      cuda->self_comm=1;
      if (comm_x_only)
      {
            if (sendnum[iswap])
                {
          n = Cuda_CommCuda_PackComm_Self(&cuda->shared_data,sendnum[iswap],iswap,firstrecv[iswap],pbc[iswap],pbc_flag[iswap]);
          if(n<0) error->all(FLERR," # CUDA ERRROR on PackComm_Self");
          if((sizeof(X_CFLOAT)!=sizeof(double)) && n)
            n=(n+1)*sizeof(X_CFLOAT)/sizeof(double);
                }
      }
      else if (ghost_velocity)
      {
                n = avec->pack_comm_vel(sendnum[iswap],&iswap,
                                (double*) firstrecv,pbc_flag[iswap],pbc[iswap]);
            //avec->unpack_comm_vel(recvnum[iswap],firstrecv[iswap],(double*) firstrecv);
      }
      else
      {
                n = avec->pack_comm(sendnum[iswap],&iswap,
                            (double*) firstrecv,pbc_flag[iswap],pbc[iswap]);
                //avec->unpack_comm(recvnum[iswap],firstrecv[iswap],(double*) firstrecv);
      }
      cuda->self_comm=0;
    }
  }
}

void CommCuda::forward_comm_pack_cuda()
{
        static int count=0;
        static double kerneltime=0.0;
        static double copytime=0.0;
    my_times time1,time2,time3;
  int n;  // initialize comm buffers & exchange memory

  MPI_Request request;
  MPI_Status status;
  AtomVec *avec = atom->avec;
  double **x = atom->x;

  cuda->shared_data.domain.xy=domain->xy;
  cuda->shared_data.domain.xz=domain->xz;
  cuda->shared_data.domain.yz=domain->yz;
  cuda->shared_data.domain.prd[0]=domain->prd[0];
  cuda->shared_data.domain.prd[1]=domain->prd[1];
  cuda->shared_data.domain.prd[2]=domain->prd[2];
  cuda->shared_data.domain.triclinic=domain->triclinic;
  if(not comm_x_only && not avec->cudable) cuda->downloadAll();  //if not comm_x_only the communication routine of the atom_vec style class is used

  // exchange data with another proc
  // if other proc is self, just copy
  // if comm_x_only set, exchange or copy directly to x, don't unpack

  for (int iswap = 0; iswap < nswap; iswap++) {
    if (sendproc[iswap] != me)
    {
      if (comm_x_only)
      {


my_gettime(CLOCK_REALTIME,&time1);

      //  n = Cuda_CommCuda_PackComm(&cuda->shared_data,sendnum[iswap],iswap,(void*) cuda->shared_data.comm.buf_send[iswap],pbc[iswap],pbc_flag[iswap]);
                  n = Cuda_CommCuda_PackComm(&cuda->shared_data,sendnum[iswap],iswap,(void*)buf_send,pbc[iswap],pbc_flag[iswap]);

my_gettime(CLOCK_REALTIME,&time2);

        if((sizeof(X_CFLOAT)!=sizeof(double)) && n) //some complicated way to safe some transfer size if single precision is used
          n=(n+1)*sizeof(X_CFLOAT)/sizeof(double);
                cuda->shared_data.comm.send_size[iswap]=n;
      }
      else if (ghost_velocity)
      {
my_gettime(CLOCK_REALTIME,&time1);

       // n = Cuda_CommCuda_PackComm_Vel(&cuda->shared_data,sendnum[iswap],iswap,(void*) &buf_send[iswap*maxsend],pbc[iswap],pbc_flag[iswap]);

my_gettime(CLOCK_REALTIME,&time2);

        if((sizeof(X_CFLOAT)!=sizeof(double)) && n) //some complicated way to safe some transfer size if single precision is used
          n=(n+1)*sizeof(X_CFLOAT)/sizeof(double);
                cuda->shared_data.comm.send_size[iswap]=n;
       }
      else
      {
            MPI_Irecv(buf_recv,size_forward_recv[iswap],MPI_DOUBLE,
                  recvproc[iswap],0,world,&request);

        if(avec->cudable)
          n = avec->pack_comm(sendnum[iswap],&iswap,
                           cuda->shared_data.comm.buf_send[iswap],pbc_flag[iswap],pbc[iswap]);
        else
              n = avec->pack_comm(sendnum[iswap],sendlist[iswap],
                            cuda->shared_data.comm.buf_send[iswap],pbc_flag[iswap],pbc[iswap]);

            MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap],0,world);
            MPI_Wait(&request,&status);
            avec->unpack_comm(recvnum[iswap],firstrecv[iswap],buf_recv);
      }

    }
    else  //sendproc == me
    {
      if (comm_x_only)
      {
            if (sendnum[iswap])
                {
          n = Cuda_CommCuda_PackComm_Self(&cuda->shared_data,sendnum[iswap],iswap,firstrecv[iswap],pbc[iswap],pbc_flag[iswap]);
          if(n<0) error->all(FLERR," # CUDA ERRROR on PackComm_Self");
          if((sizeof(X_CFLOAT)!=sizeof(double)) && n)
            n=(n+1)*sizeof(X_CFLOAT)/sizeof(double);
                }
      }
      else if (ghost_velocity)
      {
                n = avec->pack_comm_vel(sendnum[iswap],sendlist[iswap],
                                buf_send,pbc_flag[iswap],pbc[iswap]);
            avec->unpack_comm_vel(recvnum[iswap],firstrecv[iswap],buf_send);
      }
      else
      {
                n = avec->pack_comm(sendnum[iswap],sendlist[iswap],
                            buf_send,pbc_flag[iswap],pbc[iswap]);
                avec->unpack_comm(recvnum[iswap],firstrecv[iswap],buf_send);
      }
    }
  }
  if(not comm_x_only && not avec->cudable) cuda->uploadAll();
}

void CommCuda::forward_comm_transfer_cuda()
{
        static int count=0;
        static double kerneltime=0.0;
        static double copytime=0.0;
    my_times time1,time2,time3;
  int n;
  MPI_Request request;
  MPI_Status status;
  AtomVec *avec = atom->avec;
  double **x = atom->x;
  cuda->shared_data.domain.xy=domain->xy;
  cuda->shared_data.domain.xz=domain->xz;
  cuda->shared_data.domain.yz=domain->yz;
  cuda->shared_data.domain.prd[0]=domain->prd[0];
  cuda->shared_data.domain.prd[1]=domain->prd[1];
  cuda->shared_data.domain.prd[2]=domain->prd[2];
  cuda->shared_data.domain.triclinic=domain->triclinic;
  if(not comm_x_only && not avec->cudable) cuda->downloadAll();  //if not comm_x_only the communication routine of the atom_vec style class is used
//printf("A\n");
  // exchange data with another proc
  // if other proc is self, just copy
  // if comm_x_only set, exchange or copy directly to x, don't unpack

  for (int iswap = 0; iswap < nswap; iswap++) {
    if (sendproc[iswap] != me)
    {
      if (comm_x_only)
      {

        int size_forward_recv_now=0;

        if((sizeof(X_CFLOAT)!=sizeof(double)) && size_forward_recv[iswap]) //some complicated way to safe some transfer size if single precision is used
          size_forward_recv_now=(size_forward_recv[iswap]+1)*sizeof(X_CFLOAT)/sizeof(double);
        else
          size_forward_recv_now=size_forward_recv[iswap];

        //printf("A: %i \n",size_forward_recv_now/1024*4);
        //MPI_Irecv(cuda->shared_data.comm.buf_recv[iswap],size_forward_recv_now,MPI_DOUBLE,
        //         recvproc[iswap],0,world,&request);
        MPI_Irecv(buf_recv,size_forward_recv_now,MPI_DOUBLE,
                 recvproc[iswap],0,world,&request);
                //printf("%p %p %i\n",buf_send, cuda->shared_data.comm.buf_send_dev[iswap], cuda->shared_data.comm.send_size[iswap]*sizeof(double));
        //memcpy(buf_send,cuda->shared_data.comm.buf_send[iswap],cuda->shared_data.comm.send_size[iswap]*sizeof(double));
        //        CudaWrapper_SyncStream(1);
        //printf("B: %i \n",cuda->shared_data.comm.send_size[iswap]/1024*4);
                CudaWrapper_DownloadCudaDataAsync((void*) buf_send, cuda->shared_data.comm.buf_send_dev[iswap], cuda->shared_data.comm.send_size[iswap]*sizeof(double),2);
            //MPI_Send(cuda->shared_data.comm.buf_send[iswap],cuda->shared_data.comm.send_size[iswap],MPI_DOUBLE,sendproc[iswap],0,world);
my_gettime(CLOCK_REALTIME,&time1);
        CudaWrapper_SyncStream(2);
        //printf("C: %i \n",cuda->shared_data.comm.send_size[iswap]/1024*4);
my_gettime(CLOCK_REALTIME,&time2);
cuda->shared_data.cuda_timings.comm_forward_download+=
      time2.tv_sec-time1.tv_sec+1.0*(time2.tv_nsec-time1.tv_nsec)/1000000000;
            MPI_Send(buf_send,cuda->shared_data.comm.send_size[iswap],MPI_DOUBLE,sendproc[iswap],0,world);
            MPI_Wait(&request,&status);
        //printf("D: %i \n",cuda->shared_data.comm.send_size[iswap]/1024*4);
                CudaWrapper_UploadCudaDataAsync((void*) buf_recv,cuda->shared_data.comm.buf_recv_dev[iswap], size_forward_recv_now*sizeof(double),2);
my_gettime(CLOCK_REALTIME,&time1);
        CudaWrapper_SyncStream(2);
        //printf("E: %i \n",cuda->shared_data.comm.send_size[iswap]/1024*4);
        //memcpy(cuda->shared_data.comm.buf_recv[iswap],buf_recv,size_forward_recv_now*sizeof(double));
                 //printf("RecvSize: %i SendSize: %i\n",size_forward_recv_now*sizeof(double),cuda->shared_data.comm.send_size[iswap]*sizeof(double));
my_gettime(CLOCK_REALTIME,&time3);
cuda->shared_data.cuda_timings.comm_forward_upload+=
      time3.tv_sec-time1.tv_sec+1.0*(time3.tv_nsec-time1.tv_nsec)/1000000000;
cuda->shared_data.cuda_timings.comm_forward_mpi_lower+=
      time3.tv_sec-time2.tv_sec+1.0*(time3.tv_nsec-time2.tv_nsec)/1000000000;
my_gettime(CLOCK_REALTIME,&time3);
cuda->shared_data.cuda_timings.comm_forward_mpi_upper+=
      time3.tv_sec-time1.tv_sec+1.0*(time3.tv_nsec-time1.tv_nsec)/1000000000;
      }
      else if (ghost_velocity)
      {
 /*       int size_forward_recv_now=0;

        if((sizeof(X_CFLOAT)!=sizeof(double)) && size_forward_recv[iswap]) //some complicated way to safe some transfer size if single precision is used
          size_forward_recv_now=(size_forward_recv[iswap]+1)*sizeof(X_CFLOAT)/sizeof(double);
        else
          size_forward_recv_now=size_forward_recv[iswap];

my_gettime(CLOCK_REALTIME,&time1);

        MPI_Irecv(cuda->shared_data.comm.buf_recv[iswap],size_forward_recv_now,MPI_DOUBLE,
                 recvproc[iswap],0,world,&request);

my_gettime(CLOCK_REALTIME,&time2);

            MPI_Send(cuda->shared_data.comm.buf_send[iswap],cuda->shared_data.comm.send_size[iswap],MPI_DOUBLE,sendproc[iswap],0,world);
            MPI_Wait(&request,&status);

my_gettime(CLOCK_REALTIME,&time3);
cuda->shared_data.cuda_timings.comm_forward_mpi_upper+=
      time3.tv_sec-time1.tv_sec+1.0*(time3.tv_nsec-time1.tv_nsec)/1000000000;
cuda->shared_data.cuda_timings.comm_forward_mpi_lower+=
      time3.tv_sec-time2.tv_sec+1.0*(time3.tv_nsec-time2.tv_nsec)/1000000000;*/

       }
      else
      {
            MPI_Irecv(buf_recv,size_forward_recv[iswap],MPI_DOUBLE,
                  recvproc[iswap],0,world,&request);

        if(avec->cudable)
          n = avec->pack_comm(sendnum[iswap],&iswap,
                           buf_send,pbc_flag[iswap],pbc[iswap]);
        else
              n = avec->pack_comm(sendnum[iswap],sendlist[iswap],
                            buf_send,pbc_flag[iswap],pbc[iswap]);

            MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap],0,world);
            MPI_Wait(&request,&status);
            avec->unpack_comm(recvnum[iswap],firstrecv[iswap],buf_recv);
      }

    }
    else  //sendproc == me
    {
      if (comm_x_only)
      {
            if (sendnum[iswap])
                {
                }
      }
      else if (ghost_velocity)
      {
      }
      else
      {
                n = avec->pack_comm(sendnum[iswap],sendlist[iswap],
                            buf_send,pbc_flag[iswap],pbc[iswap]);
                avec->unpack_comm(recvnum[iswap],firstrecv[iswap],buf_send);
      }
    }
  }
  if(not comm_x_only && not avec->cudable) cuda->uploadAll();
}

void CommCuda::forward_comm_unpack_cuda()
{
        static int count=0;
        static double kerneltime=0.0;
        static double copytime=0.0;
    my_times time1,time2,time3;
  int n;
  MPI_Request request;
  MPI_Status status;
  AtomVec *avec = atom->avec;
  double **x = atom->x;

  cuda->shared_data.domain.xy=domain->xy;
  cuda->shared_data.domain.xz=domain->xz;
  cuda->shared_data.domain.yz=domain->yz;
  cuda->shared_data.domain.prd[0]=domain->prd[0];
  cuda->shared_data.domain.prd[1]=domain->prd[1];
  cuda->shared_data.domain.prd[2]=domain->prd[2];
  cuda->shared_data.domain.triclinic=domain->triclinic;
  if(not comm_x_only && not avec->cudable) cuda->downloadAll();  //if not comm_x_only the communication routine of the atom_vec style class is used

  // exchange data with another proc
  // if other proc is self, just copy
  // if comm_x_only set, exchange or copy directly to x, don't unpack

  for (int iswap = 0; iswap < nswap; iswap++) {
    if (sendproc[iswap] != me)
    {
      if (comm_x_only)
      {

        //Cuda_CommCuda_UnpackComm(&cuda->shared_data,recvnum[iswap],firstrecv[iswap],cuda->shared_data.comm.buf_recv[iswap],iswap); //Unpack for cpu exchange happens implicitely since buf==x[firstrecv]
        Cuda_CommCuda_UnpackComm(&cuda->shared_data,recvnum[iswap],firstrecv[iswap],buf_recv,iswap); //Unpack for cpu exchange happens implicitely since buf==x[firstrecv]

      }
      else if (ghost_velocity)
      {
        //Cuda_CommCuda_UnpackComm_Vel(&cuda->shared_data,recvnum[iswap],firstrecv[iswap],(void*)&buf_recv[iswap*maxrecv]); //Unpack for cpu exchange happens implicitely since buf==x[firstrecv]
      }
      else
      {
            MPI_Irecv(buf_recv,size_forward_recv[iswap],MPI_DOUBLE,
                  recvproc[iswap],0,world,&request);

        if(avec->cudable)
          n = avec->pack_comm(sendnum[iswap],&iswap,
                           buf_send,pbc_flag[iswap],pbc[iswap]);
        else
              n = avec->pack_comm(sendnum[iswap],sendlist[iswap],
                            buf_send,pbc_flag[iswap],pbc[iswap]);

            MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap],0,world);
            MPI_Wait(&request,&status);
            avec->unpack_comm(recvnum[iswap],firstrecv[iswap],buf_recv);
      }

    }
    else  //sendproc == me
    {
      if (comm_x_only)
      {
            if (sendnum[iswap])
                {
                }
      }
      else if (ghost_velocity)
      {
      }
      else
      {
                n = avec->pack_comm(sendnum[iswap],sendlist[iswap],
                            buf_send,pbc_flag[iswap],pbc[iswap]);
                avec->unpack_comm(recvnum[iswap],firstrecv[iswap],buf_send);
      }
    }
  }
  if(not comm_x_only && not avec->cudable) cuda->uploadAll();
}

void CommCuda::forward_comm_pair(Pair *pair)
{
  if(not cuda->shared_data.pair.cudable_force)
  {
          return CommBrick::forward_comm_pair(pair);
  }

  int iswap,n;
  double *buf;
  MPI_Request request;
  MPI_Status status;

  int nsize = pair->comm_forward;

  for (iswap = 0; iswap < nswap; iswap++) {

    // pack buffer

    n = pair->pack_forward_comm(sendnum[iswap],&iswap,
                                buf_send,pbc_flag[iswap],pbc[iswap]);
        int nrecv = recvnum[iswap]*nsize;
        if(nrecv<0) nrecv=-(nrecv+1)/2;
        int nsend = n;
        if(nsend<0) nsend=-(nsend+1)/2;

    // exchange with another proc
    // if self, set recv buffer to send buffer

    if (sendproc[iswap] != me) {
      MPI_Irecv(buf_recv,nrecv,MPI_DOUBLE,recvproc[iswap],0,
                world,&request);
      MPI_Send(buf_send,nsend,MPI_DOUBLE,sendproc[iswap],0,world);
      MPI_Wait(&request,&status);
      buf = buf_recv;
    } else buf = buf_send;

    // unpack buffer

    pair->unpack_forward_comm(recvnum[iswap],firstrecv[iswap],buf);
  }
}

/* ----------------------------------------------------------------------
   reverse communication of forces on atoms every timestep
   other per-atom attributes may also be sent via pack/unpack routines
------------------------------------------------------------------------- */

void CommCuda::reverse_comm()
{
  int n;
  MPI_Request request;
  MPI_Status status;
  AtomVec *avec = atom->avec;
  double **f = atom->f;
  double *buf;

  if(not comm_f_only && not avec->cudable) cuda->downloadAll();  //not yet implemented in CUDA but only needed for non standard atom styles

  // exchange data with another proc
  // if other proc is self, just copy
  // if comm_f_only set, exchange or copy directly from f, don't pack

  for (int iswap = nswap-1; iswap >= 0; iswap--) {
    if (sendproc[iswap] != me) {
      if (comm_f_only) {

    int size_recv_now=size_reverse_recv[iswap];
        if((sizeof(F_CFLOAT)!=sizeof(double))&& size_reverse_recv[iswap])
          size_recv_now=(size_recv_now+1)*sizeof(F_CFLOAT)/sizeof(double);
        MPI_Irecv(buf_recv,size_recv_now,MPI_DOUBLE,
                  sendproc[iswap],0,world,&request);

    buf=buf_send;
    if (size_reverse_send[iswap])
    {
      Cuda_CommCuda_PackReverse(&cuda->shared_data,size_reverse_send[iswap]/3,firstrecv[iswap],buf);
    }
    else buf=NULL;
    int size_reverse_send_now=size_reverse_send[iswap];
        if((sizeof(F_CFLOAT)!=sizeof(double))&& size_reverse_send[iswap])
          size_reverse_send_now=(size_reverse_send_now+1)*sizeof(F_CFLOAT)/sizeof(double);
        MPI_Send(buf,size_reverse_send_now,MPI_DOUBLE,
                 recvproc[iswap],0,world);
        MPI_Wait(&request,&status);
        Cuda_CommCuda_UnpackReverse(&cuda->shared_data,sendnum[iswap],iswap,buf_recv);

      } else {
        MPI_Irecv(buf_recv,size_reverse_recv[iswap],MPI_DOUBLE,
                  sendproc[iswap],0,world,&request);
        n = avec->pack_reverse(recvnum[iswap],firstrecv[iswap],buf_send);
        MPI_Send(buf_send,n,MPI_DOUBLE,recvproc[iswap],0,world);
        MPI_Wait(&request,&status);

      avec->unpack_reverse(sendnum[iswap],sendlist[iswap],buf_recv);
      }

    } else {
      if (comm_f_only) {
        if (sendnum[iswap])
              Cuda_CommCuda_UnpackReverse_Self(&cuda->shared_data,sendnum[iswap],iswap,firstrecv[iswap]);
      } else {
        n = avec->pack_reverse(recvnum[iswap],firstrecv[iswap],buf_send);
        avec->unpack_reverse(sendnum[iswap],sendlist[iswap],buf_send);
      }
    }
  }
  if(not comm_f_only && not avec->cudable) cuda->uploadAll();  //not yet implemented in CUDA but only needed for non standard atom styles
}

/* ----------------------------------------------------------------------
   exchange: move atoms to correct processors
   atoms exchanged with all 6 stencil neighbors
   send out atoms that have left my box, receive ones entering my box
   atoms will be lost if not inside some proc's box
     can happen if atom moves outside of non-periodic bounary
     or if atom moves more than one proc away
   this routine called before every reneighboring
   for triclinic, atoms must be in lamda coords (0-1) before exchange is called
------------------------------------------------------------------------- */

void CommCuda::exchange()
{
  AtomVec *avec = atom->avec;

  if(not cuda->oncpu && avec->cudable)
            return exchange_cuda();

  if(not cuda->oncpu) cuda->downloadAll();

  CommBrick::exchange();
}


void CommCuda::exchange_cuda()
{
  int i,m,nsend,nrecv,nrecv1,nrecv2,nlocal;
  double lo,hi,value;
  double **x;
  double *sublo,*subhi,*buf;
  MPI_Request request;
  MPI_Status status;
  AtomVec *avec = atom->avec;
    my_times time1,time2,time3;

  // clear global->local map for owned and ghost atoms
  // b/c atoms migrate to new procs in exchange() and
  // new ghosts are created in borders()
  // map_set() is done at end of borders()


  if(map_style) cuda->cu_tag->download();

  if (map_style) atom->map_clear();

  // subbox bounds for orthogonal or triclinic

  if (triclinic == 0) {
    sublo = domain->sublo;
    subhi = domain->subhi;
  } else {
    sublo = domain->sublo_lamda;
    subhi = domain->subhi_lamda;
  }

  // loop over dimensions

  for (int dim = 0; dim < 3; dim++) {
    // fill buffer with atoms leaving my box, using < and >=
    // when atom is deleted, fill it in with last atom

          cuda->shared_data.exchange_dim=dim;

    nlocal = atom->nlocal;
    avec->maxsend=&maxsend;
    nsend=avec->pack_exchange(dim,(double*) &buf_send);
    nlocal = atom->nlocal;


    atom->nlocal = nlocal;

    // send/recv atoms in both directions
    // if 1 proc in dimension, no send/recv, set recv buf to send buf
    // if 2 procs in dimension, single send/recv
    // if more than 2 procs in dimension, send/recv to both neighbors

 my_gettime(CLOCK_REALTIME,&time1);

    if (procgrid[dim] == 1) {
      nrecv = nsend;
      buf = buf_send;

    } else {
      MPI_Sendrecv(&nsend,1,MPI_INT,procneigh[dim][0],0,
                   &nrecv1,1,MPI_INT,procneigh[dim][1],0,world,&status);
      nrecv = nrecv1;
      if (procgrid[dim] > 2) {
        MPI_Sendrecv(&nsend,1,MPI_INT,procneigh[dim][1],0,
                     &nrecv2,1,MPI_INT,procneigh[dim][0],0,world,&status);
        nrecv += nrecv2;
      }
      if (nrecv+1 > maxrecv) grow_recv(nrecv+1);

      MPI_Irecv(buf_recv,nrecv1,MPI_DOUBLE,procneigh[dim][1],0,
                world,&request);
      MPI_Send(buf_send,nsend,MPI_DOUBLE,procneigh[dim][0],0,world);
      MPI_Wait(&request,&status);

      if (procgrid[dim] > 2) {
        MPI_Irecv(&buf_recv[nrecv1],nrecv2,MPI_DOUBLE,procneigh[dim][0],0,
                  world,&request);
        MPI_Send(buf_send,nsend,MPI_DOUBLE,procneigh[dim][1],0,world);
        MPI_Wait(&request,&status);

            if((nrecv1==0)||(nrecv2==0)) buf_recv[nrecv]=0;
      }

      buf = buf_recv;
    }
        //printf("nsend: %i nrecv: %i\n",nsend,nrecv);
    // check incoming atoms to see if they are in my box
    // if so, add to my list
my_gettime(CLOCK_REALTIME,&time2);
cuda->shared_data.cuda_timings.comm_exchange_mpi+=
      time2.tv_sec-time1.tv_sec+1.0*(time2.tv_nsec-time1.tv_nsec)/1000000000;

    if(nrecv)
    {
      avec->maxsend=&maxsend;
      avec->unpack_exchange(buf);
    }
  }

  if(atom->firstgroupname) cuda->downloadAll();

  if(atom->firstgroupname) atom->first_reorder();

  if(atom->firstgroupname) cuda->uploadAll();
}

/* ----------------------------------------------------------------------
   borders: list nearby atoms to send to neighboring procs at every timestep
   one list is created for every swap that will be made
   as list is made, actually do swaps
   this does equivalent of a communicate (so don't need to explicitly
     call communicate routine on reneighboring timestep)
   this routine is called before every reneighboring
   for triclinic, atoms must be in lamda coords (0-1) before borders is called
------------------------------------------------------------------------- */


void CommCuda::borders()
{
  AtomVec *avec = atom->avec;
  if(not cuda->oncpu && avec->cudable)
  {
          if(cuda->shared_data.overlap_comm&&cuda->finished_setup)
             borders_cuda_overlap_forward_comm();
           else
             borders_cuda();

           return;
  }

  CommBrick::borders();

  cuda->setSystemParams();
  if(cuda->finished_setup) {cuda->checkResize(); cuda->uploadAll();}
  cuda->shared_data.atom.nghost=atom->nghost;
  cu_sendlist->upload();
}

void CommCuda::borders_cuda()
{
  int i,n,itype,iswap,dim,ineed,twoneed,smax,rmax;
  int nsend,nrecv,nfirst,nlast,ngroup;
  double lo,hi;
  int *type;
  double **x;
  double *buf,*mlo,*mhi;
  MPI_Request request;
  MPI_Status status;
  AtomVec *avec = atom->avec;
    my_times time1,time2,time3;

  // clear old ghosts

  atom->nghost = 0;

  // do swaps over all 3 dimensions

  iswap = 0;
  smax = rmax = 0;

  cuda->shared_data.comm.nsend=0;
  for (dim = 0; dim < 3; dim++) {
    nlast = 0;
    twoneed = 2*maxneed[dim];
    for (ineed = 0; ineed < twoneed; ineed++) {

      // find atoms within slab boundaries lo/hi using <= and >=
      // check atoms between nfirst and nlast
      //   for first swaps in a dim, check owned and ghost
      //   for later swaps in a dim, only check newly arrived ghosts
      // store sent atom indices in list for use in future timesteps

      x = atom->x;
      if (style == SINGLE) {
        lo = slablo[iswap];
        hi = slabhi[iswap];
      } else {
        type = atom->type;
        mlo = multilo[iswap];
        mhi = multihi[iswap];
      }
      if (ineed % 2 == 0) {
        nfirst = nlast;
        nlast = atom->nlocal + atom->nghost;
      }

      nsend = 0;

      // find send atoms according to SINGLE vs MULTI
      // all atoms eligible versus atoms in bordergroup
      // only need to limit loop to bordergroup for first sends (ineed < 2)
      // on these sends, break loop in two: owned (in group) and ghost
     do
     {
       if(nsend>=maxsendlist[iswap]) grow_list(iswap,static_cast <int> (nsend*1.05));
               nsend=Cuda_CommCuda_BuildSendlist(&cuda->shared_data,bordergroup,ineed,style==SINGLE?1:0,atom->nfirst,nfirst,nlast,dim,iswap);
     }while(nsend>=maxsendlist[iswap]);
      // pack up list of border atoms

      if (nsend*size_border > maxsend)
        grow_send(nsend*size_border,0);

      if (ghost_velocity)
        n = avec->pack_border_vel(nsend,&iswap,buf_send,
                           pbc_flag[iswap],pbc[iswap]);
      else
        n = avec->pack_border(nsend,&iswap,buf_send,
                           pbc_flag[iswap],pbc[iswap]);

      // swap atoms with other proc
      // put incoming ghosts at end of my atom arrays
      // if swapping with self, simply copy, no messages

my_gettime(CLOCK_REALTIME,&time1);
      if (sendproc[iswap] != me) {
        MPI_Sendrecv(&nsend,1,MPI_INT,sendproc[iswap],0,
                     &nrecv,1,MPI_INT,recvproc[iswap],0,world,&status);
        if (nrecv*size_border > maxrecv)
          grow_recv(nrecv*size_border);
        MPI_Irecv(buf_recv,nrecv*size_border,MPI_DOUBLE,
                  recvproc[iswap],0,world,&request);
        MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap],0,world);
        MPI_Wait(&request,&status);
        buf = buf_recv;
      } else {
        nrecv = nsend;
        buf = buf_send;
      }

my_gettime(CLOCK_REALTIME,&time2);
cuda->shared_data.cuda_timings.comm_border_mpi+=
      time2.tv_sec-time1.tv_sec+1.0*(time2.tv_nsec-time1.tv_nsec)/1000000000;

      // unpack buffer

      if (ghost_velocity)
        avec->unpack_border_vel(nrecv,atom->nlocal+atom->nghost,buf);
      else
        avec->unpack_border(nrecv,atom->nlocal+atom->nghost,buf);

      // set all pointers & counters

      smax = MAX(smax,nsend);
      rmax = MAX(rmax,nrecv);
      sendnum[iswap] = nsend;
      recvnum[iswap] = nrecv;
      size_forward_recv[iswap] = nrecv*size_forward;
      size_reverse_send[iswap] = nrecv*size_reverse;
      size_reverse_recv[iswap] = nsend*size_reverse;
      firstrecv[iswap] = atom->nlocal + atom->nghost;
      atom->nghost += nrecv;
      iswap++;
    }
  }

  // insure send/recv buffers are long enough for all forward & reverse comm

  int max = MAX(maxforward*smax,maxreverse*rmax);
  if (max > maxsend) grow_send(max,0);
  max = MAX(maxforward*rmax,maxreverse*smax);
  if (max > maxrecv) grow_recv(max);

  // reset global->local map
  if(map_style)
  {
          cuda->cu_tag->download();
         atom->map_set();
  }

  cuda->setSystemParams();
  cuda->shared_data.atom.nghost+=n;
}

void CommCuda::borders_cuda_overlap_forward_comm()
{
  int i,n,itype,iswap,dim,ineed,twoneed,smax,rmax;
  int nsend,nrecv,nfirst,nlast,ngroup;
  double lo,hi;
  int *type;
  double **x;
  double *buf,*mlo,*mhi;
  MPI_Request request;
  MPI_Status status;
  AtomVec *avec = atom->avec;
    my_times time1,time2,time3;

  // clear old ghosts

  atom->nghost = 0;

  // do swaps over all 3 dimensions

  iswap = 0;
  smax = rmax = 0;

  cuda->shared_data.comm.nsend=0;
  for (dim = 0; dim < 3; dim++) {
    nlast = 0;
    twoneed = 2*maxneed[dim];
    for (ineed = 0; ineed < twoneed; ineed++) {

      // find atoms within slab boundaries lo/hi using <= and >=
      // check atoms between nfirst and nlast
      //   for first swaps in a dim, check owned and ghost
      //   for later swaps in a dim, only check newly arrived ghosts
      // store sent atom indices in list for use in future timesteps

      x = atom->x;
      if (style == SINGLE) {
        lo = slablo[iswap];
        hi = slabhi[iswap];
      } else {
        type = atom->type;
        mlo = multilo[iswap];
        mhi = multihi[iswap];
      }
      if (ineed % 2 == 0) {
        nfirst = nlast;
        nlast = atom->nlocal + atom->nghost;
      }

      nsend = 0;

      // find send atoms according to SINGLE vs MULTI
      // all atoms eligible versus atoms in bordergroup
      // only need to limit loop to bordergroup for first sends (ineed < 2)
      // on these sends, break loop in two: owned (in group) and ghost
     do
     {
       if(nsend>=maxsendlist[iswap]) grow_list(iswap,static_cast <int> (nsend*1.05));
               nsend=Cuda_CommCuda_BuildSendlist(&cuda->shared_data,bordergroup,ineed,style==SINGLE?1:0,atom->nfirst,nfirst,nlast,dim,iswap);
     }while(nsend>=maxsendlist[iswap]);
         cuda->shared_data.comm.nsend_swap[iswap]=nsend;
          // pack up list of border atoms

      if (nsend*size_border > maxsend)
        grow_send(nsend*size_border,0);

      if (ghost_velocity)
        n = avec->pack_border_vel(nsend,&iswap,buf_send,
                           pbc_flag[iswap],pbc[iswap]);
      else
        n = avec->pack_border(nsend,&iswap,buf_send,
                           pbc_flag[iswap],pbc[iswap]);

      // swap atoms with other proc
      // put incoming ghosts at end of my atom arrays
      // if swapping with self, simply copy, no messages

my_gettime(CLOCK_REALTIME,&time1);
      if (sendproc[iswap] != me) {
        MPI_Sendrecv(&nsend,1,MPI_INT,sendproc[iswap],0,
                     &nrecv,1,MPI_INT,recvproc[iswap],0,world,&status);
        if (nrecv*size_border > maxrecv)
          grow_recv(nrecv*size_border);
        MPI_Irecv(buf_recv,nrecv*size_border,MPI_DOUBLE,
                  recvproc[iswap],0,world,&request);
        MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap],0,world);
        MPI_Wait(&request,&status);
        buf = buf_recv;
      } else {
        nrecv = nsend;
        buf = buf_send;
      }

my_gettime(CLOCK_REALTIME,&time2);
cuda->shared_data.cuda_timings.comm_border_mpi+=
      time2.tv_sec-time1.tv_sec+1.0*(time2.tv_nsec-time1.tv_nsec)/1000000000;

      // unpack buffer

      if (ghost_velocity)
        avec->unpack_border_vel(nrecv,atom->nlocal+atom->nghost,buf);
      else
        avec->unpack_border(nrecv,atom->nlocal+atom->nghost,buf);

      // set all pointers & counters

      smax = MAX(smax,nsend);
      rmax = MAX(rmax,nrecv);
      sendnum[iswap] = nsend;
      recvnum[iswap] = nrecv;
      size_forward_recv[iswap] = nrecv*size_forward;
      size_reverse_send[iswap] = nrecv*size_reverse;
      size_reverse_recv[iswap] = nsend*size_reverse;
      firstrecv[iswap] = atom->nlocal + atom->nghost;
      atom->nghost += nrecv;
      iswap++;
    }
  }

  // insure send/recv buffers are long enough for all forward & reverse comm

  int max = MAX(maxforward*smax,maxreverse*rmax);
  if (max > maxsend) grow_send(max,0);
  max = MAX(maxforward*rmax,maxreverse*smax);
  if (max > maxrecv) grow_recv(max);

  // reset global->local map
  if(map_style)
  {
          cuda->cu_tag->download();
         atom->map_set();
  }

  cuda->setSystemParams();
  cuda->shared_data.atom.nghost+=n;
}




void CommCuda::forward_comm_fix(Fix *fix)
{
  int iswap,n;
  double *buf;
  MPI_Request request;
  MPI_Status status;

  int nsize = fix->comm_forward;

  for (iswap = 0; iswap < nswap; iswap++) {
    // pack buffer
    if(fix->cudable_comm&&cuda->finished_setup)
    {
            int swap=iswap;
        if(sendproc[iswap] == me) {swap=-iswap-1; buf=(double*)&(firstrecv[iswap]);}
        else buf=buf_send;

        n = fix->pack_forward_comm(sendnum[iswap],&swap,
                                   buf,pbc_flag[iswap],pbc[iswap]);
        if(sendproc[iswap] == me)
        {
                continue;
        }
    }
    else
      n = fix->pack_forward_comm(sendnum[iswap],sendlist[iswap],
                                 buf_send,pbc_flag[iswap],pbc[iswap]);

     // exchange with another proc
    // if self, set recv buffer to send buffer

    if (sendproc[iswap] != me) {
      MPI_Irecv(buf_recv,nsize*recvnum[iswap],MPI_DOUBLE,recvproc[iswap],0,
                world,&request);
      MPI_Send(buf_send,n,MPI_DOUBLE,sendproc[iswap],0,world);
      MPI_Wait(&request,&status);
      buf = buf_recv;
    } else buf = buf_send;

    // unpack buffer

    fix->unpack_forward_comm(recvnum[iswap],firstrecv[iswap],buf);
  }
}


void CommCuda::grow_send(int n, int flag)
{
  int oldmaxsend = (maxsend+BUFEXTRA)*sizeof(double);
  maxsend = static_cast<int> (BUFFACTOR * n);
  if (flag){
    if(cuda->pinned)
    {
      double* tmp = new double[oldmaxsend];
      memcpy((void*) tmp,(void*) buf_send,oldmaxsend*sizeof(double));
      if(buf_send) CudaWrapper_FreePinnedHostData((void*) (buf_send));
      buf_send = (double*) CudaWrapper_AllocPinnedHostData((maxsend+BUFEXTRA)*sizeof(double),false);
      memcpy(buf_send,tmp,oldmaxsend*sizeof(double));
      delete [] tmp;
    }
    else
    {
    buf_send = (double *)
      memory->srealloc(buf_send,(maxsend+BUFEXTRA)*sizeof(double),
                       "comm:buf_send");printf("srealloc\n");
    }
  }
  else {
    if(cuda->pinned)
    {
      if(buf_send) CudaWrapper_FreePinnedHostData((void*) buf_send);
      buf_send = (double*) CudaWrapper_AllocPinnedHostData((maxsend+BUFEXTRA)*sizeof(double),false);
    }
    else
    {
      memory->sfree(buf_send);
      buf_send = (double *) memory->smalloc((maxsend+BUFEXTRA)*sizeof(double),
                                          "comm:buf_send");
    }
    for(int i=0;i<maxswap;i++)
    {
      if(cuda->shared_data.comm.buf_send_dev[i]) CudaWrapper_FreeCudaData(cuda->shared_data.comm.buf_send_dev[i],oldmaxsend);
      cuda->shared_data.comm.buf_send_dev[i]=CudaWrapper_AllocCudaData((maxsend+BUFEXTRA)*sizeof(double));
    }
  }
}
/* ----------------------------------------------------------------------
   free/malloc the size of the recv buffer as needed with BUFFACTOR
------------------------------------------------------------------------- */


void CommCuda::grow_recv(int n)
{
  int oldmaxrecv = maxrecv*sizeof(double);
  maxrecv = static_cast<int> (BUFFACTOR * n);
  if(cuda->pinned)
  {
    if(buf_recv) CudaWrapper_FreePinnedHostData((void*)buf_recv);
    buf_recv = (double*) CudaWrapper_AllocPinnedHostData(maxrecv*sizeof(double), false,true);
  }
  else
  {
    memory->sfree(buf_recv);
    buf_recv = (double *) memory->smalloc(maxrecv*sizeof(double),
                                        "comm:buf_recv");
  }
  for(int i=0;i<maxswap;i++)
  {
    if(cuda->shared_data.comm.buf_recv_dev[i]) CudaWrapper_FreeCudaData(cuda->shared_data.comm.buf_recv_dev[i],oldmaxrecv);
    cuda->shared_data.comm.buf_recv_dev[i]=CudaWrapper_AllocCudaData((maxrecv)*sizeof(double));
  }
}

/* ----------------------------------------------------------------------
   realloc the size of the iswap sendlist as needed with BUFFACTOR
------------------------------------------------------------------------- */

void CommCuda::grow_list(int iswap, int n)
{

  MYDBG(printf(" # CUDA CommCuda::grow_list\n");)
  if(cuda->finished_setup&&cu_sendlist) cu_sendlist->download();
  if(!cu_sendlist||n*BUFFACTOR>cu_sendlist->get_dim()[1]||n*BUFFACTOR>maxsendlist[iswap])
  {
          for(int i=0;i<maxswap;i++)
          {
            maxsendlist[i] = static_cast<int> (BUFFACTOR * n);
            sendlist[i] = (int *)
                    memory->srealloc(sendlist[i],maxsendlist[i]*sizeof(int),
                                     "comm:sendlist[iswap]");
          }
          delete cu_sendlist;
          cu_sendlist=new cCudaData<int, int, xy> ((int*)sendlist,maxswap,maxsendlist[iswap]);
          cuda->shared_data.comm.sendlist.dev_data=cu_sendlist->dev_data();
    cuda->shared_data.comm.maxlistlength=maxsendlist[iswap];
    cu_sendlist->upload();
  }
 }

/* ----------------------------------------------------------------------
   realloc the buffers needed for swaps
------------------------------------------------------------------------- */

void CommCuda::grow_swap(int n)
{
  int oldmaxswap=maxswap;
  CommBrick::grow_swap(n);
  if(n>cu_sendlist->get_dim()[0])
  {
   MYDBG(printf(" # CUDA CommCuda::grow_swap\n");)

          delete cu_sendlist;
          cu_sendlist=new cCudaData<int, int, xy> ((int*)sendlist,n,BUFMIN);
          cuda->shared_data.comm.sendlist.dev_data=cu_sendlist->dev_data();
    cuda->shared_data.comm.maxlistlength=BUFMIN;
    cuda->shared_data.comm.maxswap=n;
    cuda->shared_data.comm.nsend_swap=new int[n];
    cuda->shared_data.comm.send_size=new int[n];
    cuda->shared_data.comm.recv_size=new int[n];
  }
  for(int i=0;i<oldmaxswap;i++)
  {
    if(cuda->shared_data.comm.buf_recv_dev[i]) CudaWrapper_FreeCudaData(cuda->shared_data.comm.buf_recv_dev[i],maxrecv*sizeof(double));
    if(cuda->shared_data.comm.buf_send_dev[i]) CudaWrapper_FreeCudaData(cuda->shared_data.comm.buf_send_dev[i],maxsend*sizeof(double));
    cuda->shared_data.comm.buf_recv_dev[i]=NULL;
    cuda->shared_data.comm.buf_send_dev[i]=NULL;
  }
  cuda->shared_data.comm.buf_send= new double*[n];
  cuda->shared_data.comm.buf_recv= new double*[n];
  cuda->shared_data.comm.buf_send_dev= new void*[n];
  cuda->shared_data.comm.buf_recv_dev= new void*[n];
  for(int i=0;i<n;i++)
  {
    cuda->shared_data.comm.buf_recv[i]=NULL;
    cuda->shared_data.comm.buf_send[i]=NULL;
    cuda->shared_data.comm.buf_recv_dev[i]=NULL;
    cuda->shared_data.comm.buf_send_dev[i]=NULL;
  }
  grow_send(maxsend,0);
  grow_recv(maxrecv);

  maxswap=n;
}

/* ----------------------------------------------------------------------
   allocation of swap info
------------------------------------------------------------------------- */

void CommCuda::allocate_swap(int n)
{
   CommBrick::allocate_swap(n);

          delete cu_pbc;
          delete cu_slablo;
          delete cu_slabhi;

    cuda->shared_data.comm.maxswap=n;
          if(cu_sendlist)
          {
            cu_pbc=new cCudaData<int, int, xy> ((int*)pbc,n,6);
            cu_slablo = new cCudaData<double, X_CFLOAT,x>(slablo,n);
            cu_slabhi = new cCudaData<double, X_CFLOAT,x>(slabhi,n);

            cuda->shared_data.comm.pbc.dev_data=cu_pbc->dev_data();
            cuda->shared_data.comm.slablo.dev_data=cu_slablo->dev_data();
            cuda->shared_data.comm.slabhi.dev_data=cu_slabhi->dev_data();
          }
    cuda->shared_data.comm.nsend_swap=new int[n];
    cuda->shared_data.comm.send_size=new int[n];
    cuda->shared_data.comm.recv_size=new int[n];
    cuda->shared_data.comm.buf_send= new double*[n];
    cuda->shared_data.comm.buf_recv= new double*[n];
    cuda->shared_data.comm.buf_send_dev= new void*[n];
    cuda->shared_data.comm.buf_recv_dev= new void*[n];
    for(int i=0;i<n;i++) cuda->shared_data.comm.buf_send_dev[i]=NULL;
    for(int i=0;i<n;i++) cuda->shared_data.comm.buf_recv_dev[i]=NULL;
}


/* ----------------------------------------------------------------------
   allocation of multi-type swap info
------------------------------------------------------------------------- */

void CommCuda::allocate_multi(int n)
{
  CommBrick::allocate_multi(n);

          delete cu_multilo;
          delete cu_multihi;
          cu_multilo = new cCudaData<double, X_CFLOAT,xy>(slablo,n,atom->ntypes+1);
          cu_multihi = new cCudaData<double, X_CFLOAT,xy>(slabhi,n,atom->ntypes+1);

          cuda->shared_data.comm.multilo.dev_data=cu_multilo->dev_data();
          cuda->shared_data.comm.multihi.dev_data=cu_multihi->dev_data();
}

/* ----------------------------------------------------------------------
   free memory for swaps
------------------------------------------------------------------------- */

void CommCuda::free_swap()
{

  CommBrick::free_swap();

  delete cuda->shared_data.comm.nsend_swap; cuda->shared_data.comm.nsend_swap=NULL;
  delete cu_pbc; cu_pbc = NULL;
  delete cu_slablo; cu_slablo = NULL;
  delete cu_slabhi; cu_slabhi = NULL;
  for(int i=0;i<maxswap;i++)
  {
    if(cuda->shared_data.comm.buf_recv_dev[i]) CudaWrapper_FreeCudaData(cuda->shared_data.comm.buf_recv_dev[i],maxrecv*sizeof(double));
    if(cuda->shared_data.comm.buf_send_dev[i]) CudaWrapper_FreeCudaData(cuda->shared_data.comm.buf_send_dev[i],maxsend*sizeof(double));
  }


}

/* ----------------------------------------------------------------------
   free memory for multi-type swaps
------------------------------------------------------------------------- */

void CommCuda::free_multi()
{
  CommBrick::free_multi();
  delete cu_multilo; cu_multilo = NULL;
  delete cu_multihi; cu_multihi = NULL;
}
