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

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "remap_kokkos.h"
#include "pack_kokkos.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MIN(A,B) ((A) < (B) ? (A) : (B))
#define MAX(A,B) ((A) > (B) ? (A) : (B))

/* ---------------------------------------------------------------------- */

template<class DeviceType>
RemapKokkos<DeviceType>::RemapKokkos(LAMMPS *lmp) : Pointers(lmp)
{
  plan = NULL;
}

template<class DeviceType>
RemapKokkos<DeviceType>::RemapKokkos(LAMMPS *lmp, MPI_Comm comm,
             int in_ilo, int in_ihi, int in_jlo, int in_jhi,
             int in_klo, int in_khi,
             int out_ilo, int out_ihi, int out_jlo, int out_jhi,
             int out_klo, int out_khi,
             int nqty, int permute, int memory,
             int precision, int usecollective) : Pointers(lmp)
{
  plan = remap_3d_create_plan_kokkos(comm,
                              in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,
                              out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi,
                              nqty,permute,memory,precision,usecollective);
  if (plan == NULL) error->one(FLERR,"Could not create 3d remap plan");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
RemapKokkos<DeviceType>::~RemapKokkos()
{
  remap_3d_destroy_plan_kokkos(plan);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void RemapKokkos<DeviceType>::perform(typename AT::t_FFT_SCALAR_1d d_in, typename AT::t_FFT_SCALAR_1d d_out, typename AT::t_FFT_SCALAR_1d d_buf)
{
  remap_3d_kokkos(d_in,d_out,d_buf,plan);
}


/* ----------------------------------------------------------------------
   Data layout for 3d remaps:

   data set of Nfast x Nmid x Nslow elements is owned by P procs
   each element = nqty contiguous datums
   on input, each proc owns a subsection of the elements
   on output, each proc will own a (presumably different) subsection
   my subsection must not overlap with any other proc's subsection,
     i.e. the union of all proc's input (or output) subsections must
     exactly tile the global Nfast x Nmid x Nslow data set
   when called from C, all subsection indices are
     C-style from 0 to N-1 where N = Nfast or Nmid or Nslow
   when called from F77, all subsection indices are
     F77-style from 1 to N where N = Nfast or Nmid or Nslow
   a proc can own 0 elements on input or output
     by specifying hi index < lo index
   on both input and output, data is stored contiguously on a processor
     with a fast-varying, mid-varying, and slow-varying index
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Perform 3d remap

   Arguments:
   in           starting address of input data on this proc
   out          starting address of where output data for this proc
                  will be placed (can be same as in)
   buf          extra memory required for remap
                if memory=0 was used in call to remap_3d_create_plan
                  then buf must be big enough to hold output result
                  i.e. nqty * (out_ihi-out_ilo+1) * (out_jhi-out_jlo+1) *
                              (out_khi-out_klo+1)
                if memory=1 was used in call to remap_3d_create_plan
                  then buf is not used, can just be a dummy pointer
   plan         plan returned by previous call to remap_3d_create_plan
------------------------------------------------------------------------- */

template<class DeviceType>
void RemapKokkos<DeviceType>::remap_3d_kokkos(typename AT::t_FFT_SCALAR_1d d_in, typename AT::t_FFT_SCALAR_1d d_out, typename AT::t_FFT_SCALAR_1d d_buf,
              struct remap_plan_3d_kokkos<DeviceType> *plan)
{
  // collective flag not yet supported

  // use point-to-point communication

  int i,isend,irecv;
  typename AT::t_FFT_SCALAR_1d d_scratch;

  if (plan->memory == 0)
    d_scratch = d_buf;
  else
    d_scratch = plan->d_scratch;

  // post all recvs into scratch space

  for (irecv = 0; irecv < plan->nrecv; irecv++) {
    FFT_SCALAR* scratch = d_scratch.ptr_on_device() + plan->recv_bufloc[irecv];
    MPI_Irecv(scratch,plan->recv_size[irecv],
              MPI_FFT_SCALAR,plan->recv_proc[irecv],0,
              plan->comm,&plan->request[irecv]);
  }

  // send all messages to other procs

  for (isend = 0; isend < plan->nsend; isend++) {
    int in_offset = plan->send_offset[isend];
    plan->pack(d_in,in_offset,
               plan->d_sendbuf,0,&plan->packplan[isend]);
    MPI_Send(plan->d_sendbuf.ptr_on_device(),plan->send_size[isend],MPI_FFT_SCALAR,
             plan->send_proc[isend],0,plan->comm);
  }

  // copy in -> scratch -> out for self data

  if (plan->self) {
    isend = plan->nsend;
    irecv = plan->nrecv;

    int in_offset = plan->send_offset[isend];
    int scratch_offset = plan->recv_bufloc[irecv];
    int out_offset = plan->recv_offset[irecv];

    plan->pack(d_in,in_offset,
               d_scratch,scratch_offset,
               &plan->packplan[isend]);
    plan->unpack(d_scratch,scratch_offset,
                 d_out,out_offset,&plan->unpackplan[irecv]);
  }

  // unpack all messages from scratch -> out

  for (i = 0; i < plan->nrecv; i++) {
    MPI_Waitany(plan->nrecv,plan->request,&irecv,MPI_STATUS_IGNORE);

    int scratch_offset = plan->recv_bufloc[irecv];
    int out_offset = plan->recv_offset[irecv];

    plan->unpack(d_scratch,scratch_offset,
                 d_out,out_offset,&plan->unpackplan[irecv]);
  }
}

/* ----------------------------------------------------------------------
   Create plan for performing a 3d remap

   Arguments:
   comm                 MPI communicator for the P procs which own the data
   in_ilo,in_ihi        input bounds of data I own in fast index
   in_jlo,in_jhi        input bounds of data I own in mid index
   in_klo,in_khi        input bounds of data I own in slow index
   out_ilo,out_ihi      output bounds of data I own in fast index
   out_jlo,out_jhi      output bounds of data I own in mid index
   out_klo,out_khi      output bounds of data I own in slow index
   nqty                 # of datums per element
   permute              permutation in storage order of indices on output
                          0 = no permutation
                          1 = permute once = mid->fast, slow->mid, fast->slow
                          2 = permute twice = slow->fast, fast->mid, mid->slow
   memory               user provides buffer memory for remap or system does
                          0 = user provides memory
                          1 = system provides memory
   precision            precision of data
                          1 = single precision (4 bytes per datum)
                          2 = double precision (8 bytes per datum)
   usecollective        whether to use collective MPI or point-to-point
------------------------------------------------------------------------- */

template<class DeviceType>
struct remap_plan_3d_kokkos<DeviceType>* RemapKokkos<DeviceType>::remap_3d_create_plan_kokkos(
  MPI_Comm comm,
  int in_ilo, int in_ihi, int in_jlo, int in_jhi,
  int in_klo, int in_khi,
  int out_ilo, int out_ihi, int out_jlo, int out_jhi,
  int out_klo, int out_khi,
  int nqty, int permute, int memory, int precision, int usecollective)
{

  struct remap_plan_3d_kokkos<DeviceType> *plan;
  struct extent_3d *inarray, *outarray;
  struct extent_3d in,out,overlap;
  int i,iproc,nsend,nrecv,ibuf,size,me,nprocs;

  // query MPI info

  MPI_Comm_rank(comm,&me);
  MPI_Comm_size(comm,&nprocs);

  // allocate memory for plan data struct

  plan = new struct remap_plan_3d_kokkos<DeviceType>;
  if (plan == NULL) return NULL;
  plan->usecollective = usecollective;

  // store parameters in local data structs

  in.ilo = in_ilo;
  in.ihi = in_ihi;
  in.isize = in.ihi - in.ilo + 1;

  in.jlo = in_jlo;
  in.jhi = in_jhi;
  in.jsize = in.jhi - in.jlo + 1;

  in.klo = in_klo;
  in.khi = in_khi;
  in.ksize = in.khi - in.klo + 1;

  out.ilo = out_ilo;
  out.ihi = out_ihi;
  out.isize = out.ihi - out.ilo + 1;

  out.jlo = out_jlo;
  out.jhi = out_jhi;
  out.jsize = out.jhi - out.jlo + 1;

  out.klo = out_klo;
  out.khi = out_khi;
  out.ksize = out.khi - out.klo + 1;

  // combine output extents across all procs

  inarray = (struct extent_3d *) malloc(nprocs*sizeof(struct extent_3d));
  if (inarray == NULL) return NULL;

  outarray = (struct extent_3d *) malloc(nprocs*sizeof(struct extent_3d));
  if (outarray == NULL) return NULL;

  MPI_Allgather(&out,sizeof(struct extent_3d),MPI_BYTE,
                outarray,sizeof(struct extent_3d),MPI_BYTE,comm);

  // count send collides, including self

  nsend = 0;
  iproc = me;
  for (i = 0; i < nprocs; i++) {
    iproc++;
    if (iproc == nprocs) iproc = 0;
    nsend += remap_3d_collide(&in,&outarray[iproc],&overlap);
  }

  // malloc space for send info

  if (nsend) {
    plan->pack = PackKokkos<DeviceType>::pack_3d;

    plan->send_offset = (int *) malloc(nsend*sizeof(int));
    plan->send_size = (int *) malloc(nsend*sizeof(int));
    plan->send_proc = (int *) malloc(nsend*sizeof(int));
    plan->packplan = (struct pack_plan_3d *)
      malloc(nsend*sizeof(struct pack_plan_3d));

    if (plan->send_offset == NULL || plan->send_size == NULL ||
        plan->send_proc == NULL || plan->packplan == NULL) return NULL;
  }

  // store send info, with self as last entry

  nsend = 0;
  iproc = me;
  for (i = 0; i < nprocs; i++) {
    iproc++;
    if (iproc == nprocs) iproc = 0;
    if (remap_3d_collide(&in,&outarray[iproc],&overlap)) {
      plan->send_proc[nsend] = iproc;
      plan->send_offset[nsend] = nqty *
        ((overlap.klo-in.klo)*in.jsize*in.isize +
         ((overlap.jlo-in.jlo)*in.isize + overlap.ilo-in.ilo));
      plan->packplan[nsend].nfast = nqty*overlap.isize;
      plan->packplan[nsend].nmid = overlap.jsize;
      plan->packplan[nsend].nslow = overlap.ksize;
      plan->packplan[nsend].nstride_line = nqty*in.isize;
      plan->packplan[nsend].nstride_plane = nqty*in.jsize*in.isize;
      plan->packplan[nsend].nqty = nqty;
      plan->send_size[nsend] = nqty*overlap.isize*overlap.jsize*overlap.ksize;
      nsend++;
    }
  }

  // plan->nsend = # of sends not including self

  if (nsend && plan->send_proc[nsend-1] == me) {
    if (plan->usecollective) // for collectives include self in nsend list
      plan->nsend = nsend;
    else
      plan->nsend = nsend - 1;
  } else
    plan->nsend = nsend;

  // combine input extents across all procs

  MPI_Allgather(&in,sizeof(struct extent_3d),MPI_BYTE,
                inarray,sizeof(struct extent_3d),MPI_BYTE,comm);

  // count recv collides, including self

  nrecv = 0;
  iproc = me;
  for (i = 0; i < nprocs; i++) {
    iproc++;
    if (iproc == nprocs) iproc = 0;
    nrecv += remap_3d_collide(&out,&inarray[iproc],&overlap);
  }

  // malloc space for recv info

  if (nrecv) {
    if (permute == 0)
      plan->unpack = PackKokkos<DeviceType>::unpack_3d;
    else if (permute == 1) {
      if (nqty == 1)
        plan->unpack = PackKokkos<DeviceType>::unpack_3d_permute1_1;
      else if (nqty == 2)
        plan->unpack = PackKokkos<DeviceType>::unpack_3d_permute1_2;
      else
        plan->unpack = PackKokkos<DeviceType>::unpack_3d_permute1_n;
    }
    else if (permute == 2) {
      if (nqty == 1)
        plan->unpack = PackKokkos<DeviceType>::unpack_3d_permute2_1;
      else if (nqty == 2)
        plan->unpack = PackKokkos<DeviceType>::unpack_3d_permute2_2;
      else
        plan->unpack = PackKokkos<DeviceType>::unpack_3d_permute2_n;
    }

    plan->recv_offset = (int *) malloc(nrecv*sizeof(int));
    plan->recv_size = (int *) malloc(nrecv*sizeof(int));
    plan->recv_proc = (int *) malloc(nrecv*sizeof(int));
    plan->recv_bufloc = (int *) malloc(nrecv*sizeof(int));
    plan->request = (MPI_Request *) malloc(nrecv*sizeof(MPI_Request));
    plan->unpackplan = (struct pack_plan_3d *)
      malloc(nrecv*sizeof(struct pack_plan_3d));

    if (plan->recv_offset == NULL || plan->recv_size == NULL ||
        plan->recv_proc == NULL || plan->recv_bufloc == NULL ||
        plan->request == NULL || plan->unpackplan == NULL) return NULL;
  }

  // store recv info, with self as last entry

  ibuf = 0;
  nrecv = 0;
  iproc = me;

  for (i = 0; i < nprocs; i++) {
    iproc++;
    if (iproc == nprocs) iproc = 0;
    if (remap_3d_collide(&out,&inarray[iproc],&overlap)) {
      plan->recv_proc[nrecv] = iproc;
      plan->recv_bufloc[nrecv] = ibuf;

      if (permute == 0) {
        plan->recv_offset[nrecv] = nqty *
          ((overlap.klo-out.klo)*out.jsize*out.isize +
           (overlap.jlo-out.jlo)*out.isize + (overlap.ilo-out.ilo));
        plan->unpackplan[nrecv].nfast = nqty*overlap.isize;
        plan->unpackplan[nrecv].nmid = overlap.jsize;
        plan->unpackplan[nrecv].nslow = overlap.ksize;
        plan->unpackplan[nrecv].nstride_line = nqty*out.isize;
        plan->unpackplan[nrecv].nstride_plane = nqty*out.jsize*out.isize;
        plan->unpackplan[nrecv].nqty = nqty;
      }
      else if (permute == 1) {
        plan->recv_offset[nrecv] = nqty *
          ((overlap.ilo-out.ilo)*out.ksize*out.jsize +
           (overlap.klo-out.klo)*out.jsize + (overlap.jlo-out.jlo));
        plan->unpackplan[nrecv].nfast = overlap.isize;
        plan->unpackplan[nrecv].nmid = overlap.jsize;
        plan->unpackplan[nrecv].nslow = overlap.ksize;
        plan->unpackplan[nrecv].nstride_line = nqty*out.jsize;
        plan->unpackplan[nrecv].nstride_plane = nqty*out.ksize*out.jsize;
        plan->unpackplan[nrecv].nqty = nqty;
      }
      else {
        plan->recv_offset[nrecv] = nqty *
          ((overlap.jlo-out.jlo)*out.isize*out.ksize +
           (overlap.ilo-out.ilo)*out.ksize + (overlap.klo-out.klo));
        plan->unpackplan[nrecv].nfast = overlap.isize;
        plan->unpackplan[nrecv].nmid = overlap.jsize;
        plan->unpackplan[nrecv].nslow = overlap.ksize;
        plan->unpackplan[nrecv].nstride_line = nqty*out.ksize;
        plan->unpackplan[nrecv].nstride_plane = nqty*out.isize*out.ksize;
        plan->unpackplan[nrecv].nqty = nqty;
      }

      plan->recv_size[nrecv] = nqty*overlap.isize*overlap.jsize*overlap.ksize;
      ibuf += plan->recv_size[nrecv];
      nrecv++;
    }
  }

  // plan->nrecv = # of recvs not including self
  // for collectives include self in the nsend list

  if (nrecv && plan->recv_proc[nrecv-1] == me) {
    if (plan->usecollective) plan->nrecv = nrecv;
    else plan->nrecv = nrecv - 1;
  } else plan->nrecv = nrecv;

  // init remaining fields in remap plan

  plan->memory = memory;

  if (nrecv == plan->nrecv) plan->self = 0;
  else plan->self = 1;

  // free locally malloced space

  free(inarray);
  free(outarray);

  // find biggest send message (not including self) and malloc space for it

  size = 0;
  for (nsend = 0; nsend < plan->nsend; nsend++)
    size = MAX(size,plan->send_size[nsend]);

  if (size) {
    plan->d_sendbuf = typename AT::t_FFT_SCALAR_1d("remap3d:sendbuf",size);
    if (!plan->d_sendbuf.data()) return NULL;
  }

  // if requested, allocate internal scratch space for recvs,
  // only need it if I will receive any data (including self)

  if (memory == 1) {
    if (nrecv > 0) {
      plan->d_scratch =
        typename AT::t_FFT_SCALAR_1d("remap3d:scratch",nqty*out.isize*out.jsize*out.ksize);
      if (!plan->d_scratch.data()) return NULL;
    }
  }

  // not using collective - dup comm

  MPI_Comm_dup(comm,&plan->comm);

  // return pointer to plan

  return plan;
}

/* ----------------------------------------------------------------------
   Destroy a 3d remap plan
------------------------------------------------------------------------- */

template<class DeviceType>
void RemapKokkos<DeviceType>::remap_3d_destroy_plan_kokkos(struct remap_plan_3d_kokkos<DeviceType> *plan)
{
  if (plan == NULL) return;

  // free MPI communicator

  if (!((plan->usecollective) && (plan->commringlen == 0)))
    MPI_Comm_free(&plan->comm);

  // free internal arrays

  if (plan->nsend || plan->self) {
    free(plan->send_offset);
    free(plan->send_size);
    free(plan->send_proc);
    free(plan->packplan);
  }

  if (plan->nrecv || plan->self) {
    free(plan->recv_offset);
    free(plan->recv_size);
    free(plan->recv_proc);
    free(plan->recv_bufloc);
    free(plan->request);
    free(plan->unpackplan);
  }

  // free plan itself

  delete plan;
}

namespace LAMMPS_NS {
template class RemapKokkos<LMPDeviceType>;
#ifdef KOKKOS_HAVE_CUDA
template class RemapKokkos<LMPHostType>;
#endif
}
