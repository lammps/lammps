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

#include "stdio.h"
#include "stdlib.h"
#include "remap.h"

#define PACK_DATA FFT_SCALAR

#include "pack.h"

#define MIN(A,B) ((A) < (B) ? (A) : (B))
#define MAX(A,B) ((A) > (B) ? (A) : (B))

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

void remap_3d(FFT_SCALAR *in, FFT_SCALAR *out, FFT_SCALAR *buf,
              struct remap_plan_3d *plan)
{
  // use point-to-point communication

  if (!plan->usecollective) { 
    MPI_Status status;
    int i,isend,irecv;
    FFT_SCALAR *scratch;

    if (plan->memory == 0)
      scratch = buf;
    else
      scratch = plan->scratch;

    // post all recvs into scratch space

    for (irecv = 0; irecv < plan->nrecv; irecv++)
      MPI_Irecv(&scratch[plan->recv_bufloc[irecv]],plan->recv_size[irecv],
                MPI_FFT_SCALAR,plan->recv_proc[irecv],0,
                plan->comm,&plan->request[irecv]);

    // send all messages to other procs

    for (isend = 0; isend < plan->nsend; isend++) {
      plan->pack(&in[plan->send_offset[isend]],
                 plan->sendbuf,&plan->packplan[isend]);
      MPI_Send(plan->sendbuf,plan->send_size[isend],MPI_FFT_SCALAR,
               plan->send_proc[isend],0,plan->comm);
    }

    // copy in -> scratch -> out for self data

    if (plan->self) {
      isend = plan->nsend;
      irecv = plan->nrecv;
      plan->pack(&in[plan->send_offset[isend]],
                 &scratch[plan->recv_bufloc[irecv]],
                 &plan->packplan[isend]);
      plan->unpack(&scratch[plan->recv_bufloc[irecv]],
                   &out[plan->recv_offset[irecv]],&plan->unpackplan[irecv]);
    }

    // unpack all messages from scratch -> out

    for (i = 0; i < plan->nrecv; i++) {
      MPI_Waitany(plan->nrecv,plan->request,&irecv,&status);
      plan->unpack(&scratch[plan->recv_bufloc[irecv]],
                   &out[plan->recv_offset[irecv]],&plan->unpackplan[irecv]);
    }

  // use All2Allv collective for remap communication

  } else { 
    if (plan->commringlen > 0) {
      MPI_Status status;
      int i,isend,irecv;
      FFT_SCALAR *scratch;

      if (plan->memory == 0) scratch = buf;
      else scratch = plan->scratch;

      // create send and recv buffers for alltoallv collective

      int sendBufferSize = 0;
      int recvBufferSize = 0;
      for (int i=0;i<plan->nsend;i++)
        sendBufferSize += plan->send_size[i];
      for (int i=0;i<plan->nrecv;i++)
        recvBufferSize += plan->recv_size[i];

      FFT_SCALAR *packedSendBuffer 
        = (FFT_SCALAR *) malloc(sizeof(FFT_SCALAR) * sendBufferSize);
      FFT_SCALAR *packedRecvBuffer 
        = (FFT_SCALAR *) malloc(sizeof(FFT_SCALAR) * recvBufferSize);

      int *sendcnts = (int *) malloc(sizeof(int) * plan->commringlen);
      int *rcvcnts = (int *) malloc(sizeof(int) * plan->commringlen);
      int *sdispls = (int *) malloc(sizeof(int) * plan->commringlen);
      int *rdispls = (int *) malloc(sizeof(int) * plan->commringlen);
      int *nrecvmap = (int *) malloc(sizeof(int) * plan->commringlen);

      // create and populate send data, count and displacement buffers

      int currentSendBufferOffset = 0;
      for (isend = 0; isend < plan->commringlen; isend++) {
        sendcnts[isend] = 0;
        sdispls[isend] = 0;
        int foundentry = 0;
        for (int i=0;(i<plan->nsend && !foundentry); i++) {
          if (plan->send_proc[i] == plan->commringlist[isend]) {
            foundentry = 1;
            sendcnts[isend] = plan->send_size[i];
            sdispls[isend] = currentSendBufferOffset;
            plan->pack(&in[plan->send_offset[i]],
                       &packedSendBuffer[currentSendBufferOffset],
                       &plan->packplan[i]);
            currentSendBufferOffset += plan->send_size[i];
          }
        }
      }

      // create and populate recv count and displacement buffers

      int currentRecvBufferOffset = 0;
      for (irecv = 0; irecv < plan->commringlen; irecv++) {
        rcvcnts[irecv] = 0;
        rdispls[irecv] = 0;
        nrecvmap[irecv] = -1;
        int foundentry = 0;
        for (int i=0;(i<plan->nrecv && !foundentry); i++) {
          if (plan->recv_proc[i] == plan->commringlist[irecv]) {
            foundentry = 1;
            rcvcnts[irecv] = plan->recv_size[i];
            rdispls[irecv] = currentRecvBufferOffset;
            currentRecvBufferOffset += plan->recv_size[i];
            nrecvmap[irecv] = i;
          }
        }
      }

      int mpirc = MPI_Alltoallv(packedSendBuffer, sendcnts, sdispls,
                                MPI_FFT_SCALAR, packedRecvBuffer, rcvcnts,
                                rdispls, MPI_FFT_SCALAR, plan->comm);

      // unpack the data from the recv buffer into out

      currentRecvBufferOffset = 0;
      for (irecv = 0; irecv < plan->commringlen; irecv++) {
        if (nrecvmap[irecv] > -1) {
          plan->unpack(&packedRecvBuffer[currentRecvBufferOffset],
                       &out[plan->recv_offset[nrecvmap[irecv]]],
                       &plan->unpackplan[nrecvmap[irecv]]);
          currentRecvBufferOffset += plan->recv_size[nrecvmap[irecv]];
        }
      }

      // free temporary data structures

      free(sendcnts);
      free(rcvcnts);
      free(sdispls);
      free(rdispls);
      free(nrecvmap);
      free(packedSendBuffer);
      free(packedRecvBuffer);
    }
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

struct remap_plan_3d *remap_3d_create_plan(
  MPI_Comm comm,
  int in_ilo, int in_ihi, int in_jlo, int in_jhi,
  int in_klo, int in_khi,
  int out_ilo, int out_ihi, int out_jlo, int out_jhi,
  int out_klo, int out_khi,
  int nqty, int permute, int memory, int precision, int usecollective)

{

  struct remap_plan_3d *plan;
  struct extent_3d *inarray, *outarray;
  struct extent_3d in,out,overlap;
  int i,iproc,nsend,nrecv,ibuf,size,me,nprocs;

  // query MPI info

  MPI_Comm_rank(comm,&me);
  MPI_Comm_size(comm,&nprocs);

  // allocate memory for plan data struct

  plan = (struct remap_plan_3d *) malloc(sizeof(struct remap_plan_3d));
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
    plan->pack = pack_3d;

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
      plan->unpack = unpack_3d;
    else if (permute == 1) {
      if (nqty == 1)
        plan->unpack = unpack_3d_permute1_1;
      else if (nqty == 2)
        plan->unpack = unpack_3d_permute1_2;
      else
        plan->unpack = unpack_3d_permute1_n;
    }
    else if (permute == 2) {
      if (nqty == 1)
        plan->unpack = unpack_3d_permute2_1;
      else if (nqty == 2)
        plan->unpack = unpack_3d_permute2_2;
      else
        plan->unpack = unpack_3d_permute2_n;
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

  // create sub-comm rank list

  if (plan->usecollective) { 
    plan->commringlist = NULL;

    // merge recv and send rank lists
    // ask Steve Plimpton about method to more accurately determine 
    // maximum number of procs contributing to pencil

    int maxcommsize = nprocs;  
    int *commringlist = (int *) malloc(maxcommsize*sizeof(int));
    int commringlen = 0;

    for (int i = 0; i < nrecv; i++) {
      commringlist[i] = plan->recv_proc[i];
      commringlen++;
    }

    for (int i = 0; i < nsend; i++) {
      int foundentry = 0;
      for (int j=0;j<commringlen;j++)
        if (commringlist[j] == plan->send_proc[i]) foundentry = 1;
      if (!foundentry) {
        commringlist[commringlen] = plan->send_proc[i];
        commringlen++;
      }
    }

    // sort initial commringlist

    int swap = 0;
    for (int c = 0 ; c < (commringlen - 1); c++) {
      for (int d = 0 ; d < commringlen - c - 1; d++) {
        if (commringlist[d] > commringlist[d+1]) {
          swap = commringlist[d];
          commringlist[d]   = commringlist[d+1];
          commringlist[d+1] = swap;
        }
      }
    }

    // collide all inarray extents for the comm ring with all output
    // extents and all outarray extents for the comm ring with all input
    // extents - if there is a collison add the rank to the comm ring,
    // keep iterating until nothing is added to commring

    int commringappend = 1;
    while (commringappend) {
      int newcommringlen = commringlen;
      commringappend = 0;
      for (int i=0;i<commringlen;i++) {
        for (int j=0;j<nprocs;j++) {
          if (remap_3d_collide(&inarray[commringlist[i]],
                               &outarray[j],&overlap)) {
            int alreadyinlist = 0;
            for (int k=0;k<newcommringlen;k++) {
              if (commringlist[k] == j) {
                alreadyinlist = 1;
              }
            }
            if (!alreadyinlist) {
              commringlist[newcommringlen++] = j;
              commringappend = 1;
            }
          }
          if (remap_3d_collide(&outarray[commringlist[i]],
                               &inarray[j],&overlap)) {
            int alreadyinlist = 0;
            for (int k=0;k<newcommringlen;k++) {
              if (commringlist[k] == j) alreadyinlist = 1;
            }
            if (!alreadyinlist) {
              commringlist[newcommringlen++] = j;
              commringappend = 1;
            }
          }
        }
      }
      commringlen = newcommringlen;
    }

    // sort the final commringlist

    for (int c = 0 ; c < ( commringlen - 1 ); c++) {
      for (int d = 0 ; d < commringlen - c - 1; d++) {
        if (commringlist[d] > commringlist[d+1]) {
          swap = commringlist[d];
          commringlist[d]   = commringlist[d+1];
          commringlist[d+1] = swap;
        }
      }
    }

    // resize commringlist to final size

    commringlist = (int *) realloc(commringlist, commringlen*sizeof(int));

    // set the plan->commringlist

    plan->commringlen = commringlen;
    plan->commringlist = commringlist;
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

  plan->sendbuf = NULL;

  size = 0;
  for (nsend = 0; nsend < plan->nsend; nsend++)
    size = MAX(size,plan->send_size[nsend]);

  if (size) {
    plan->sendbuf = (FFT_SCALAR *) malloc(size*sizeof(FFT_SCALAR));
    if (plan->sendbuf == NULL) return NULL;
  }

  // if requested, allocate internal scratch space for recvs,
  // only need it if I will receive any data (including self)

  plan->scratch = NULL;

  if (memory == 1) {
    if (nrecv > 0) {
      plan->scratch =
        (FFT_SCALAR *) malloc(nqty*out.isize*out.jsize*out.ksize *
                              sizeof(FFT_SCALAR));
      if (plan->scratch == NULL) return NULL;
    }
  }

  // if using collective and the commringlist is NOT empty create a
  // communicator for the plan based off an MPI_Group created with
  // ranks from the commringlist

  if ((plan->usecollective && (plan->commringlen > 0))) {
    MPI_Group orig_group, new_group;
    MPI_Comm_group(comm, &orig_group);
    MPI_Group_incl(orig_group, plan->commringlen,
                   plan->commringlist, &new_group);
    MPI_Comm_create(comm, new_group, &plan->comm);
  }

  // if using collective and the comm ring list is empty create
  // a communicator for the plan with an empty group

  else if ((plan->usecollective) && (plan->commringlen == 0)) {
    MPI_Comm_create(comm, MPI_GROUP_EMPTY, &plan->comm);
  }

  // not using collective - dup comm

  else MPI_Comm_dup(comm,&plan->comm);

  // return pointer to plan

  return plan;
}

/* ----------------------------------------------------------------------
   Destroy a 3d remap plan
------------------------------------------------------------------------- */

void remap_3d_destroy_plan(struct remap_plan_3d *plan)
{
  // free MPI communicator

  if (!((plan->usecollective) && (plan->commringlen == 0)))
    MPI_Comm_free(&plan->comm);

  if (plan->usecollective) {
    if (plan->commringlist != NULL)
      free(plan->commringlist);
  }

  // free internal arrays

  if (plan->nsend || plan->self) {
    free(plan->send_offset);
    free(plan->send_size);
    free(plan->send_proc);
    free(plan->packplan);
    if (plan->sendbuf) free(plan->sendbuf);
  }

  if (plan->nrecv || plan->self) {
    free(plan->recv_offset);
    free(plan->recv_size);
    free(plan->recv_proc);
    free(plan->recv_bufloc);
    free(plan->request);
    free(plan->unpackplan);
    if (plan->scratch) free(plan->scratch);
  }

  // free plan itself

  free(plan);
}

/* ----------------------------------------------------------------------
   collide 2 sets of indices to determine overlap
   compare bounds of block1 with block2 to see if they overlap
   return 1 if they do and put bounds of overlapping section in overlap
   return 0 if they do not overlap
------------------------------------------------------------------------- */

int remap_3d_collide(struct extent_3d *block1, struct extent_3d *block2,
                     struct extent_3d *overlap)

{
  overlap->ilo = MAX(block1->ilo,block2->ilo);
  overlap->ihi = MIN(block1->ihi,block2->ihi);
  overlap->jlo = MAX(block1->jlo,block2->jlo);
  overlap->jhi = MIN(block1->jhi,block2->jhi);
  overlap->klo = MAX(block1->klo,block2->klo);
  overlap->khi = MIN(block1->khi,block2->khi);

  if (overlap->ilo > overlap->ihi ||
      overlap->jlo > overlap->jhi ||
      overlap->klo > overlap->khi) return 0;

  overlap->isize = overlap->ihi - overlap->ilo + 1;
  overlap->jsize = overlap->jhi - overlap->jlo + 1;
  overlap->ksize = overlap->khi - overlap->klo + 1;

  return 1;
}
