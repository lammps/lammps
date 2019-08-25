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
#include <cstdlib>
#include <cstring>
#include "irregular.h"
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "comm.h"
#include "memory.h"

using namespace LAMMPS_NS;

#if defined(LMP_QSORT)
// allocate space for static class variable
// prototype for non-class function
int *Irregular::proc_recv_copy;
static int compare_standalone(const void *, const void *);
#else
#include "mergesort.h"
// prototype for non-class function
static int compare_standalone(const int, const int, void *);
#endif

#define BUFFACTOR 1.5
#define BUFMIN 1000
#define BUFEXTRA 1000

/* ---------------------------------------------------------------------- */

Irregular::Irregular(LAMMPS *lmp) : Pointers(lmp)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  triclinic = domain->triclinic;
  map_style = atom->map_style;

  // migrate work vectors

  maxlocal = 0;
  mproclist = NULL;
  msizes = NULL;

  // send buffers

  maxdbuf = 0;
  dbuf = NULL;
  maxbuf = 0;
  buf = NULL;

  // universal work vectors

  memory->create(work1,nprocs,"irregular:work1");
  memory->create(work2,nprocs,"irregular:work2");

  // initialize buffers for migrate atoms, not used for datum comm
  // these can persist for multiple irregular operations

  maxsend = BUFMIN;
  memory->create(buf_send,maxsend+BUFEXTRA,"comm:buf_send");
  maxrecv = BUFMIN;
  memory->create(buf_recv,maxrecv,"comm:buf_recv");
}

/* ---------------------------------------------------------------------- */

Irregular::~Irregular()
{
  memory->destroy(mproclist);
  memory->destroy(msizes);
  memory->destroy(dbuf);
  memory->destroy(buf);
  memory->destroy(work1);
  memory->destroy(work2);
  memory->destroy(buf_send);
  memory->destroy(buf_recv);
}

/* ----------------------------------------------------------------------
   communicate atoms to new owning procs via irregular communication
   can be used in place of comm->exchange()
   unlike exchange(), allows atoms to have moved arbitrarily long distances
   sets up irregular plan, invokes it, destroys it
   sortflag = flag for sorting order of received messages by proc ID
   preassign = 1 if already know procs that atoms are assigned to via RCB
   procassign = list of proc assignments for each owned atom
   atoms MUST be remapped to be inside simulation box before this is called
   for triclinic: atoms must be in lamda coords (0-1) before this is called
------------------------------------------------------------------------- */

void Irregular::migrate_atoms(int sortflag, int preassign, int *procassign)
{
  // clear global->local map since atoms move to new procs
  // clear old ghosts so map_set() at end will operate only on local atoms
  // exchange() doesn't need to clear ghosts b/c borders()
  //   is called right after and it clears ghosts and calls map_set()

  if (map_style) atom->map_clear();
  atom->nghost = 0;
  atom->avec->clear_bonus();

  // subbox bounds for orthogonal or triclinic box

  double *sublo,*subhi;
  if (triclinic == 0) {
    sublo = domain->sublo;
    subhi = domain->subhi;
  } else {
    sublo = domain->sublo_lamda;
    subhi = domain->subhi_lamda;
  }

  // if Comm will be called to assign new atom coords to procs,
  // may need to setup RCB info

  if (!preassign) comm->coord2proc_setup();

  // loop over atoms, flag any that are not in my sub-box
  // fill buffer with atoms leaving my box, using < and >=
  // assign which proc it belongs to via Comm::coord2proc()
  // if coord2proc() returns me, due to round-off
  //   in triclinic x2lamda(), then keep atom and don't send
  // when atom is deleted, fill it in with last atom

  AtomVec *avec = atom->avec;
  double **x = atom->x;
  int nlocal = atom->nlocal;

  if (nlocal > maxlocal) {
    maxlocal = nlocal;
    memory->destroy(mproclist);
    memory->destroy(msizes);
    memory->create(mproclist,maxlocal,"irregular:mproclist");
    memory->create(msizes,maxlocal,"irregular:msizes");
  }

  int igx,igy,igz;
  int nsend = 0;
  int nsendatom = 0;
  int i = 0;

  if (preassign) {
    while (i < nlocal) {
      if (procassign[i] == me) i++;
      else {
        mproclist[nsendatom] = procassign[i];
        if (nsend > maxsend) grow_send(nsend,1);
        msizes[nsendatom] = avec->pack_exchange(i,&buf_send[nsend]);
        nsend += msizes[nsendatom];
        nsendatom++;
        avec->copy(nlocal-1,i,1);
        procassign[i] = procassign[nlocal-1];
        nlocal--;
      }
    }

  } else {
    while (i < nlocal) {
      if (x[i][0] < sublo[0] || x[i][0] >= subhi[0] ||
          x[i][1] < sublo[1] || x[i][1] >= subhi[1] ||
          x[i][2] < sublo[2] || x[i][2] >= subhi[2]) {
        mproclist[nsendatom] = comm->coord2proc(x[i],igx,igy,igz);
        if (mproclist[nsendatom] == me) i++;
        else {
          if (nsend > maxsend) grow_send(nsend,1);
          msizes[nsendatom] = avec->pack_exchange(i,&buf_send[nsend]);
          nsend += msizes[nsendatom];
          nsendatom++;
          avec->copy(nlocal-1,i,1);
          nlocal--;
        }
      } else i++;
    }
  }

  atom->nlocal = nlocal;

  // create irregular communication plan, perform comm, destroy plan
  // returned nrecv = size of buffer needed for incoming atoms

  int nrecv = create_atom(nsendatom,msizes,mproclist,sortflag);
  if (nrecv > maxrecv) grow_recv(nrecv);
  exchange_atom(buf_send,msizes,buf_recv);
  destroy_atom();

  // add received atoms to my list

  int m = 0;
  while (m < nrecv) m += avec->unpack_exchange(&buf_recv[m]);

  // reset global->local map

  if (map_style) atom->map_set();
}

/* ----------------------------------------------------------------------
   check if any atoms need to migrate further than one proc away in any dim
   if not, caller can decide to use comm->exchange() instead
   should not be called for layout = TILED
   atoms must be remapped to be inside simulation box before this is called
   for triclinic: atoms must be in lamda coords (0-1) before this is called
   return 1 if migrate required, 0 if not
------------------------------------------------------------------------- */

int Irregular::migrate_check()
{
  // migrate required if comm layout is tiled
  // cannot use myloc[] logic below

  if (comm->layout == Comm::LAYOUT_TILED) return 1;

  // subbox bounds for orthogonal or triclinic box

  double *sublo,*subhi;
  if (triclinic == 0) {
    sublo = domain->sublo;
    subhi = domain->subhi;
  } else {
    sublo = domain->sublo_lamda;
    subhi = domain->subhi_lamda;
  }

  // loop over atoms, check for any that are not in my sub-box
  // assign which proc it belongs to via Comm::coord2proc()
  // if logical igx,igy,igz of newproc > one away from myloc, set flag = 1
  // this check needs to observe PBC
  // cannot check via comm->procneigh since it ignores PBC

  double **x = atom->x;
  int nlocal = atom->nlocal;

  int *periodicity = domain->periodicity;
  int *myloc = comm->myloc;
  int *procgrid = comm->procgrid;
  int igx,igy,igz,glo,ghi;

  int flag = 0;
  for (int i = 0; i < nlocal; i++) {
    if (x[i][0] < sublo[0] || x[i][0] >= subhi[0] ||
        x[i][1] < sublo[1] || x[i][1] >= subhi[1] ||
        x[i][2] < sublo[2] || x[i][2] >= subhi[2]) {
      comm->coord2proc(x[i],igx,igy,igz);

      glo = myloc[0] - 1;
      ghi = myloc[0] + 1;
      if (periodicity[0]) {
        if (glo < 0) glo = procgrid[0] - 1;
        if (ghi >= procgrid[0]) ghi = 0;
      }
      if (igx != myloc[0] && igx != glo && igx != ghi) flag = 1;

      glo = myloc[1] - 1;
      ghi = myloc[1] + 1;
      if (periodicity[1]) {
        if (glo < 0) glo = procgrid[1] - 1;
        if (ghi >= procgrid[1]) ghi = 0;
      }
      if (igy != myloc[1] && igy != glo && igy != ghi) flag = 1;

      glo = myloc[2] - 1;
      ghi = myloc[2] + 1;
      if (periodicity[2]) {
        if (glo < 0) glo = procgrid[2] - 1;
        if (ghi >= procgrid[2]) ghi = 0;
      }
      if (igz != myloc[2] && igz != glo && igz != ghi) flag = 1;
    }
  }

  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_MAX,world);
  return flagall;
}

/* ----------------------------------------------------------------------
   create a communication plan for atoms
   n = # of atoms to send
   sizes = # of doubles for each atom
   proclist = proc to send each atom to (not including self)
   sortflag = flag for sorting order of received messages by proc ID
   return total # of doubles I will recv (not including self)
------------------------------------------------------------------------- */

int Irregular::create_atom(int n, int *sizes, int *proclist, int sortflag)
{
  int i;

  // setup for collective comm
  // work1 = 1 for procs I send a message to, not including self
  // work2 = 1 for all procs, used for ReduceScatter

  for (i = 0; i < nprocs; i++) {
    work1[i] = 0;
    work2[i] = 1;
  }
  for (i = 0; i < n; i++) work1[proclist[i]] = 1;
  work1[me] = 0;

  // nrecv_proc = # of procs I receive messages from, not including self
  // options for performing ReduceScatter operation
  // some are more efficient on some machines at big sizes

#ifdef LAMMPS_RS_ALLREDUCE_INPLACE
  MPI_Allreduce(MPI_IN_PLACE,work1,nprocs,MPI_INT,MPI_SUM,world);
  nrecv_proc = work1[me];
#else
#ifdef LAMMPS_RS_ALLREDUCE
  MPI_Allreduce(work1,work2,nprocs,MPI_INT,MPI_SUM,world);
  nrecv_proc = work2[me];
#else
  MPI_Reduce_scatter(work1,&nrecv_proc,work2,MPI_INT,MPI_SUM,world);
#endif
#endif

  // allocate receive arrays

  proc_recv = new int[nrecv_proc];
  length_recv = new int[nrecv_proc];
  request = new MPI_Request[nrecv_proc];
  status = new MPI_Status[nrecv_proc];

  // nsend_proc = # of messages I send

  for (i = 0; i < nprocs; i++) work1[i] = 0;
  for (i = 0; i < n; i++) work1[proclist[i]] += sizes[i];

  nsend_proc = 0;
  for (i = 0; i < nprocs; i++)
    if (work1[i]) nsend_proc++;

  // allocate send arrays

  proc_send = new int[nsend_proc];
  length_send = new int[nsend_proc];
  num_send = new int[nsend_proc];
  index_send = new int[n];
  offset_send = new int[n];

  // list still stores size of message for procs I send to
  // proc_send = procs I send to
  // length_send = # of doubles I send to each proc
  // to balance pattern of send messages:
  //   each proc begins with iproc > me, continues until iproc = me
  // reset list to store which send message each proc corresponds to

  int iproc = me;
  int isend = 0;
  for (i = 0; i < nprocs; i++) {
    iproc++;
    if (iproc == nprocs) iproc = 0;
    if (work1[iproc] > 0) {
      proc_send[isend] = iproc;
      length_send[isend] = work1[iproc];
      work1[iproc] = isend;
      isend++;
    }
  }

  // num_send = # of atoms I send to each proc

  for (i = 0; i < nsend_proc; i++) num_send[i] = 0;
  for (i = 0; i < n; i++) {
    isend = work1[proclist[i]];
    num_send[isend]++;
  }

  // work2 = offsets into index_send for each proc I send to
  // index_send = list of which atoms to send to each proc
  //   1st N1 values are atom indices for 1st proc,
  //   next N2 values are atom indices for 2nd proc, etc
  // offset_send = where each atom starts in send buffer

  work2[0] = 0;
  for (i = 1; i < nsend_proc; i++) work2[i] = work2[i-1] + num_send[i-1];

  for (i = 0; i < n; i++) {
    isend = work1[proclist[i]];
    index_send[work2[isend]++] = i;
    if (i) offset_send[i] = offset_send[i-1] + sizes[i-1];
    else offset_send[i] = 0;
  }

  // tell receivers how much data I send
  // sendmax_proc = # of doubles I send in largest single message

  sendmax_proc = 0;
  for (i = 0; i < nsend_proc; i++) {
    MPI_Request tmpReq; // Use non-blocking send to avoid possible deadlock
    MPI_Isend(&length_send[i],1,MPI_INT,proc_send[i],0,world,&tmpReq);
    MPI_Request_free(&tmpReq); // the MPI_Barrier below marks completion
    sendmax_proc = MAX(sendmax_proc,length_send[i]);
  }

  // receive incoming messages
  // proc_recv = procs I recv from
  // length_recv = # of doubles each proc sends me
  // nrecvsize = total size of atom data I recv

  int nrecvsize = 0;
  for (i = 0; i < nrecv_proc; i++) {
    MPI_Recv(&length_recv[i],1,MPI_INT,MPI_ANY_SOURCE,0,world,status);
    proc_recv[i] = status->MPI_SOURCE;
    nrecvsize += length_recv[i];
  }

  // sort proc_recv and length_recv by proc ID if requested
  // useful for debugging to insure reproducible ordering of received atoms
  // invoke by adding final arg = 1 to create_atom() call in migrate_atoms()

  if (sortflag) {
    int *order = new int[nrecv_proc];
    int *proc_recv_ordered = new int[nrecv_proc];
    int *length_recv_ordered = new int[nrecv_proc];

    for (i = 0; i < nrecv_proc; i++) order[i] = i;

#if defined(LMP_QSORT)
    proc_recv_copy = proc_recv;
    qsort(order,nrecv_proc,sizeof(int),compare_standalone);
#else
    merge_sort(order,nrecv_proc,(void *)proc_recv,compare_standalone);
#endif

    int j;
    for (i = 0; i < nrecv_proc; i++) {
      j = order[i];
      proc_recv_ordered[i] = proc_recv[j];
      length_recv_ordered[i] = length_recv[j];
    }

    memcpy(proc_recv,proc_recv_ordered,nrecv_proc*sizeof(int));
    memcpy(length_recv,length_recv_ordered,nrecv_proc*sizeof(int));
    delete [] order;
    delete [] proc_recv_ordered;
    delete [] length_recv_ordered;
  }

  // barrier to insure all MPI_ANY_SOURCE messages are received
  // else another proc could proceed to exchange_atom() and send to me

  MPI_Barrier(world);

  // return size of atom data I will receive

  return nrecvsize;
}

#if defined(LMP_QSORT)

/* ----------------------------------------------------------------------
   comparison function invoked by qsort()
   accesses static class member proc_recv_copy, set before call to qsort()
------------------------------------------------------------------------- */

int compare_standalone(const void *iptr, const void *jptr)
{
  int i = *((int *) iptr);
  int j = *((int *) jptr);
  int *proc_recv = Irregular::proc_recv_copy;
  if (proc_recv[i] < proc_recv[j]) return -1;
  if (proc_recv[i] > proc_recv[j]) return 1;
  return 0;
}

#else

/* ----------------------------------------------------------------------
   comparison function invoked by merge_sort()
   void pointer contains proc_recv list;
------------------------------------------------------------------------- */

int compare_standalone(const int i, const int j, void *ptr)
{
  int *proc_recv = (int *) ptr;
  if (proc_recv[i] < proc_recv[j]) return -1;
  if (proc_recv[i] > proc_recv[j]) return 1;
  return 0;
}

#endif

/* ----------------------------------------------------------------------
   communicate atoms via PlanAtom
   sendbuf = list of atoms to send
   sizes = # of doubles for each atom
   recvbuf = received atoms
------------------------------------------------------------------------- */

void Irregular::exchange_atom(double *sendbuf, int *sizes, double *recvbuf)
{
  int i,m,n,count;
  bigint offset;

  // post all receives

  offset = 0;
  for (int irecv = 0; irecv < nrecv_proc; irecv++) {
    MPI_Irecv(&recvbuf[offset],length_recv[irecv],MPI_DOUBLE,
              proc_recv[irecv],0,world,&request[irecv]);
    offset += length_recv[irecv];
  }

  // reallocate buf for largest send if necessary

  if (sendmax_proc > maxdbuf) {
    memory->destroy(dbuf);
    maxdbuf = sendmax_proc;
    memory->create(dbuf,maxdbuf,"irregular:dbuf");
  }

  // send each message
  // pack buf with list of atoms
  // m = index of atom in sendbuf

  n = 0;
  for (int isend = 0; isend < nsend_proc; isend++) {
    offset = 0;
    count = num_send[isend];
    for (i = 0; i < count; i++) {
      m = index_send[n++];
      memcpy(&dbuf[offset],&sendbuf[offset_send[m]],sizes[m]*sizeof(double));
      offset += sizes[m];
    }
    MPI_Send(dbuf,length_send[isend],MPI_DOUBLE,proc_send[isend],0,world);
  }

  // wait on all incoming messages

  if (nrecv_proc) MPI_Waitall(nrecv_proc,request,status);
}

/* ----------------------------------------------------------------------
   destroy vectors in communication plan for atoms
------------------------------------------------------------------------- */

void Irregular::destroy_atom()
{
  delete [] proc_send;
  delete [] length_send;
  delete [] num_send;
  delete [] index_send;
  delete [] offset_send;
  delete [] proc_recv;
  delete [] length_recv;
  delete [] request;
  delete [] status;
}

/* ----------------------------------------------------------------------
   create communication plan based on list of datums of uniform size
   n = # of datums to send
   proclist = proc to send each datum to, can include self
   sortflag = flag for sorting order of received messages by proc ID
   return total # of datums I will recv, including any to self
------------------------------------------------------------------------- */

int Irregular::create_data(int n, int *proclist, int sortflag)
{
  int i,m;

  // setup for collective comm
  // work1 = 1 for procs I send a message to, not including self
  // work2 = 1 for all procs, used for ReduceScatter

  for (i = 0; i < nprocs; i++) {
    work1[i] = 0;
    work2[i] = 1;
  }
  for (i = 0; i < n; i++) work1[proclist[i]] = 1;
  work1[me] = 0;

  // nrecv_proc = # of procs I receive messages from, not including self
  // options for performing ReduceScatter operation
  // some are more efficient on some machines at big sizes

#ifdef LAMMPS_RS_ALLREDUCE_INPLACE
  MPI_Allreduce(MPI_IN_PLACE,work1,nprocs,MPI_INT,MPI_SUM,world);
  nrecv_proc = work1[me];
#else
#ifdef LAMMPS_RS_ALLREDUCE
  MPI_Allreduce(work1,work2,nprocs,MPI_INT,MPI_SUM,world);
  nrecv_proc = work2[me];
#else
  MPI_Reduce_scatter(work1,&nrecv_proc,work2,MPI_INT,MPI_SUM,world);
#endif
#endif

  // allocate receive arrays

  proc_recv = new int[nrecv_proc];
  num_recv = new int[nrecv_proc];
  request = new MPI_Request[nrecv_proc];
  status = new MPI_Status[nrecv_proc];

  // work1 = # of datums I send to each proc, including self
  // nsend_proc = # of procs I send messages to, not including self

  for (i = 0; i < nprocs; i++) work1[i] = 0;
  for (i = 0; i < n; i++) work1[proclist[i]]++;

  nsend_proc = 0;
  for (i = 0; i < nprocs; i++)
    if (work1[i]) nsend_proc++;
  if (work1[me]) nsend_proc--;

  // allocate send and self arrays

  proc_send = new int[nsend_proc];
  num_send = new int[nsend_proc];
  index_send = new int[n-work1[me]];
  index_self = new int[work1[me]];
  maxindex = n;

  // proc_send = procs I send to
  // num_send = # of datums I send to each proc
  // num_self = # of datums I copy to self
  // to balance pattern of send messages:
  //   each proc begins with iproc > me, continues until iproc = me
  // reset work1 to store which send message each proc corresponds to

  int iproc = me;
  int isend = 0;
  for (i = 0; i < nprocs; i++) {
    iproc++;
    if (iproc == nprocs) iproc = 0;
    if (iproc == me) {
      num_self = work1[iproc];
      work1[iproc] = 0;
    } else if (work1[iproc] > 0) {
      proc_send[isend] = iproc;
      num_send[isend] = work1[iproc];
      work1[iproc] = isend;
      isend++;
    }
  }

  // work2 = offsets into index_send for each proc I send to
  // m = ptr into index_self
  // index_send = list of which datums to send to each proc
  //   1st N1 values are datum indices for 1st proc,
  //   next N2 values are datum indices for 2nd proc, etc
  // index_self = list of which datums to copy to self

  work2[0] = 0;
  for (i = 1; i < nsend_proc; i++) work2[i] = work2[i-1] + num_send[i-1];

  m = 0;
  for (i = 0; i < n; i++) {
    iproc = proclist[i];
    if (iproc == me) index_self[m++] = i;
    else {
      isend = work1[iproc];
      index_send[work2[isend]++] = i;
    }
  }

  // tell receivers how much data I send
  // sendmax_proc = largest # of datums I send in a single message

  sendmax_proc = 0;
  for (i = 0; i < nsend_proc; i++) {
    MPI_Request tmpReq; // Use non-blocking send to avoid possible deadlock
    MPI_Isend(&num_send[i],1,MPI_INT,proc_send[i],0,world,&tmpReq);
    MPI_Request_free(&tmpReq); // the MPI_Barrier below marks completion
    sendmax_proc = MAX(sendmax_proc,num_send[i]);
  }

  // receive incoming messages
  // proc_recv = procs I recv from
  // num_recv = # of datums each proc sends me
  // nrecvdatum = total # of datums I recv

  int nrecvdatum = 0;
  for (i = 0; i < nrecv_proc; i++) {
    MPI_Recv(&num_recv[i],1,MPI_INT,MPI_ANY_SOURCE,0,world,status);
    proc_recv[i] = status->MPI_SOURCE;
    nrecvdatum += num_recv[i];
  }
  nrecvdatum += num_self;

  // sort proc_recv and num_recv by proc ID if requested
  // useful for debugging to insure reproducible ordering of received datums

  if (sortflag) {
    int *order = new int[nrecv_proc];
    int *proc_recv_ordered = new int[nrecv_proc];
    int *num_recv_ordered = new int[nrecv_proc];

    for (i = 0; i < nrecv_proc; i++) order[i] = i;

#if defined(LMP_QSORT)
    proc_recv_copy = proc_recv;
    qsort(order,nrecv_proc,sizeof(int),compare_standalone);
#else
    merge_sort(order,nrecv_proc,(void *)proc_recv,compare_standalone);
#endif

    int j;
    for (i = 0; i < nrecv_proc; i++) {
      j = order[i];
      proc_recv_ordered[i] = proc_recv[j];
      num_recv_ordered[i] = num_recv[j];
    }

    memcpy(proc_recv,proc_recv_ordered,nrecv_proc*sizeof(int));
    memcpy(num_recv,num_recv_ordered,nrecv_proc*sizeof(int));
    delete [] order;
    delete [] proc_recv_ordered;
    delete [] num_recv_ordered;
  }

  // barrier to insure all MPI_ANY_SOURCE messages are received
  // else another proc could proceed to exchange_data() and send to me

  MPI_Barrier(world);

  // return # of datums I will receive

  return nrecvdatum;
}

/* ----------------------------------------------------------------------
   create communication plan based on list of datums of uniform size
   n = # of datums to send
   procs = how many datums to send to each proc, must include self
   sort = flag for sorting order of received messages by proc ID
   return total # of datums I will recv, including any to self
------------------------------------------------------------------------- */

int Irregular::create_data_grouped(int n, int *procs, int sortflag)
{
  int i,j,k,m;

  // setup for collective comm
  // work1 = # of datums I send to each proc, set self to 0
  // work2 = 1 for all procs, used for ReduceScatter

  for (i = 0; i < nprocs; i++) {
    work1[i] = procs[i];
    work2[i] = 1;
  }
  work1[me] = 0;

  // nrecv_proc = # of procs I receive messages from, not including self
  // options for performing ReduceScatter operation
  // some are more efficient on some machines at big sizes

#ifdef LAMMPS_RS_ALLREDUCE_INPLACE
  MPI_Allreduce(MPI_IN_PLACE,work1,nprocs,MPI_INT,MPI_SUM,world);
  nrecv_proc = work1[me];
#else
#ifdef LAMMPS_RS_ALLREDUCE
  MPI_Allreduce(work1,work2,nprocs,MPI_INT,MPI_SUM,world);
  nrecv_proc = work2[me];
#else
  MPI_Reduce_scatter(work1,&nrecv_proc,work2,MPI_INT,MPI_SUM,world);
#endif
#endif

  // allocate receive arrays

  proc_recv = new int[nrecv_proc];
  num_recv = new int[nrecv_proc];
  request = new MPI_Request[nrecv_proc];
  status = new MPI_Status[nrecv_proc];

  // work1 = # of datums I send to each proc, including self
  // nsend_proc = # of procs I send messages to, not including self

  for (i = 0; i < nprocs; i++) work1[i] = procs[i];

  nsend_proc = 0;
  for (i = 0; i < nprocs; i++)
    if (work1[i]) nsend_proc++;
  if (work1[me]) nsend_proc--;

  // allocate send and self arrays

  proc_send = new int[nsend_proc];
  num_send = new int[nsend_proc];
  index_send = new int[n-work1[me]];
  index_self = new int[work1[me]];
  maxindex = n;

  // proc_send = procs I send to
  // num_send = # of datums I send to each proc
  // num_self = # of datums I copy to self
  // to balance pattern of send messages:
  //   each proc begins with iproc > me, continues until iproc = me
  // reset work1 to store which send message each proc corresponds to

  int iproc = me;
  int isend = 0;
  for (i = 0; i < nprocs; i++) {
    iproc++;
    if (iproc == nprocs) iproc = 0;
    if (iproc == me) {
      num_self = work1[iproc];
      work1[iproc] = 0;
    } else if (work1[iproc] > 0) {
      proc_send[isend] = iproc;
      num_send[isend] = work1[iproc];
      work1[iproc] = isend;
      isend++;
    }
  }

  // work2 = offsets into index_send for each proc I send to
  // m = ptr into index_self
  // index_send = list of which datums to send to each proc
  //   1st N1 values are datum indices for 1st proc,
  //   next N2 values are datum indices for 2nd proc, etc
  // index_self = list of which datums to copy to self

  work2[0] = 0;
  for (i = 1; i < nsend_proc; i++) work2[i] = work2[i-1] + num_send[i-1];

  m = 0;
  i = 0;
  for (iproc = 0; iproc < nprocs; iproc++) {
    k = procs[iproc];
    for (j = 0; j < k; j++) {
      if (iproc == me) index_self[m++] = i++;
      else {
        isend = work1[iproc];
        index_send[work2[isend]++] = i++;
      }
    }
  }

  // tell receivers how much data I send
  // sendmax_proc = largest # of datums I send in a single message

  sendmax_proc = 0;
  for (i = 0; i < nsend_proc; i++) {
    MPI_Request tmpReq; // Use non-blocking send to avoid possible deadlock
    MPI_Isend(&num_send[i],1,MPI_INT,proc_send[i],0,world,&tmpReq);
    MPI_Request_free(&tmpReq); // the MPI_Barrier below marks completion
    sendmax_proc = MAX(sendmax_proc,num_send[i]);
  }

  // receive incoming messages
  // proc_recv = procs I recv from
  // num_recv = # of datums each proc sends me
  // nrecvdatum = total # of datums I recv

  int nrecvdatum = 0;
  for (i = 0; i < nrecv_proc; i++) {
    MPI_Recv(&num_recv[i],1,MPI_INT,MPI_ANY_SOURCE,0,world,status);
    proc_recv[i] = status->MPI_SOURCE;
    nrecvdatum += num_recv[i];
  }
  nrecvdatum += num_self;

  // sort proc_recv and num_recv by proc ID if requested
  // useful for debugging to insure reproducible ordering of received datums

  if (sortflag) {
    int *order = new int[nrecv_proc];
    int *proc_recv_ordered = new int[nrecv_proc];
    int *num_recv_ordered = new int[nrecv_proc];

    for (i = 0; i < nrecv_proc; i++) order[i] = i;

#if defined(LMP_QSORT)
    proc_recv_copy = proc_recv;
    qsort(order,nrecv_proc,sizeof(int),compare_standalone);
#else
    merge_sort(order,nrecv_proc,(void *)proc_recv,compare_standalone);
#endif

    int j;
    for (i = 0; i < nrecv_proc; i++) {
      j = order[i];
      proc_recv_ordered[i] = proc_recv[j];
      num_recv_ordered[i] = num_recv[j];
    }

    memcpy(proc_recv,proc_recv_ordered,nrecv_proc*sizeof(int));
    memcpy(num_recv,num_recv_ordered,nrecv_proc*sizeof(int));
    delete [] order;
    delete [] proc_recv_ordered;
    delete [] num_recv_ordered;
  }

  // barrier to insure all MPI_ANY_SOURCE messages are received
  // else another proc could proceed to exchange_data() and send to me

  MPI_Barrier(world);

  // return # of datums I will receive

  return nrecvdatum;
}

/* ----------------------------------------------------------------------
   communicate datums via PlanData
   sendbuf = list of datums to send
   nbytes = size of each datum
   recvbuf = received datums (including copied from me)
------------------------------------------------------------------------- */

void Irregular::exchange_data(char *sendbuf, int nbytes, char *recvbuf)
{
  int i,n,count;
  bigint m;       // these 2 lines enable send/recv buf to be larger than 2 GB
  char *dest;

  // post all receives, starting after self copies

  bigint offset = num_self*(bigint)nbytes;
  for (int irecv = 0; irecv < nrecv_proc; irecv++) {
    MPI_Irecv(&recvbuf[offset],num_recv[irecv]*nbytes,MPI_CHAR,
              proc_recv[irecv],0,world,&request[irecv]);
    offset += num_recv[irecv]*nbytes;
  }

  // reallocate buf for largest send if necessary

  if (sendmax_proc*nbytes > maxbuf) {
    memory->destroy(buf);
    maxbuf = sendmax_proc*nbytes;
    memory->create(buf,maxbuf,"irregular:buf");
  }

  // send each message
  // pack buf with list of datums
  // m = index of datum in sendbuf

  n = 0;
  for (int isend = 0; isend < nsend_proc; isend++) {
    count = num_send[isend];
    dest = buf;
    for (i = 0; i < count; i++) {
      m = index_send[n++];
      memcpy(dest,&sendbuf[m*nbytes],nbytes);
      dest += nbytes;
    }
    MPI_Send(buf,count*nbytes,MPI_CHAR,proc_send[isend],0,world);
  }

  // copy datums to self, put at beginning of recvbuf

  dest = recvbuf;
  for (i = 0; i < num_self; i++) {
    m = index_self[i];
    memcpy(dest,&sendbuf[m*nbytes],nbytes);
    dest += nbytes;
  }

  // wait on all incoming messages

  if (nrecv_proc) MPI_Waitall(nrecv_proc,request,status);
}

/* ----------------------------------------------------------------------
   destroy vectors in communication plan for datums
------------------------------------------------------------------------- */

void Irregular::destroy_data()
{
  delete [] proc_send;
  delete [] num_send;
  delete [] index_send;
  delete [] proc_recv;
  delete [] num_recv;
  delete [] index_self;
  delete [] request;
  delete [] status;
}

/* ----------------------------------------------------------------------
   realloc the size of the send buffer as needed with BUFFACTOR & BUFEXTRA
   if flag = 1, realloc
   if flag = 0, don't need to realloc with copy, just free/malloc
------------------------------------------------------------------------- */

void Irregular::grow_send(int n, int flag)
{
  maxsend = static_cast<int> (BUFFACTOR * n);
  if (flag)
    memory->grow(buf_send,maxsend+BUFEXTRA,"comm:buf_send");
  else {
    memory->destroy(buf_send);
    memory->create(buf_send,maxsend+BUFEXTRA,"comm:buf_send");
  }
}

/* ----------------------------------------------------------------------
   free/malloc the size of the recv buffer as needed with BUFFACTOR
------------------------------------------------------------------------- */

void Irregular::grow_recv(int n)
{
  maxrecv = static_cast<int> (BUFFACTOR * n);
  memory->destroy(buf_recv);
  memory->create(buf_recv,maxrecv,"comm:buf_recv");
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

bigint Irregular::memory_usage()
{
  bigint bytes = 0;
  bytes += maxsend*sizeof(double);   // buf_send
  bytes += maxrecv*sizeof(double);   // buf_recv
  bytes += maxdbuf*sizeof(double);   // dbuf
  bytes += maxbuf;                   // buf
  bytes += 2*maxlocal*sizeof(int);   // mproclist,msizes
  bytes += 2*nprocs*sizeof(int);     // work1,work2
  return bytes;
}
