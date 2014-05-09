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

#include "lmptype.h"
#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "irregular.h"
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "comm.h"
#include "memory.h"

using namespace LAMMPS_NS;

// allocate space for static class variable
// prototype for non-class function

int *Irregular::proc_recv_copy;
int compare_standalone(const void *, const void *);

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
  procgrid = comm->procgrid;
  grid2proc = comm->grid2proc;

  aplan = NULL;
  dplan = NULL;

  // initialize buffers for atom comm, not used for datum comm
  // these can persist for multiple irregular operations

  maxsend = BUFMIN;
  memory->create(buf_send,maxsend+BUFEXTRA,"comm:buf_send");
  maxrecv = BUFMIN;
  memory->create(buf_recv,maxrecv,"comm:buf_recv");
}

/* ---------------------------------------------------------------------- */

Irregular::~Irregular()
{
  if (aplan) destroy_atom();
  if (dplan) destroy_data();

  memory->destroy(buf_send);
  memory->destroy(buf_recv);
}

/* ----------------------------------------------------------------------
   communicate atoms to new owning procs via irregular communication
   can be used in place of comm->exchange()
   unlike exchange(), allows atoms to have moved arbitrarily long distances
   sets up irregular plan, invokes it, destroys it
   atoms must be remapped to be inside simulation box before this is called
   for triclinic: atoms must be in lamda coords (0-1) before this is called
------------------------------------------------------------------------- */

void Irregular::migrate_atoms(int sortflag)
{
  // clear global->local map since atoms move to new procs
  // clear old ghosts so map_set() at end will operate only on local atoms
  // exchange() doesn't need to clear ghosts b/c borders()
  //   is called right after and it clears ghosts and calls map_set()

  if (map_style) atom->map_clear();
  atom->nghost = 0;
  atom->avec->clear_bonus();

  // subbox bounds for orthogonal or triclinic box
  // other comm/domain data used by coord2proc()

  double *sublo,*subhi;
  if (triclinic == 0) {
    sublo = domain->sublo;
    subhi = domain->subhi;
  } else {
    sublo = domain->sublo_lamda;
    subhi = domain->subhi_lamda;
  }

  uniform = comm->uniform;
  xsplit = comm->xsplit;
  ysplit = comm->ysplit;
  zsplit = comm->zsplit;
  boxlo = domain->boxlo;
  prd = domain->prd;

  // loop over atoms, flag any that are not in my sub-box
  // fill buffer with atoms leaving my box, using < and >=
  // assign which proc it belongs to via coord2proc()
  // if coord2proc() returns me, due to round-off
  //   in triclinic x2lamda(), then keep atom and don't send
  // when atom is deleted, fill it in with last atom

  AtomVec *avec = atom->avec;
  double **x = atom->x;
  int nlocal = atom->nlocal;

  int nsend = 0;
  int nsendatom = 0;
  int *sizes = new int[nlocal];
  int *proclist = new int[nlocal];
  int igx,igy,igz;

  int i = 0;
  while (i < nlocal) {
    if (x[i][0] < sublo[0] || x[i][0] >= subhi[0] ||
        x[i][1] < sublo[1] || x[i][1] >= subhi[1] ||
        x[i][2] < sublo[2] || x[i][2] >= subhi[2]) {
      proclist[nsendatom] = coord2proc(x[i],igx,igy,igz);
      if (proclist[nsendatom] != me) {
        if (nsend > maxsend) grow_send(nsend,1);
        sizes[nsendatom] = avec->pack_exchange(i,&buf_send[nsend]);
        nsend += sizes[nsendatom];
        nsendatom++;
        avec->copy(nlocal-1,i,1);
        nlocal--;
      } else i++;
    } else i++;
  }
  atom->nlocal = nlocal;

  // create irregular communication plan, perform comm, destroy plan
  // returned nrecv = size of buffer needed for incoming atoms

  int nrecv = create_atom(nsendatom,sizes,proclist,sortflag);
  if (nrecv > maxrecv) grow_recv(nrecv);
  exchange_atom(buf_send,sizes,buf_recv);
  destroy_atom();

  delete [] sizes;
  delete [] proclist;

  // add received atoms to my list

  int m = 0;
  while (m < nrecv) m += avec->unpack_exchange(&buf_recv[m]);

  // reset global->local map

  if (map_style) atom->map_set();
}

/* ----------------------------------------------------------------------
   check if any atoms need to migrate further than one proc away in any dim
   if not, caller can decide to use comm->exchange() instead
   atoms must be remapped to be inside simulation box before this is called
   for triclinic: atoms must be in lamda coords (0-1) before this is called
   return 1 if migrate required, 0 if not
------------------------------------------------------------------------- */

int Irregular::migrate_check()
{
  // subbox bounds for orthogonal or triclinic box
  // other comm/domain data used by coord2proc()

  double *sublo,*subhi;
  if (triclinic == 0) {
    sublo = domain->sublo;
    subhi = domain->subhi;
  } else {
    sublo = domain->sublo_lamda;
    subhi = domain->subhi_lamda;
  }

  uniform = comm->uniform;
  xsplit = comm->xsplit;
  ysplit = comm->ysplit;
  zsplit = comm->zsplit;
  boxlo = domain->boxlo;
  prd = domain->prd;

  // loop over atoms, check for any that are not in my sub-box
  // assign which proc it belongs to via coord2proc()
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
      coord2proc(x[i],igx,igy,igz);

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
   return total # of doubles I will recv (not including self)
------------------------------------------------------------------------- */

int Irregular::create_atom(int n, int *sizes, int *proclist, int sortflag)
{
  int i;

  // allocate plan and work vectors

  if (aplan) destroy_atom();
  aplan = (PlanAtom *) memory->smalloc(sizeof(PlanAtom),"irregular:aplan");
  int *list = new int[nprocs];
  int *count = new int[nprocs];

  // nrecv = # of messages I receive

  for (i = 0; i < nprocs; i++) {
    list[i] = 0;
    count[i] = 1;
  }
  for (i = 0; i < n; i++) list[proclist[i]] = 1;

  int nrecv;
  MPI_Reduce_scatter(list,&nrecv,count,MPI_INT,MPI_SUM,world);

  // allocate receive arrays

  int *proc_recv = new int[nrecv];
  int *length_recv = new int[nrecv];
  MPI_Request *request = new MPI_Request[nrecv];
  MPI_Status *status = new MPI_Status[nrecv];

  // nsend = # of messages I send

  for (i = 0; i < nprocs; i++) list[i] = 0;
  for (i = 0; i < n; i++) list[proclist[i]] += sizes[i];

  int nsend = 0;
  for (i = 0; i < nprocs; i++)
    if (list[i]) nsend++;

  // allocate send arrays

  int *proc_send = new int[nsend];
  int *length_send = new int[nsend];
  int *num_send = new int[nsend];
  int *index_send = new int[n];
  int *offset_send = new int[n];

  // list still stores size of message for procs I send to
  // proc_send = procs I send to
  // length_send = total size of message I send to each proc
  // to balance pattern of send messages:
  //   each proc begins with iproc > me, continues until iproc = me
  // reset list to store which send message each proc corresponds to

  int iproc = me;
  int isend = 0;
  for (i = 0; i < nprocs; i++) {
    iproc++;
    if (iproc == nprocs) iproc = 0;
    if (list[iproc] > 0) {
      proc_send[isend] = iproc;
      length_send[isend] = list[iproc];
      list[iproc] = isend;
      isend++;
    }
  }

  // num_send = # of atoms I send to each proc

  for (i = 0; i < nsend; i++) num_send[i] = 0;
  for (i = 0; i < n; i++) {
    isend = list[proclist[i]];
    num_send[isend]++;
  }

  // count = offsets into index_send for each proc I send to
  // index_send = list of which atoms to send to each proc
  //   1st N1 values are atom indices for 1st proc,
  //   next N2 values are atom indices for 2nd proc, etc
  // offset_send = where each atom starts in send buffer

  count[0] = 0;
  for (i = 1; i < nsend; i++) count[i] = count[i-1] + num_send[i-1];

  for (i = 0; i < n; i++) {
    isend = list[proclist[i]];
    index_send[count[isend]++] = i;
    if (i) offset_send[i] = offset_send[i-1] + sizes[i-1];
    else offset_send[i] = 0;
  }

  // tell receivers how much data I send
  // sendmax = largest # of doubles I send in a single message

  int sendmax = 0;
  for (i = 0; i < nsend; i++) {
    MPI_Send(&length_send[i],1,MPI_INT,proc_send[i],0,world);
    sendmax = MAX(sendmax,length_send[i]);
  }

  // receive incoming messages
  // proc_recv = procs I recv from
  // length_recv = total size of message each proc sends me
  // nrecvsize = total size of data I recv

  int nrecvsize = 0;
  for (i = 0; i < nrecv; i++) {
    MPI_Recv(&length_recv[i],1,MPI_INT,MPI_ANY_SOURCE,0,world,status);
    proc_recv[i] = status->MPI_SOURCE;
    nrecvsize += length_recv[i];
  }

  // sort proc_recv and num_recv by proc ID if requested
  // useful for debugging to insure reproducible ordering of received atoms
  // invoke by adding final arg = 1 to create_atom() call in migrate_atoms()

  if (sortflag) {
    int *order = new int[nrecv];
    int *proc_recv_ordered = new int[nrecv];
    int *length_recv_ordered = new int[nrecv];

    for (i = 0; i < nrecv; i++) order[i] = i;
    proc_recv_copy = proc_recv;
    qsort(order,nrecv,sizeof(int),compare_standalone);

    int j;
    for (i = 0; i < nrecv; i++) {
      j = order[i];
      proc_recv_ordered[i] = proc_recv[j];
      length_recv_ordered[i] = length_recv[j];
    }

    memcpy(proc_recv,proc_recv_ordered,nrecv*sizeof(int));
    memcpy(length_recv,length_recv_ordered,nrecv*sizeof(int));
    delete [] order;
    delete [] proc_recv_ordered;
    delete [] length_recv_ordered;
  }

  // barrier to insure all MPI_ANY_SOURCE messages are received
  // else another proc could proceed to exchange_atom() and send to me

  MPI_Barrier(world);

  // free work vectors

  delete [] count;
  delete [] list;

  // initialize plan

  aplan->nsend = nsend;
  aplan->nrecv = nrecv;
  aplan->sendmax = sendmax;

  aplan->proc_send = proc_send;
  aplan->length_send = length_send;
  aplan->num_send = num_send;
  aplan->index_send = index_send;
  aplan->offset_send = offset_send;
  aplan->proc_recv = proc_recv;
  aplan->length_recv = length_recv;

  aplan->request = request;
  aplan->status = status;

  return nrecvsize;
}

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

/* ----------------------------------------------------------------------
   communicate atoms via PlanAtom
   sendbuf = list of atoms to send
   sizes = # of doubles for each atom
   recvbuf = received atoms
------------------------------------------------------------------------- */

void Irregular::exchange_atom(double *sendbuf, int *sizes, double *recvbuf)
{
  int i,m,n,offset,num_send;

  // post all receives

  offset = 0;
  for (int irecv = 0; irecv < aplan->nrecv; irecv++) {
    MPI_Irecv(&recvbuf[offset],aplan->length_recv[irecv],MPI_DOUBLE,
              aplan->proc_recv[irecv],0,world,&aplan->request[irecv]);
    offset += aplan->length_recv[irecv];
  }

  // allocate buf for largest send

  double *buf;
  memory->create(buf,aplan->sendmax,"irregular:buf");

  // send each message
  // pack buf with list of atoms
  // m = index of atom in sendbuf

  int *index_send = aplan->index_send;
  int nsend = aplan->nsend;
  n = 0;

  for (int isend = 0; isend < nsend; isend++) {
    offset = 0;
    num_send = aplan->num_send[isend];
    for (i = 0; i < num_send; i++) {
      m = index_send[n++];
      memcpy(&buf[offset],&sendbuf[aplan->offset_send[m]],
             sizes[m]*sizeof(double));
      offset += sizes[m];
    }
    MPI_Send(buf,aplan->length_send[isend],MPI_DOUBLE,
             aplan->proc_send[isend],0,world);
  }

  // free temporary send buffer

  memory->destroy(buf);

  // wait on all incoming messages

  if (aplan->nrecv) MPI_Waitall(aplan->nrecv,aplan->request,aplan->status);
}

/* ----------------------------------------------------------------------
   destroy communication plan for atoms
------------------------------------------------------------------------- */

void Irregular::destroy_atom()
{
  delete [] aplan->proc_send;
  delete [] aplan->length_send;
  delete [] aplan->num_send;
  delete [] aplan->index_send;
  delete [] aplan->offset_send;
  delete [] aplan->proc_recv;
  delete [] aplan->length_recv;
  delete [] aplan->request;
  delete [] aplan->status;
  memory->sfree(aplan);
  aplan = NULL;
}

/* ----------------------------------------------------------------------
   create a communication plan for datums
   n = # of datums to send
   proclist = proc to send each datum to (including self)
   return total # of datums I will recv (including self)
------------------------------------------------------------------------- */

int Irregular::create_data(int n, int *proclist)
{
  int i,m;

  // allocate plan and work vectors

  dplan = (PlanData *) memory->smalloc(sizeof(PlanData),"irregular:dplan");
  int *list = new int[nprocs];
  int *count = new int[nprocs];

  // nrecv = # of messages I receive

  for (i = 0; i < nprocs; i++) {
    list[i] = 0;
    count[i] = 1;
  }
  for (i = 0; i < n; i++) list[proclist[i]] = 1;

  int nrecv;
  MPI_Reduce_scatter(list,&nrecv,count,MPI_INT,MPI_SUM,world);
  if (list[me]) nrecv--;

  // allocate receive arrays

  int *proc_recv = new int[nrecv];
  int *num_recv = new int[nrecv];
  MPI_Request *request = new MPI_Request[nrecv];
  MPI_Status *status = new MPI_Status[nrecv];

  // nsend = # of messages I send

  for (i = 0; i < nprocs; i++) list[i] = 0;
  for (i = 0; i < n; i++) list[proclist[i]]++;

  int nsend = 0;
  for (i = 0; i < nprocs; i++)
    if (list[i]) nsend++;
  if (list[me]) nsend--;

  // allocate send and self arrays

  int *proc_send = new int[nsend];
  int *num_send = new int[nsend];
  int *index_send = new int[n-list[me]];
  int *index_self = new int[list[me]];

  // proc_send = procs I send to
  // num_send = # of datums I send to each proc
  // num_self = # of datums I copy to self
  // to balance pattern of send messages:
  //   each proc begins with iproc > me, continues until iproc = me
  // reset list to store which send message each proc corresponds to

  int num_self;

  int iproc = me;
  int isend = 0;
  for (i = 0; i < nprocs; i++) {
    iproc++;
    if (iproc == nprocs) iproc = 0;
    if (iproc == me) num_self = list[iproc];
    else if (list[iproc] > 0) {
      proc_send[isend] = iproc;
      num_send[isend] = list[iproc];
      list[iproc] = isend;
      isend++;
    }
  }
  list[me] = 0;

  // count = offsets into index_send for each proc I send to
  // m = ptr into index_self
  // index_send = list of which datums to send to each proc
  //   1st N1 values are datum indices for 1st proc,
  //   next N2 values are datum indices for 2nd proc, etc

  count[0] = 0;
  for (i = 1; i < nsend; i++) count[i] = count[i-1] + num_send[i-1];

  m = 0;
  for (i = 0; i < n; i++) {
    iproc = proclist[i];
    if (iproc == me) index_self[m++] = i;
    else {
      isend = list[iproc];
      index_send[count[isend]++] = i;
    }
  }

  // tell receivers how much data I send
  // sendmax = largest # of datums I send in a single message

  int sendmax = 0;
  for (i = 0; i < nsend; i++) {
    MPI_Send(&num_send[i],1,MPI_INT,proc_send[i],0,world);
    sendmax = MAX(sendmax,num_send[i]);
  }

  // receive incoming messages
  // proc_recv = procs I recv from
  // num_recv = total size of message each proc sends me
  // nrecvsize = total size of data I recv

  int nrecvsize = 0;
  for (i = 0; i < nrecv; i++) {
    MPI_Recv(&num_recv[i],1,MPI_INT,MPI_ANY_SOURCE,0,world,status);
    proc_recv[i] = status->MPI_SOURCE;
    nrecvsize += num_recv[i];
  }
  nrecvsize += num_self;

  // barrier to insure all MPI_ANY_SOURCE messages are received
  // else another proc could proceed to exchange_data() and send to me

  MPI_Barrier(world);

  // free work vectors

  delete [] count;
  delete [] list;

  // initialize plan and return it

  dplan->nsend = nsend;
  dplan->nrecv = nrecv;
  dplan->sendmax = sendmax;

  dplan->proc_send = proc_send;
  dplan->num_send = num_send;
  dplan->index_send = index_send;
  dplan->proc_recv = proc_recv;
  dplan->num_recv = num_recv;
  dplan->num_self = num_self;
  dplan->index_self = index_self;

  dplan->request = request;
  dplan->status = status;

  return nrecvsize;
}

/* ----------------------------------------------------------------------
   communicate datums via PlanData
   sendbuf = list of datums to send
   nbytes = size of each datum
   recvbuf = received datums (including copied from me)
------------------------------------------------------------------------- */

void Irregular::exchange_data(char *sendbuf, int nbytes, char *recvbuf)
{
  int i,m,n,offset,num_send;

  // post all receives, starting after self copies

  offset = dplan->num_self*nbytes;
  for (int irecv = 0; irecv < dplan->nrecv; irecv++) {
    MPI_Irecv(&recvbuf[offset],dplan->num_recv[irecv]*nbytes,MPI_CHAR,
              dplan->proc_recv[irecv],0,world,&dplan->request[irecv]);
    offset += dplan->num_recv[irecv]*nbytes;
  }

  // allocate buf for largest send

  char *buf;
  memory->create(buf,dplan->sendmax*nbytes,"irregular:buf");

  // send each message
  // pack buf with list of datums
  // m = index of datum in sendbuf

  int *index_send = dplan->index_send;
  int nsend = dplan->nsend;
  n = 0;

  for (int isend = 0; isend < nsend; isend++) {
    num_send = dplan->num_send[isend];
    for (i = 0; i < num_send; i++) {
      m = index_send[n++];
      memcpy(&buf[i*nbytes],&sendbuf[m*nbytes],nbytes);
    }
    MPI_Send(buf,dplan->num_send[isend]*nbytes,MPI_CHAR,
             dplan->proc_send[isend],0,world);
  }

  // free temporary send buffer

  memory->destroy(buf);

  // copy datums to self, put at beginning of recvbuf

  int *index_self = dplan->index_self;
  int num_self = dplan->num_self;

  for (i = 0; i < num_self; i++) {
    m = index_self[i];
    memcpy(&recvbuf[i*nbytes],&sendbuf[m*nbytes],nbytes);
  }

  // wait on all incoming messages

  if (dplan->nrecv) MPI_Waitall(dplan->nrecv,dplan->request,dplan->status);
}

/* ----------------------------------------------------------------------
   destroy communication plan for datums
------------------------------------------------------------------------- */

void Irregular::destroy_data()
{
  delete [] dplan->proc_send;
  delete [] dplan->num_send;
  delete [] dplan->index_send;
  delete [] dplan->proc_recv;
  delete [] dplan->num_recv;
  delete [] dplan->index_self;
  delete [] dplan->request;
  delete [] dplan->status;
  memory->sfree(dplan);
  dplan = NULL;
}

/* ----------------------------------------------------------------------
   determine which proc owns atom with coord x[3]
   x will be in box (orthogonal) or lamda coords (triclinic)
   for uniform = 1, directly calculate owning proc
   for non-uniform, iteratively find owning proc via binary search
   return owning proc ID via grid2proc
   return igx,igy,igz = logical grid loc of owing proc within 3d grid of procs
------------------------------------------------------------------------- */

int Irregular::coord2proc(double *x, int &igx, int &igy, int &igz)
{
  if (uniform) {
    if (triclinic == 0) {
      igx = static_cast<int> (procgrid[0] * (x[0]-boxlo[0]) / prd[0]);
      igy = static_cast<int> (procgrid[1] * (x[1]-boxlo[1]) / prd[1]);
      igz = static_cast<int> (procgrid[2] * (x[2]-boxlo[2]) / prd[2]);
    } else {
      igx = static_cast<int> (procgrid[0] * x[0]);
      igy = static_cast<int> (procgrid[1] * x[1]);
      igz = static_cast<int> (procgrid[2] * x[2]);
    }

  } else {
    if (triclinic == 0) {
      igx = binary((x[0]-boxlo[0])/prd[0],procgrid[0],xsplit);
      igy = binary((x[1]-boxlo[1])/prd[1],procgrid[1],ysplit);
      igz = binary((x[2]-boxlo[2])/prd[2],procgrid[2],zsplit);
    } else {
      igx = binary(x[0],procgrid[0],xsplit);
      igy = binary(x[1],procgrid[1],ysplit);
      igz = binary(x[2],procgrid[2],zsplit);
    }
  }

  if (igx < 0) igx = 0;
  if (igx >= procgrid[0]) igx = procgrid[0] - 1;
  if (igy < 0) igy = 0;
  if (igy >= procgrid[1]) igy = procgrid[1] - 1;
  if (igz < 0) igz = 0;
  if (igz >= procgrid[2]) igz = procgrid[2] - 1;

  return grid2proc[igx][igy][igz];
}

/* ----------------------------------------------------------------------
   binary search for value in N-length ascending vec
   value may be outside range of vec limits
   always return index from 0 to N-1 inclusive
   return 0 if value < vec[0]
   reutrn N-1 if value >= vec[N-1]
   return index = 1 to N-2 if vec[index] <= value < vec[index+1]
------------------------------------------------------------------------- */

int Irregular::binary(double value, int n, double *vec)
{
  int lo = 0;
  int hi = n-1;

  if (value < vec[lo]) return lo;
  if (value >= vec[hi]) return hi;

  // insure vec[lo] <= value < vec[hi] at every iteration
  // done when lo,hi are adjacent

  int index = (lo+hi)/2;
  while (lo < hi-1) {
    if (value < vec[index]) hi = index;
    else if (value >= vec[index]) lo = index;
    index = (lo+hi)/2;
  }

  return index;
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
  bigint bytes = memory->usage(buf_send,maxsend);
  bytes += memory->usage(buf_recv,maxrecv);
  return bytes;
}
