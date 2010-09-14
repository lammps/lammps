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

#define BUFFACTOR 1.5
#define BUFMIN 1000
#define BUFEXTRA 1000

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

Irregular::Irregular(LAMMPS *lmp) : Pointers(lmp)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  triclinic = domain->triclinic;
  map_style = atom->map_style;
  procgrid = comm->procgrid;
  grid2proc = comm->grid2proc;

  // initialize comm buffers

  maxsend = BUFMIN;
  buf_send = (double *) 
    memory->smalloc((maxsend+BUFEXTRA)*sizeof(double),"comm:buf_send");
  maxrecv = BUFMIN;
  buf_recv = (double *) 
    memory->smalloc(maxrecv*sizeof(double),"comm:buf_recv");
}

/* ---------------------------------------------------------------------- */

Irregular::~Irregular()
{
  memory->sfree(buf_send);
  memory->sfree(buf_recv);
}

/* ----------------------------------------------------------------------
   communicate atoms to new owning procs via irregular communication
   can be used in place of comm->exchange()
   unlike exchange(), allows atoms to have moved arbitrarily long distances
   first setup irregular comm pattern, then invoke it
   atoms must be remapped to be inside simulation box before this is called
   for triclinic: atoms must be in lamda coords (0-1) before this is called
------------------------------------------------------------------------- */

void Irregular::migrate_atoms()
{
  // clear global->local map since atoms move to new procs
  // zero out ghosts so map_set() at end will operate only on local atoms
  // exchange() doesn't need to zero ghosts b/c borders()
  //   is called right after and it zeroes ghosts and calls map_set()

  if (map_style) atom->map_clear();
  atom->nghost = 0;

  // subbox bounds for orthogonal or triclinic

  double *sublo,*subhi;
  if (triclinic == 0) {
    sublo = domain->sublo;
    subhi = domain->subhi;
  } else {
    sublo = domain->sublo_lamda;
    subhi = domain->subhi_lamda;
  }

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

  int i = 0;
  while (i < nlocal) {
    if (x[i][0] < sublo[0] || x[i][0] >= subhi[0] ||
	x[i][1] < sublo[1] || x[i][1] >= subhi[1] ||
	x[i][2] < sublo[2] || x[i][2] >= subhi[2]) {
      proclist[nsendatom] = coord2proc(x[i]);
      if (proclist[nsendatom] != me) {
	if (nsend > maxsend) grow_send(nsend,1);
	sizes[nsendatom] = avec->pack_exchange(i,&buf_send[nsend]);
	nsend += sizes[nsendatom];
	nsendatom++;
	avec->copy(nlocal-1,i);
	nlocal--;
      } else i++;
    } else i++;
  }
  atom->nlocal = nlocal;

  // create irregular communication plan, perform comm, destroy plan
  // returned nrecv = size of buffer needed for incoming atoms

  int nrecv;
  PlanAtom *plan = create_atom(nsendatom,sizes,proclist,&nrecv);
  if (nrecv > maxrecv) grow_recv(nrecv);
  exchange_atom(plan,buf_send,sizes,buf_recv);
  destroy_atom(plan);

  delete [] sizes;
  delete [] proclist;

  // add received atoms to my list

  int m = 0;
  while (m < nrecv) m += avec->unpack_exchange(&buf_recv[m]);

  // reset global->local map

  if (map_style) atom->map_set();
}

/* ----------------------------------------------------------------------
   create an irregular communication plan for atoms
   n = # of atoms to send
   sizes = # of doubles for each atom
   proclist = proc to send each atom to (none to me)
   return nrecvsize = total # of doubles I will recv
------------------------------------------------------------------------- */

Irregular::PlanAtom *Irregular::create_atom(int n, int *sizes, 
					    int *proclist, int *nrecvsize)
{
  int i;

  // allocate plan and work vectors

  PlanAtom *plan = (struct PlanAtom *)
    memory->smalloc(sizeof(PlanAtom),"irregular:plan");
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

  // num_send = # of datums I send to each proc

  for (i = 0; i < nsend; i++) num_send[i] = 0;
  for (i = 0; i < n; i++) {
    isend = list[proclist[i]];
    num_send[isend]++;
  }

  // count = offsets into n-length index_send for each proc I send to
  // index_send = list of which datums to send to each proc
  //   1st N1 values are datum indices for 1st proc,
  //   next N2 values are datum indices for 2nd proc, etc
  // offset_send = where each datum starts in send buffer

  count[0] = 0;
  for (i = 1; i < nsend; i++) count[i] = count[i-1] + num_send[i-1];

  for (i = 0; i < n; i++) {
    isend = list[proclist[i]];
    index_send[count[isend]++] = i;
    if (i) offset_send[i] = offset_send[i-1] + sizes[i-1];
    else offset_send[i] = 0;
  }

  // tell receivers how much data I send
  // sendmax = largest # of datums I send in a single message

  int sendmax = 0;
  for (i = 0; i < nsend; i++) {
    MPI_Send(&length_send[i],1,MPI_INT,proc_send[i],0,world);
    sendmax = MAX(sendmax,length_send[i]);
  }

  // receive incoming messages
  // proc_recv = procs I recv from
  // length_recv = total size of message each proc sends me
  // nrecvsize = total size of data I recv

  *nrecvsize = 0;
  for (i = 0; i < nrecv; i++) {
    MPI_Recv(&length_recv[i],1,MPI_INT,MPI_ANY_SOURCE,0,world,status);
    proc_recv[i] = status->MPI_SOURCE;
    *nrecvsize += length_recv[i];
  }

  // barrier to insure all MPI_ANY_SOURCE messages are received
  // else another proc could proceed to irregular_perform() and send to me

  MPI_Barrier(world);

  // free work vectors

  delete [] count;
  delete [] list;
    
  // initialize plan and return it

  plan->nsend = nsend;
  plan->nrecv = nrecv;
  plan->sendmax = sendmax;

  plan->proc_send = proc_send;
  plan->length_send = length_send;
  plan->num_send = num_send;
  plan->index_send = index_send;
  plan->offset_send = offset_send;
  plan->proc_recv = proc_recv;
  plan->length_recv = length_recv;
  plan->request = request;
  plan->status = status;

  return plan;
}

/* ----------------------------------------------------------------------
   create an irregular communication plan for datums
   n = # of datums to send
   proclist = proc to send each datum to (none to me)
   return nrecvsize = total # of datums I will recv
------------------------------------------------------------------------- */

Irregular::PlanData *Irregular::create_data(int n, int *proclist,
					    int *nrecvsize)
{
  int i;

  // allocate plan and work vectors

  PlanData *plan = (struct PlanData *) 
    memory->smalloc(sizeof(PlanData),"irregular:plan");
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
  int *num_recv = new int[nrecv];
  MPI_Request *request = new MPI_Request[nrecv];
  MPI_Status *status = new MPI_Status[nrecv];

  // nsend = # of messages I send

  for (i = 0; i < nprocs; i++) list[i] = 0;
  for (i = 0; i < n; i++) list[proclist[i]]++;

  int nsend = 0;
  for (i = 0; i < nprocs; i++)
    if (list[i]) nsend++;

  // allocate send arrays

  int *proc_send = new int[nsend];
  int *num_send = new int[nsend];
  int *index_send = new int[n];

  // proc_send = procs I send to
  // num_send = # of datums I send to each proc
  // to balance pattern of send messages:
  //   each proc begins with iproc > me, continues until iproc = me

  int iproc = me;
  int isend = 0;
  for (i = 0; i < nprocs; i++) {
    iproc++;
    if (iproc == nprocs) iproc = 0;
    if (list[iproc] > 0) {
      proc_send[isend] = iproc;
      num_send[isend] = list[iproc];
      isend++;
    }
  }

  // count = offsets into n-length index_send for each proc I send to
  // index_send = list of which datums to send to each proc
  //   1st N1 values are datum indices for 1st proc,
  //   next N2 values are datum indices for 2nd proc, etc
  // offset_send = where each datum starts in send buffer

  count[0] = 0;
  for (i = 1; i < nsend; i++) count[i] = count[i-1] + num_send[i-1];

  for (i = 0; i < n; i++) {
    isend = list[proclist[i]];
    index_send[count[isend]++] = i;
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

  *nrecvsize = 0;
  for (i = 0; i < nrecv; i++) {
    MPI_Recv(&num_recv[i],1,MPI_INT,MPI_ANY_SOURCE,0,world,status);
    proc_recv[i] = status->MPI_SOURCE;
    *nrecvsize += num_recv[i];
  }

  // barrier to insure all MPI_ANY_SOURCE messages are received
  // else another proc could proceed to irregular_perform() and send to me

  MPI_Barrier(world);

  // free work vectors

  delete [] count;
  delete [] list;
    
  // initialize plan and return it

  plan->nsend = nsend;
  plan->nrecv = nrecv;
  plan->sendmax = sendmax;

  plan->proc_send = proc_send;
  plan->num_send = num_send;
  plan->index_send = index_send;
  plan->proc_recv = proc_recv;
  plan->num_recv = num_recv;
  plan->request = request;
  plan->status = status;

  return plan;
}

/* ----------------------------------------------------------------------
   perform irregular communication of atoms
   plan = previouly computed PlanAtom via create_atom()
   sendbuf = list of atoms to send
   sizes = # of doubles for each atom
   recvbuf = received atoms
------------------------------------------------------------------------- */

void Irregular::exchange_atom(PlanAtom *plan, double *sendbuf, int *sizes,
			      double *recvbuf)
{
  int i,m,n,offset;

  // post all receives

  offset = 0;
  for (int irecv = 0; irecv < plan->nrecv; irecv++) {
    MPI_Irecv(&recvbuf[offset],plan->length_recv[irecv],MPI_DOUBLE,
	      plan->proc_recv[irecv],0,world,&plan->request[irecv]);
    offset += plan->length_recv[irecv];
  }

  // allocate buf for largest send

  double *buf = (double *) memory->smalloc(plan->sendmax*sizeof(double),
					   "irregular:buf");

  // send each message
  // pack buf with list of datums (datum = one atom)
  // m = index of datum in sendbuf

  n = 0;
  for (int isend = 0; isend < plan->nsend; isend++) {
    offset = 0;
    for (i = 0; i < plan->num_send[isend]; i++) {
      m = plan->index_send[n++];
      memcpy(&buf[offset],&sendbuf[plan->offset_send[m]],
	     sizes[m]*sizeof(double));
      offset += sizes[m];
    }
    MPI_Send(buf,plan->length_send[isend],MPI_DOUBLE,
	     plan->proc_send[isend],0,world);
  }       

  // free temporary send buffer

  memory->sfree(buf);

  // wait on all incoming messages

  if (plan->nrecv) MPI_Waitall(plan->nrecv,plan->request,plan->status);
}

/* ----------------------------------------------------------------------
   perform irregular communication of datums
   plan = previouly computed PlanData via create_data()
   sendbuf = list of datums to send
   nbytes = size of each datum
   recvbuf = received datums
------------------------------------------------------------------------- */

void Irregular::exchange_data(PlanData *plan, char *sendbuf, int nbytes,
			      char *recvbuf)
{
  int i,m,n,offset;

  // post all receives

  offset = 0;
  for (int irecv = 0; irecv < plan->nrecv; irecv++) {
    MPI_Irecv(&recvbuf[offset],nbytes*plan->num_recv[irecv],MPI_CHAR,
	      plan->proc_recv[irecv],0,world,&plan->request[irecv]);
    offset += nbytes*plan->num_recv[irecv];
  }

  // allocate buf for largest send

  char *buf = (char *) memory->smalloc(plan->sendmax*nbytes,"irregular:buf");

  // send each message
  // pack buf with list of datums (datum = one atom)
  // m = index of datum in sendbuf

  n = 0;
  for (int isend = 0; isend < plan->nsend; isend++) {
    for (i = 0; i < plan->num_send[isend]; i++) {
      m = plan->index_send[n++];
      memcpy(&buf[i*nbytes],&sendbuf[m*nbytes],nbytes);
    }
    MPI_Send(buf,nbytes*plan->num_send[isend],MPI_CHAR,
	     plan->proc_send[isend],0,world);
  }       

  // free temporary send buffer

  memory->sfree(buf);

  // copy datums owned by self

  for (i = 0; i < plan->num_self; i++) {
    m = plan->index_self[i++];
    memcpy(&recvbuf[i*nbytes],&sendbuf[m*nbytes],nbytes);
  }

  // wait on all incoming messages

  if (plan->nrecv) MPI_Waitall(plan->nrecv,plan->request,plan->status);
}

/* ----------------------------------------------------------------------
   destroy an irregular communication plan for atoms
------------------------------------------------------------------------- */

void Irregular::destroy_atom(PlanAtom *plan)
{
  delete [] plan->proc_send;
  delete [] plan->length_send;
  delete [] plan->num_send;
  delete [] plan->index_send;
  delete [] plan->offset_send;
  delete [] plan->proc_recv;
  delete [] plan->length_recv;
  delete [] plan->request;
  delete [] plan->status;
  memory->sfree(plan);
}

/* ----------------------------------------------------------------------
   destroy an irregular communication plan for datums
------------------------------------------------------------------------- */

void Irregular::destroy_data(PlanData *plan)
{
  delete [] plan->proc_send;
  delete [] plan->num_send;
  delete [] plan->index_send;
  delete [] plan->proc_recv;
  delete [] plan->num_recv;
  delete [] plan->request;
  delete [] plan->status;
  memory->sfree(plan);
}

/* ----------------------------------------------------------------------
   determine which proc owns atom with x coord
   x will be in box (orthogonal) or lamda coords (triclinic)
------------------------------------------------------------------------- */

int Irregular::coord2proc(double *x)
{
  int loc[3];
  if (triclinic == 0) {
    double *boxlo = domain->boxlo;
    double *boxhi = domain->boxhi;
    loc[0] = static_cast<int>
      (procgrid[0] * (x[0]-boxlo[0]) / (boxhi[0]-boxlo[0]));
    loc[1] = static_cast<int>
      (procgrid[1] * (x[1]-boxlo[1]) / (boxhi[1]-boxlo[1]));
    loc[2] = static_cast<int>
      (procgrid[2] * (x[2]-boxlo[2]) / (boxhi[2]-boxlo[2]));
  } else {
    loc[0] = static_cast<int> (procgrid[0] * x[0]);
    loc[1] = static_cast<int> (procgrid[1] * x[1]);
    loc[2] = static_cast<int> (procgrid[2] * x[2]);
  }

  if (loc[0] < 0) loc[0] = 0;
  if (loc[0] >= procgrid[0]) loc[0] = procgrid[0] - 1;
  if (loc[1] < 0) loc[1] = 0;
  if (loc[1] >= procgrid[1]) loc[1] = procgrid[1] - 1;
  if (loc[2] < 0) loc[2] = 0;
  if (loc[2] >= procgrid[2]) loc[2] = procgrid[2] - 1;

  return grid2proc[loc[0]][loc[1]][loc[2]];
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
    buf_send = (double *) 
      memory->srealloc(buf_send,(maxsend+BUFEXTRA)*sizeof(double),
		       "comm:buf_send");
  else {
    memory->sfree(buf_send);
    buf_send = (double *) memory->smalloc((maxsend+BUFEXTRA)*sizeof(double),
					  "comm:buf_send");
  }
}

/* ----------------------------------------------------------------------
   free/malloc the size of the recv buffer as needed with BUFFACTOR 
------------------------------------------------------------------------- */

void Irregular::grow_recv(int n)
{
  maxrecv = static_cast<int> (BUFFACTOR * n);
  memory->sfree(buf_recv);
  buf_recv = (double *) memory->smalloc(maxrecv*sizeof(double),
					"comm:buf_recv");
}

