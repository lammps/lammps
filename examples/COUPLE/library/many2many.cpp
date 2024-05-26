#include <mpi.h>
#include <cstdlib>
#include <cstdio>
#include "many2many.h"
#include "irregular.h"
#include "memorylib.h"
#include "errorlib.h"

#include <map>

#define MAX(A,B) ((A) > (B)) ? (A) : (B)

/* ---------------------------------------------------------------------- */

Many2Many::Many2Many(MPI_Comm caller)
{
  comm = caller;
  MPI_Comm_rank(comm,&me);
  MPI_Comm_size(comm,&nprocs);

  memory = new MemoryLib(comm);
  error = new ErrorLib(comm);

  src_own = dest_own = NULL;
  src_off = dest_off = NULL;
  src_iwork = dest_iwork = NULL;
  src_dwork = dest_dwork = NULL;
  irregular = NULL;
}

/* ---------------------------------------------------------------------- */

Many2Many::~Many2Many()
{
  delete memory;
  delete error;
  deallocate();
}

/* ----------------------------------------------------------------------
   create a many2many pattern, deallocating any previous pattern
   each proc will contribute nsrc items with IDs listed in id_src
   each proc will receive ndest items with IDs listed in id_dest
   only sets up communication via rendezvous algorithm and Irregular class
     if id_src does not match id_dest on all procs
------------------------------------------------------------------------- */

void Many2Many::setup(int nsrc, int *id_src, int ndest, int *id_dest)
{
  int i,j,isrc,idest,nsend,nrecv;
  int *proclist,*work;
  std::map<int,int> hash;
  std::map<int,int>::iterator loc;

  // free any previous many2many info

  deallocate();

  // allocate on-proc and off-proc index lists

  src_own = 
    (int *) memory->smalloc(nsrc*sizeof(int),"many2many:src_own");
  dest_own = 
    (int *) memory->smalloc(ndest*sizeof(int),"many2many:dest_own");
  src_off = 
    (int *) memory->smalloc(nsrc*sizeof(int),"many2many:src_off");
  dest_off = 
    (int *) memory->smalloc(ndest*sizeof(int),"many2many:dest_off");

  // store destination IDs in hash

  for (int i = 0; i < ndest; i++)
    hash.insert(std::pair<int,int> (id_dest[i],i));

  // src_own, dest_own = list of IDs in both src and dest
  // nsrc_off = list of IDs in src owned by other procs

  nown = nsrc_off = 0;
  nsrc_off = 0;
  for (i = 0; i < nsrc; i++) {
    loc = hash.find(id_src[i]);
    if (loc != hash.end()) {
      src_own[nown] = i;
      dest_own[nown] = loc->second;
      nown++;
    } else src_off[nsrc_off++] = i;
  }

  // all done if all procs have one-to-one mapping of src and dest IDs
  // else figure out irregular comm needed

  int flag = 0;
  if (nown == nsrc && nown == ndest) flag = 1;
  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_MIN,comm);
  if (flagall) return;

  // ndest_off = list of IDs in dest owned by other procs

  work = (int *) memory->smalloc(ndest*sizeof(int),"many2many:work");

  for (i = 0; i < ndest; i++) work[i] = 0;
  for (i = 0; i < nown; i++) work[dest_own[i]] = 1;

  ndest_off = 0;
  for (i = 0; i < ndest; i++)
    if (work[i] == 0) dest_off[ndest_off++] = i;

  memory->sfree(work);

  // realloc off-proc arrays to smaller size

  src_off = (int *) 
    memory->srealloc(src_off,nsrc_off*sizeof(int),"many2many:src_off");
  dest_off = (int *) 
    memory->srealloc(dest_off,ndest_off*sizeof(int),"many2many:dest_off");

  // send off-proc src and dest Datums to 3rd-party proc via irregular comm
  // proc = ID % nprocs

  nsend = nsrc_off + ndest_off;
  proclist = new int[nsend];
  Datum1 *send1 = new Datum1[nsend];

  for (i = 0; i < nsrc_off; i++) {
    proclist[i] = id_src[src_off[i]] % nprocs;
    send1[i].id = id_src[src_off[i]];
    send1[i].proc = me;
    send1[i].index = src_off[i];
  }
  for (i = 0, j = nsrc_off; i < ndest_off; i++, j++) {
    proclist[j] = id_dest[dest_off[i]] % nprocs;
    send1[j].id = -id_dest[dest_off[i]];
    send1[j].proc = me;
    send1[j].index = dest_off[i];
  }

  irregular = new Irregular(comm);
  irregular->pattern(nsend,proclist);
  nrecv = irregular->size(sizeof(Datum1)) / sizeof(Datum1);
  Datum1 *recv1 = new Datum1[nrecv];
  irregular->exchange((char *) send1, (char *) recv1);
  delete irregular;
  delete [] proclist;

  // as 3rd-party proc, now have matching pairs of off-proc IDs
  // store src IDs (which are positive) in hash
  // loop over dest IDs (which are negative) to find matches
  // send match info back to src procs via a 2nd irregular comm

  nsend = nrecv/2;
  proclist = new int[nsend];
  Datum2 *send2 = new Datum2[nsend];
  nsend = 0;

  hash.clear();
  for (isrc = 0; isrc < nrecv; isrc++)
    if (recv1[isrc].id > 0)
      hash.insert(std::pair<int,int> (recv1[isrc].id,isrc));

  for (idest = 0; idest < nrecv; idest++)
    if (recv1[idest].id < 0) {
      loc = hash.find(-recv1[idest].id);
      if (loc != hash.end()) {
	isrc = loc->second;
	proclist[nsend] = recv1[isrc].proc;
	send2[nsend].slocal = recv1[isrc].index;
	send2[nsend].dlocal = recv1[idest].index;
	send2[nsend].dproc = recv1[idest].proc;
	nsend++;
      } else error->one("Did not receive matching src/dest ID");
    }

  irregular = new Irregular(comm);
  irregular->pattern(nsend,proclist);
  nrecv = irregular->size(sizeof(Datum2)) / sizeof(Datum2);
  Datum2 *recv2 = new Datum2[nrecv];
  irregular->exchange((char *) send2, (char *) recv2);
  delete irregular;
  delete [] proclist;

  // use list of received src->dest Datums to build final irregular commm
  // irregular comm will communicate off-proc info from src to dest directly
  // work = local indices of dest IDs to send initially

  nsend = nrecv;
  proclist = new int[nsend];
  work = new int[nsend];

  for (i = 0; i < nrecv; i++) {
    src_off[i] = recv2[i].slocal;
    work[i] = recv2[i].dlocal;
    proclist[i] = recv2[i].dproc;
  }

  irregular = new Irregular(comm);
  irregular->pattern(nsend,proclist);

  // send receiver's local indices
  // receiver stores them as indirection list in dest_off

  nrecv = irregular->size(sizeof(int)) / sizeof(int);
  irregular->exchange((char *) work, (char *) dest_off);

  delete [] proclist;
  delete [] work;

  // create work arrays for data exchange of int/double data

  src_iwork = 
    (int *) memory->smalloc(nsrc_off*sizeof(int),"many2many:src_iwork");
  dest_iwork = 
    (int *) memory->smalloc(ndest_off*sizeof(int),"many2many:dest_iwork");
  src_dwork = 
    (double *) memory->smalloc(nsrc_off*sizeof(double),"many2many:src_dwork");
  dest_dwork = 
    (double *) memory->smalloc(ndest_off*sizeof(double),
			       "many2many:dest_dwork");

  // clean up

  delete [] send1;
  delete [] recv1;
  delete [] send2;
  delete [] recv2;

  // debug checks for full coverage of srd/dest - can delete eventually

  work = new int[MAX(nsrc,ndest)];

  for (i = 0; i < nsrc; i++) work[i] = 0;
  for (i = 0; i < nown; i++) work[src_own[i]]++;
  for (i = 0; i < nsrc_off; i++) work[src_off[i]]++;
  for (i = 0; i < nsrc; i++)
    if (work[i] != 1) {
      char str[128];
      sprintf(str,"BAD SRC %d: %d %d\n",me,i,work[i]);
      error->one(str);
    }

  for (i = 0; i < ndest; i++) work[i] = 0;
  for (i = 0; i < nown; i++) work[dest_own[i]]++;
  for (i = 0; i < ndest_off; i++) work[dest_off[i]]++;
  for (i = 0; i < ndest; i++)
    if (work[i] != 1) {
      char str[128];
      sprintf(str,"BAD DEST %d: %d %d\n",me,i,work[i]);
      error->one(str);
    }

  delete [] work;
}

/* ----------------------------------------------------------------------
   transfer one src entity to dest entity, matched by IDs in create()
   operates on an int vector
------------------------------------------------------------------------- */

void Many2Many::exchange(int *src, int *dest)
{
  int i;

  // copy on-proc info

  for (i = 0; i < nown; i++)
    dest[dest_own[i]] = src[src_own[i]];

  // communicate off-proc info
  // user src_off and dest_off to pack/unpack data

  if (irregular) {
    int nrecv = irregular->size(sizeof(int)) / sizeof(int);
    for (i = 0; i < nsrc_off; i++) src_iwork[i] = src[src_off[i]];
    irregular->exchange((char *) src_iwork, (char *) dest_iwork);
    for (i = 0; i < ndest_off; i++) dest[dest_off[i]] = dest_iwork[i];
  }
}

/* ----------------------------------------------------------------------
   transfer one src entity to dest entity, matched by IDs in create()
   operates on a double vector
------------------------------------------------------------------------- */

void Many2Many::exchange(double *src, double *dest)
{
  int i;

  // copy on-proc info

  for (int i = 0; i < nown; i++)
    dest[dest_own[i]] = src[src_own[i]];

  // communicate off-proc info
  // user src_off and dest_off to pack/unpack data

  if (irregular) {
    int nrecv = irregular->size(sizeof(double)) / sizeof(double);
    for (i = 0; i < nsrc_off; i++) src_dwork[i] = src[src_off[i]];
    irregular->exchange((char *) src_dwork, (char *) dest_dwork);
    for (i = 0; i < ndest_off; i++) dest[dest_off[i]] = dest_dwork[i];
  }
}

/* ---------------------------------------------------------------------- */

void Many2Many::deallocate()
{
  memory->sfree(src_own);
  memory->sfree(dest_own);
  memory->sfree(src_off);
  memory->sfree(dest_off);
  memory->sfree(src_iwork);
  memory->sfree(dest_iwork);
  memory->sfree(src_dwork);
  memory->sfree(dest_dwork);

  delete irregular;
}
