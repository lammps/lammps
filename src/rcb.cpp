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
#include "string.h"
#include "rcb.h"
#include "irregular.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MYHUGE 1.0e30
#define TINY 1.0e-6
#define DELTA 16384

// prototypes for non-class functions

void box_merge(void *, void *, int *, MPI_Datatype *);
void median_merge(void *, void *, int *, MPI_Datatype *);

// NOTE: if want to have reuse flag, need to sum Tree across procs

/* ---------------------------------------------------------------------- */

RCB::RCB(LAMMPS *lmp) : Pointers(lmp)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  ndot = maxdot = 0;
  dots = NULL;

  nlist = maxlist = 0;
  dotlist = dotmark = NULL;

  maxbuf = 0;
  buf = NULL;

  maxrecv = maxsend = 0;
  recvproc = recvindex = sendproc = sendindex = NULL;

  tree = NULL;
  irregular = NULL;

  // create MPI data and function types for box and median AllReduce ops

  MPI_Type_contiguous(6,MPI_DOUBLE,&box_type);
  MPI_Type_commit(&box_type);
  MPI_Type_contiguous(sizeof(Median),MPI_CHAR,&med_type);
  MPI_Type_commit(&med_type);

  MPI_Op_create(box_merge,1,&box_op);
  MPI_Op_create(median_merge,1,&med_op);

  reuse = 0;
}

/* ---------------------------------------------------------------------- */

RCB::~RCB()
{
  memory->sfree(dots);
  memory->destroy(dotlist);
  memory->destroy(dotmark);
  memory->sfree(buf);

  memory->destroy(recvproc);
  memory->destroy(recvindex);
  memory->destroy(sendproc);
  memory->destroy(sendindex);

  memory->sfree(tree);
  delete irregular;

  MPI_Type_free(&med_type);
  MPI_Type_free(&box_type);
  MPI_Op_free(&box_op);
  MPI_Op_free(&med_op);
}

/* ----------------------------------------------------------------------
   perform RCB balancing of N particles at coords X in bounding box LO/HI
   if wt = NULL, ignore per-particle weights
   if wt defined, per-particle weights > 0.0
   dimension = 2 or 3
   as documented in rcb.h:
     sets noriginal,nfinal,nkeep,recvproc,recvindex,lo,hi
   all proc particles will be inside or on surface of 3-d box
     defined by final lo/hi
   // NOTE: worry about re-use of data structs for fix balance?
------------------------------------------------------------------------- */

void RCB::compute(int dimension, int n, double **x, double *wt,
                  double *bboxlo, double *bboxhi)
{
  int i,j,k;
  int keep,outgoing,incoming,incoming2;
  int dim,markactive;
  int indexlo,indexhi;
  int first_iteration,breakflag;
  double wttot,wtlo,wthi,wtsum,wtok,wtupto,wtmax;
  double targetlo,targethi;
  double valuemin,valuemax,valuehalf;
  double tolerance;
  MPI_Comm comm,comm_half;
  MPI_Request request,request2;
  Median med,medme;

  // create list of my Dots

  ndot = nkeep = noriginal = n;

  if (ndot > maxdot) {
    maxdot = ndot;
    memory->sfree(dots);
    dots = (Dot *) memory->smalloc(ndot*sizeof(Dot),"RCB:dots");
  }

  for (i = 0; i < ndot; i++) {
    dots[i].x[0] = x[i][0];
    dots[i].x[1] = x[i][1];
    dots[i].x[2] = x[i][2];
    dots[i].proc = me;
    dots[i].index = i;
  }

  if (wt)
    for (i = 0; i < ndot; i++) dots[i].wt = wt[i];
  else
    for (i = 0; i < ndot; i++) dots[i].wt = 1.0;

  // initial bounding box = simulation box
  // includes periodic or shrink-wrapped boundaries
  
  lo = bbox.lo;
  hi = bbox.hi;

  lo[0] = bboxlo[0];
  lo[1] = bboxlo[1];
  lo[2] = bboxlo[2];
  hi[0] = bboxhi[0];
  hi[1] = bboxhi[1];
  hi[2] = bboxhi[2];

  cut = 0.0;
  cutdim = -1;

  // initialize counters

  counters[0] = 0;
  counters[1] = 0;
  counters[2] = 0;
  counters[3] = ndot;
  counters[4] = maxdot;
  counters[5] = 0;
  counters[6] = 0;

  // create communicator for use in recursion

  MPI_Comm_dup(world,&comm);

  // recurse until partition is a single proc = me
  // proclower,procupper = lower,upper procs in partition
  // procmid = 1st proc in upper half of partition

  int procpartner,procpartner2;

  int procmid;
  int proclower = 0;
  int procupper = nprocs - 1;

  while (proclower != procupper) {

    // if odd # of procs, lower partition gets extra one

    procmid = proclower + (procupper - proclower) / 2 + 1;

    // determine communication partner(s)
    // readnumber = # of proc partners to read from

    if (me < procmid)
      procpartner = me + (procmid - proclower);
    else
      procpartner = me - (procmid - proclower);
    
    int readnumber = 1;
    if (procpartner > procupper) {
      readnumber = 0;
      procpartner--;
    }
    if (me == procupper && procpartner != procmid - 1) {
      readnumber = 2;
      procpartner2 = procpartner + 1;
    }
    
    // wttot = summed weight of entire partition
    // search tolerance = largest single weight (plus epsilon)
    // targetlo = desired weight in lower half of partition
    // targethi = desired weight in upper half of partition

    wtmax = wtsum = 0.0;

    if (wt) {
      for (i = 0; i < ndot; i++) {
        wtsum += dots[i].wt;
        if (dots[i].wt > wtmax) wtmax = dots[i].wt;
      }
    } else {
      for (i = 0; i < ndot; i++) wtsum += dots[i].wt;
      wtmax = 1.0;
    }

    MPI_Allreduce(&wtsum,&wttot,1,MPI_DOUBLE,MPI_SUM,comm);
    if (wt) MPI_Allreduce(&wtmax,&tolerance,1,MPI_DOUBLE,MPI_MAX,comm);
    else tolerance = 1.0;

    tolerance *= 1.0 + TINY;
    targetlo = wttot * (procmid - proclower) / (procupper + 1 - proclower);
    targethi = wttot - targetlo;

    // dim = dimension to bisect on
    // do not allow choice of z dimension for 2d system

    dim = 0;
    if (hi[1]-lo[1] > hi[0]-lo[0]) dim = 1;
    if (dimension == 3) {
      if (dim == 0 && hi[2]-lo[2] > hi[0]-lo[0]) dim = 2;
      if (dim == 1 && hi[2]-lo[2] > hi[1]-lo[1]) dim = 2;
    }

    // create active list and mark array for dots
    // initialize active list to all dots

    if (ndot > maxlist) {
      memory->destroy(dotlist);
      memory->destroy(dotmark);
      maxlist = maxdot;
      memory->create(dotlist,maxlist,"RCB:dotlist");
      memory->create(dotmark,maxlist,"RCB:dotmark");
    }

    nlist = ndot;
    for (i = 0; i < nlist; i++) dotlist[i] = i;

    // median iteration
    // zoom in on bisector until correct # of dots in each half of partition
    // as each iteration of median-loop begins, require:
    //   all non-active dots are marked with 0/1 in dotmark
    //   valuemin <= every active dot <= valuemax
    //   wtlo, wthi = total wt of non-active dots
    // when leave median-loop, require only:
    //   valuehalf = correct cut position
    //   all dots <= valuehalf are marked with 0 in dotmark
    //   all dots >= valuehalf are marked with 1 in dotmark
    // markactive = which side of cut is active = 0/1
    // indexlo,indexhi = indices of dot closest to median

    wtlo = wthi = 0.0;
    valuemin = lo[dim];
    valuemax = hi[dim];
    first_iteration = 1;

    while (1) {

      // choose bisector value
      // use old value on 1st iteration if old cut dimension is the same
      // on 2nd option: could push valuehalf towards geometric center 
      //   with "1.0-factor" to force overshoot

      if (first_iteration && reuse && dim == tree[procmid].dim) {
	counters[5]++;
	valuehalf = tree[procmid].cut;
	if (valuehalf < valuemin || valuehalf > valuemax)
	  valuehalf = 0.5 * (valuemin + valuemax);	  
      } else if (wt)
	valuehalf = valuemin + (targetlo - wtlo) /
	  (wttot - wtlo - wthi) * (valuemax - valuemin);
      else
	valuehalf = 0.5 * (valuemin + valuemax);

      first_iteration = 0;
      
      // initialize local median data structure

      medme.totallo = medme.totalhi = 0.0;
      medme.valuelo = -MYHUGE;
      medme.valuehi = MYHUGE;
      medme.wtlo = medme.wthi = 0.0;
      medme.countlo = medme.counthi = 0;
      medme.proclo = medme.prochi = me;

      // mark all active dots on one side or other of bisector
      // also set all fields in median data struct
      // save indices of closest dots on either side

      for (j = 0; j < nlist; j++) {
	i = dotlist[j];
	if (dots[i].x[dim] <= valuehalf) {            // in lower part
	  medme.totallo += dots[i].wt;
	  dotmark[i] = 0;
	  if (dots[i].x[dim] > medme.valuelo) {       // my closest dot
	    medme.valuelo = dots[i].x[dim];
	    medme.wtlo = dots[i].wt;
	    medme.countlo = 1;
	    indexlo = i;
	  } else if (dots[i].x[dim] == medme.valuelo) {   // tied for closest
	    medme.wtlo += dots[i].wt;
	    medme.countlo++;
	  }
	}
	else {                                        // in upper part
	  medme.totalhi += dots[i].wt;
	  dotmark[i] = 1;
	  if (dots[i].x[dim] < medme.valuehi) {       // my closest dot
	    medme.valuehi = dots[i].x[dim];
	    medme.wthi = dots[i].wt;
	    medme.counthi = 1;
	    indexhi = i;
	  } else if (dots[i].x[dim] == medme.valuehi) {   // tied for closest
	    medme.wthi += dots[i].wt;
	    medme.counthi++;
	  }
	}
      }

      // combine median data struct across current subset of procs

      counters[0]++;
      MPI_Allreduce(&medme,&med,1,med_type,med_op,comm);

      // test median guess for convergence
      // move additional dots that are next to cut across it

      if (wtlo + med.totallo < targetlo) {    // lower half TOO SMALL

	wtlo += med.totallo;
	valuehalf = med.valuehi;

	if (med.counthi == 1) {                  // only one dot to move
	  if (wtlo + med.wthi < targetlo) {  // move it, keep iterating
	    if (me == med.prochi) dotmark[indexhi] = 0;
	  }
	  else {                                 // only move if beneficial
	    if (wtlo + med.wthi - targetlo < targetlo - wtlo)
	      if (me == med.prochi) dotmark[indexhi] = 0;
	    break;                               // all done
	  }
	}
	else {                                   // multiple dots to move
	  breakflag = 0;
	  wtok = 0.0;
	  if (medme.valuehi == med.valuehi) wtok = medme.wthi;   
	  if (wtlo + med.wthi >= targetlo) {                // all done
	    MPI_Scan(&wtok,&wtupto,1,MPI_DOUBLE,MPI_SUM,comm);
	    wtmax = targetlo - wtlo;
	    if (wtupto > wtmax) wtok = wtok - (wtupto - wtmax);
	    breakflag = 1;
	  }                                      // wtok = most I can move
	  for (j = 0, wtsum = 0.0; j < nlist && wtsum < wtok; j++) {
	    i = dotlist[j];
	    if (dots[i].x[dim] == med.valuehi) { // only move if better
	      if (wtsum + dots[i].wt - wtok < wtok - wtsum)
		dotmark[i] = 0;
	      wtsum += dots[i].wt;
	    }
	  }
	  if (breakflag) break;                   // done if moved enough
	}

	wtlo += med.wthi;
	if (targetlo-wtlo <= tolerance) break;  // close enough

	valuemin = med.valuehi;                   // iterate again
	markactive = 1;
      }

      else if (wthi + med.totalhi < targethi) {  // upper half TOO SMALL

	wthi += med.totalhi;
	valuehalf = med.valuelo;

	if (med.countlo == 1) {                  // only one dot to move
	  if (wthi + med.wtlo < targethi) {  // move it, keep iterating
	    if (me == med.proclo) dotmark[indexlo] = 1;
	  }
	  else {                                 // only move if beneficial
	    if (wthi + med.wtlo - targethi < targethi - wthi)
	      if (me == med.proclo) dotmark[indexlo] = 1;
	    break;                               // all done
	  }
	}
	else {                                   // multiple dots to move
	  breakflag = 0;
	  wtok = 0.0;
	  if (medme.valuelo == med.valuelo) wtok = medme.wtlo;   
	  if (wthi + med.wtlo >= targethi) {                // all done
	    MPI_Scan(&wtok,&wtupto,1,MPI_DOUBLE,MPI_SUM,comm);
	    wtmax = targethi - wthi;
	    if (wtupto > wtmax) wtok = wtok - (wtupto - wtmax);
	    breakflag = 1;
	  }                                      // wtok = most I can move
	  for (j = 0, wtsum = 0.0; j < nlist && wtsum < wtok; j++) {
	    i = dotlist[j];
	    if (dots[i].x[dim] == med.valuelo) { // only move if better
	      if (wtsum + dots[i].wt - wtok < wtok - wtsum) 
		dotmark[i] = 1;
	      wtsum += dots[i].wt;
	    }
	  }
	  if (breakflag) break;                   // done if moved enough
	}

	wthi += med.wtlo;
	if (targethi-wthi <= tolerance) break;  // close enough

	valuemax = med.valuelo;                   // iterate again
	markactive = 0;
      }

      else                  // Goldilocks result: both partitions just right
	break;

      // shrink the active list
      
      k = 0;
      for (j = 0; j < nlist; j++) {
	i = dotlist[j];
	if (dotmark[i] == markactive) dotlist[k++] = i;
      }
      nlist = k;
    }

    // found median
    // store cut info only if I am procmid

    if (me == procmid) {
      cut = valuehalf;
      cutdim = dim;
    }

    // use cut to shrink my RCB bounding box

    if (me < procmid) hi[dim] = valuehalf;
    else lo[dim] = valuehalf;

    // outgoing = number of dots to ship to partner
    // nkeep = number of dots that have never migrated

    markactive = (me < procpartner);
    for (i = 0, keep = 0, outgoing = 0; i < ndot; i++)
      if (dotmark[i] == markactive) outgoing++;
      else if (i < nkeep) keep++;
    nkeep = keep;
    
    // alert partner how many dots I'll send, read how many I'll recv

    MPI_Send(&outgoing,1,MPI_INT,procpartner,0,world);
    incoming = 0;
    if (readnumber) {
      MPI_Recv(&incoming,1,MPI_INT,procpartner,0,world,MPI_STATUS_IGNORE);
      if (readnumber == 2) {
	MPI_Recv(&incoming2,1,MPI_INT,procpartner2,0,world,MPI_STATUS_IGNORE);
	incoming += incoming2;
      }
    }

    // check if need to alloc more space

    int ndotnew = ndot - outgoing + incoming;
    if (ndotnew > maxdot) {
      while (maxdot < ndotnew) maxdot += DELTA;
      dots = (Dot *) memory->srealloc(dots,maxdot*sizeof(Dot),"RCB::dots");
      counters[6]++;
    }

    counters[1] += outgoing;
    counters[2] += incoming;
    if (ndotnew > counters[3]) counters[3] = ndotnew;
    if (maxdot > counters[4]) counters[4] = maxdot;
    
    // malloc comm send buffer

    if (outgoing > maxbuf) {
      memory->sfree(buf);
      maxbuf = outgoing;
      buf = (Dot *) memory->smalloc(maxbuf*sizeof(Dot),"RCB:buf");
    }

    // fill buffer with dots that are marked for sending
    // pack down the unmarked ones
    
    keep = outgoing = 0;
    for (i = 0; i < ndot; i++) {
      if (dotmark[i] == markactive)
	memcpy(&buf[outgoing++],&dots[i],sizeof(Dot));
      else
	memcpy(&dots[keep++],&dots[i],sizeof(Dot));
    }

    // post receives for dots

    if (readnumber > 0) {
      MPI_Irecv(&dots[keep],incoming*sizeof(Dot),MPI_CHAR,
                procpartner,1,world,&request);
      if (readnumber == 2) {
	keep += incoming - incoming2;
	MPI_Irecv(&dots[keep],incoming2*sizeof(Dot),MPI_CHAR,
                  procpartner2,1,world,&request2);
      }
    }
    
    // handshake before sending dots to insure recvs have been posted
    
    if (readnumber > 0) {
      MPI_Send(NULL,0,MPI_INT,procpartner,0,world);
      if (readnumber == 2) MPI_Send(NULL,0,MPI_INT,procpartner2,0,world);
    }
    MPI_Recv(NULL,0,MPI_INT,procpartner,0,world,MPI_STATUS_IGNORE);

    // send dots to partner

    MPI_Rsend(buf,outgoing*sizeof(Dot),MPI_CHAR,procpartner,1,world);
    
    // wait until all dots are received

    if (readnumber > 0) {
      MPI_Wait(&request,MPI_STATUS_IGNORE);
      if (readnumber == 2) MPI_Wait(&request2,MPI_STATUS_IGNORE);
    }

    ndot = ndotnew;

    // cut partition in half, create new communicators of 1/2 size

    int split;
    if (me < procmid) {
      procupper = procmid - 1;
      split = 0;
    } else {
      proclower = procmid;
      split = 1;
    }

    MPI_Comm_split(comm,split,me,&comm_half);
    MPI_Comm_free(&comm);
    comm = comm_half;
  }

  // clean up

  MPI_Comm_free(&comm);

  // set public variables with results of rebalance

  nfinal = ndot;

  if (nfinal > maxrecv) {
    memory->destroy(recvproc);
    memory->destroy(recvindex);
    maxrecv = nfinal;
    memory->create(recvproc,maxrecv,"RCB:recvproc");
    memory->create(recvindex,maxrecv,"RCB:recvindex");
  }

  for (i = 0; i < nfinal; i++) {
    recvproc[i] = dots[i].proc;
    recvindex[i] = dots[i].index;
  }
}

/* ----------------------------------------------------------------------
   custom MPI reduce operation
   merge of each component of an RCB bounding box
------------------------------------------------------------------------- */

void box_merge(void *in, void *inout, int *len, MPI_Datatype *dptr)

{
  RCB::BBox *box1 = (RCB::BBox *) in;
  RCB::BBox *box2 = (RCB::BBox *) inout;

  for (int i = 0; i < 3; i++) {
    if (box1->lo[i] < box2->lo[i]) box2->lo[i] = box1->lo[i];
    if (box1->hi[i] > box2->hi[i]) box2->hi[i] = box1->hi[i];
  }
}

/* ----------------------------------------------------------------------
   custom MPI reduce operation
   merge median data structure
   on input:
     in,inout->totallo, totalhi = weight in both partitions on this proc
               valuelo, valuehi = pos of nearest dot(s) to cut on this proc
               wtlo, wthi       = total wt of dot(s) at that pos on this proc
               countlo, counthi = # of dot(s) nearest to cut on this proc
               proclo, prochi   = not used
   on exit:
     inout->   totallo, totalhi = total # of active dots in both partitions
               valuelo, valuehi = pos of nearest dot(s) to cut
               wtlo, wthi       = total wt of dot(s) at that position
               countlo, counthi = total # of dot(s) nearest to cut
               proclo, prochi   = one unique proc who owns a nearest dot
	                          all procs must get same proclo,prochi
------------------------------------------------------------------------- */

void median_merge(void *in, void *inout, int *len, MPI_Datatype *dptr)

{
  RCB::Median *med1 = (RCB::Median *) in;
  RCB::Median *med2 = (RCB::Median *) inout;
  
  med2->totallo += med1->totallo;
  if (med1->valuelo > med2->valuelo) {
    med2->valuelo = med1->valuelo;
    med2->wtlo = med1->wtlo;
    med2->countlo = med1->countlo;
    med2->proclo = med1->proclo;
  }
  else if (med1->valuelo == med2->valuelo) {
    med2->wtlo += med1->wtlo;
    med2->countlo += med1->countlo;
    if (med1->proclo < med2->proclo) med2->proclo = med1->proclo;
  }

  med2->totalhi += med1->totalhi;
  if (med1->valuehi < med2->valuehi) {
    med2->valuehi = med1->valuehi;
    med2->wthi = med1->wthi;
    med2->counthi = med1->counthi;
    med2->prochi = med1->prochi;
  }
  else if (med1->valuehi == med2->valuehi) {
    med2->wthi += med1->wthi;
    med2->counthi += med1->counthi;
    if (med1->prochi < med2->prochi) med2->prochi = med1->prochi;
  }
}

/* ----------------------------------------------------------------------
   invert the RCB rebalance result to convert receive info into send info
   sortflag = flag for sorting order of received messages by proc ID
------------------------------------------------------------------------- */

void RCB::invert(int sortflag)
{

  // only create Irregular if not previously created
  // allows Irregular to persist for multiple RCB calls by fix balance

  if (!irregular) irregular = new Irregular(lmp);

  // nsend = # of dots to request from other procs

  int nsend = nfinal-nkeep;

  int *proclist;
  memory->create(proclist,nsend,"RCB:proclist");

  Invert *sinvert = 
    (Invert *) memory->smalloc(nsend*sizeof(Invert),"RCB:sinvert");

  int m = 0;
  for (int i = nkeep; i < nfinal; i++) {
    proclist[m] = recvproc[i];
    sinvert[m].rindex = recvindex[i];
    sinvert[m].sproc = me;
    sinvert[m].sindex = i;
    m++;
  }
  
  // perform inversion via irregular comm
  // nrecv = # of my dots to send to other procs

  int nrecv = irregular->create_data(nsend,proclist,sortflag);
  Invert *rinvert = 
    (Invert *) memory->smalloc(nrecv*sizeof(Invert),"RCB:rinvert");
  irregular->exchange_data((char *) sinvert,sizeof(Invert),(char *) rinvert);
  irregular->destroy_data();

  // set public variables from requests to send my dots

  if (noriginal > maxsend) {
    memory->destroy(sendproc);
    memory->destroy(sendindex);
    maxsend = noriginal;
    memory->create(sendproc,maxsend,"RCB:sendproc");
    memory->create(sendindex,maxsend,"RCB:sendindex");
  }

  for (int i = 0; i < nkeep; i++) {
    sendproc[recvindex[i]] = me;
    sendindex[recvindex[i]] = i;
  }

  for (int i = 0; i < nrecv; i++) {
    m = rinvert[i].rindex;
    sendproc[m] = rinvert[i].sproc;
    sendindex[m] = rinvert[i].sindex;
  }

  // clean-up

  memory->destroy(proclist);
  memory->destroy(sinvert);
  memory->destroy(rinvert);
}

/* ----------------------------------------------------------------------
   memory use of Irregular
------------------------------------------------------------------------- */

bigint RCB::memory_usage()
{
  bigint bytes = 0;
  if (irregular) bytes += irregular->memory_usage();
  return bytes;
}

// -----------------------------------------------------------------------
// DEBUG methods
// -----------------------------------------------------------------------
/*
// consistency checks on RCB results

void RCB::check()
{
  int i,iflag,total1,total2;
  double weight,wtmax,wtmin,wtone,tolerance;

  // check that total # of dots remained the same

  MPI_Allreduce(&ndotorig,&total1,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&ndot,&total2,1,MPI_INT,MPI_SUM,world);
  if (total1 != total2) {
    if (me == 0) 
      printf("ERROR: Points before RCB = %d, Points after RCB = %d\n",
	     total1,total2);
  }
  
  // check that result is load-balanced within log2(P)*max-wt

  weight = wtone = 0.0;
  for (i = 0; i < ndot; i++) {
    weight += dots[i].wt;
    if (dots[i].wt > wtone) wtone = dots[i].wt;
  }

  MPI_Allreduce(&weight,&wtmin,1,MPI_DOUBLE,MPI_MIN,world);
  MPI_Allreduce(&weight,&wtmax,1,MPI_DOUBLE,MPI_MAX,world);
  MPI_Allreduce(&wtone,&tolerance,1,MPI_DOUBLE,MPI_MAX,world);

  // i = smallest power-of-2 >= nprocs
  // tolerance = largest-single-weight*log2(nprocs)

  for (i = 0; (nprocs >> i) != 0; i++);
  tolerance = tolerance * i * (1.0 + TINY);

  if (wtmax - wtmin > tolerance) {
    if (me == 0) 
      printf("ERROR: Load-imbalance > tolerance of %g\n",tolerance);
    MPI_Barrier(world);
    if (weight == wtmin) printf("  Proc %d has weight = %g\n",me,weight);
    if (weight == wtmax) printf("  Proc %d has weight = %g\n",me,weight);
  }
  
  MPI_Barrier(world);
  
  // check that final set of points is inside RCB box of each proc
  
  iflag = 0;
  for (i = 0; i < ndot; i++) {
    if (dots[i].x[0] < lo[0] || dots[i].x[0] > hi[0] ||
	dots[i].x[1] < lo[1] || dots[i].x[1] > hi[1] ||
	dots[i].x[2] < lo[2] || dots[i].x[2] > hi[2])
      iflag++;
  }
  if (iflag > 0) 
    printf("ERROR: %d points are out-of-box on proc %d\n",iflag,me);
}

// stats for RCB decomposition

void RCB::stats(int flag)
{
  int i,iflag,sum,min,max;
  double ave,rsum,rmin,rmax;
  double weight,wttot,wtmin,wtmax;

  if (me == 0) printf("RCB Statistics:\n");

  // distribution info
  
  for (i = 0, weight = 0.0; i < ndot; i++) weight += dots[i].wt;
  MPI_Allreduce(&weight,&wttot,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&weight,&wtmin,1,MPI_DOUBLE,MPI_MIN,world);
  MPI_Allreduce(&weight,&wtmax,1,MPI_DOUBLE,MPI_MAX,world);

  if (me == 0) {
    printf(" Total weight of dots = %g\n",wttot);
    printf(" Weight on each proc: ave = %g, max = %g, min = %g\n",
	   wttot/nprocs,wtmax,wtmin);
  }
  if (flag) {
    MPI_Barrier(world);
    printf("    Proc %d has weight = %g\n",me,weight);
  }

  for (i = 0, weight = 0.0; i < ndot; i++) 
    if (dots[i].wt > weight) weight = dots[i].wt;
  MPI_Allreduce(&weight,&wtmax,1,MPI_DOUBLE,MPI_MAX,world);
  
  if (me == 0) printf(" Maximum weight of single dot = %g\n",wtmax);
  if (flag) {
    MPI_Barrier(world);
    printf("    Proc %d max weight = %g\n",me,weight);
  }

  // counter info

  MPI_Allreduce(&counters[0],&sum,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&counters[0],&min,1,MPI_INT,MPI_MIN,world);
  MPI_Allreduce(&counters[0],&max,1,MPI_INT,MPI_MAX,world);
  ave = ((double) sum)/nprocs;
  if (me == 0) 
    printf(" Median iter: ave = %g, min = %d, max = %d\n",ave,min,max);
  if (flag) {
    MPI_Barrier(world);
    printf("    Proc %d median count = %d\n",me,counters[0]);
  }

  MPI_Allreduce(&counters[1],&sum,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&counters[1],&min,1,MPI_INT,MPI_MIN,world);
  MPI_Allreduce(&counters[1],&max,1,MPI_INT,MPI_MAX,world);
  ave = ((double) sum)/nprocs;
  if (me == 0) 
    printf(" Send count: ave = %g, min = %d, max = %d\n",ave,min,max);
  if (flag) {
    MPI_Barrier(world);
    printf("    Proc %d send count = %d\n",me,counters[1]);
  }

  MPI_Allreduce(&counters[2],&sum,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&counters[2],&min,1,MPI_INT,MPI_MIN,world);
  MPI_Allreduce(&counters[2],&max,1,MPI_INT,MPI_MAX,world);
  ave = ((double) sum)/nprocs;
  if (me == 0) 
    printf(" Recv count: ave = %g, min = %d, max = %d\n",ave,min,max);
  if (flag) {
    MPI_Barrier(world);
    printf("    Proc %d recv count = %d\n",me,counters[2]);
  }

  MPI_Allreduce(&counters[3],&sum,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&counters[3],&min,1,MPI_INT,MPI_MIN,world);
  MPI_Allreduce(&counters[3],&max,1,MPI_INT,MPI_MAX,world);
  ave = ((double) sum)/nprocs;
  if (me == 0) 
    printf(" Max dots: ave = %g, min = %d, max = %d\n",ave,min,max);
  if (flag) {
    MPI_Barrier(world);
    printf("    Proc %d max dots = %d\n",me,counters[3]);
  }

  MPI_Allreduce(&counters[4],&sum,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&counters[4],&min,1,MPI_INT,MPI_MIN,world);
  MPI_Allreduce(&counters[4],&max,1,MPI_INT,MPI_MAX,world);
  ave = ((double) sum)/nprocs;
  if (me == 0) 
    printf(" Max memory: ave = %g, min = %d, max = %d\n",ave,min,max);
  if (flag) {
    MPI_Barrier(world);
    printf("    Proc %d max memory = %d\n",me,counters[4]);
  }

  if (reuse) {
    MPI_Allreduce(&counters[5],&sum,1,MPI_INT,MPI_SUM,world);
    MPI_Allreduce(&counters[5],&min,1,MPI_INT,MPI_MIN,world);
    MPI_Allreduce(&counters[5],&max,1,MPI_INT,MPI_MAX,world);
    ave = ((double) sum)/nprocs;
    if (me == 0) 
      printf(" # of Reuse: ave = %g, min = %d, max = %d\n",ave,min,max);
    if (flag) {
      MPI_Barrier(world);
      printf("    Proc %d # of Reuse = %d\n",me,counters[5]);
    }
  }
  
  MPI_Allreduce(&counters[6],&sum,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&counters[6],&min,1,MPI_INT,MPI_MIN,world);
  MPI_Allreduce(&counters[6],&max,1,MPI_INT,MPI_MAX,world);
  ave = ((double) sum)/nprocs;
  if (me == 0) 
    printf(" # of OverAlloc: ave = %g, min = %d, max = %d\n",ave,min,max);
  if (flag) {
    MPI_Barrier(world);
    printf("    Proc %d # of OverAlloc = %d\n",me,counters[6]);
  }

  // RCB boxes for each proc
  
  if (flag) {
    if (me == 0) printf(" RCB sub-domain boxes:\n");
    for (i = 0; i < 3; i++) {
      MPI_Barrier(world);
      if (me == 0) printf("    Dimension %d\n",i+1);
      MPI_Barrier(world);
      printf("      Proc = %d: Box = %g %g\n",me,lo[i],hi[i]);
    }
  }
}
*/
