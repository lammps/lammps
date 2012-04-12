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
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "balance.h"
#include "atom.h"
#include "comm.h"
#include "irregular.h"
#include "domain.h"
#include "force.h"
#include "update.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{NONE,UNIFORM,USER,DYNAMIC};
enum{X,Y,Z};
enum{EXPAND,CONTRACT};

//#define BALANCE_DEBUG 1

/* ---------------------------------------------------------------------- */

Balance::Balance(LAMMPS *lmp) : Pointers(lmp)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  memory->create(pcount,nprocs,"balance:pcount");
  memory->create(allcount,nprocs,"balance:allcount");

  user_xsplit = user_ysplit = user_zsplit = NULL;
  dflag = 0;

  fp = NULL;
}

/* ---------------------------------------------------------------------- */

Balance::~Balance()
{
  memory->destroy(pcount);
  memory->destroy(allcount);

  delete [] user_xsplit;
  delete [] user_ysplit;
  delete [] user_zsplit;

  if (dflag) {
    delete [] bstr;
    delete [] ops;
    delete [] counts[0];
    delete [] counts[1];
    delete [] counts[2];
    delete [] cuts;
    delete [] onecount;
    //MPI_Comm_free(&commslice[0]);
    //MPI_Comm_free(&commslice[1]);
    //MPI_Comm_free(&commslice[2]);
  }

  if (fp) fclose(fp);
}

/* ----------------------------------------------------------------------
   called as balance command in input script
------------------------------------------------------------------------- */

void Balance::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Balance command before simulation box is defined");

  if (comm->me == 0 && screen) fprintf(screen,"Balancing ...\n");

  // parse arguments

  if (narg < 1) error->all(FLERR,"Illegal balance command");

  int dimension = domain->dimension;
  int *procgrid = comm->procgrid;
  xflag = yflag = zflag = NONE;
  dflag = 0;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"x") == 0) {
      if (xflag == DYNAMIC) error->all(FLERR,"Illegal balance command");
      if (strcmp(arg[iarg+1],"uniform") == 0) {
	if (iarg+2 > narg) error->all(FLERR,"Illegal balance command");
	xflag = UNIFORM;
	iarg += 2;
      } else {
	if (1 + procgrid[0]-1 > narg)
	  error->all(FLERR,"Illegal balance command");
	xflag = USER;
	delete [] user_xsplit;
	user_xsplit = new double[procgrid[0]+1];
	user_xsplit[0] = 0.0;
	iarg++;
	for (int i = 1; i < procgrid[0]; i++)
	  user_xsplit[i] = force->numeric(arg[iarg++]);
	user_xsplit[procgrid[0]] = 1.0;
      }
    } else if (strcmp(arg[iarg],"y") == 0) {
      if (yflag == DYNAMIC) error->all(FLERR,"Illegal balance command");
      if (strcmp(arg[iarg+1],"uniform") == 0) {
	if (iarg+2 > narg) error->all(FLERR,"Illegal balance command");
	yflag = UNIFORM;
	iarg += 2;
      } else {
	if (1 + procgrid[1]-1 > narg)
	  error->all(FLERR,"Illegal balance command");
	yflag = USER;
	delete [] user_ysplit;
	user_ysplit = new double[procgrid[1]+1];
	user_ysplit[0] = 0.0;
	iarg++;
	for (int i = 1; i < procgrid[1]; i++)
	  user_ysplit[i] = force->numeric(arg[iarg++]);
	user_ysplit[procgrid[1]] = 1.0;
      }
    } else if (strcmp(arg[iarg],"z") == 0) {
      if (zflag == DYNAMIC) error->all(FLERR,"Illegal balance command");
      if (strcmp(arg[iarg+1],"uniform") == 0) {
	if (iarg+2 > narg) error->all(FLERR,"Illegal balance command");
	zflag = UNIFORM;
	iarg += 2;
      } else {
	if (1 + procgrid[2]-1 > narg)
	  error->all(FLERR,"Illegal balance command");
	zflag = USER;
	delete [] user_zsplit;
	user_zsplit = new double[procgrid[2]+1];
	user_zsplit[0] = 0.0;
	iarg++;
	for (int i = 1; i < procgrid[2]; i++)
	  user_zsplit[i] = force->numeric(arg[iarg++]);
	user_zsplit[procgrid[2]] = 1.0;
      }

    } else if (strcmp(arg[iarg],"dynamic") == 0) {
      if (xflag != NONE || yflag != NONE || zflag != NONE)
	error->all(FLERR,"Illegal balance command");
      if (iarg+5 > narg) error->all(FLERR,"Illegal balance command");
      dflag = 1;
      xflag = yflag = DYNAMIC;
      if (dimension == 3) zflag = DYNAMIC;
      nrepeat = atoi(arg[iarg+1]);
      niter = atoi(arg[iarg+2]);
      if (nrepeat <= 0 || niter <= 0)
	error->all(FLERR,"Illegal balance command");
      int n = strlen(arg[iarg+3]) + 1;
      bstr = new char[n];
      strcpy(bstr,arg[iarg+3]);
      thresh = atof(arg[iarg+4]);
      if (thresh < 1.0) error->all(FLERR,"Illegal balance command");
      iarg += 5;

    } else error->all(FLERR,"Illegal balance command");
  }

  // error check

  if (zflag && dimension == 2)
    error->all(FLERR,"Cannot balance in z dimension for 2d simulation");

  if (xflag == USER)
    for (int i = 1; i <= procgrid[0]; i++)
      if (user_xsplit[i-1] >= user_xsplit[i]) 
	error->all(FLERR,"Illegal balance command");
  if (yflag == USER)
    for (int i = 1; i <= procgrid[1]; i++)
      if (user_ysplit[i-1] >= user_ysplit[i]) 
	error->all(FLERR,"Illegal balance command");
  if (zflag == USER)
    for (int i = 1; i <= procgrid[2]; i++)
      if (user_zsplit[i-1] >= user_zsplit[i]) 
	error->all(FLERR,"Illegal balance command");

  if (dflag) {
    for (int i = 0; i < strlen(bstr); i++) {
      if (bstr[i] != 'x' && bstr[i] != 'y' && bstr[i] != 'z') 
	error->all(FLERR,"Balance dynamic string is invalid");
      if (bstr[i] == 'z' && dimension == 2) 
	error->all(FLERR,"Balance dynamic string is invalid for 2d simulation");
    }
  }

  // insure atoms are in current box & update box via shrink-wrap
  // no exchange() since doesn't matter if atoms are assigned to correct procs

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  domain->reset_box();
  if (domain->triclinic) domain->lamda2x(atom->nlocal);

  // imbinit = initial imbalance
  // use current splits instead of nlocal since atoms may not be in sub-box

  domain->x2lamda(atom->nlocal);
  int maxinit;
  double imbinit = imbalance_splits(maxinit);
  domain->lamda2x(atom->nlocal);

  // debug output of initial state

#ifdef BALANCE_DEBUG
  dumpout(update->ntimestep);
#endif

  // explicit setting of sub-domain sizes

  if (xflag == UNIFORM) {
    for (int i = 0; i < procgrid[0]; i++) 
      comm->xsplit[i] = i * 1.0/procgrid[0];
    comm->xsplit[procgrid[0]] = 1.0;
  }

  if (yflag == UNIFORM) {
    for (int i = 0; i < procgrid[1]; i++) 
      comm->ysplit[i] = i * 1.0/procgrid[1];
    comm->ysplit[procgrid[1]] = 1.0;
  }

  if (zflag == UNIFORM) {
    for (int i = 0; i < procgrid[2]; i++) 
      comm->zsplit[i] = i * 1.0/procgrid[2];
    comm->zsplit[procgrid[2]] = 1.0;
  }

  if (xflag == USER)
    for (int i = 0; i <= procgrid[0]; i++) comm->xsplit[i] = user_xsplit[i];

  if (yflag == USER)
    for (int i = 0; i <= procgrid[1]; i++) comm->ysplit[i] = user_ysplit[i];

  if (zflag == USER)
    for (int i = 0; i <= procgrid[2]; i++) comm->zsplit[i] = user_zsplit[i];

  // static load-balance of sub-domain sizes

  int count = 0;
  if (dflag) {
    dynamic_setup(bstr);
    count = dynamic_once();
  }

  // debug output of final result

#ifdef BALANCE_DEBUG
  dumpout(-1);
#endif

  // reset comm->uniform flag if necessary

  if (comm->uniform) {
    if (xflag == USER || xflag == DYNAMIC) comm->uniform = 0;
    if (yflag == USER || yflag == DYNAMIC) comm->uniform = 0;
    if (zflag == USER || zflag == DYNAMIC) comm->uniform = 0;
  } else {
    if (dimension == 3) {
      if (xflag == UNIFORM && yflag == UNIFORM && zflag == UNIFORM)
	comm->uniform = 1;
    } else {
      if (xflag == UNIFORM && yflag == UNIFORM) comm->uniform = 1;
    }
  }

  // reset proc sub-domains

  if (domain->triclinic) domain->set_lamda_box();
  domain->set_local_box();

  // move atoms to new processors via irregular()

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  Irregular *irregular = new Irregular(lmp);
  irregular->migrate_atoms();
  delete irregular;
  if (domain->triclinic) domain->lamda2x(atom->nlocal);

  // check if any atoms were lost

  bigint natoms;
  bigint nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal,&natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);
  if (natoms != atom->natoms) {
    char str[128];
    sprintf(str,"Lost atoms via balance: original " BIGINT_FORMAT 
	    " current " BIGINT_FORMAT,atom->natoms,natoms);
    error->all(FLERR,str);
  }

  // imbfinal = final imbalance based on final nlocal

  int maxfinal;
  double imbfinal = imbalance_nlocal(maxfinal);

  if (me == 0) {
    if (screen) {
      fprintf(screen,"  iteration count = %d\n",count);
      fprintf(screen,"  initial/final max atoms/proc = %d %d\n",
	      maxinit,maxfinal);
      fprintf(screen,"  initial/final imbalance factor = %g %g\n",
	      imbinit,imbfinal);
    }
    if (logfile) {
      fprintf(logfile,"  iteration count = %d\n",count);
      fprintf(logfile,"  initial/final max atoms/proc = %d %d\n",
	      maxinit,maxfinal);
      fprintf(logfile,"  initial/final imbalance factor = %g %g\n",
	      imbinit,imbfinal);
    }
  }

  if (me == 0) {
    if (screen) {
      fprintf(screen,"  x cuts:");
      for (int i = 0; i <= comm->procgrid[0]; i++)
	fprintf(screen," %g",comm->xsplit[i]);
      fprintf(screen,"\n");
      fprintf(screen,"  y cuts:");
      for (int i = 0; i <= comm->procgrid[1]; i++)
	fprintf(screen," %g",comm->ysplit[i]);
      fprintf(screen,"\n");
      fprintf(screen,"  z cuts:");
      for (int i = 0; i <= comm->procgrid[2]; i++)
	fprintf(screen," %g",comm->zsplit[i]);
      fprintf(screen,"\n");
    }
    if (logfile) {
      fprintf(logfile,"  x cuts:");
      for (int i = 0; i <= comm->procgrid[0]; i++)
	fprintf(logfile," %g",comm->xsplit[i]);
      fprintf(logfile,"\n");
      fprintf(logfile,"  y cuts:");
      for (int i = 0; i <= comm->procgrid[1]; i++)
	fprintf(logfile," %g",comm->ysplit[i]);
      fprintf(logfile,"\n");
      fprintf(logfile,"  z cuts:");
      for (int i = 0; i <= comm->procgrid[2]; i++)
	fprintf(logfile," %g",comm->zsplit[i]);
      fprintf(logfile,"\n");
    }
  }
}

/* ----------------------------------------------------------------------
   calculate imbalance based on nlocal
   return max = max atom per proc
   return imbalance factor = max atom per proc / ave atom per proc
------------------------------------------------------------------------- */

double Balance::imbalance_nlocal(int &max)
{
  MPI_Allreduce(&atom->nlocal,&max,1,MPI_INT,MPI_MAX,world);
  double imbalance = max / (1.0 * atom->natoms / nprocs);
  return imbalance;
}

/* ----------------------------------------------------------------------
   calculate imbalance based on processor splits in 3 dims
   atoms must be in lamda coords (0-1) before called
   map atoms to 3d grid of procs
   return max = max atom per proc
   return imbalance factor = max atom per proc / ave atom per proc
------------------------------------------------------------------------- */

double Balance::imbalance_splits(int &max)
{
  double *xsplit = comm->xsplit;
  double *ysplit = comm->ysplit;
  double *zsplit = comm->zsplit;

  int nx = comm->procgrid[0];
  int ny = comm->procgrid[1];
  int nz = comm->procgrid[2];

  for (int i = 0; i < nprocs; i++) pcount[i] = 0;

  double **x = atom->x;
  int nlocal = atom->nlocal;
  int ix,iy,iz;

  for (int i = 0; i < nlocal; i++) {
    ix = binary(x[i][0],nx,xsplit);
    iy = binary(x[i][1],ny,ysplit);
    iz = binary(x[i][2],nz,zsplit);
    pcount[iz*nx*ny + iy*nx + ix]++;
  }

  MPI_Allreduce(pcount,allcount,nprocs,MPI_INT,MPI_SUM,world);
  max = 0;
  for (int i = 0; i < nprocs; i++) max = MAX(max,allcount[i]);
  double imbalance = max / (1.0 * atom->natoms / nprocs);
  return imbalance;
}

/* ----------------------------------------------------------------------
   setup static load balance operations
   called from command
------------------------------------------------------------------------- */

void Balance::dynamic_setup(char *str)
{
  nops = strlen(str);
  ops = new int[nops];

  for (int i = 0; i < strlen(str); i++) {
    if (str[i] == 'x') ops[i] = X;
    if (str[i] == 'y') ops[i] = Y;
    if (str[i] == 'z') ops[i] = Z;
  }
  
  splits[0] = comm->xsplit;
  splits[1] = comm->ysplit;
  splits[2] = comm->zsplit;

  counts[0] = new bigint[comm->procgrid[0]];
  counts[1] = new bigint[comm->procgrid[1]];
  counts[2] = new bigint[comm->procgrid[2]];

  int max = MAX(comm->procgrid[0],comm->procgrid[1]);
  max = MAX(max,comm->procgrid[2]);
  cuts = new double[max+1];
  onecount = new bigint[max];

  //MPI_Comm_split(world,comm->myloc[0],0,&commslice[0]);
  //MPI_Comm_split(world,comm->myloc[1],0,&commslice[1]);
  //MPI_Comm_split(world,comm->myloc[2],0,&commslice[2]);
}

/* ----------------------------------------------------------------------
   setup dynamic load balance operations
   called from fix balance
------------------------------------------------------------------------- */

void Balance::dynamic_setup(int nrepeat_in, int niter_in,
			    char *str, double thresh_in)
{
  nrepeat = nrepeat_in;
  niter = niter_in;
  thresh = thresh_in;

  dynamic_setup(str);
}

/* ----------------------------------------------------------------------
   perform static load balance by changing xyz split proc boundaries in Comm
   called from command
   return actual iteration count
------------------------------------------------------------------------- */

int Balance::dynamic_once()
{
  int i,m,max;
  double imbfactor;

  int *procgrid = comm->procgrid;

  domain->x2lamda(atom->nlocal);

  int count = 0;
  for (int irepeat = 0; irepeat < nrepeat; irepeat++) {
    for (i = 0; i < nops; i++) {
      for (m = 0; m < niter; m++) {
	stats(ops[i],procgrid[ops[i]],splits[ops[i]],counts[ops[i]]);
	adjust(procgrid[ops[i]],counts[ops[i]],splits[ops[i]]);
	count++;

#ifdef BALANCE_DEBUG
	dumpout(-1);
#endif
      }
      imbfactor = imbalance_splits(max);
      if (comm->me == 0) printf("AAA %d %d %g\n",irepeat,i,imbfactor);
      if (imbfactor <= thresh) break;
    }
    if (i < nops) break;
  }

  domain->lamda2x(atom->nlocal);

  return count;
}

/* ----------------------------------------------------------------------
   perform dynamic load balance by changing xyz split proc boundaries in Comm
   called from fix balance
   return actual iteration count
------------------------------------------------------------------------- */

int Balance::dynamic()
{
  int i,m,max;
  double imbfactor;

#ifdef BALANCE_DEBUG
  dumpout(update->ntimestep);
#endif

  int *procgrid = comm->procgrid;

  domain->x2lamda(atom->nlocal);

  int count = 0;
  for (int irepeat = 0; irepeat < nrepeat; irepeat++) {
    for (i = 0; i < nops; i++) {
      for (m = 0; m < niter; m++) {
	stats(ops[i],procgrid[ops[i]],splits[ops[i]],counts[ops[i]]);
	adjust(procgrid[ops[i]],counts[ops[i]],splits[ops[i]]);
	count++;

#ifdef BALANCE_DEBUG
	dumpout(-1);
#endif
      }
      imbfactor = imbalance_splits(max);
      if (imbfactor <= thresh) break;
    }
    if (i < nops) break;
  }

  domain->lamda2x(atom->nlocal);

#ifdef BALANCE_DEBUG
  dumpout(-1);
#endif

  return count;
}

/* ----------------------------------------------------------------------
   count atoms in each slice
   current cuts may be very different than original cuts,
   so use binary search to find which slice each atom is in
------------------------------------------------------------------------- */

void Balance::stats(int dim, int n, double *split, bigint *count)
{
  for (int i = 0; i < n; i++) onecount[i] = 0;

  double **x = atom->x;
  int nlocal = atom->nlocal;
  int index;

  for (int i = 0; i < nlocal; i++) {
    index = binary(x[i][dim],n,split);
    onecount[index]++;
  }

  MPI_Allreduce(onecount,count,n,MPI_LMP_BIGINT,MPI_SUM,world);
}

/* ----------------------------------------------------------------------
   adjust cuts between N slices in a dim via diffusive method
   count = atoms per slice
   split = current N+1 cuts, with 0.0 and 1.0 at end points
   overwrite split with new cuts
   diffusion means slices with more atoms than their neighbors "send" atoms,
     by moving cut closer to sender, further from receiver
------------------------------------------------------------------------- */

void Balance::adjust(int n, bigint *count, double *split)
{
  // damping factor

  double damp = 0.5;

  // loop over slices
  // cut I is between 2 slices (I-1 and I) with counts
  // cut I+1 is between 2 slices (I and I+1) with counts
  // for a cut between 2 slices, only slice with larger count adjusts it
  // special treatment of end slices with only 1 neighbor

  bigint leftcount,mycount,rightcount;
  double rho,target,targetleft,targetright;

  for (int i = 0; i < n; i++) {
    if (i == 0) leftcount = MAXBIGINT;
    else leftcount = count[i-1];
    mycount = count[i];
    if (i == n-1) rightcount = MAXBIGINT;
    else rightcount = count[i+1];

    // middle slice is <= both left and right, so do nothing
    // special case if 2 slices both have count = 0 -> no change in cut

    if (mycount <= leftcount && mycount <= rightcount) {
      if (leftcount == 0) cuts[i] = split[i];
      if (rightcount == 0) cuts[i+1] = split[i+1];
      continue;
    }

    // rho = density of atoms in the slice

    rho = mycount / (split[i+1] - split[i]);

    // middle slice has more atoms than left or right slice
    // send atoms in that dir

    if (mycount > leftcount) {
      target = damp * 0.5*(mycount-leftcount);
      cuts[i] = split[i] + target/rho;
    }
    if (mycount > rightcount) {
      target = damp * 0.5*(mycount-rightcount);
      cuts[i+1] = split[i+1] - target/rho;
    }

    /*
    // middle slice has more atoms then left or right slice
    // if delta from middle to top slice > delta between top and bottom slice
    //   then send atoms both dirs to bring all 3 slices to same count
    // else bottom slice is very low, so send atoms only in that dir

    if (mycount > leftcount && mycount > rightcount) {
      if (mycount-MAX(leftcount,rightcount) >= fabs(leftcount-rightcount)) {
	if (leftcount <= rightcount) {
	  targetleft = damp * 
	    (rightcount-leftcount + (mycount-rightcount)/3.0);
	  targetright = damp * (mycount-rightcount)/3.0;
	  cuts[i] = split[i] + targetleft/rho;
	  cuts[i+1] = split[i+1] - targetright/rho;
	} else {
	  targetleft = damp * (mycount-leftcount)/3.0;
	  targetright = damp * 
	    (leftcount-rightcount + (mycount-leftcount)/3.0);
	  cuts[i] = split[i] + targetleft/rho;
	  cuts[i+1] = split[i+1] - targetright/rho;
	}
      } else if (leftcount < rightcount) {
	target = damp * 0.5*(mycount-leftcount);
	cuts[i] = split[i] + target/rho;
	cuts[i+1] = split[i+1];
      } else if (rightcount < leftcount) {
	target = damp * 0.5*(mycount-rightcount);
	cuts[i+1] = split[i+1] - target/rho;
	cuts[i] = split[i];
      }

    // middle slice has more atoms than only left or right slice
    // send atoms only in that dir

    } else if (mycount > leftcount) {
      target = damp * 0.5*(mycount-leftcount);
      cuts[i] = split[i] + target/rho;
    } else if (mycount > rightcount) {
      target = damp * 0.5*(mycount-rightcount);
      cuts[i+1] = split[i+1] - target/rho;
    }
    */
  }

  // overwrite adjustable splits with new cuts

  for (int i = 1; i < n; i++) split[i] = cuts[i];
}

/* ----------------------------------------------------------------------
   binary search for value in N-length ascending vec
   value may be outside range of vec limits
   always return index from 0 to N-1 inclusive
   return 0 if value < vec[0]
   reutrn N-1 if value >= vec[N-1]
   return index = 1 to N-2 if vec[index] <= value < vec[index+1]
------------------------------------------------------------------------- */

int Balance::binary(double value, int n, double *vec)
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
   create dump file of line segments in Pizza.py mdump mesh format
   just for debug/viz purposes
   write to tmp.mdump.*, open first time if necessary
   write xy lines around each proc's sub-domain for 2d
   write xyz cubes around each proc's sub-domain for 3d
   all procs but 0 just return
------------------------------------------------------------------------- */

void Balance::dumpout(bigint tstamp)
{
  if (me) return;

  bigint tstep;
  if (tstamp >= 0) tstep = tstamp;
  else tstep = laststep + 1;
  laststep = tstep;

  if (fp == NULL) {
    char file[32];
    sprintf(file,"tmp.mdump.%ld",tstep);
    fp = fopen(file,"w");
    if (!fp) error->one(FLERR,"Cannot open balance output file");

    // write out one square/cute per processor for 2d/3d
    // only write once since topology is static

    fprintf(fp,"ITEM: TIMESTEP\n");
    fprintf(fp,"%ld\n",tstep);
    fprintf(fp,"ITEM: NUMBER OF SQUARES\n");
    fprintf(fp,"%d\n",nprocs);
    fprintf(fp,"ITEM: SQUARES\n");

    int nx = comm->procgrid[0] + 1;
    int ny = comm->procgrid[1] + 1;
    int nz = comm->procgrid[2] + 1;

    if (domain->dimension == 2) {
      int m = 0;
      for (int j = 0; j < comm->procgrid[1]; j++)
	for (int i = 0; i < comm->procgrid[0]; i++) {
	  int c1 = j*nx + i + 1;
	  int c2 = c1 + 1;
	  int c3 = c2 + nx;
	  int c4 = c3 - 1;
	  fprintf(fp,"%d %d %d %d %d %d\n",m+1,m+1,c1,c2,c3,c4);
	  m++;
	}

    } else {
      int m = 0;
      for (int k = 0; k < comm->procgrid[2]; k++)
	for (int j = 0; j < comm->procgrid[1]; j++)
	  for (int i = 0; i < comm->procgrid[0]; i++) {
	    int c1 = k*ny*nx + j*nx + i + 1;
	    int c2 = c1 + 1;
	    int c3 = c2 + nx;
	    int c4 = c3 - 1;
	    int c5 = c1 + ny*nx;
	    int c6 = c2 + ny*nx;
	    int c7 = c3 + ny*nx;
	    int c8 = c4 + ny*nx;
	    fprintf(fp,"%d %d %d %d %d %d %d %d %d %d\n",
		    m+1,m+1,c1,c2,c3,c4,c5,c6,c7,c8);
	    m++;
	  }
    }
  }

  // write out nodal coords, can be different every call
  // scale xsplit,ysplit,zsplit values to full box
  // only implmented for orthogonal boxes, not triclinic

  int nx = comm->procgrid[0] + 1;
  int ny = comm->procgrid[1] + 1;
  int nz = comm->procgrid[2] + 1;

  double *boxlo = domain->boxlo;
  double *boxhi = domain->boxhi;
  double *prd = domain->prd;

  fprintf(fp,"ITEM: TIMESTEP\n");
  fprintf(fp,"%ld\n",tstep);
  fprintf(fp,"ITEM: NUMBER OF NODES\n");
  if (domain->dimension == 2) 
    fprintf(fp,"%d\n",nx*ny);
  else 
    fprintf(fp,"%d\n",nx*ny*nz);
  fprintf(fp,"ITEM: BOX BOUNDS\n");
  fprintf(fp,"%g %g\n",boxlo[0],boxhi[0]);
  fprintf(fp,"%g %g\n",boxlo[1],boxhi[1]);
  fprintf(fp,"%g %g\n",boxlo[2],boxhi[2]);
  fprintf(fp,"ITEM: NODES\n");

  if (domain->dimension == 2) {
    int m = 0;
    for (int j = 0; j < ny; j++)
      for (int i = 0; i < nx; i++) {
	fprintf(fp,"%d %d %g %g %g\n",m+1,1,
		boxlo[0] + prd[0]*comm->xsplit[i],
		boxlo[1] + prd[1]*comm->ysplit[j],
		0.0);
	m++;
      }
  } else {
    int m = 0;
    for (int k = 0; k < nz; k++)
      for (int j = 0; j < ny; j++)
	for (int i = 0; i < nx; i++) {
	  fprintf(fp,"%d %d %g %g %g\n",m+1,1,
		  boxlo[0] + prd[0]*comm->xsplit[i],
		  boxlo[1] + prd[1]*comm->ysplit[j],
		  boxlo[2] + prd[2]*comm->zsplit[j]);
	  m++;
      }
  }
}
