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
   Contributing author: Tzu-Ray Shan (U Florida, rayshan@ufl.edu)
------------------------------------------------------------------------- */

#include "mpi.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_qeq.h"
#include "atom.h"
#include "force.h"
#include "group.h"
#include "respa.h"
#include "pair_comb.h"
#include "update.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

/* ---------------------------------------------------------------------- */

FixQEQ::FixQEQ(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  if (narg < 5) error->all("Illegal fix qeq command");

  peratom_flag = 1;
  size_peratom_cols = 0;
  peratom_freq = 1;

  nevery = force->inumeric(arg[3]);
  precision = force->numeric(arg[4]);

  if (nevery <= 0 || precision <= 0.0) error->all("Illegal fix qeq command");

  MPI_Comm_rank(world,&me);

  // optional args

  fp = NULL;

  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"file") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix qeq command");
      if (me == 0) {
	fp = fopen(arg[iarg+1],"w");
	if (fp == NULL) {
	  char str[128];
	  sprintf(str,"Cannot open fix qeq file %s",arg[iarg+1]);
	  error->one(str);
	}
      }
      iarg += 2;
    } else error->all("Illegal fix qeq command");
  }
  
  nmax = atom->nmax;
  qf = (double *) memory->smalloc(nmax*sizeof(double),"qeq:qf");
  q1 = (double *) memory->smalloc(nmax*sizeof(double),"qeq:q1");
  q2 = (double *) memory->smalloc(nmax*sizeof(double),"qeq:q2");
  vector_atom = qf;

  // zero the vector since dump may access it on timestep 0
  // zero the vector since a variable may access it before first run

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) qf[i] = 0.0;
}

/* ---------------------------------------------------------------------- */

FixQEQ::~FixQEQ()
{
  if (me == 0 && fp) fclose(fp);
  memory->sfree(qf);
  memory->sfree(q1);
  memory->sfree(q2);
}

/* ---------------------------------------------------------------------- */

int FixQEQ::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixQEQ::init()
{
  if (atom->q == NULL) error->all("Must use atom style charge with fix qeq");

  comb = (PairComb *) force->pair_match("comb",1);
  if (comb == NULL) error->all("Fix qeq must be used with pair_style comb");

  if (strcmp(update->integrate_style,"respa") == 0)
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

  ngroup = group->count(igroup);
  if (ngroup == 0.0) error->all("Fix qeq group has no atoms");
}

/* ---------------------------------------------------------------------- */

void FixQEQ::setup(int vflag)
{
  firstflag = 1;
  if (strcmp(update->integrate_style,"verlet") == 0)
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
  firstflag = 0;
}

/* ---------------------------------------------------------------------- */

void FixQEQ::post_force(int vflag)
{
  int i,iloop,loopmax;
  double heatpq,qmass,dtq,dtq2;
  double enegchkall,enegmaxall;

  if (update->ntimestep % nevery) return;

  // reallocate work arrays if necessary
  // qf = charge force
  // q1 = charge displacement
  // q2 = tmp storage of charge force for next iteration

  if (atom->nmax > nmax) {
    memory->sfree(qf);
    memory->sfree(q1);
    memory->sfree(q2);
    nmax = atom->nmax;
    qf = (double *) memory->smalloc(nmax*sizeof(double),"qeq:qf");
    q1 = (double *) memory->smalloc(nmax*sizeof(double),"qeq:q1");
    q2 = (double *) memory->smalloc(nmax*sizeof(double),"qeq:q2");
    vector_atom = qf;
  }
  
  // more loops for first-time charge equilibrium

  iloop = 0; 
  if (firstflag) loopmax = 1000;
  else loopmax = 500;

  // charge-equilibration loop
  
  heatpq = 0.01;
  qmass = 0.06;
  dtq = 0.040;
  dtq2 = 0.5*dtq*dtq/qmass;

  double enegchk = 0.0;
  double enegtot = 0.0; 
  double enegmax = 0.0;

  double *q = atom->q;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++)
    q1[i] = q2[i] = qf[i] = 0.0;

  for (iloop = 0; iloop < loopmax; iloop ++ ) {
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	q1[i] = qf[i]*dtq2 - heatpq*q1[i];
	q[i] += q1[i]; 
      }

    enegtot = comb->yasu_char(qf,igroup);
    enegtot /= ngroup;
    enegchk = enegmax = 0.0;

    for (i = 0; i < nlocal ; i++)
      if (mask[i] & groupbit) {
	q2[i] = enegtot-qf[i];
	enegmax = MAX(enegmax,fabs(q2[i]));
	enegchk += fabs(q2[i]);
	qf[i] = q2[i];
      }

    MPI_Allreduce(&enegchk,&enegchkall,1,MPI_DOUBLE,MPI_SUM,world);
    enegchk = enegchkall/ngroup;
    MPI_Allreduce(&enegmax,&enegmaxall,1,MPI_DOUBLE,MPI_MAX,world);
    enegmax = enegmaxall;
    
    if (enegchk <= precision && enegmax <= 10.0*precision) break;

    if (me == 0 && fp)
      fprintf(fp,"Step: %d, loop: %d, enegtot %.6g, "
	      "enegmax %.6g, fq deviation: %.6g\n",
	      update->ntimestep,iloop,enegtot,enegmax,enegchk); 
    
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
	q1[i] += qf[i]*dtq2 - heatpq*q1[i]; 
  } 

  if (me == 0 && fp) {
    if (iloop == loopmax)
      fprintf(fp,"Charges did not converge in %d iterations\n",iloop);
    else
      fprintf(fp,"Charges converged in %d iterations to %.10f\n",
	      iloop,enegchk);
  }
}

/* ---------------------------------------------------------------------- */

void FixQEQ::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixQEQ::memory_usage()
{
  double bytes = atom->nmax*3 * sizeof(double);
  return bytes;
}
