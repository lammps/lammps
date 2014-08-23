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
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "pair_list.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "force.h"
#include "memory.h"

#include "error.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

using namespace LAMMPS_NS;

static const char * const stylename[] = {
  "none", "harmonic", "morse", "lj126", NULL
};

// fast power function for integer exponent > 0
static double mypow(double x, int n) {
  double yy;

  if (x == 0.0) return 0.0;

  for (yy = 1.0; n != 0; n >>= 1, x *=x)
    if (n & 1) yy *= x;

  return yy;
}

typedef struct { double x,y,z; } dbl3_t;

/* ---------------------------------------------------------------------- */

PairList::PairList(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  respa_enable = 0;
  cut_global = 0.0;
  style = NULL;
  params = NULL;
  check_flag = 1;
}

/* ---------------------------------------------------------------------- */

PairList::~PairList()
{
  memory->destroy(setflag);
  memory->destroy(cutsq);
  memory->destroy(style);
  memory->destroy(params);
}

/* ----------------------------------------------------------------------
   in this pair style we don't use a neighbor list, but loop through
   a list of pairwise interactions, determines the corresponding local
   atom indices and compute those forces.
------------------------------------------------------------------------- */

void PairList::compute(int eflag, int vflag)
{
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = eflag_global = 
    vflag_global = eflag_atom = vflag_atom = 0;

  const int nlocal = atom->nlocal;
  const int newton_pair = force->newton_pair;
  const dbl3_t * _noalias const x = (dbl3_t *) atom->x[0];
  dbl3_t * _noalias const f = (dbl3_t *) atom->f[0];

  double fpair,epair;
  int i,j;

  int pc = 0;
  for (int n=0; n < npairs; ++n) {
    const list_parm_t &par = params[n];
    i = atom->map(par.id1);
    j = atom->map(par.id2);

    // if one of the two atoms is missing on the node skip
    if ((i < 0) || (j < 0)) continue;

    // both atoms are ghosts -> skip
    if ((i >= nlocal) && (j >= nlocal)) continue;

    // with newton pair and one ghost we have to skip half the cases.
    // if id1 is a ghost, we skip if the sum of both ids is even.
    // if id2 is a ghost, we skip if the sum of both ids is odd.
    if (newton_pair) {
      if ((i >= nlocal) && ((par.id1+par.id2) & 1) == 0) continue;
      if ((j >= nlocal) && ((par.id1+par.id2) & 1) == 1) continue;
    }

    const double dx = x[i].x - x[j].x;
    const double dy = x[i].y - x[j].y;
    const double dz = x[i].z - x[j].z;
    const double rsq = dx*dx + dy*dy + dz*dz;

    fpair = epair = 0.0;
    if (check_flag) {
      if (newton_pair || i < nlocal) ++pc;
      if (newton_pair || j < nlocal) ++pc;
    }

    if (rsq < par.cutsq) {
      const double r2inv = 1.0/rsq;

      if (style[n] == HARM) {
	const double r = sqrt(rsq);
	const double dr = par.parm.harm.r0 - r;
	fpair = 2.0*par.parm.harm.k*dr/r;

	if (eflag_either)
	  epair = par.parm.harm.k*dr*dr - par.offset;

      } else if (style[n] == MORSE) {

	const double r = sqrt(rsq);
	const double dr = par.parm.morse.r0 - r;
	const double dexp = exp(par.parm.morse.alpha * dr);
	fpair = 2.0*par.parm.morse.d0*par.parm.morse.alpha 
                * (dexp*dexp - dexp) / r;

	if (eflag_either)
	  epair = par.parm.morse.d0 * (dexp*dexp - 2.0*dexp) - par.offset;

      } else if (style[n] == LJ126) {

	const double r6inv = r2inv*r2inv*r2inv;
	const double sig6  = mypow(par.parm.lj126.sigma,6);
	fpair =  24.0*par.parm.lj126.epsilon*r6inv 
                 * (2.0*sig6*sig6*r6inv - sig6) * r2inv;

	if (eflag_either)
	  epair = 4.0*par.parm.lj126.epsilon*r6inv 
                  * (sig6*sig6*r6inv - sig6) - par.offset;
      }

      if (newton_pair || i < nlocal) {
	f[i].x += dx*fpair;
	f[i].y += dy*fpair;
	f[i].z += dz*fpair;
      }

      if (newton_pair || j < nlocal) {
	f[j].x -= dx*fpair;
	f[j].y -= dy*fpair;
	f[j].z -= dz*fpair;
      }

      if (evflag) ev_tally(i,j,nlocal,newton_pair,epair,0.0,fpair,dx,dy,dz);
    } 
  }
  if (vflag_fdotr) virial_fdotr_compute();

  if (check_flag) {
     int tmp;
     MPI_Allreduce(&pc,&tmp,1,MPI_INT,MPI_SUM,world);
     if (tmp != 2*npairs)
       error->all(FLERR,"Not all pairs processed in pair_style list");
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairList::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");
}

/* ----------------------------------------------------------------------
   create one pair style for each arg in list
------------------------------------------------------------------------- */

void PairList::settings(int narg, char **arg)
{
  if (narg < 2)
    error->all(FLERR,"Illegal pair_style command");

  cut_global = force->numeric(FLERR,arg[1]);
  if (narg > 2) {
    if (strcmp(arg[2],"nocheck") == 0) check_flag = 0;
    if (strcmp(arg[2],"check") == 0) check_flag = 1;
  }

  FILE *fp = force->open_potential(arg[0]);
  char line[1024];
  if (fp == NULL)
    error->all(FLERR,"Cannot open pair list file");

  // count lines in file for upper limit of storage needed
  int num = 1;
  while(fgets(line,1024,fp)) ++num;
  rewind(fp);
  memory->create(style,num,"pair_list:style");
  memory->create(params,num,"pair_list:params");

  // now read and parse pair list file for real
  npairs = 0;
  char *ptr;
  tagint id1, id2;
  int nharm=0, nmorse=0, nlj126=0;

  while(fgets(line,1024,fp)) {
    ptr = strtok(line," \t\n\r\f");

    // skip empty lines
    if (!ptr) continue;
    
    // skip comment lines starting with #
    if (*ptr == '#') continue;

    // get atom ids of pair
    id1 = ATOTAGINT(ptr);
    ptr = strtok(NULL," \t\n\r\f");
    if (!ptr)
      error->all(FLERR,"Incorrectly formatted pair list file");
    id2 = ATOTAGINT(ptr);

    // get potential type
    ptr = strtok(NULL," \t\n\r\f");
    if (!ptr)
      error->all(FLERR,"Incorrectly formatted pair list file");

    style[npairs] = NONE;
    list_parm_t &par = params[npairs];
    par.id1 = id1;
    par.id2 = id2;

    // harmonic potential
    if (strcmp(ptr,stylename[HARM]) == 0) {
      style[npairs] = HARM;

      ptr = strtok(NULL," \t\n\r\f");
      if ((ptr == NULL) || (*ptr == '#'))
	error->all(FLERR,"Incorrectly formatted harmonic pair parameters");
      par.parm.harm.k = force->numeric(FLERR,ptr);

      ptr = strtok(NULL," \t\n\r\f");
      if ((ptr == NULL) || (*ptr == '#'))
	error->all(FLERR,"Incorrectly formatted harmonic pair parameters");
      par.parm.harm.r0 = force->numeric(FLERR,ptr);

      ++nharm;

      // morse potential
    } else if (strcmp(ptr,stylename[MORSE]) == 0) {
      style[npairs] = MORSE;

      ptr = strtok(NULL," \t\n\r\f");
      if (!ptr)
	error->all(FLERR,"Incorrectly formatted morse pair parameters");
      par.parm.morse.d0 = force->numeric(FLERR,ptr);

      ptr = strtok(NULL," \t\n\r\f");
      if (!ptr)
	error->all(FLERR,"Incorrectly formatted morse pair parameters");
      par.parm.morse.alpha = force->numeric(FLERR,ptr);

      ptr = strtok(NULL," \t\n\r\f");
      if (!ptr)
	error->all(FLERR,"Incorrectly formatted morse pair parameters");
      par.parm.morse.r0 = force->numeric(FLERR,ptr);

      ++nmorse;

    } else if (strcmp(ptr,stylename[LJ126]) == 0) {
      // 12-6 lj potential
      style[npairs] = LJ126;

      ptr = strtok(NULL," \t\n\r\f");
      if (!ptr)
	error->all(FLERR,"Incorrectly formatted 12-6 LJ pair parameters");
      par.parm.lj126.epsilon = force->numeric(FLERR,ptr);

      ptr = strtok(NULL," \t\n\r\f");
      if (!ptr)
	error->all(FLERR,"Incorrectly formatted 12-6 LJ pair parameters");
      par.parm.lj126.sigma = force->numeric(FLERR,ptr);

      ++nlj126;

    } else {
      error->all(FLERR,"Unknown pair list potential style");
    }

    // optional cutoff parameter. if not specified use global value
    ptr = strtok(NULL," \t\n\r\f");
    if ((ptr != NULL) && (*ptr != '#')) {
      double cut = force->numeric(FLERR,ptr);
      par.cutsq = cut*cut;
    } else {
      par.cutsq = cut_global*cut_global;
    }

    // count complete entry
    ++npairs;
  }
  fclose(fp);

  // informative output
  if (comm->me == 0) {
    if (screen)
      fprintf(screen,"Read %d (%d/%d/%d) interacting pairs from %s\n",
	      npairs, nharm, nmorse, nlj126, arg[0]);
    if (logfile)
      fprintf(logfile,"Read %d (%d/%d/%d) interacting pairs from %s\n",
	      npairs, nharm, nmorse, nlj126, arg[0]);
  }
}

/* ----------------------------------------------------------------------
   there are no coeffs to be set, but we need to update setflag and pretend
------------------------------------------------------------------------- */

void PairList::coeff(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style: compute energy offset at cutoff
------------------------------------------------------------------------- */

void PairList::init_style()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style list requires atom IDs");

  if (atom->map_style == 0)
    error->all(FLERR,"Pair style list requires an atom map");

  if (offset_flag) {
    for (int n=0; n < npairs; ++n) {
      list_parm_t &par = params[n];

      if (style[n] == HARM) {
	const double dr = sqrt(par.cutsq) - par.parm.harm.r0;
	par.offset = par.parm.harm.k*dr*dr;

      } else if (style[n] == MORSE) {
	const double dr = par.parm.morse.r0 - sqrt(par.cutsq);
	const double dexp = exp(par.parm.morse.alpha * dr);
	par.offset = par.parm.morse.d0 * (dexp*dexp - 2.0*dexp);

      } else if (style[n] == LJ126) {
	const double r6inv = par.cutsq*par.cutsq*par.cutsq;
	const double sig6  = mypow(par.parm.lj126.sigma,6);
	par.offset = 4.0*par.parm.lj126.epsilon*r6inv * (sig6*sig6*r6inv - sig6);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
   since we don't use atom types or neighbor lists, this is a NOP.
------------------------------------------------------------------------- */

double PairList::init_one(int, int)
{
  return cut_global;
}

/* ----------------------------------------------------------------------
   memory usage of each sub-style
------------------------------------------------------------------------- */

double PairList::memory_usage()
{
  double bytes = npairs * sizeof(int);
  bytes += npairs * sizeof(list_parm_t);
  const int n = atom->ntypes+1;
  bytes += n*(n*sizeof(int) + sizeof(int *));
  bytes += n*(n*sizeof(double) + sizeof(double *));
  return bytes;
}
