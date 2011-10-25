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

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_line_lj.h"
#include "atom.h"
#include "atom_vec_line.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define DELTA 10000

/* ---------------------------------------------------------------------- */

PairLineLJ::PairLineLJ(LAMMPS *lmp) : Pair(lmp)
{
  avec = (AtomVecLine *) atom->style_match("line");
  if (!avec) error->all(FLERR,"Pair line/lj requires atom style line");

  dmax = nmax = 0;
  discrete = NULL;
  dnum = dfirst = NULL;

  single_enable = 0;
  restartinfo = 0;
}

/* ---------------------------------------------------------------------- */

PairLineLJ::~PairLineLJ()
{
  memory->sfree(discrete);
  memory->destroy(dnum);
  memory->destroy(dfirst);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(epsilon);
    memory->destroy(sigma);
    memory->destroy(lj1);
    memory->destroy(lj2);
    memory->destroy(lj3);
    memory->destroy(lj4);
  }
}

/* ---------------------------------------------------------------------- */

void PairLineLJ::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  int ni,nj,npi,npj,ifirst,jfirst;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r2inv,r6inv,term1,term2,sig,sig3,forcelj;
  double xi[2],xj[2],fi[2],fj[2],dxi,dxj,dyi,dyj,ti,tj;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  double **torque = atom->torque;
  int *line = atom->line;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int newton_pair = force->newton_pair;
  
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  
  // grow discrete list if necessary and initialize

  if (nall > nmax) {
    nmax = nall;
    memory->destroy(dnum);
    memory->destroy(dfirst);
    memory->create(dnum,nall,"pair:dnum");
    memory->create(dfirst,nall,"pair:dfirst");
  }
  for (i = 0; i < nall; i++) dnum[i] = 0;
  ndiscrete = 0;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];
      
      if (rsq >= cutsq[itype][jtype]) continue;
      
      // line/line interactions = NxN particles
      
      evdwl = 0.0;
      if (line[i] >= 0 && line[j] >= 0) {
	if (dnum[i] == 0) discretize(i,sigma[itype][itype]);
	npi = dnum[i];
	ifirst = dfirst[i];
	if (dnum[j] == 0) discretize(j,sigma[jtype][jtype]);
	npj = dnum[j];
	jfirst = dfirst[j];

	for (ni = 0; ni < npi; ni++) {
	  dxi = discrete[ifirst+ni].dx;
	  dyi = discrete[ifirst+ni].dy;

	  for (nj = 0; nj < npj; nj++) {
	    dxj = discrete[jfirst+nj].dx;
	    dyj = discrete[jfirst+nj].dy;

	    xi[0] = x[i][0] + dxi;
	    xi[1] = x[i][1] + dyi;
	    xj[0] = x[j][0] + dxj;
	    xj[1] = x[j][1] + dyj;

	    delx = xi[0] - xj[0];
	    dely = xi[1] - xj[1];
	    rsq = delx*delx + dely*dely;

	    sig = 0.5 * (discrete[ifirst+ni].sigma+discrete[jfirst+nj].sigma);
	    sig3 = sig*sig*sig;
	    term2 = 24.0*epsilon[itype][jtype] * sig3*sig3;
	    term1 = 2.0 * term2 * sig3*sig3;
	    r2inv = 1.0/rsq;
	    r6inv = r2inv*r2inv*r2inv;
	    forcelj = r6inv * (term1*r6inv - term2);
	    fpair = forcelj*r2inv;

	    if (eflag) evdwl += r6inv*(term1/12.0*r6inv-term2/6.0);

	    fi[0] = delx*fpair;
	    fi[1] = dely*fpair;
	    f[i][0] += fi[0];
	    f[i][1] += fi[1];
	    torque[i][2] += dxi*fi[1] - dyi*fi[0];

	    if (newton_pair || j < nlocal) {
	      fj[0] = -delx*fpair;
	      fj[1] = -dely*fpair;
	      f[j][0] += fj[0];
	      f[j][1] += fj[1];
	      torque[j][2] += dxj*fj[1] - dyj*fj[0];
	    }
	  }
	}

      // line/particle interaction = Nx1 particles
      // convert line into Np particles based on sigma and line length

      } else if (line[i] >= 0) {
	if (dnum[i] == 0) discretize(i,sigma[itype][itype]);
	npi = dnum[i];
	ifirst = dfirst[i];

	for (ni = 0; ni < npi; ni++) {
	  dxi = discrete[ifirst+ni].dx;
	  dyi = discrete[ifirst+ni].dy;

	  xi[0] = x[i][0] + dxi;
	  xi[1] = x[i][1] + dyi;
	  xj[0] = x[j][0];
	  xj[1] = x[j][1];

	  delx = xi[0] - xj[0];
	  dely = xi[1] - xj[1];
	  rsq = delx*delx + dely*dely;

	  sig = 0.5 * (discrete[ifirst+ni].sigma+sigma[jtype][jtype]);
	  sig3 = sig*sig*sig;
	  term2 = 24.0*epsilon[itype][jtype] * sig3*sig3;
	  term1 = 2.0 * term2 * sig3*sig3;
	  r2inv = 1.0/rsq;
	  r6inv = r2inv*r2inv*r2inv;
	  forcelj = r6inv * (term1*r6inv - term2);
	  fpair = forcelj*r2inv;
	  
	  if (eflag) evdwl += r6inv*(term1/12.0*r6inv-term2/6.0);

	  fi[0] = delx*fpair;
	  fi[1] = dely*fpair;
	  f[i][0] += fi[0];
	  f[i][1] += fi[1];
	  torque[i][2] += dxi*fi[1] - dyi*fi[0];
	  
	  if (newton_pair || j < nlocal) {
	    fj[0] = -delx*fpair;
	    fj[1] = -dely*fpair;
	    f[j][0] += fj[0];
	    f[j][1] += fj[1];
	  }
	}

      // particle/line interaction = Nx1 particles
      // convert line into Np particles based on sigma and line length

      } else if (line[j] >= 0) {
	if (dnum[j] == 0) discretize(j,sigma[jtype][jtype]);
	npj = dnum[j];
	jfirst = dfirst[j];

	for (nj = 0; nj < npj; nj++) {
	  dxj = discrete[jfirst+nj].dx;
	  dyj = discrete[jfirst+nj].dy;

	  xi[0] = x[i][0];
	  xi[1] = x[i][1];
	  xj[0] = x[j][0] + dxj;
	  xj[1] = x[j][1] + dyj;

	  delx = xi[0] - xj[0];
	  dely = xi[1] - xj[1];
	  rsq = delx*delx + dely*dely;

	  sig = 0.5 * (sigma[itype][itype]+discrete[jfirst+nj].sigma);
	  sig3 = sig*sig*sig;
	  term2 = 24.0*epsilon[itype][jtype] * sig3*sig3;
	  term1 = 2.0 * term2 * sig3*sig3;
	  r2inv = 1.0/rsq;
	  r6inv = r2inv*r2inv*r2inv;
	  forcelj = r6inv * (term1*r6inv - term2);
	  fpair = forcelj*r2inv;
	  
	  if (eflag) evdwl += r6inv*(term1/12.0*r6inv-term2/6.0);

	  fi[0] = delx*fpair;
	  fi[1] = dely*fpair;
	  f[i][0] += fi[0];
	  f[i][1] += fi[1];
	  
	  if (newton_pair || j < nlocal) {
	    f[j][0] += fj[0];
	    f[j][1] += fj[1];
	    fj[0] = -delx*fpair;
	    fj[1] = -dely*fpair;
	    torque[j][2] += dxj*fj[1] - dyj*fj[0];
	  }
	}
	
      // particle/particle interaction = 1x1 particles

      } else {
	r2inv = 1.0/rsq;
	r6inv = r2inv*r2inv*r2inv;
	forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
	fpair = forcelj*r2inv;

	if (eflag)
	  evdwl += r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]);

	f[i][0] += delx*fpair;
	f[i][1] += dely*fpair;
	f[i][2] += delz*fpair;
	if (newton_pair || j < nlocal) {
	  f[j][0] -= delx*fpair;
	  f[j][1] -= dely*fpair;
	  f[j][2] -= delz*fpair;
	}
      }

      if (evflag) ev_tally(i,j,nlocal,newton_pair,
			   evdwl,0.0,fpair,delx,dely,delz);
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays 
------------------------------------------------------------------------- */

void PairLineLJ::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(epsilon,n+1,n+1,"pair:epsilon");
  memory->create(sigma,n+1,n+1,"pair:sigma");
  memory->create(lj1,n+1,n+1,"pair:lj1");
  memory->create(lj2,n+1,n+1,"pair:lj2");
  memory->create(lj3,n+1,n+1,"pair:lj3");
  memory->create(lj4,n+1,n+1,"pair:lj4");
}

/* ----------------------------------------------------------------------
   global settings 
------------------------------------------------------------------------- */

void PairLineLJ::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal pair_style command");

  cut_global = force->numeric(arg[0]);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
	if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairLineLJ::coeff(int narg, char **arg)
{
  if (narg < 4 || narg > 5) 
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  double epsilon_one = force->numeric(arg[2]);
  double sigma_one = force->numeric(arg[3]);

  double cut_one = cut_global;
  if (narg == 5) cut_one = force->numeric(arg[4]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLineLJ::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    epsilon[i][j] = mix_energy(epsilon[i][i],epsilon[j][j],
			       sigma[i][i],sigma[j][j]);
    sigma[i][j] = mix_distance(sigma[i][i],sigma[j][j]);
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
  }

  lj1[i][j] = 48.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj2[i][j] = 24.0 * epsilon[i][j] * pow(sigma[i][j],6.0);
  lj3[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj4[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],6.0);

  epsilon[j][i] = epsilon[i][j];
  sigma[j][i] = sigma[i][j];
  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];
  lj3[j][i] = lj3[i][j];
  lj4[j][i] = lj4[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   discretize line segment I into N sub-segments no more than sigma in length
   store new discrete particles in Discrete list
------------------------------------------------------------------------- */

void PairLineLJ::discretize(int i, double sigma)
{
  AtomVecLine::Bonus *bonus = avec->bonus;
  double length = bonus[atom->line[i]].length;
  double theta = bonus[atom->line[i]].theta;
  int n = static_cast<int> (length/sigma) + 1;
  dnum[i] = n;
  dfirst[i] = ndiscrete;

  if (ndiscrete + n > dmax) {
    dmax += DELTA;
    discrete = (Discrete *)
      memory->srealloc(discrete,dmax*sizeof(Discrete),"pair:discrete");
  }

  double *x = atom->x[i];
  sigma = length/n;
  double delta;

  for (int m = 0; m < n; m++) {
    delta = -0.5 + (2*m+1)/(2.0*n);
    discrete[ndiscrete].dx = delta*length*cos(theta);
    discrete[ndiscrete].dy = delta*length*sin(theta);
    discrete[ndiscrete].sigma = sigma;
    ndiscrete++;
  }
}
