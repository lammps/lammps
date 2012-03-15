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
   Contributing author: Pieter J. in 't Veld (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math_vector.h"
#include "pair_lj_coul.h"
#include "atom.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "kspace.h"
#include "update.h"
#include "integrate.h"
#include "respa.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define EWALD_F   1.12837917
#define EWALD_P   0.3275911
#define A1        0.254829592
#define A2       -0.284496736
#define A3        1.421413741
#define A4       -1.453152027
#define A5        1.061405429

/* ---------------------------------------------------------------------- */

PairLJCoul::PairLJCoul(LAMMPS *lmp) : Pair(lmp)
{
  respa_enable = 1;
  ftable = NULL;
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

#define PAIR_ILLEGAL	"Illegal pair_style lj/coul command"
#define PAIR_CUTOFF	"Only one cut-off allowed when requesting all long"
#define PAIR_MISSING	"Cut-offs missing in pair_style lj/coul"
#define PAIR_COUL_CUT	"Coulombic cut not supported in pair_style lj/coul"
#define PAIR_LARGEST	"Using largest cut-off for lj/coul long long"
#define PAIR_MIX	"Mixing forced for lj coefficients"

void PairLJCoul::options(char **arg, int order)
{
  char *option[] = {"long", "cut", "off", NULL};
  int i;

  if (!*arg) error->all(FLERR,PAIR_ILLEGAL);
  for (i=0; option[i]&&strcmp(arg[0], option[i]); ++i);
  switch (i) {
    default: error->all(FLERR,PAIR_ILLEGAL);
    case 0: ewald_order |= 1<<order; break;		// set kspace r^-order
    case 2: ewald_off |= 1<<order;			// turn r^-order off
    case 1: break;
  }
}


void PairLJCoul::settings(int narg, char **arg)
{
  if (narg != 3 && narg != 4) error->all(FLERR,"Illegal pair_style command");

  ewald_off = 0;
  ewald_order = 0;
  options(arg, 6);
  options(++arg, 1);
  if (!comm->me && ewald_order&(1<<6)) error->warning(FLERR,PAIR_MIX);
  if (!comm->me && ewald_order==((1<<1)|(1<<6))) error->warning(FLERR,PAIR_LARGEST);
  if (!*(++arg)) error->all(FLERR,PAIR_MISSING);
  if (!((ewald_order^ewald_off)&(1<<1))) error->all(FLERR,PAIR_COUL_CUT);
  cut_lj_global = force->numeric(*(arg++));
  if (*arg&&(ewald_order&0x42==0x42)) error->all(FLERR,PAIR_CUTOFF);
  if (narg == 4) cut_coul = force->numeric(*arg);
  else cut_coul = cut_lj_global;

  if (allocated) {					// reset explicit cuts
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
	if (setflag[i][j]) cut_lj[i][j] = cut_lj_global;
  }
}

/* ----------------------------------------------------------------------
   free all arrays
------------------------------------------------------------------------- */

PairLJCoul::~PairLJCoul()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut_lj_read);
    memory->destroy(cut_lj);
    memory->destroy(cut_ljsq);
    memory->destroy(epsilon_read);
    memory->destroy(epsilon);
    memory->destroy(sigma_read);
    memory->destroy(sigma);
    memory->destroy(lj1);
    memory->destroy(lj2);
    memory->destroy(lj3);
    memory->destroy(lj4);
    memory->destroy(offset);
  }
  if (ftable) free_tables();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairLJCoul::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut_lj_read,n+1,n+1,"pair:cut_lj_read");
  memory->create(cut_lj,n+1,n+1,"pair:cut_lj");
  memory->create(cut_ljsq,n+1,n+1,"pair:cut_ljsq");
  memory->create(epsilon_read,n+1,n+1,"pair:epsilon_read");
  memory->create(epsilon,n+1,n+1,"pair:epsilon");
  memory->create(sigma_read,n+1,n+1,"pair:sigma_read");
  memory->create(sigma,n+1,n+1,"pair:sigma");
  memory->create(lj1,n+1,n+1,"pair:lj1");
  memory->create(lj2,n+1,n+1,"pair:lj2");
  memory->create(lj3,n+1,n+1,"pair:lj3");
  memory->create(lj4,n+1,n+1,"pair:lj4");
  memory->create(offset,n+1,n+1,"pair:offset");
}

/* ----------------------------------------------------------------------
   extract protected data from object
------------------------------------------------------------------------- */

void *PairLJCoul::extract(const char *id, int &dim)
{
  char *ids[] = {
    "B", "sigma", "epsilon", "ewald_order", "ewald_cut", "ewald_mix",
    "cut_coul", "cut_LJ", NULL};
  void *ptrs[] = {
    lj4, sigma, epsilon, &ewald_order, &cut_coul, &mix_flag, &cut_coul, &cut_lj_global, NULL};
  int i;

  for (i=0; ids[i]&&strcmp(ids[i], id); ++i);
  if (i <= 2) dim = 2;
  else dim = 0;
  return ptrs[i];
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairLJCoul::coeff(int narg, char **arg)
{
  if (narg < 4 || narg > 5) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  double epsilon_one = force->numeric(arg[2]);
  double sigma_one = force->numeric(arg[3]);

  double cut_lj_one = cut_lj_global;
  if (narg == 5) cut_lj_one = force->numeric(arg[4]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon_read[i][j] = epsilon_one;
      sigma_read[i][j] = sigma_one;
      cut_lj_read[i][j] = cut_lj_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLJCoul::init_style()
{
  char *style1[] = {"ewald", "ewald/n", "pppm", NULL};
  char *style6[] = {"ewald/n", NULL};
  int i;

  // require an atom style with charge defined

  if (!atom->q_flag && (ewald_order&(1<<1)))
    error->all(FLERR,
	"Invoking coulombic in pair style lj/coul requires atom attribute q");

  // request regular or rRESPA neighbor lists

  int irequest;

  if (update->whichflag == 0 && strstr(update->integrate_style,"respa")) {
    int respa = 0;
    if (((Respa *) update->integrate)->level_inner >= 0) respa = 1;
    if (((Respa *) update->integrate)->level_middle >= 0) respa = 2;

    if (respa == 0) irequest = neighbor->request(this);
    else if (respa == 1) {
      irequest = neighbor->request(this);
      neighbor->requests[irequest]->id = 1;
      neighbor->requests[irequest]->half = 0;
      neighbor->requests[irequest]->respainner = 1;
      irequest = neighbor->request(this);
      neighbor->requests[irequest]->id = 3;
      neighbor->requests[irequest]->half = 0;
      neighbor->requests[irequest]->respaouter = 1;
    } else {
      irequest = neighbor->request(this);
      neighbor->requests[irequest]->id = 1;
      neighbor->requests[irequest]->half = 0;
      neighbor->requests[irequest]->respainner = 1;
      irequest = neighbor->request(this);
      neighbor->requests[irequest]->id = 2;
      neighbor->requests[irequest]->half = 0;
      neighbor->requests[irequest]->respamiddle = 1;
      irequest = neighbor->request(this);
      neighbor->requests[irequest]->id = 3;
      neighbor->requests[irequest]->half = 0;
      neighbor->requests[irequest]->respaouter = 1;
    }

  } else irequest = neighbor->request(this);

  cut_coulsq = cut_coul * cut_coul;

  // set rRESPA cutoffs

  if (strstr(update->integrate_style,"respa") &&
      ((Respa *) update->integrate)->level_inner >= 0)
    cut_respa = ((Respa *) update->integrate)->cutoff;
  else cut_respa = NULL;

  // ensure use of KSpace long-range solver, set g_ewald

  if (ewald_order&(1<<1)) {				// r^-1 kspace
    if (force->kspace == NULL) 
      error->all(FLERR,"Pair style is incompatible with KSpace style");
    for (i=0; style1[i]&&strcmp(force->kspace_style, style1[i]); ++i);
    if (!style1[i]) error->all(FLERR,"Pair style is incompatible with KSpace style");
  }
  if (ewald_order&(1<<6)) {				// r^-6 kspace
    if (force->kspace == NULL) 
      error->all(FLERR,"Pair style is incompatible with KSpace style");
    for (i=0; style6[i]&&strcmp(force->kspace_style, style6[i]); ++i);
    if (!style6[i]) error->all(FLERR,"Pair style is incompatible with KSpace style");
  }
  if (force->kspace) g_ewald = force->kspace->g_ewald;

  // setup force tables

  if (ncoultablebits) init_tables();
}

/* ----------------------------------------------------------------------
   neighbor callback to inform pair style of neighbor list to use
   regular or rRESPA
------------------------------------------------------------------------- */

void PairLJCoul::init_list(int id, NeighList *ptr)
{
  if (id == 0) list = ptr;
  else if (id == 1) listinner = ptr;
  else if (id == 2) listmiddle = ptr;
  else if (id == 3) listouter = ptr;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLJCoul::init_one(int i, int j)
{
  if ((ewald_order&(1<<6))||(setflag[i][j] == 0)) {
    epsilon[i][j] = mix_energy(epsilon_read[i][i],epsilon_read[j][j],
			       sigma_read[i][i],sigma_read[j][j]);
    sigma[i][j] = mix_distance(sigma_read[i][i],sigma_read[j][j]);
    if (ewald_order&(1<<6))
      cut_lj[i][j] = cut_lj_global;
    else
      cut_lj[i][j] = mix_distance(cut_lj_read[i][i],cut_lj_read[j][j]);
  }
  else {
    sigma[i][j] = sigma_read[i][j];
    epsilon[i][j] = epsilon_read[i][j];
    cut_lj[i][j] = cut_lj_read[i][j];
  }

  double cut = MAX(cut_lj[i][j], cut_coul);
  cutsq[i][j] = cut*cut;
  cut_ljsq[i][j] = cut_lj[i][j] * cut_lj[i][j];

  lj1[i][j] = 48.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj2[i][j] = 24.0 * epsilon[i][j] * pow(sigma[i][j],6.0);
  lj3[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj4[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],6.0);

  // check interior rRESPA cutoff

  if (cut_respa && MIN(cut_lj[i][j],cut_coul) < cut_respa[3])
    error->all(FLERR,"Pair cutoff < Respa interior cutoff");
 
  if (offset_flag) {
    double ratio = sigma[i][j] / cut_lj[i][j];
    offset[i][j] = 4.0 * epsilon[i][j] * (pow(ratio,12.0) - pow(ratio,6.0));
  } else offset[i][j] = 0.0;

  cutsq[j][i] = cutsq[i][j];
  cut_ljsq[j][i] = cut_ljsq[i][j];
  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];
  lj3[j][i] = lj3[i][j];
  lj4[j][i] = lj4[i][j];
  offset[j][i] = offset[i][j];

  return cut;
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJCoul::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
	fwrite(&epsilon_read[i][j],sizeof(double),1,fp);
	fwrite(&sigma_read[i][j],sizeof(double),1,fp);
	fwrite(&cut_lj_read[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJCoul::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
	if (me == 0) {
	  fread(&epsilon_read[i][j],sizeof(double),1,fp);
	  fread(&sigma_read[i][j],sizeof(double),1,fp);
	  fread(&cut_lj_read[i][j],sizeof(double),1,fp);
	}
	MPI_Bcast(&epsilon_read[i][j],1,MPI_DOUBLE,0,world);
	MPI_Bcast(&sigma_read[i][j],1,MPI_DOUBLE,0,world);
	MPI_Bcast(&cut_lj_read[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJCoul::write_restart_settings(FILE *fp)
{
  fwrite(&cut_lj_global,sizeof(double),1,fp);
  fwrite(&cut_coul,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&ewald_order,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJCoul::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&cut_lj_global,sizeof(double),1,fp);
    fread(&cut_coul,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
    fread(&ewald_order,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_lj_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_coul,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&ewald_order,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   compute pair interactions
------------------------------------------------------------------------- */

void PairLJCoul::compute(int eflag, int vflag)
{
  double evdwl,ecoul,fpair;
  evdwl = ecoul = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;
  
  double **x = atom->x, *x0 = x[0];
  double **f = atom->f, *f0 = f[0], *fi = f0;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  int i, j, order1 = ewald_order&(1<<1), order6 = ewald_order&(1<<6);
  int *ineigh, *ineighn, *jneigh, *jneighn, typei, typej, ni;
  double qi = 0.0, qri = 0.0;
  double *cutsqi, *cut_ljsqi, *lj1i, *lj2i, *lj3i, *lj4i, *offseti;
  double rsq, r2inv, force_coul, force_lj;
  double g2 = g_ewald*g_ewald, g6 = g2*g2*g2, g8 = g6*g2;
  vector xi, d;

  ineighn = (ineigh = list->ilist)+list->inum;

  for (; ineigh<ineighn; ++ineigh) {			// loop over my atoms
    i = *ineigh; fi = f0+3*i;
    if (order1) qri = (qi = q[i])*qqrd2e;		// initialize constants
    offseti = offset[typei = type[i]];
    lj1i = lj1[typei]; lj2i = lj2[typei]; lj3i = lj3[typei]; lj4i = lj4[typei];
    cutsqi = cutsq[typei]; cut_ljsqi = cut_ljsq[typei];
    memcpy(xi, x0+(i+(i<<1)), sizeof(vector));
    jneighn = (jneigh = list->firstneigh[i])+list->numneigh[i];

    for (; jneigh<jneighn; ++jneigh) {			// loop over neighbors
      j = *jneigh;
      ni = sbmask(j);
      j &= NEIGHMASK;
      
      { register double *xj = x0+(j+(j<<1));
	d[0] = xi[0] - xj[0];				// pair vector
	d[1] = xi[1] - xj[1];
	d[2] = xi[2] - xj[2]; }

      if ((rsq = vec_dot(d, d)) >= cutsqi[typej = type[j]]) continue;
      r2inv = 1.0/rsq;

      if (order1 && (rsq < cut_coulsq)) {		// coulombic
	if (!ncoultablebits || rsq <= tabinnersq) {	// series real space
	  register double r = sqrt(rsq), x = g_ewald*r;
	  register double s = qri*q[j], t = 1.0/(1.0+EWALD_P*x);
	  if (ni == 0) {
	    s *= g_ewald*exp(-x*x);
	    force_coul = (t *= ((((t*A5+A4)*t+A3)*t+A2)*t+A1)*s/x)+EWALD_F*s;
	    if (eflag) ecoul = t;
	  }
	  else {					// special case
	    r = s*(1.0-special_coul[ni])/r; s *= g_ewald*exp(-x*x);
	    force_coul = (t *= ((((t*A5+A4)*t+A3)*t+A2)*t+A1)*s/x)+EWALD_F*s-r;
	    if (eflag) ecoul = t-r;
	  }
	}						// table real space
	else {
	  register union_int_float_t t;
	  t.f = rsq;
	  register const int k = (t.i & ncoulmask)>>ncoulshiftbits;
	  register double f = (rsq-rtable[k])*drtable[k], qiqj = qi*q[j];
	  if (ni == 0) {
	    force_coul = qiqj*(ftable[k]+f*dftable[k]);
	    if (eflag) ecoul = qiqj*(etable[k]+f*detable[k]);
	  }
	  else {					// special case
	    t.f = (1.0-special_coul[ni])*(ctable[k]+f*dctable[k]);
	    force_coul = qiqj*(ftable[k]+f*dftable[k]-t.f);
	    if (eflag) ecoul = qiqj*(etable[k]+f*detable[k]-t.f);
	  }
	}
      }
      else force_coul = ecoul = 0.0;

      if (rsq < cut_ljsqi[typej]) {			// lj
       	if (order6) {					// long-range lj
	  register double rn = r2inv*r2inv*r2inv;
	  register double x2 = g2*rsq, a2 = 1.0/x2;
	  x2 = a2*exp(-x2)*lj4i[typej];
	  if (ni == 0) {
	    force_lj =
	      (rn*=rn)*lj1i[typej]-g8*(((6.0*a2+6.0)*a2+3.0)*a2+1.0)*x2*rsq;
	    if (eflag)
	      evdwl = rn*lj3i[typej]-g6*((a2+1.0)*a2+0.5)*x2;
	  }
	  else {					// special case
	    register double f = special_lj[ni], t = rn*(1.0-f);
	    force_lj = f*(rn *= rn)*lj1i[typej]-
	      g8*(((6.0*a2+6.0)*a2+3.0)*a2+1.0)*x2*rsq+t*lj2i[typej];
	    if (eflag) 
	      evdwl = f*rn*lj3i[typej]-g6*((a2+1.0)*a2+0.5)*x2+t*lj4i[typej];
	  }
	}
	else {						// cut lj
	  register double rn = r2inv*r2inv*r2inv;
	  if (ni == 0) {
	    force_lj = rn*(rn*lj1i[typej]-lj2i[typej]);
	    if (eflag) evdwl = rn*(rn*lj3i[typej]-lj4i[typej])-offseti[typej];
	  }
	  else {					// special case
	    register double f = special_lj[ni];
	    force_lj = f*rn*(rn*lj1i[typej]-lj2i[typej]);
	    if (eflag)
	      evdwl = f * (rn*(rn*lj3i[typej]-lj4i[typej])-offseti[typej]);
	  }
	}
      }
      else force_lj = evdwl = 0.0;

      fpair = (force_coul+force_lj)*r2inv;

      if (newton_pair || j < nlocal) {
	register double *fj = f0+(j+(j<<1)), f;
	fi[0] += f = d[0]*fpair; fj[0] -= f;
	fi[1] += f = d[1]*fpair; fj[1] -= f;
	fi[2] += f = d[2]*fpair; fj[2] -= f;
      }
      else {
	fi[0] += d[0]*fpair;
	fi[1] += d[1]*fpair;
	fi[2] += d[2]*fpair;
      }
      
      if (evflag) ev_tally(i,j,nlocal,newton_pair,
			   evdwl,ecoul,fpair,d[0],d[1],d[2]);
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

void PairLJCoul::compute_inner()
{
  double rsq, r2inv, force_coul = 0.0, force_lj, fpair;

  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *x0 = atom->x[0], *f0 = atom->f[0], *fi = f0, *q = atom->q;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  double cut_out_on = cut_respa[0];
  double cut_out_off = cut_respa[1];
  
  double cut_out_diff = cut_out_off - cut_out_on;
  double cut_out_on_sq = cut_out_on*cut_out_on;
  double cut_out_off_sq = cut_out_off*cut_out_off;

  int *ineigh, *ineighn, *jneigh, *jneighn, typei, typej, ni;
  int i, j, order1 = (ewald_order|(ewald_off^-1))&(1<<1);
  double qri, *cut_ljsqi, *lj1i, *lj2i;
  vector xi, d;

  ineighn = (ineigh = list->ilist)+list->inum;

  for (; ineigh<ineighn; ++ineigh) {			// loop over my atoms
    i = *ineigh; fi = f0+3*i;
    qri = qqrd2e*q[i];
    memcpy(xi, x0+(i+(i<<1)), sizeof(vector));
    cut_ljsqi = cut_ljsq[typei = type[i]];
    lj1i = lj1[typei]; lj2i = lj2[typei];
    jneighn = (jneigh = list->firstneigh[i])+list->numneigh[i];

    for (; jneigh<jneighn; ++jneigh) {			// loop over neighbors
      j = *jneigh;
      ni = sbmask(j);
      j &= NEIGHMASK;
      
      { register double *xj = x0+(j+(j<<1));
	d[0] = xi[0] - xj[0];				// pair vector
	d[1] = xi[1] - xj[1];
	d[2] = xi[2] - xj[2]; }

      if ((rsq = vec_dot(d, d)) >= cut_out_off_sq) continue;
      r2inv = 1.0/rsq;

      if (order1 && (rsq < cut_coulsq))			// coulombic
	force_coul = ni == 0 ?
	  qri*q[j]*sqrt(r2inv) : qri*q[j]*sqrt(r2inv)*special_coul[ni];

      if (rsq < cut_ljsqi[typej = type[j]]) {		// lennard-jones
	register double rn = r2inv*r2inv*r2inv;
	force_lj = ni == 0 ?
	  rn*(rn*lj1i[typej]-lj2i[typej]) :
	  rn*(rn*lj1i[typej]-lj2i[typej])*special_lj[ni];
      }
      else force_lj = 0.0;

      fpair = (force_coul + force_lj) * r2inv;
      
      if (rsq > cut_out_on_sq) {			// switching
        register double rsw = (sqrt(rsq) - cut_out_on)/cut_out_diff; 
	fpair  *= 1.0 + rsw*rsw*(2.0*rsw-3.0);
      }

      if (newton_pair || j < nlocal) {			// force update
	register double *fj = f0+(j+(j<<1)), f;
	fi[0] += f = d[0]*fpair; fj[0] -= f;
	fi[1] += f = d[1]*fpair; fj[1] -= f;
	fi[2] += f = d[2]*fpair; fj[2] -= f;
      }
      else {
	fi[0] += d[0]*fpair;
	fi[1] += d[1]*fpair;
	fi[2] += d[2]*fpair;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairLJCoul::compute_middle()
{
  double rsq, r2inv, force_coul = 0.0, force_lj, fpair;

  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *x0 = atom->x[0], *f0 = atom->f[0], *fi = f0, *q = atom->q;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  double cut_in_off = cut_respa[0];
  double cut_in_on = cut_respa[1];
  double cut_out_on = cut_respa[2];
  double cut_out_off = cut_respa[3];

  double cut_in_diff = cut_in_on - cut_in_off;
  double cut_out_diff = cut_out_off - cut_out_on;
  double cut_in_off_sq = cut_in_off*cut_in_off;
  double cut_in_on_sq = cut_in_on*cut_in_on;
  double cut_out_on_sq = cut_out_on*cut_out_on;
  double cut_out_off_sq = cut_out_off*cut_out_off;

  int *ineigh, *ineighn, *jneigh, *jneighn, typei, typej, ni;
  int i, j, order1 = (ewald_order|(ewald_off^-1))&(1<<1);
  double qri, *cut_ljsqi, *lj1i, *lj2i;
  vector xi, d;

  ineighn = (ineigh = list->ilist)+list->inum;

  for (; ineigh<ineighn; ++ineigh) {			// loop over my atoms
    i = *ineigh; fi = f0+3*i;
    qri = qqrd2e*q[i];
    memcpy(xi, x0+(i+(i<<1)), sizeof(vector));
    cut_ljsqi = cut_ljsq[typei = type[i]];
    lj1i = lj1[typei]; lj2i = lj2[typei];
    jneighn = (jneigh = list->firstneigh[i])+list->numneigh[i];

    for (; jneigh<jneighn; ++jneigh) {
      j = *jneigh;
      ni = sbmask(j);
      j &= NEIGHMASK;
      
      { register double *xj = x0+(j+(j<<1));
	d[0] = xi[0] - xj[0];				// pair vector
	d[1] = xi[1] - xj[1];
	d[2] = xi[2] - xj[2]; }

      if ((rsq = vec_dot(d, d)) >= cut_out_off_sq) continue;
      if (rsq <= cut_in_off_sq) continue;
      r2inv = 1.0/rsq;

      if (order1 && (rsq < cut_coulsq))			// coulombic
	force_coul = ni == 0 ?
	  qri*q[j]*sqrt(r2inv) : qri*q[j]*sqrt(r2inv)*special_coul[ni];

      if (rsq < cut_ljsqi[typej = type[j]]) {		// lennard-jones
	register double rn = r2inv*r2inv*r2inv;
	force_lj = ni == 0 ?
	  rn*(rn*lj1i[typej]-lj2i[typej]) :
	  rn*(rn*lj1i[typej]-lj2i[typej])*special_lj[ni];
      }
      else force_lj = 0.0;

      fpair = (force_coul + force_lj) * r2inv;
      
      if (rsq < cut_in_on_sq) {				// switching
        register double rsw = (sqrt(rsq) - cut_in_off)/cut_in_diff; 
	fpair  *= rsw*rsw*(3.0 - 2.0*rsw);
      }
      if (rsq > cut_out_on_sq) {
        register double rsw = (sqrt(rsq) - cut_out_on)/cut_out_diff; 
	fpair  *= 1.0 + rsw*rsw*(2.0*rsw-3.0);
      }

      if (newton_pair || j < nlocal) {			// force update
	register double *fj = f0+(j+(j<<1)), f;
	fi[0] += f = d[0]*fpair; fj[0] -= f;
	fi[1] += f = d[1]*fpair; fj[1] -= f;
	fi[2] += f = d[2]*fpair; fj[2] -= f;
      }
      else {
	fi[0] += d[0]*fpair;
	fi[1] += d[1]*fpair;
	fi[2] += d[2]*fpair;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairLJCoul::compute_outer(int eflag, int vflag)
{
  double evdwl,ecoul,fpair;
  evdwl = ecoul = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;

  double **x = atom->x, *x0 = x[0];
  double **f = atom->f, *f0 = f[0], *fi = f0;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  int i, j, order1 = ewald_order&(1<<1), order6 = ewald_order&(1<<6);
  int *ineigh, *ineighn, *jneigh, *jneighn, typei, typej, ni, respa_flag;
  double qi = 0.0, qri = 0.0; 
  double *cutsqi, *cut_ljsqi, *lj1i, *lj2i, *lj3i, *lj4i, *offseti;
  double rsq, r2inv, force_coul, force_lj;
  double g2 = g_ewald*g_ewald, g6 = g2*g2*g2, g8 = g6*g2;
  double respa_lj = 0.0, respa_coul = 0.0, frespa = 0.0;
  vector xi, d;

  double cut_in_off = cut_respa[2];
  double cut_in_on = cut_respa[3];

  double cut_in_diff = cut_in_on - cut_in_off;
  double cut_in_off_sq = cut_in_off*cut_in_off;
  double cut_in_on_sq = cut_in_on*cut_in_on;

  ineighn = (ineigh = list->ilist)+list->inum;

  for (; ineigh<ineighn; ++ineigh) {			// loop over my atoms
    i = *ineigh; fi = f0+3*i;
    if (order1) qri = (qi = q[i])*qqrd2e;		// initialize constants
    offseti = offset[typei = type[i]];
    lj1i = lj1[typei]; lj2i = lj2[typei]; lj3i = lj3[typei]; lj4i = lj4[typei];
    cutsqi = cutsq[typei]; cut_ljsqi = cut_ljsq[typei];
    memcpy(xi, x0+(i+(i<<1)), sizeof(vector));
    jneighn = (jneigh = list->firstneigh[i])+list->numneigh[i];

    for (; jneigh<jneighn; ++jneigh) {			// loop over neighbors
      j = *jneigh;
      ni = sbmask(j);
      j &= NEIGHMASK;

      { register double *xj = x0+(j+(j<<1));
	d[0] = xi[0] - xj[0];				// pair vector
	d[1] = xi[1] - xj[1];
	d[2] = xi[2] - xj[2]; }

      if ((rsq = vec_dot(d, d)) >= cutsqi[typej = type[j]]) continue;
      r2inv = 1.0/rsq;

      if ((respa_flag = (rsq>cut_in_off_sq)&&(rsq<cut_in_on_sq))) {
	register double rsw = (sqrt(rsq)-cut_in_off)/cut_in_diff;
	frespa = rsw*rsw*(3.0-2.0*rsw);
      }

      if (order1 && (rsq < cut_coulsq)) {		// coulombic
	if (!ncoultablebits || rsq <= tabinnersq) {	// series real space
	  register double r = sqrt(rsq), s = qri*q[j];
	  if (respa_flag)				// correct for respa
	    respa_coul = ni == 0 ? frespa*s/r : frespa*s/r*special_coul[ni];
	  register double x = g_ewald*r, t = 1.0/(1.0+EWALD_P*x);
	  if (ni == 0) {
	    s *= g_ewald*exp(-x*x);
	    force_coul = (t *= ((((t*A5+A4)*t+A3)*t+A2)*t+A1)*s/x)+EWALD_F*s;
	    if (eflag) ecoul = t;
	  }
	  else {					// correct for special
	    r = s*(1.0-special_coul[ni])/r; s *= g_ewald*exp(-x*x);
	    force_coul = (t *= ((((t*A5+A4)*t+A3)*t+A2)*t+A1)*s/x)+EWALD_F*s-r;
	    if (eflag) ecoul = t-r;
	  }
	}						// table real space
	else {
	  if (respa_flag) respa_coul = ni == 0 ?	// correct for respa
	      frespa*qri*q[j]/sqrt(rsq) :
	      frespa*qri*q[j]/sqrt(rsq)*special_coul[ni];
	  register union_int_float_t t;
	  t.f = rsq;
	  register const int k = (t.i & ncoulmask) >> ncoulshiftbits;
	  register double f = (rsq-rtable[k])*drtable[k], qiqj = qi*q[j];
	  if (ni == 0) {
	    force_coul = qiqj*(ftable[k]+f*dftable[k]);
	    if (eflag) ecoul = qiqj*(etable[k]+f*detable[k]);
	  }
	  else {					// correct for special
	    t.f = (1.0-special_coul[ni])*(ctable[k]+f*dctable[k]);
	    force_coul = qiqj*(ftable[k]+f*dftable[k]-t.f);
	    if (eflag) ecoul = qiqj*(etable[k]+f*detable[k]-t.f);
	  }
	}
      }
      else force_coul = respa_coul = ecoul = 0.0;

      if (rsq < cut_ljsqi[typej]) {			// lennard-jones
	register double rn = r2inv*r2inv*r2inv;
	if (respa_flag) respa_lj = ni == 0 ? 		// correct for respa
	    frespa*rn*(rn*lj1i[typej]-lj2i[typej]) :
	    frespa*rn*(rn*lj1i[typej]-lj2i[typej])*special_lj[ni];
	if (order6) {					// long-range form
	  register double x2 = g2*rsq, a2 = 1.0/x2;
	  x2 = a2*exp(-x2)*lj4i[typej];
	  if (ni == 0) {
	    force_lj =
	      (rn*=rn)*lj1i[typej]-g8*(((6.0*a2+6.0)*a2+3.0)*a2+1.0)*x2*rsq;
	    if (eflag) evdwl = rn*lj3i[typej]-g6*((a2+1.0)*a2+0.5)*x2;
	  }
	  else {					// correct for special
	    register double f = special_lj[ni], t = rn*(1.0-f);
	    force_lj = f*(rn *= rn)*lj1i[typej]-
	      g8*(((6.0*a2+6.0)*a2+3.0)*a2+1.0)*x2*rsq+t*lj2i[typej];
	    if (eflag)
	      evdwl = f*rn*lj3i[typej]-g6*((a2+1.0)*a2+0.5)*x2+t*lj4i[typej];
	  }
	}
	else {						// cut form
	  if (ni == 0) {
	    force_lj = rn*(rn*lj1i[typej]-lj2i[typej]);
	    if (eflag) evdwl = rn*(rn*lj3i[typej]-lj4i[typej])-offseti[typej];
	  }
	  else {					// correct for special
	    register double f = special_lj[ni];
	    force_lj = f*rn*(rn*lj1i[typej]-lj2i[typej]);
	    if (eflag)
	      evdwl = f*(rn*(rn*lj3i[typej]-lj4i[typej])-offseti[typej]);
	  }
	}
      }
      else force_lj = respa_lj = evdwl = 0.0;

      fpair = (force_coul+force_lj)*r2inv;
      frespa = fpair-(respa_coul+respa_lj)*r2inv;

      if (newton_pair || j < nlocal) {
	register double *fj = f0+(j+(j<<1)), f;
	fi[0] += f = d[0]*frespa; fj[0] -= f;
	fi[1] += f = d[1]*frespa; fj[1] -= f;
	fi[2] += f = d[2]*frespa; fj[2] -= f;
      }
      else {
	fi[0] += d[0]*frespa;
	fi[1] += d[1]*frespa;
	fi[2] += d[2]*frespa;
      }

      if (evflag) ev_tally(i,j,nlocal,newton_pair,
			   evdwl,ecoul,fpair,d[0],d[1],d[2]);
    }
  }
}

/* ----------------------------------------------------------------------
   setup force tables used in compute routines
------------------------------------------------------------------------- */

void PairLJCoul::init_tables()
{
  int masklo,maskhi;
  double r,grij,expm2,derfc,rsw;
  double qqrd2e = force->qqrd2e;

  tabinnersq = tabinner*tabinner;
  init_bitmap(tabinner,cut_coul,ncoultablebits,
	      masklo,maskhi,ncoulmask,ncoulshiftbits);
  
  int ntable = 1;
  for (int i = 0; i < ncoultablebits; i++) ntable *= 2;
  
  // linear lookup tables of length N = 2^ncoultablebits
  // stored value = value at lower edge of bin
  // d values = delta from lower edge to upper edge of bin

  if (ftable) free_tables();
  
  memory->create(rtable,ntable,"pair:rtable");
  memory->create(ftable,ntable,"pair:ftable");
  memory->create(ctable,ntable,"pair:ctable");
  memory->create(etable,ntable,"pair:etable");
  memory->create(drtable,ntable,"pair:drtable");
  memory->create(dftable,ntable,"pair:dftable");
  memory->create(dctable,ntable,"pair:dctable");
  memory->create(detable,ntable,"pair:detable");

  if (cut_respa == NULL) {
    vtable = ptable = dvtable = dptable = NULL;
  } else {
    memory->create(vtable,ntable,"pair:vtable");
    memory->create(ptable,ntable,"pair:ptable");
    memory->create(dvtable,ntable,"pair:dvtable");
    memory->create(dptable,ntable,"pair:dptable");
  }

  union_int_float_t rsq_lookup;
  union_int_float_t minrsq_lookup;
  int itablemin;
  minrsq_lookup.i = 0 << ncoulshiftbits;
  minrsq_lookup.i |= maskhi;
    
  for (int i = 0; i < ntable; i++) {
    rsq_lookup.i = i << ncoulshiftbits;
    rsq_lookup.i |= masklo;
    if (rsq_lookup.f < tabinnersq) {
      rsq_lookup.i = i << ncoulshiftbits;
      rsq_lookup.i |= maskhi;
    }
    r = sqrtf(rsq_lookup.f);
    grij = g_ewald * r;
    expm2 = exp(-grij*grij);
    derfc = erfc(grij);
    if (cut_respa == NULL) {
      rtable[i] = rsq_lookup.f;
      ftable[i] = qqrd2e/r * (derfc + EWALD_F*grij*expm2);
      ctable[i] = qqrd2e/r;
      etable[i] = qqrd2e/r * derfc;
    } else {
      rtable[i] = rsq_lookup.f;
      ftable[i] = qqrd2e/r * (derfc + EWALD_F*grij*expm2 - 1.0);
      ctable[i] = 0.0;
      etable[i] = qqrd2e/r * derfc;
      ptable[i] = qqrd2e/r;
      vtable[i] = qqrd2e/r * (derfc + EWALD_F*grij*expm2);
      if (rsq_lookup.f > cut_respa[2]*cut_respa[2]) {
	if (rsq_lookup.f < cut_respa[3]*cut_respa[3]) {
	  rsw = (r - cut_respa[2])/(cut_respa[3] - cut_respa[2]); 
	  ftable[i] += qqrd2e/r * rsw*rsw*(3.0 - 2.0*rsw);
	  ctable[i] = qqrd2e/r * rsw*rsw*(3.0 - 2.0*rsw);
	} else {
	  ftable[i] = qqrd2e/r * (derfc + EWALD_F*grij*expm2);
	  ctable[i] = qqrd2e/r;
	}
      }
    }
    minrsq_lookup.f = MIN(minrsq_lookup.f,rsq_lookup.f);
  }

  tabinnersq = minrsq_lookup.f;
  
  int ntablem1 = ntable - 1;
  
  for (int i = 0; i < ntablem1; i++) {
    drtable[i] = 1.0/(rtable[i+1] - rtable[i]);
    dftable[i] = ftable[i+1] - ftable[i];
    dctable[i] = ctable[i+1] - ctable[i];
    detable[i] = etable[i+1] - etable[i];
  }

  if (cut_respa) {
    for (int i = 0; i < ntablem1; i++) {
      dvtable[i] = vtable[i+1] - vtable[i];
      dptable[i] = ptable[i+1] - ptable[i];
    }
  }
  
  // get the delta values for the last table entries 
  // tables are connected periodically between 0 and ntablem1
    
  drtable[ntablem1] = 1.0/(rtable[0] - rtable[ntablem1]);
  dftable[ntablem1] = ftable[0] - ftable[ntablem1];
  dctable[ntablem1] = ctable[0] - ctable[ntablem1];
  detable[ntablem1] = etable[0] - etable[ntablem1];
  if (cut_respa) {
    dvtable[ntablem1] = vtable[0] - vtable[ntablem1];
    dptable[ntablem1] = ptable[0] - ptable[ntablem1];
  }

  // get the correct delta values at itablemax    
  // smallest r is in bin itablemin
  // largest r is in bin itablemax, which is itablemin-1,
  //   or ntablem1 if itablemin=0
  // deltas at itablemax only needed if corresponding rsq < cut*cut
  // if so, compute deltas between rsq and cut*cut 
	
  double f_tmp,c_tmp,e_tmp,p_tmp = 0.0,v_tmp = 0.0;
  itablemin = minrsq_lookup.i & ncoulmask;
  itablemin >>= ncoulshiftbits;  
  int itablemax = itablemin - 1; 
  if (itablemin == 0) itablemax = ntablem1;     
  rsq_lookup.i = itablemax << ncoulshiftbits;
  rsq_lookup.i |= maskhi;

  if (rsq_lookup.f < cut_coulsq) {
    rsq_lookup.f = cut_coulsq;  
    r = sqrtf(rsq_lookup.f);
    grij = g_ewald * r;
    expm2 = exp(-grij*grij);
    derfc = erfc(grij);

    if (cut_respa == NULL) {
      f_tmp = qqrd2e/r * (derfc + EWALD_F*grij*expm2);
      c_tmp = qqrd2e/r;
      e_tmp = qqrd2e/r * derfc;
    } else {
      f_tmp = qqrd2e/r * (derfc + EWALD_F*grij*expm2 - 1.0);
      c_tmp = 0.0;
      e_tmp = qqrd2e/r * derfc;
      p_tmp = qqrd2e/r;
      v_tmp = qqrd2e/r * (derfc + EWALD_F*grij*expm2);
      if (rsq_lookup.f > cut_respa[2]*cut_respa[2]) {
        if (rsq_lookup.f < cut_respa[3]*cut_respa[3]) {
          rsw = (r - cut_respa[2])/(cut_respa[3] - cut_respa[2]); 
          f_tmp += qqrd2e/r * rsw*rsw*(3.0 - 2.0*rsw);
          c_tmp = qqrd2e/r * rsw*rsw*(3.0 - 2.0*rsw);
        } else {
          f_tmp = qqrd2e/r * (derfc + EWALD_F*grij*expm2);
          c_tmp = qqrd2e/r;
        }
      }
    }

    drtable[itablemax] = 1.0/(rsq_lookup.f - rtable[itablemax]);   
    dftable[itablemax] = f_tmp - ftable[itablemax];
    dctable[itablemax] = c_tmp - ctable[itablemax];
    detable[itablemax] = e_tmp - etable[itablemax];
    if (cut_respa) {
      dvtable[itablemax] = v_tmp - vtable[itablemax];
      dptable[itablemax] = p_tmp - ptable[itablemax];
    }   
  }
}

/* ----------------------------------------------------------------------
   free memory for tables used in pair computations
------------------------------------------------------------------------- */

void PairLJCoul::free_tables()
{
  memory->destroy(rtable);
  memory->destroy(drtable);
  memory->destroy(ftable);
  memory->destroy(dftable);
  memory->destroy(ctable);
  memory->destroy(dctable);
  memory->destroy(etable);
  memory->destroy(detable);
  memory->destroy(vtable);
  memory->destroy(dvtable);
  memory->destroy(ptable);
  memory->destroy(dptable);
}

/* ---------------------------------------------------------------------- */

double PairLJCoul::single(int i, int j, int itype, int jtype,
			  double rsq, double factor_coul, double factor_lj,
			  double &fforce)
{
  double r2inv, r6inv, force_coul, force_lj;
  double g2 = g_ewald*g_ewald, g6 = g2*g2*g2, g8 = g6*g2, *q = atom->q;

  double eng = 0.0;

  r2inv = 1.0/rsq;
  if ((ewald_order&2) && (rsq < cut_coulsq)) {		// coulombic
    if (!ncoultablebits || rsq <= tabinnersq) {		// series real space
      register double r = sqrt(rsq), x = g_ewald*r;
      register double s = force->qqrd2e*q[i]*q[j], t = 1.0/(1.0+EWALD_P*x);
      r = s*(1.0-factor_coul)/r; s *= g_ewald*exp(-x*x);
      force_coul = (t *= ((((t*A5+A4)*t+A3)*t+A2)*t+A1)*s/x)+EWALD_F*s-r;
      eng += t-r;
    }
    else {						// table real space
      register union_int_float_t t;
      t.f = rsq;
      register const int k = (t.i & ncoulmask) >> ncoulshiftbits;
      register double f = (rsq-rtable[k])*drtable[k], qiqj = q[i]*q[j];
      t.f = (1.0-factor_coul)*(ctable[k]+f*dctable[k]);
      force_coul = qiqj*(ftable[k]+f*dftable[k]-t.f);
      eng += qiqj*(etable[k]+f*detable[k]-t.f);
    }
  } else force_coul = 0.0;
  
  if (rsq < cut_ljsq[itype][jtype]) {			// lennard-jones
    r6inv = r2inv*r2inv*r2inv;
    if (ewald_order&64) {				// long-range
      register double x2 = g2*rsq, a2 = 1.0/x2, t = r6inv*(1.0-factor_lj);
      x2 = a2*exp(-x2)*lj4[itype][jtype];
      force_lj = factor_lj*(r6inv *= r6inv)*lj1[itype][jtype]-
       	g8*(((6.0*a2+6.0)*a2+3.0)*a2+a2)*x2*rsq+t*lj2[itype][jtype];
      eng += factor_lj*r6inv*lj3[itype][jtype]-
	g6*((a2+1.0)*a2+0.5)*x2+t*lj4[itype][jtype];
    }
    else {						// cut
      force_lj = factor_lj*r6inv*(lj1[itype][jtype]*r6inv-lj2[itype][jtype]);
      eng += factor_lj*(r6inv*(r6inv*lj3[itype][jtype]-
			       lj4[itype][jtype])-offset[itype][jtype]);
    }
  } else force_lj = 0.0;

  fforce = (force_coul+force_lj)*r2inv;
  return eng;
}
