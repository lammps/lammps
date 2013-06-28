/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
-------------------------------------------------------------------------

/* ----------------------------------------------------------------------
   Contributing author: Pieter J. in 't Veld and Stan Moore (Sandia)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math_const.h"
#include "math_vector.h"
#include "pair_lj_long_dipole_long.h"
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
using namespace MathConst;

#define EWALD_F   1.12837917
#define EWALD_P   0.3275911
#define A1        0.254829592
#define A2       -0.284496736
#define A3        1.421413741
#define A4       -1.453152027
#define A5        1.061405429

// ----------------------------------------------------------------------

PairLJLongDipoleLong::PairLJLongDipoleLong(LAMMPS *lmp) : Pair(lmp)
{
  dispersionflag = ewaldflag = dipoleflag = 1;
  respa_enable = 0;
  single_enable = 0;
}

// ----------------------------------------------------------------------
// global settings
// ----------------------------------------------------------------------

#define PAIR_ILLEGAL	"Illegal pair_style lj/long/dipole/long command"
#define PAIR_CUTOFF	"Only one cut-off allowed when requesting all long"
#define PAIR_MISSING	"Cut-offs missing in pair_style lj/long/dipole/long"
#define PAIR_COUL_CUT	"Coulombic cut not supported in pair_style lj/long/dipole/long"
#define PAIR_LARGEST	"Using largest cut-off for lj/long/dipole/long long long"
#define PAIR_MIX	"Geometric mixing assumed for 1/r^6 coefficients"

void PairLJLongDipoleLong::options(char **arg, int order)
{
  const char *option[] = {"long", "cut", "off", NULL};
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

void PairLJLongDipoleLong::settings(int narg, char **arg)
{
  if (narg != 3 && narg != 4) error->all(FLERR,"Illegal pair_style command");

  ewald_off = 0;
  ewald_order = 0;
  options(arg, 6);
  options(++arg, 3);
  options(arg, 1);
  if (!comm->me && ewald_order&(1<<6))
    error->warning(FLERR,PAIR_MIX);
  if (!comm->me && ewald_order==((1<<3)|(1<<6)))
    error->warning(FLERR,PAIR_LARGEST);
  if (!*(++arg))
    error->all(FLERR,PAIR_MISSING);
  if (!((ewald_order^ewald_off)&(1<<3)))
    error->all(FLERR,PAIR_COUL_CUT);
  cut_lj_global = force->numeric(FLERR,*(arg++));
  if (narg == 4 && (ewald_order==74))
    error->all(FLERR,PAIR_CUTOFF);
  if (narg == 4) cut_coul = force->numeric(FLERR,*(arg++));
  else cut_coul = cut_lj_global;

  if (allocated) {					// reset explicit cuts
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
	if (setflag[i][j]) cut_lj[i][j] = cut_lj_global;
  }
}

// ----------------------------------------------------------------------
// free all arrays
// ----------------------------------------------------------------------

PairLJLongDipoleLong::~PairLJLongDipoleLong()
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
  //if (ftable) free_tables();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairLJLongDipoleLong::allocate()
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

void *PairLJLongDipoleLong::extract(const char *id, int &dim)
{
  const char *ids[] = {
    "B", "sigma", "epsilon", "ewald_order", "ewald_cut", "ewald_mix",
    "cut_coul", "cut_vdwl", NULL};
  void *ptrs[] = {
    lj4, sigma, epsilon, &ewald_order, &cut_coul, &mix_flag, &cut_coul, 
    &cut_lj_global, NULL};
  int i;

  for (i=0; ids[i]&&strcmp(ids[i], id); ++i);
  if (i <= 2) dim = 2;
  else dim = 0;
  return ptrs[i];
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairLJLongDipoleLong::coeff(int narg, char **arg)
{
  if (narg < 4 || narg > 5) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  double epsilon_one = force->numeric(FLERR,arg[2]);
  double sigma_one = force->numeric(FLERR,arg[3]);

  double cut_lj_one = cut_lj_global;
  if (narg == 5) cut_lj_one = force->numeric(FLERR,arg[4]);

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

void PairLJLongDipoleLong::init_style()
{
  const char *style3[] = {"ewald/disp", NULL};
  const char *style6[] = {"ewald/disp", NULL};
  int i;

  if (strcmp(update->unit_style,"electron") == 0)
    error->all(FLERR,"Cannot (yet) use 'electron' units with dipoles");

  // require an atom style with charge defined

  if (!atom->q_flag && (ewald_order&(1<<1)))
    error->all(FLERR,
	"Invoking coulombic in pair style lj/long/dipole/long requires atom attribute q");
  if (!atom->mu && (ewald_order&(1<<3)))
    error->all(FLERR,"Pair lj/long/dipole/long requires atom attributes mu, torque");
  if (!atom->torque && (ewald_order&(1<<3)))
    error->all(FLERR,"Pair lj/long/dipole/long requires atom attributes mu, torque");

  neighbor->request(this);

  cut_coulsq = cut_coul * cut_coul;

  // ensure use of KSpace long-range solver, set g_ewald

  if (ewald_order&(1<<3)) {				// r^-1 kspace
    if (force->kspace == NULL) 
      error->all(FLERR,"Pair style is incompatible with KSpace style");
    for (i=0; style3[i]&&strcmp(force->kspace_style, style3[i]); ++i);
    if (!style3[i])
      error->all(FLERR,"Pair style is incompatible with KSpace style");
  }
  if (ewald_order&(1<<6)) {				// r^-6 kspace
    if (force->kspace == NULL) 
      error->all(FLERR,"Pair style is incompatible with KSpace style");
    for (i=0; style6[i]&&strcmp(force->kspace_style, style6[i]); ++i);
    if (!style6[i])
      error->all(FLERR,"Pair style is incompatible with KSpace style");
  }
  if (force->kspace) g_ewald = force->kspace->g_ewald;
}

/* ----------------------------------------------------------------------
   neighbor callback to inform pair style of neighbor list to use
   regular or rRESPA
------------------------------------------------------------------------- */

void PairLJLongDipoleLong::init_list(int id, NeighList *ptr)
{
  if (id == 0) list = ptr;
  else if (id == 1) listinner = ptr;
  else if (id == 2) listmiddle = ptr;
  else if (id == 3) listouter = ptr;

  if (id)
    error->all(FLERR,"Pair style lj/long/dipole/long does not currently support respa");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLJLongDipoleLong::init_one(int i, int j)
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

  //if (cut_respa && MIN(cut_lj[i][j],cut_coul) < cut_respa[3])
    //error->all(FLERR,"Pair cutoff < Respa interior cutoff");
 
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

void PairLJLongDipoleLong::write_restart(FILE *fp)
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

void PairLJLongDipoleLong::read_restart(FILE *fp)
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

void PairLJLongDipoleLong::write_restart_settings(FILE *fp)
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

void PairLJLongDipoleLong::read_restart_settings(FILE *fp)
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

void PairLJLongDipoleLong::compute(int eflag, int vflag)
{
  double evdwl,ecoul,fpair;
  evdwl = ecoul = 0.0;

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;
  
  double **x = atom->x, *x0 = x[0];
  double **mu = atom->mu, *mu0 = mu[0], *imu, *jmu;
  double **tq = atom->torque, *tq0 = tq[0], *tqi;
  double **f = atom->f, *f0 = f[0], *fi = f0, fx, fy, fz;
  double *q = atom->q, qi = 0, qj;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  int i, j;
  int order1 = ewald_order&(1<<1), order3 = ewald_order&(1<<3),
      order6 = ewald_order&(1<<6);
  int *ineigh, *ineighn, *jneigh, *jneighn, typei, typej, ni;
  double *cutsqi, *cut_ljsqi, *lj1i, *lj2i, *lj3i, *lj4i, *offseti;
  double rsq, r2inv, force_coul, force_lj;
  double g2 = g_ewald*g_ewald, g6 = g2*g2*g2, g8 = g6*g2;
  double B0, B1, B2, B3, G0, G1, G2, mudi, mudj, muij;
  vector force_d = VECTOR_NULL, ti = VECTOR_NULL, tj = VECTOR_NULL;
  vector mui, muj, xi, d;
  
  double C1 = 2.0 * g_ewald / MY_PIS;
  double C2 = 2.0 * g2 * C1;
  double C3 = 2.0 * g2 * C2;

  ineighn = (ineigh = list->ilist)+list->inum;

  for (; ineigh<ineighn; ++ineigh) {			// loop over all neighs
    i = *ineigh; fi = f0+3*i; tqi = tq0+3*i;
    qi = q[i];				// initialize constants
    offseti = offset[typei = type[i]];
    lj1i = lj1[typei]; lj2i = lj2[typei]; lj3i = lj3[typei]; lj4i = lj4[typei];
    cutsqi = cutsq[typei]; cut_ljsqi = cut_ljsq[typei];
    memcpy(xi, x0+(i+(i<<1)), sizeof(vector));
    memcpy(mui, imu = mu0+(i<<2), sizeof(vector));
    
    jneighn = (jneigh = list->firstneigh[i])+list->numneigh[i];

    for (; jneigh<jneighn; ++jneigh) {			// loop over neighbors
      j = *jneigh;
      ni = sbmask(j);					// special index
      j &= NEIGHMASK;
      
      { register double *xj = x0+(j+(j<<1));
	d[0] = xi[0] - xj[0];				// pair vector
	d[1] = xi[1] - xj[1];
	d[2] = xi[2] - xj[2]; }

      if ((rsq = vec_dot(d, d)) >= cutsqi[typej = type[j]]) continue;
      r2inv = 1.0/rsq;

      if (order3 && (rsq < cut_coulsq)) {		// dipole
	memcpy(muj, jmu = mu0+(j<<2), sizeof(vector));
	{						// series real space
	  register double r = sqrt(rsq);
	  register double x = g_ewald*r;
	  register double f = exp(-x*x)*qqrd2e;

	  B0 = 1.0/(1.0+EWALD_P*x);			// eqn 2.8
	  B0 *= ((((A5*B0+A4)*B0+A3)*B0+A2)*B0+A1)*f/r;
	  B1 = (B0 + C1 * f) * r2inv;
	  B2 = (3.0*B1 + C2 * f) * r2inv;
	  B3 = (5.0*B2 + C3 * f) * r2inv;

	  mudi = mui[0]*d[0]+mui[1]*d[1]+mui[2]*d[2];
	  mudj = muj[0]*d[0]+muj[1]*d[1]+muj[2]*d[2];
	  muij = mui[0]*muj[0]+mui[1]*muj[1]+mui[2]*muj[2];
	  G0 = qi*(qj = q[j]);				// eqn 2.10
	  G1 = qi*mudj-qj*mudi+muij;
	  G2 = -mudi*mudj;
	  force_coul = G0*B1+G1*B2+G2*B3;
	  
	  mudi *= B2; mudj *= B2;			// torque contribs
	  ti[0] = mudj*d[0]+(qj*d[0]-muj[0])*B1;
	  ti[1] = mudj*d[1]+(qj*d[1]-muj[1])*B1;
	  ti[2] = mudj*d[2]+(qj*d[2]-muj[2])*B1;

	  if (newton_pair || j < nlocal) {
	    tj[0] = mudi*d[0]-(qi*d[0]+mui[0])*B1;
	    tj[1] = mudi*d[1]-(qi*d[1]+mui[1])*B1;
	    tj[2] = mudi*d[2]-(qi*d[2]+mui[2])*B1;
	  }

	  if (eflag) ecoul = G0*B0+G1*B1+G2*B2;
	  if (ni > 0) {					// adj part, eqn 2.13
	    force_coul -= (f = qqrd2e*(1.0-special_coul[ni])/r)*(
	       	(3.0*G1+15.0*G2*r2inv)*r2inv+G0)*r2inv;
	    if (eflag)
	      ecoul -= f*((G1+3.0*G2*r2inv)*r2inv+G0);
	    B1 -= f*r2inv;
	  }
	  B0 = mudj+qj*B1; B3 = -qi*B1+mudi;		// position independent
      if (ni > 0) B0 -= f*3.0*mudj*r2inv*r2inv/B2;
      if (ni > 0) B3 -= f*3.0*mudi*r2inv*r2inv/B2;
	  force_d[0] = B0*mui[0]+B3*muj[0];		// force contribs
	  force_d[1] = B0*mui[1]+B3*muj[1];
	  force_d[2] = B0*mui[2]+B3*muj[2];
      if (ni > 0) {
	    ti[0] -= f*(3.0*mudj*r2inv*r2inv*d[0]/B2+(qj*r2inv*d[0]-muj[0]*r2inv));
	    ti[1] -= f*(3.0*mudj*r2inv*r2inv*d[1]/B2+(qj*r2inv*d[1]-muj[1]*r2inv));
	    ti[2] -= f*(3.0*mudj*r2inv*r2inv*d[2]/B2+(qj*r2inv*d[2]-muj[2]*r2inv));
	    if (newton_pair || j < nlocal) {
	      tj[0] -= f*(3.0*mudi*r2inv*r2inv*d[0]/B2-(qi*r2inv*d[0]+mui[0]*r2inv));
	      tj[1] -= f*(3.0*mudi*r2inv*r2inv*d[1]/B2-(qi*r2inv*d[1]+mui[1]*r2inv));
	      tj[2] -= f*(3.0*mudi*r2inv*r2inv*d[2]/B2-(qi*r2inv*d[2]+mui[2]*r2inv));
	    }
      }
	}						// table real space
      } else {
	force_coul = ecoul = 0.0;
	memset(force_d, 0, 3*sizeof(double));
      }

      if (rsq < cut_ljsqi[typej]) {			// lj
       	if (order6) {					// long-range lj
	  register double rn = r2inv*r2inv*r2inv;
	  register double x2 = g2*rsq, a2 = 1.0/x2;
	  x2 = a2*exp(-x2)*lj4i[typej];
	  if (ni < 0) {
	    force_lj =
	      (rn*=rn)*lj1i[typej]-g8*(((6.0*a2+6.0)*a2+3.0)*a2+1.0)*x2*rsq;
	    if (eflag) evdwl = rn*lj3i[typej]-g6*((a2+1.0)*a2+0.5)*x2;
	  }
	  else {					// special case
	    register double f = special_lj[ni], t = rn*(1.0-f);
	    force_lj = f*(rn *= rn)*lj1i[typej]-
	      g8*(((6.0*a2+6.0)*a2+3.0)*a2+1.0)*x2*rsq+t*lj2i[typej];
	    if (eflag) evdwl = 
		f*rn*lj3i[typej]-g6*((a2+1.0)*a2+0.5)*x2+t*lj4i[typej];
	  }
	}
	else {						// cut lj
	  register double rn = r2inv*r2inv*r2inv;
	  if (ni < 0) {
	    force_lj = rn*(rn*lj1i[typej]-lj2i[typej]);
	    if (eflag) evdwl = rn*(rn*lj3i[typej]-lj4i[typej])-offseti[typej];
	  }
	  else {					// special case
	    register double f = special_lj[ni];
	    force_lj = f*rn*(rn*lj1i[typej]-lj2i[typej]);
	    if (eflag) evdwl = f*(
		rn*(rn*lj3i[typej]-lj4i[typej])-offseti[typej]);
	  }
	}
	force_lj *= r2inv;
      }
      else force_lj = evdwl = 0.0;

      fpair = force_coul+force_lj;			// force
      if (newton_pair || j < nlocal) {
	register double *fj = f0+(j+(j<<1));
	fi[0] += fx = d[0]*fpair+force_d[0]; fj[0] -= fx;
	fi[1] += fy = d[1]*fpair+force_d[1]; fj[1] -= fy;
	fi[2] += fz = d[2]*fpair+force_d[2]; fj[2] -= fz;
	tqi[0] += mui[1]*ti[2]-mui[2]*ti[1];		// torque
	tqi[1] += mui[2]*ti[0]-mui[0]*ti[2];
	tqi[2] += mui[0]*ti[1]-mui[1]*ti[0];
	register double *tqj = tq0+(j+(j<<1));
	tqj[0] += muj[1]*tj[2]-muj[2]*tj[1];
	tqj[1] += muj[2]*tj[0]-muj[0]*tj[2];
	tqj[2] += muj[0]*tj[1]-muj[1]*tj[0];
      }
      else {
	fi[0] += fx = d[0]*fpair+force_d[0];		// force
	fi[1] += fy = d[1]*fpair+force_d[1];
	fi[2] += fz = d[2]*fpair+force_d[2];
	tqi[0] += mui[1]*ti[2]-mui[2]*ti[1];		// torque
	tqi[1] += mui[2]*ti[0]-mui[0]*ti[2];
	tqi[2] += mui[0]*ti[1]-mui[1]*ti[0];
      }

      if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,
			   evdwl,ecoul,fx,fy,fz,d[0],d[1],d[2]);
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

/*
double PairLJLongDipoleLong::single(int i, int j, int itype, int jtype,
			    double rsq, double factor_coul, double factor_lj,
			    double &fforce)
{
  double r6inv, force_coul, force_lj;
  double g2 = g_ewald*g_ewald, g6 = g2*g2*g2, g8 = g6*g2, *q = atom->q;

  double eng = 0.0;
  double r2inv = 1.0/rsq;

  if ((ewald_order&(1<<3)) && (rsq < cut_coulsq)) {	// coulombic
    double *mui = atom->mu[i], *muj = atom->mu[j];
    double *xi = atom->x[i], *xj = atom->x[j];
    double qi = q[i], qj = q[j];
    double G0, G1, G2, B0, B1, B2, B3, mudi, mudj, muij;
    vector d = {xi[0]-xj[0], xi[1]-xj[1], xi[2]-xj[2]};
    {							// series real space
      register double r = sqrt(rsq);
      register double x = g_ewald*r;
      register double f = exp(-x*x)*qqrd2e;

      B0 = 1.0/(1.0+EWALD_P*x);			// eqn 2.8
      B0 *= ((((A5*B0+A4)*B0+A3)*B0+A2)*B0+A1)*f/r;
      B1 = (B0 + C1 * f) * r2inv;
      B2 = (3.0*B1 + C2 * f) * r2inv;
      B3 = (5.0*B2 + C3 * f) * r2inv;

      mudi = mui[0]*d[0]+mui[1]*d[1]+mui[2]*d[2];
      mudj = muj[0]*d[0]+muj[1]*d[1]+muj[2]*d[2];
      muij = mui[0]*muj[0]+mui[1]*muj[1]+mui[2]*muj[2];
      G0 = qi*(qj = q[j]);				// eqn 2.10
      G1 = qi*mudj-qj*mudi+muij;
      G2 = -mudi*mudj;
      force_coul = G0*B1+G1*B2+G2*B3;
	  
      eng += G0*B0+G1*B1+G2*B2;	
      if (factor_coul < 1.0) {			      	// adj part, eqn 2.13
	force_coul -= (f = force->qqrd2e*(1.0-factor_coul)/r)*(
	    (3.0*G1+6.0*muij+15.0*G2*r2inv)*r2inv+G0);
	eng -= f*((G1+3.0*G2*r2inv)*r2inv+G0);
	B1 -= f*r2inv;
      }
      B0 = mudj*B2-qj*B1; B3 = qi*B1+mudi*B2;		// position independent
      //force_d[0] = B0*mui[0]+B3*muj[0];		// force contributions
      //force_d[1] = B0*mui[1]+B3*muj[1];
      //force_d[2] = B0*mui[2]+B3*muj[2];
    }							// table real space
  }
  else force_coul = 0.0;

  if (rsq < cut_ljsq[itype][jtype]) {			// lennard-jones
    r6inv = r2inv*r2inv*r2inv;
    if (ewald_order&0x40) {				// long-range
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
  } 
  else force_lj = 0.0;

  fforce = (force_coul+force_lj)*r2inv;
  return eng;
}
*/

