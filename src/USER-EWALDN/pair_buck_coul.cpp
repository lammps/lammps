/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

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
#include "pair_buck_coul.h"
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

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#define EWALD_F   1.12837917
#define EWALD_P   0.3275911
#define A1        0.254829592
#define A2       -0.284496736
#define A3        1.421413741
#define A4       -1.453152027
#define A5        1.061405429

/* ---------------------------------------------------------------------- */

PairBuckCoul::PairBuckCoul(LAMMPS *lmp) : Pair(lmp)
{
  respa_enable = 1;
  ftable = NULL;
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

#define PAIR_ILLEGAL	"Illegal pair_style buck/coul command"
#define PAIR_CUTOFF	"Only one cut-off allowed when requesting all long"
#define PAIR_MISSING	"Cut-offs missing in pair_style buck/coul"
#define PAIR_LJ_OFF	"LJ6 off not supported in pair_style buck/coul"
#define PAIR_COUL_CUT	"Coulombic cut not supported in pair_style buck/coul"
#define PAIR_MANY	"Too many pair_style buck/coul commands"
#define PAIR_LARGEST	"Using largest cut-off for buck/coul long long"
#define PAIR_MIX	"Geometric mixing assumed for 1/r^6 coefficients"

void PairBuckCoul::options(char **arg, int order)
{
  char *option[] = {"long", "cut", "off", NULL};
  int i;

  if (!*arg) error->all(PAIR_ILLEGAL);
  for (i=0; option[i]&&strcmp(arg[0], option[i]); ++i);
  switch (i) {
    default: error->all(PAIR_ILLEGAL);
    case 0: ewald_order |= 1<<order; break;		// set kspace r^-order
    case 2: ewald_off |= 1<<order;			// turn r^-order off
    case 1: break;
  }
}


void PairBuckCoul::settings(int narg, char **arg)
{
  ewald_order = 0;
  ewald_off = 0;
  options(arg, 6);
  options(++arg, 1);
  if (!comm->me && ewald_order&(1<<6)) error->warning(PAIR_MIX);
  if (!comm->me && ewald_order==((1<<1)|(1<<6))) error->warning(PAIR_LARGEST);
  if (!*(++arg)) error->all(PAIR_MISSING);
  if (ewald_off&(1<<6)) error->all(PAIR_LJ_OFF);
  if (!((ewald_order^ewald_off)&(1<<1))) error->all(PAIR_COUL_CUT);
  cut_buck_global = force->numeric(*(arg++));
  if (*arg&&(ewald_order&0x42==0x42)) error->all(PAIR_CUTOFF);
  cut_coul = *arg ? force->numeric(*(arg++)) : cut_buck_global;
  if (*arg) error->all(PAIR_MANY);

  if (allocated) {					// reset explicit cuts
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
	if (setflag[i][j]) cut_buck[i][j] = cut_buck_global;
  }
}

/* ----------------------------------------------------------------------
   free all arrays
------------------------------------------------------------------------- */

PairBuckCoul::~PairBuckCoul()
{
  if (allocated) {
    memory->destroy_2d_int_array(setflag);
    memory->destroy_2d_double_array(cutsq);

    memory->destroy_2d_double_array(cut_buck_read);
    memory->destroy_2d_double_array(cut_buck);
    memory->destroy_2d_double_array(cut_bucksq);
    memory->destroy_2d_double_array(buck_a_read);
    memory->destroy_2d_double_array(buck_a);
    memory->destroy_2d_double_array(buck_c_read);
    memory->destroy_2d_double_array(buck_c);
    memory->destroy_2d_double_array(buck_rho_read);
    memory->destroy_2d_double_array(buck_rho);
    memory->destroy_2d_double_array(buck1);
    memory->destroy_2d_double_array(buck2);
    memory->destroy_2d_double_array(rhoinv);
    memory->destroy_2d_double_array(offset);
  }
  if (ftable) free_tables();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairBuckCoul::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  setflag = memory->create_2d_int_array(n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  cutsq = memory->create_2d_double_array(n+1,n+1,"pair:cutsq");

  cut_buck_read = memory->create_2d_double_array(n+1,n+1,"pair:cut_buck_read");
  cut_buck = memory->create_2d_double_array(n+1,n+1,"pair:cut_buck");
  cut_bucksq = memory->create_2d_double_array(n+1,n+1,"pair:cut_bucksq");
  buck_a_read = memory->create_2d_double_array(n+1,n+1,"pair:buck_a_read");
  buck_a = memory->create_2d_double_array(n+1,n+1,"pair:buck_a");
  buck_c_read = memory->create_2d_double_array(n+1,n+1,"pair:buck_c_read");
  buck_c = memory->create_2d_double_array(n+1,n+1,"pair:buck_c");
  buck_rho_read = memory->create_2d_double_array(n+1,n+1,"pair:buck_rho_read");
  buck_rho = memory->create_2d_double_array(n+1,n+1,"pair:buck_rho");
  buck1 = memory->create_2d_double_array(n+1,n+1,"pair:buck1");
  buck2 = memory->create_2d_double_array(n+1,n+1,"pair:buck2");
  rhoinv = memory->create_2d_double_array(n+1,n+1,"pair:rhoinv");
  offset = memory->create_2d_double_array(n+1,n+1,"pair:offset");
}

/* ----------------------------------------------------------------------
   extract protected data from object
------------------------------------------------------------------------- */

void *PairBuckCoul::extract(char *id)
{
  char *ids[] = {
    "B", "sigma", "epsilon", "ewald_order", "ewald_cut", "ewald_mix",
    "cut_coul", NULL};
  void *ptrs[] = {
    buck_c, NULL, NULL, &ewald_order, &cut_coul, &mix_flag, &cut_coul, NULL};
  int i;

  for (i=0; ids[i]&&strcmp(ids[i], id); ++i);
  return ptrs[i];
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairBuckCoul::coeff(int narg, char **arg)
{
  if (narg < 5 || narg > 6) error->all("Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(*(arg++),atom->ntypes,ilo,ihi);
  force->bounds(*(arg++),atom->ntypes,jlo,jhi);

  double buck_a_one = force->numeric(*(arg++));
  double buck_rho_one = force->numeric(*(arg++));
  double buck_c_one = force->numeric(*(arg++));

  double cut_buck_one = cut_buck_global;
  if (narg == 6) cut_buck_one = force->numeric(*(arg++));

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      buck_a_read[i][j] = buck_a_one;
      buck_c_read[i][j] = buck_c_one;
      buck_rho_read[i][j] = buck_rho_one;
      cut_buck_read[i][j] = cut_buck_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all("Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairBuckCoul::init_style()
{

  // require an atom style with charge defined

  if (!atom->q_flag && (ewald_order&(1<<1)))
    error->all(
	"Invoking coulombic in pair style lj/coul requires atom attribute q");

  // request regular or rRESPA neighbor lists

  int irequest;

  if (update->whichflag == 0 && strcmp(update->integrate_style,"respa") == 0) {
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

  if (strcmp(update->integrate_style,"respa") == 0 &&
      ((Respa *) update->integrate)->level_inner >= 0)
    cut_respa = ((Respa *) update->integrate)->cutoff;
  else cut_respa = NULL;

  // ensure use of KSpace long-range solver, set g_ewald

  if (ewald_order&(1<<1)) {				// r^-1 kspace
    if (force->kspace == NULL) 
      error->all("Pair style is incompatible with KSpace style");
    g_ewald = force->kspace->g_ewald;
  }
  if (ewald_order&(1<<6)) {				// r^-6 kspace
    if (!force->kspace && strcmp(force->kspace_style,"ewald/n"))
      error->all("Pair style is incompatible with KSpace style");
    g_ewald = force->kspace->g_ewald;
  }

  // setup force tables

  if (ncoultablebits) init_tables();
}

/* ----------------------------------------------------------------------
   neighbor callback to inform pair style of neighbor list to use
   regular or rRESPA
------------------------------------------------------------------------- */

void PairBuckCoul::init_list(int id, NeighList *ptr)
{
  if (id == 0) list = ptr;
  else if (id == 1) listinner = ptr;
  else if (id == 2) listmiddle = ptr;
  else if (id == 3) listouter = ptr;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairBuckCoul::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all("All pair coeffs are not set");

  cut_buck[i][j] = cut_buck_read[i][j];
  buck_a[i][j] = buck_a_read[i][j];
  buck_c[i][j] = buck_c_read[i][j];
  buck_rho[i][j] = buck_rho_read[i][j];

  double cut = MAX(cut_buck[i][j],cut_coul);
  cutsq[i][j] = cut*cut;
  cut_bucksq[i][j] = cut_buck[i][j] * cut_buck[i][j];

  buck1[i][j] = buck_a[i][j]/buck_rho[i][j];
  buck2[i][j] = 6.0*buck_c[i][j];
  rhoinv[i][j] = 1.0/buck_rho[i][j];

  // check interior rRESPA cutoff

  if (cut_respa && MIN(cut_buck[i][j],cut_coul) < cut_respa[3])
    error->all("Pair cutoff < Respa interior cutoff");
     
  if (offset_flag) {
    double rexp = exp(-cut_buck[i][j]/buck_rho[i][j]);
    offset[i][j] = buck_a[i][j]*rexp - buck_c[i][j]/pow(cut_buck[i][j],6.0);
  } else offset[i][j] = 0.0;

  cutsq[j][i] = cutsq[i][j];
  cut_bucksq[j][i] = cut_bucksq[i][j];
  buck_a[j][i] = buck_a[i][j];
  buck_c[j][i] = buck_c[i][j];
  rhoinv[j][i] = rhoinv[i][j];
  buck1[j][i] = buck1[i][j];
  buck2[j][i] = buck2[i][j];
  offset[j][i] = offset[i][j];

  return cut;
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairBuckCoul::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
	fwrite(&buck_a_read[i][j],sizeof(double),1,fp);
	fwrite(&buck_rho_read[i][j],sizeof(double),1,fp);
	fwrite(&buck_c_read[i][j],sizeof(double),1,fp);
	fwrite(&cut_buck_read[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairBuckCoul::read_restart(FILE *fp)
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
	  fread(&buck_a_read[i][j],sizeof(double),1,fp);
	  fread(&buck_rho_read[i][j],sizeof(double),1,fp);
	  fread(&buck_c_read[i][j],sizeof(double),1,fp);
	  fread(&cut_buck_read[i][j],sizeof(double),1,fp);
	}
	MPI_Bcast(&buck_a_read[i][j],1,MPI_DOUBLE,0,world);
	MPI_Bcast(&buck_rho_read[i][j],1,MPI_DOUBLE,0,world);
	MPI_Bcast(&buck_c_read[i][j],1,MPI_DOUBLE,0,world);
	MPI_Bcast(&cut_buck_read[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairBuckCoul::write_restart_settings(FILE *fp)
{
  fwrite(&cut_buck_global,sizeof(double),1,fp);
  fwrite(&cut_coul,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&ewald_order,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairBuckCoul::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&cut_buck_global,sizeof(double),1,fp);
    fread(&cut_coul,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
    fread(&ewald_order,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_buck_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_coul,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&ewald_order,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   compute pair interactions
------------------------------------------------------------------------- */

void PairBuckCoul::compute(int eflag, int vflag)
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
  int nall = nlocal + atom->nghost;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  int i, j, order1 = ewald_order&(1<<1), order6 = ewald_order&(1<<6);
  int *ineigh, *ineighn, *jneigh, *jneighn, typei, typej, ni;
  double qi, qri, *cutsqi, *cut_bucksqi,
	 *buck1i, *buck2i, *buckai, *buckci, *rhoinvi, *offseti;
  double r, rsq, r2inv, force_coul, force_buck;
  double g2 = g_ewald*g_ewald, g6 = g2*g2*g2, g8 = g6*g2;
  vector xi, d;

  ineighn = (ineigh = list->ilist)+list->inum;

  for (; ineigh<ineighn; ++ineigh) {			// loop over my atoms
    i = *ineigh; fi = f0+3*i;
    if (order1) qri = (qi = q[i])*qqrd2e;		// initialize constants
    offseti = offset[typei = type[i]];
    buck1i = buck1[typei]; buck2i = buck2[typei];
    buckai = buck_a[typei]; buckci = buck_c[typei], rhoinvi = rhoinv[typei];
    cutsqi = cutsq[typei]; cut_bucksqi = cut_bucksq[typei];
    memcpy(xi, x0+(i+(i<<1)), sizeof(vector));
    jneighn = (jneigh = list->firstneigh[i])+list->numneigh[i];

    for (; jneigh<jneighn; ++jneigh) {			// loop over neighbors
      if ((j = *jneigh) < nall) ni = -1;
      else { ni = j/nall; j %= nall; }			// special index
      
      { register double *xj = x0+(j+(j<<1));
	d[0] = xi[0] - xj[0];				// pair vector
	d[1] = xi[1] - xj[1];
	d[2] = xi[2] - xj[2]; }

      if ((rsq = vec_dot(d, d)) >= cutsqi[typej = type[j]]) continue;
      r2inv = 1.0/rsq;
      r = sqrt(rsq);

      if (order1 && (rsq < cut_coulsq)) {		// coulombic
	if (!ncoultablebits || rsq <= tabinnersq) {	// series real space
	  register double x = g_ewald*r;
	  register double s = qri*q[j], t = 1.0/(1.0+EWALD_P*x);
	  if (ni < 0) {
	    s *= g_ewald*exp(-x*x);
	    force_coul = (t *= ((((t*A5+A4)*t+A3)*t+A2)*t+A1)*s/x)+EWALD_F*s;
	    if (eflag) ecoul = t;
	  }
	  else {					// special case
	    register double f = s*(1.0-special_coul[ni])/r;
	    s *= g_ewald*exp(-x*x);
	    force_coul = (t *= ((((t*A5+A4)*t+A3)*t+A2)*t+A1)*s/x)+EWALD_F*s-f;
	    if (eflag) ecoul = t-f;
	  }
	}						// table real space
	else {
	  register union_int_float_t t;
	  t.f = rsq;
	  register const int k = (t.i & ncoulmask) >> ncoulshiftbits;
	  register double f = (rsq-rtable[k])*drtable[k], qiqj = qi*q[j];
	  if (ni < 0) {
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

      if (rsq < cut_bucksqi[typej]) {			// buckingham
	register double rn = r2inv*r2inv*r2inv, 
			expr = exp(-r*rhoinvi[typej]);
	if (order6) {					// long-range
	  register double x2 = g2*rsq, a2 = 1.0/x2;
	  x2 = a2*exp(-x2)*buckci[typej];
	  if (ni < 0) {
	    force_buck =
	      r*expr*buck1i[typej]-g8*(((6.0*a2+6.0)*a2+3.0)*a2+1.0)*x2*rsq;
	    if (eflag) evdwl = expr*buckai[typej]-g6*((a2+1.0)*a2+0.5)*x2;
	  }
	  else {					// special case
	    register double f = special_lj[ni], t = rn*(1.0-f);
	    force_buck = f*r*expr*buck1i[typej]-
	      g8*(((6.0*a2+6.0)*a2+3.0)*a2+1.0)*x2*rsq+t*buck2i[typej];
	    if (eflag) evdwl = f*expr*buckai[typej] - 
			 g6*((a2+1.0)*a2+0.5)*x2+t*buckci[typej];
	  }
	}
	else {						// cut
	  if (ni < 0) {
	    force_buck = r*expr*buck1i[typej]-rn*buck2i[typej];
	    if (eflag) evdwl = expr*buckai[typej] - 
			 rn*buckci[typej]-offseti[typej];
	  }
	  else {					// special case
	    register double f = special_lj[ni];
	    force_buck = f*(r*expr*buck1i[typej]-rn*buck2i[typej]);
	    if (eflag) 
	      evdwl = f*(expr*buckai[typej]-rn*buckci[typej]-offseti[typej]);
	  }
	}
      }
      else force_buck = evdwl = 0.0;

      fpair = (force_coul+force_buck)*r2inv;

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

  if (vflag_fdotr) virial_compute();
}

/* ---------------------------------------------------------------------- */

void PairBuckCoul::compute_inner()
{
  double r, rsq, r2inv, force_coul, force_buck, fpair;

  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
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
  double qri, *cut_bucksqi, *buck1i, *buck2i, *rhoinvi;
  vector xi, d;

  ineighn = (ineigh = listinner->ilist)+listinner->inum;

  for (; ineigh<ineighn; ++ineigh) {			// loop over my atoms
    i = *ineigh; fi = f0+3*i;
    qri = qqrd2e*q[i];
    memcpy(xi, x0+(i+(i<<1)), sizeof(vector));
    cut_bucksqi = cut_bucksq[typei];
    buck1i = buck1[typei]; buck2i = buck2[typei]; rhoinvi = rhoinv[typei];
    jneighn = (jneigh = listinner->firstneigh[i])+listinner->numneigh[i];

    for (; jneigh<jneighn; ++jneigh) {			// loop over neighbors
      if ((j = *jneigh) < nall) ni = -1;
      else { ni = j/nall; j %= nall; }
      
      { register double *xj = x0+(j+(j<<1));
	d[0] = xi[0] - xj[0];				// pair vector
	d[1] = xi[1] - xj[1];
	d[2] = xi[2] - xj[2]; }

      if ((rsq = vec_dot(d, d)) >= cut_out_off_sq) continue;
      r2inv = 1.0/rsq;
      r = sqrt(r);

      if (order1 && (rsq < cut_coulsq))			// coulombic
	force_coul = ni<0 ?
	  qri*q[j]/r : qri*q[j]/r*special_coul[ni];

      if (rsq < cut_bucksqi[typej = type[j]]) {		// buckingham
	register double rn = r2inv*r2inv*r2inv,
			expr = exp(-r*rhoinvi[typej]);
	force_buck = ni<0 ?
	  (r*expr*buck1i[typej]-rn*buck2i[typej]) :
	  (r*expr*buck1i[typej]-rn*buck2i[typej])*special_lj[ni];
      }
      else force_buck = 0.0;

      fpair = (force_coul + force_buck) * r2inv;
      
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

void PairBuckCoul::compute_middle()
{
  double r, rsq, r2inv, force_coul, force_buck, fpair;

  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
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
  double qri, *cut_bucksqi, *buck1i, *buck2i, *rhoinvi;
  vector xi, d;

  ineighn = (ineigh = listmiddle->ilist)+listmiddle->inum;

  for (; ineigh<ineighn; ++ineigh) {			// loop over my atoms
    i = *ineigh; fi = f0+3*i;
    qri = qqrd2e*q[i];
    memcpy(xi, x0+(i+(i<<1)), sizeof(vector));
    cut_bucksqi = cut_bucksq[typei];
    buck1i = buck1[typei]; buck2i = buck2[typei]; rhoinvi = rhoinv[typei];
    jneighn = (jneigh = listmiddle->firstneigh[i])+listmiddle->numneigh[i];

    for (; jneigh<jneighn; ++jneigh) {			// loop over neighbors
      if ((j = *jneigh) < nall) ni = -1;
      else { ni = j/nall; j %= nall; }
      
      { register double *xj = x0+(j+(j<<1));
	d[0] = xi[0] - xj[0];				// pair vector
	d[1] = xi[1] - xj[1];
	d[2] = xi[2] - xj[2]; }

      if ((rsq = vec_dot(d, d)) >= cut_out_off_sq) continue;
      if (rsq <= cut_in_off_sq) continue;
      r2inv = 1.0/rsq;
      r = sqrt(rsq);

      if (order1 && (rsq < cut_coulsq))			// coulombic
	force_coul = ni<0 ?
	  qri*q[j]/r : qri*q[j]/r*special_coul[ni];

      if (rsq < cut_bucksqi[typej = type[j]]) {		// buckingham
	register double rn = r2inv*r2inv*r2inv,
			expr = exp(-r*rhoinvi[typej]);
	force_buck = ni<0 ?
	  (r*expr*buck1i[typej]-rn*buck2i[typej]) :
	  (r*expr*buck1i[typej]-rn*buck2i[typej])*special_lj[ni];
      }
      else force_buck = 0.0;

      fpair = (force_coul + force_buck) * r2inv;
      
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

void PairBuckCoul::compute_outer(int eflag, int vflag)
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
  int nall = nlocal + atom->nghost;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  int i, j, order1 = ewald_order&(1<<1), order6 = ewald_order&(1<<6);
  int *ineigh, *ineighn, *jneigh, *jneighn, typei, typej, ni, respa_flag;
  double qi, qri, *cutsqi, *cut_bucksqi,
	 *buck1i, *buck2i, *buckai, *buckci, *rhoinvi, *offseti;
  double r, rsq, r2inv, force_coul, force_buck;
  double g2 = g_ewald*g_ewald, g6 = g2*g2*g2, g8 = g6*g2;
  double respa_buck, respa_coul, frespa;
  vector xi, d;

  double cut_in_off = cut_respa[2];
  double cut_in_on = cut_respa[3];

  double cut_in_diff = cut_in_on - cut_in_off;
  double cut_in_off_sq = cut_in_off*cut_in_off;
  double cut_in_on_sq = cut_in_on*cut_in_on;

  ineighn = (ineigh = listouter->ilist)+listouter->inum;

  for (; ineigh<ineighn; ++ineigh) {			// loop over my atoms
    i = *ineigh; fi = f0+3*i;
    if (order1) qri = (qi = q[i])*qqrd2e;		// initialize constants
    offseti = offset[typei = type[i]];
    buck1i = buck1[typei]; buck2i = buck2[typei];
    buckai = buck_a[typei]; buckci = buck_c[typei]; rhoinvi = rhoinv[typei];
    cutsqi = cutsq[typei]; cut_bucksqi = cut_bucksq[typei];
    memcpy(xi, x0+(i+(i<<1)), sizeof(vector));
    jneighn = (jneigh = listouter->firstneigh[i])+listouter->numneigh[i];

    for (; jneigh<jneighn; ++jneigh) {			// loop over neighbors
      if ((j = *jneigh) < nall) ni = -1;
      else { ni = j/nall; j %= nall; }			// special index
      
      { register double *xj = x0+(j+(j<<1));
	d[0] = xi[0] - xj[0];				// pair vector
	d[1] = xi[1] - xj[1];
	d[2] = xi[2] - xj[2]; }

      if ((rsq = vec_dot(d, d)) >= cutsqi[typej = type[j]]) continue;
      r2inv = 1.0/rsq;
      r = sqrt(rsq);

      if ((respa_flag = (rsq>cut_in_off_sq)&&(rsq<cut_in_on_sq))) {
	register double rsw = (r-cut_in_off)/cut_in_diff;
	frespa = rsw*rsw*(3.0-2.0*rsw);
      }

      if (order1 && (rsq < cut_coulsq)) {		// coulombic
	if (!ncoultablebits || rsq <= tabinnersq) {	// series real space
	  register double s = qri*q[j];
	  if (respa_flag)				// correct for respa
	    respa_coul = ni<0 ? frespa*s/r : frespa*s/r*special_coul[ni];
	  register double x = g_ewald*r, t = 1.0/(1.0+EWALD_P*x);
	  if (ni < 0) {
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
	  if (respa_flag) respa_coul = ni<0 ?		// correct for respa
	      frespa*qri*q[j]/r :
	      frespa*qri*q[j]/r*special_coul[ni];
	  register union_int_float_t t;
	  t.f = rsq;
	  register const int k = (t.i & ncoulmask) >> ncoulshiftbits;
	  register double f = (rsq-rtable[k])*drtable[k], qiqj = qi*q[j];
	  if (ni < 0) {
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

      if (rsq < cut_bucksqi[typej]) {			// buckingham
	register double rn = r2inv*r2inv*r2inv,
			expr = exp(-r*rhoinvi[typej]);
	if (respa_flag) respa_buck = ni<0 ? 		// correct for respa
	    frespa*(r*expr*buck1i[typej]-rn*buck2i[typej]) :
	    frespa*(r*expr*buck1i[typej]-rn*buck2i[typej])*special_lj[ni];
	if (order6) {					// long-range form
	  register double x2 = g2*rsq, a2 = 1.0/x2;
	  x2 = a2*exp(-x2)*buckci[typej];
	  if (ni < 0) {
	    force_buck =
	      r*expr*buck1i[typej]-g8*(((6.0*a2+6.0)*a2+3.0)*a2+1.0)*x2*rsq;
	    if (eflag) evdwl = expr*buckai[typej]-g6*((a2+1.0)*a2+0.5)*x2;
	  }
	  else {					// correct for special
	    register double f = special_lj[ni], t = rn*(1.0-f);
	    force_buck = f*r*expr*buck1i[typej]-
	      g8*(((6.0*a2+6.0)*a2+3.0)*a2+1.0)*x2*rsq+t*buck2i[typej];
	    if (eflag) evdwl = f*expr*buckai[typej] - 
			 g6*((a2+1.0)*a2+0.5)*x2+t*buckci[typej];
	  }
	}
	else {						// cut form
	  if (ni < 0) {
	    force_buck = r*expr*buck1i[typej]-rn*buck2i[typej];
	    if (eflag) 
	      evdwl = expr*buckai[typej]-rn*buckci[typej]-offseti[typej];
	  }
	  else {					// correct for special
	    register double f = special_lj[ni];
	    force_buck = f*(r*expr*buck1i[typej]-rn*buck2i[typej]);
	    if (eflag)
	      evdwl = f*(expr*buckai[typej]-rn*buckci[typej]-offseti[typej]);
	  }
	}
      }
      else force_buck = respa_buck = evdwl = 0.0;

      fpair = (force_coul+force_buck)*r2inv;
      frespa = fpair-(respa_coul+respa_buck)*r2inv;

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

void PairBuckCoul::init_tables()
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
  
  rtable = (double *) memory->smalloc(ntable*sizeof(double),"pair:rtable");
  ftable = (double *) memory->smalloc(ntable*sizeof(double),"pair:ftable");
  ctable = (double *) memory->smalloc(ntable*sizeof(double),"pair:ctable");
  etable = (double *) memory->smalloc(ntable*sizeof(double),"pair:etable");
  drtable = (double *) memory->smalloc(ntable*sizeof(double),"pair:drtable");
  dftable = (double *) memory->smalloc(ntable*sizeof(double),"pair:dftable");
  dctable = (double *) memory->smalloc(ntable*sizeof(double),"pair:dctable");
  detable = (double *) memory->smalloc(ntable*sizeof(double),"pair:detable");

  if (cut_respa == NULL) {
    vtable = ptable = dvtable = dptable = NULL;
  } else {
    vtable = (double *) memory->smalloc(ntable*sizeof(double),"pair:vtable");
    ptable = (double *) memory->smalloc(ntable*sizeof(double),"pair:ptable");
    dvtable = (double *) memory->smalloc(ntable*sizeof(double),"pair:dvtable");
    dptable = (double *) memory->smalloc(ntable*sizeof(double),"pair:dptable");
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
    r = sqrt(rsq_lookup.f);
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
	
  double f_tmp,c_tmp,e_tmp,p_tmp,v_tmp;
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

void PairBuckCoul::free_tables()
{
  memory->sfree(rtable);
  memory->sfree(drtable);
  memory->sfree(ftable);
  memory->sfree(dftable);
  memory->sfree(ctable);
  memory->sfree(dctable);
  memory->sfree(etable);
  memory->sfree(detable);
  memory->sfree(vtable);
  memory->sfree(dvtable);
  memory->sfree(ptable);
  memory->sfree(dptable);
}

/* ---------------------------------------------------------------------- */

double PairBuckCoul::single(int i, int j, int itype, int jtype,
			    double rsq, double factor_coul, double factor_buck,
			    double &fforce)
{
  double f, r, r2inv, r6inv, force_coul, force_buck;
  double g2 = g_ewald*g_ewald, g6 = g2*g2*g2, g8 = g6*g2, *q = atom->q;

  r = sqrt(rsq);
  r2inv = 1.0/rsq;
  double eng = 0.0;

  if ((ewald_order&2) && (rsq < cut_coulsq)) {		// coulombic
    if (!ncoultablebits || rsq <= tabinnersq) {		// series real space
      register double x = g_ewald*r;
      register double s = force->qqrd2e*q[i]*q[j], t = 1.0/(1.0+EWALD_P*x);
      f = s*(1.0-factor_coul)/r; s *= g_ewald*exp(-x*x);
      force_coul = (t *= ((((t*A5+A4)*t+A3)*t+A2)*t+A1)*s/x)+EWALD_F*s-f;
      eng += t-f;
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
  
  if (rsq < cut_bucksq[itype][jtype]) {			// buckingham
    register double expr = factor_buck*exp(-sqrt(rsq)*rhoinv[itype][jtype]);
    r6inv = r2inv*r2inv*r2inv;
    if (ewald_order&64) {				// long-range
      register double x2 = g2*rsq, a2 = 1.0/x2, t = r6inv*(1.0-factor_buck);
      x2 = a2*exp(-x2)*buck_c[itype][jtype];
      force_buck = buck1[itype][jtype]*r*expr-
       	g8*(((6.0*a2+6.0)*a2+3.0)*a2+a2)*x2*rsq+t*buck2[itype][jtype];
      eng += buck_a[itype][jtype]*expr-
	g6*((a2+1.0)*a2+0.5)*x2+t*buck_c[itype][jtype];
    }
    else {						// cut
      force_buck = 
	buck1[itype][jtype]*r*expr-factor_buck*buck_c[itype][jtype]*r6inv;
      eng += buck_a[itype][jtype]*expr-
	factor_buck*(buck_c[itype][jtype]*r6inv-offset[itype][jtype]);
    }
  } else force_buck = 0.0;

  fforce = (force_coul+force_buck)*r2inv;
  return eng;
}

