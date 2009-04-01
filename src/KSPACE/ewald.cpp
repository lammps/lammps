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
   Contributing authors: Roy Pollock (LLNL), Paul Crozier (SNL)
------------------------------------------------------------------------- */

#include "mpi.h"
#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "ewald.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "pair.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define SMALL 0.00001

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

Ewald::Ewald(LAMMPS *lmp, int narg, char **arg) : KSpace(lmp, narg, arg)
{
  if (narg != 1) error->all("Illegal kspace_style ewald command");

  precision = atof(arg[0]);
  PI = 4.0*atan(1.0);

  kmax = 0;
  kxvecs = kyvecs = kzvecs = NULL;
  ug = NULL;
  eg = vg = NULL;
  sfacrl = sfacim = sfacrl_all = sfacim_all = NULL;

  nmax = 0;
  ek = NULL;
  cs = sn = NULL;

  kcount = 0;
}

/* ----------------------------------------------------------------------
   free all memory 
------------------------------------------------------------------------- */

Ewald::~Ewald()
{
  deallocate();
  memory->destroy_2d_double_array(ek);
  memory->destroy_3d_double_array(cs,-kmax_created);
  memory->destroy_3d_double_array(sn,-kmax_created);
}

/* ---------------------------------------------------------------------- */

void Ewald::init()
{
  if (comm->me == 0) {
    if (screen) fprintf(screen,"Ewald initialization ...\n");
    if (logfile) fprintf(logfile,"Ewald initialization ...\n");
  }

  // error check

  if (domain->triclinic) error->all("Cannot use Ewald with triclinic box");
  if (domain->dimension == 2) 
    error->all("Cannot use Ewald with 2d simulation");

  if (!atom->q_flag) error->all("Kspace style requires atom attribute q");

  if (slabflag == 0 && domain->nonperiodic > 0)
    error->all("Cannot use nonperiodic boundaries with Ewald");
  if (slabflag == 1) {
    if (domain->xperiodic != 1 || domain->yperiodic != 1 || 
	domain->boundary[2][0] != 1 || domain->boundary[2][1] != 1)
      error->all("Incorrect boundaries with slab Ewald");
  }

  // extract short-range Coulombic cutoff from pair style

  qqrd2e = force->qqrd2e;

  if (force->pair == NULL)
    error->all("KSpace style is incompatible with Pair style");
  double *p_cutoff = (double *) force->pair->extract("cut_coul");
  if (p_cutoff == NULL)
    error->all("KSpace style is incompatible with Pair style");
  double cutoff = *p_cutoff;

  qsum = qsqsum = 0.0;
  for (int i = 0; i < atom->nlocal; i++) {
    qsum += atom->q[i];
    qsqsum += atom->q[i]*atom->q[i];
  }

  double tmp;
  MPI_Allreduce(&qsum,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
  qsum = tmp;
  MPI_Allreduce(&qsqsum,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
  qsqsum = tmp;

  if (qsqsum == 0.0)
    error->all("Cannot use kspace solver on system with no charge");
  if (fabs(qsum) > SMALL && comm->me == 0) {
    char str[128];
    sprintf(str,"System is not charge neutral, net charge = %g",qsum);
    error->warning(str);
  }

  // setup K-space resolution

  g_ewald = (1.35 - 0.15*log(precision))/cutoff;
  gsqmx = -4.0*g_ewald*g_ewald*log(precision);

  if (comm->me == 0) {
    if (screen) fprintf(screen,"  G vector = %g\n",g_ewald);
    if (logfile) fprintf(logfile,"  G vector = %g\n",g_ewald);
  }

  // setup Ewald coefficients so can print stats

  setup();

  if (comm->me == 0) {
    if (screen) fprintf(screen,"  vectors: actual 1d max = %d %d %d\n",
			kcount,kmax,kmax3d);
    if (logfile) fprintf(logfile,"  vectors: actual 1d max = %d %d %d\n",
			 kcount,kmax,kmax3d);
  }
}

/* ----------------------------------------------------------------------
   adjust Ewald coeffs, called initially and whenever volume has changed 
------------------------------------------------------------------------- */

void Ewald::setup()
{
  // volume-dependent factors

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;
  
  // adjustment of z dimension for 2d slab Ewald
  // 3d Ewald just uses zprd since slab_volfactor = 1.0

  double zprd_slab = zprd*slab_volfactor;
  volume = xprd * yprd * zprd_slab;

  unitk[0] = 2.0*PI/xprd;
  unitk[1] = 2.0*PI/yprd;
  unitk[2] = 2.0*PI/zprd_slab;

  // determine kmax
  // function of current box size, precision, G_ewald (short-range cutoff)

  int nkxmx = static_cast<int> ((g_ewald*xprd/PI) * sqrt(-log(precision)));
  int nkymx = static_cast<int> ((g_ewald*yprd/PI) * sqrt(-log(precision)));
  int nkzmx = 
    static_cast<int> ((g_ewald*zprd_slab/PI) * sqrt(-log(precision)));

  int kmax_old = kmax;
  kmax = MAX(nkxmx,nkymx);
  kmax = MAX(kmax,nkzmx);
  kmax3d = 4*kmax*kmax*kmax + 6*kmax*kmax + 3*kmax;

  // if size has grown, reallocate k-dependent and nlocal-dependent arrays

  if (kmax > kmax_old) {
    deallocate();
    allocate();

    memory->destroy_2d_double_array(ek);
    memory->destroy_3d_double_array(cs,-kmax_created);
    memory->destroy_3d_double_array(sn,-kmax_created);
    nmax = atom->nmax;
    ek = memory->create_2d_double_array(nmax,3,"ewald:ek");
    cs = memory->create_3d_double_array(-kmax,kmax,3,nmax,"ewald:cs");
    sn = memory->create_3d_double_array(-kmax,kmax,3,nmax,"ewald:sn");
    kmax_created = kmax;
  }

  // pre-compute Ewald coefficients

  coeffs();
}

/* ----------------------------------------------------------------------
   compute the Ewald long-range force, energy, virial 
------------------------------------------------------------------------- */

void Ewald::compute(int eflag, int vflag)
{
  int i,k,n;

  energy = 0.0;
  if (vflag) for (n = 0; n < 6; n++) virial[n] = 0.0;

  // extend size of per-atom arrays if necessary

  if (atom->nlocal > nmax) {
    memory->destroy_2d_double_array(ek);
    memory->destroy_3d_double_array(cs,-kmax_created);
    memory->destroy_3d_double_array(sn,-kmax_created);
    nmax = atom->nmax;
    ek = memory->create_2d_double_array(nmax,3,"ewald:ek");
    cs = memory->create_3d_double_array(-kmax,kmax,3,nmax,"ewald:cs");
    sn = memory->create_3d_double_array(-kmax,kmax,3,nmax,"ewald:sn");
    kmax_created = kmax;
  }

  // partial structure factors on each processor
  // total structure factor by summing over procs

  eik_dot_r();
  MPI_Allreduce(sfacrl,sfacrl_all,kcount,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(sfacim,sfacim_all,kcount,MPI_DOUBLE,MPI_SUM,world);

  // K-space portion of electric field
  // double loop over K-vectors and local atoms

  double **f = atom->f;
  double *q = atom->q;
  int nlocal = atom->nlocal;

  int kx,ky,kz;
  double cypz,sypz,exprl,expim,partial;

  for (i = 0; i < nlocal; i++) {
    ek[i][0] = 0.0;
    ek[i][1] = 0.0;
    ek[i][2] = 0.0;
  }

  for (k = 0; k < kcount; k++) {
    kx = kxvecs[k];
    ky = kyvecs[k];
    kz = kzvecs[k];

    for (i = 0; i < nlocal; i++) {
      cypz = cs[ky][1][i]*cs[kz][2][i] - sn[ky][1][i]*sn[kz][2][i];
      sypz = sn[ky][1][i]*cs[kz][2][i] + cs[ky][1][i]*sn[kz][2][i];
      exprl = cs[kx][0][i]*cypz - sn[kx][0][i]*sypz;
      expim = sn[kx][0][i]*cypz + cs[kx][0][i]*sypz;
      partial = expim*sfacrl_all[k] - exprl*sfacim_all[k];
      ek[i][0] += partial*eg[k][0];
      ek[i][1] += partial*eg[k][1];
      ek[i][2] += partial*eg[k][2];
    }
  }

  // convert E-field to force

  for (i = 0; i < nlocal; i++) {
    f[i][0] += qqrd2e*q[i]*ek[i][0];
    f[i][1] += qqrd2e*q[i]*ek[i][1];
    f[i][2] += qqrd2e*q[i]*ek[i][2];
  }
 
  // energy if requested

  if (eflag) {
    for (k = 0; k < kcount; k++)
      energy += ug[k] * (sfacrl_all[k]*sfacrl_all[k] + 
			 sfacim_all[k]*sfacim_all[k]);
    PI = 4.0*atan(1.0);
    energy -= g_ewald*qsqsum/1.772453851 + 
      0.5*PI*qsum*qsum / (g_ewald*g_ewald*volume);
    energy *= qqrd2e;
  }

  // virial if requested

  if (vflag) {
    double uk;
    for (k = 0; k < kcount; k++) {
      uk = ug[k] * (sfacrl_all[k]*sfacrl_all[k] + sfacim_all[k]*sfacim_all[k]);
      for (n = 0; n < 6; n++) virial[n] += uk*vg[k][n];
    }
    for (n = 0; n < 6; n++) virial[n] *= qqrd2e;
  }

  if (slabflag) slabcorr(eflag);
  
}

/* ---------------------------------------------------------------------- */

void Ewald::eik_dot_r()
{
  int i,k,l,m,n,ic;
  double cstr1,sstr1,cstr2,sstr2,cstr3,sstr3,cstr4,sstr4;
  double sqk,clpm,slpm;

  double **x = atom->x;
  double *q = atom->q;
  int nlocal = atom->nlocal;

  n = 0;

  // (k,0,0), (0,l,0), (0,0,m)

  for (ic = 0; ic < 3; ic++) {
    sqk = unitk[ic]*unitk[ic];
    if (sqk <= gsqmx) {
      cstr1 = 0.0;
      sstr1 = 0.0;
      for (i = 0; i < nlocal; i++) {
	cs[0][ic][i] = 1.0;
	sn[0][ic][i] = 0.0;
	cs[1][ic][i] = cos(unitk[ic]*x[i][ic]);
	sn[1][ic][i] = sin(unitk[ic]*x[i][ic]);
	cs[-1][ic][i] = cs[1][ic][i];
	sn[-1][ic][i] = -sn[1][ic][i];
	cstr1 += q[i]*cs[1][ic][i];
	sstr1 += q[i]*sn[1][ic][i];
      }
      sfacrl[n] = cstr1;
      sfacim[n++] = sstr1;
    }
  }

  for (m = 2; m <= kmax; m++) {
    for (ic = 0; ic < 3; ic++) {
      sqk = m*unitk[ic] * m*unitk[ic];
      if (sqk <= gsqmx) {
	cstr1 = 0.0;
	sstr1 = 0.0;
	for (i = 0; i < nlocal; i++) {
	  cs[m][ic][i] = cs[m-1][ic][i]*cs[1][ic][i] - 
	    sn[m-1][ic][i]*sn[1][ic][i];
	  sn[m][ic][i] = sn[m-1][ic][i]*cs[1][ic][i] + 
	    cs[m-1][ic][i]*sn[1][ic][i];
	  cs[-m][ic][i] = cs[m][ic][i];
	  sn[-m][ic][i] = -sn[m][ic][i];
	  cstr1 += q[i]*cs[m][ic][i];
	  sstr1 += q[i]*sn[m][ic][i];
	}
	sfacrl[n] = cstr1;
	sfacim[n++] = sstr1;
      }
    }
  }

  // 1 = (k,l,0), 2 = (k,-l,0)

  for (k = 1; k <= kmax; k++) {
    for (l = 1; l <= kmax; l++) {
      sqk = (k*unitk[0] * k*unitk[0]) + (l*unitk[1] * l*unitk[1]);
      if (sqk <= gsqmx) {
	cstr1 = 0.0;
	sstr1 = 0.0;
	cstr2 = 0.0;
	sstr2 = 0.0;
	for (i = 0; i < nlocal; i++) {
	  cstr1 += q[i]*(cs[k][0][i]*cs[l][1][i] - sn[k][0][i]*sn[l][1][i]);
	  sstr1 += q[i]*(sn[k][0][i]*cs[l][1][i] + cs[k][0][i]*sn[l][1][i]);
	  cstr2 += q[i]*(cs[k][0][i]*cs[l][1][i] + sn[k][0][i]*sn[l][1][i]);
	  sstr2 += q[i]*(sn[k][0][i]*cs[l][1][i] - cs[k][0][i]*sn[l][1][i]);
	}
	sfacrl[n] = cstr1;
	sfacim[n++] = sstr1;
	sfacrl[n] = cstr2;
	sfacim[n++] = sstr2;
      }
    }
  }

  // 1 = (0,l,m), 2 = (0,l,-m)

  for (l = 1; l <= kmax; l++) {
    for (m = 1; m <= kmax; m++) {
      sqk = (l*unitk[1] * l*unitk[1]) + (m*unitk[2] * m*unitk[2]);
      if (sqk <= gsqmx) {
	cstr1 = 0.0;
	sstr1 = 0.0;
	cstr2 = 0.0;
	sstr2 = 0.0;
	for (i = 0; i < nlocal; i++) {
	  cstr1 += q[i]*(cs[l][1][i]*cs[m][2][i] - sn[l][1][i]*sn[m][2][i]);
	  sstr1 += q[i]*(sn[l][1][i]*cs[m][2][i] + cs[l][1][i]*sn[m][2][i]);
	  cstr2 += q[i]*(cs[l][1][i]*cs[m][2][i] + sn[l][1][i]*sn[m][2][i]);
	  sstr2 += q[i]*(sn[l][1][i]*cs[m][2][i] - cs[l][1][i]*sn[m][2][i]);
	}
	sfacrl[n] = cstr1;
	sfacim[n++] = sstr1;
	sfacrl[n] = cstr2;
	sfacim[n++] = sstr2;
      }
    }
  }

  // 1 = (k,0,m), 2 = (k,0,-m)

  for (k = 1; k <= kmax; k++) {
    for (m = 1; m <= kmax; m++) {
      sqk = (k*unitk[0] * k*unitk[0]) + (m*unitk[2] * m*unitk[2]);
      if (sqk <= gsqmx) {
	cstr1 = 0.0;
	sstr1 = 0.0;
	cstr2 = 0.0;
	sstr2 = 0.0;
	for (i = 0; i < nlocal; i++) {
	  cstr1 += q[i]*(cs[k][0][i]*cs[m][2][i] - sn[k][0][i]*sn[m][2][i]);
	  sstr1 += q[i]*(sn[k][0][i]*cs[m][2][i] + cs[k][0][i]*sn[m][2][i]);
	  cstr2 += q[i]*(cs[k][0][i]*cs[m][2][i] + sn[k][0][i]*sn[m][2][i]);
	  sstr2 += q[i]*(sn[k][0][i]*cs[m][2][i] - cs[k][0][i]*sn[m][2][i]);
	}
	sfacrl[n] = cstr1;
	sfacim[n++] = sstr1;
	sfacrl[n] = cstr2;
	sfacim[n++] = sstr2;
      }
    }
  }

  // 1 = (k,l,m), 2 = (k,-l,m), 3 = (k,l,-m), 4 = (k,-l,-m)

  for (k = 1; k <= kmax; k++) {
    for (l = 1; l <= kmax; l++) {
      for (m = 1; m <= kmax; m++) {
	sqk = (k*unitk[0] * k*unitk[0]) + (l*unitk[1] * l*unitk[1]) +
	  (m*unitk[2] * m*unitk[2]);
	if (sqk <= gsqmx) {
	  cstr1 = 0.0;
	  sstr1 = 0.0;
	  cstr2 = 0.0;
	  sstr2 = 0.0;
	  cstr3 = 0.0;
	  sstr3 = 0.0;
	  cstr4 = 0.0;
	  sstr4 = 0.0;
	  for (i = 0; i < nlocal; i++) {
	    clpm = cs[l][1][i]*cs[m][2][i] - sn[l][1][i]*sn[m][2][i];
	    slpm = sn[l][1][i]*cs[m][2][i] + cs[l][1][i]*sn[m][2][i];
	    cstr1 += q[i]*(cs[k][0][i]*clpm - sn[k][0][i]*slpm);
	    sstr1 += q[i]*(sn[k][0][i]*clpm + cs[k][0][i]*slpm);
	    
	    clpm = cs[l][1][i]*cs[m][2][i] + sn[l][1][i]*sn[m][2][i];
	    slpm = -sn[l][1][i]*cs[m][2][i] + cs[l][1][i]*sn[m][2][i];
	    cstr2 += q[i]*(cs[k][0][i]*clpm - sn[k][0][i]*slpm);
	    sstr2 += q[i]*(sn[k][0][i]*clpm + cs[k][0][i]*slpm);
	    
	    clpm = cs[l][1][i]*cs[m][2][i] + sn[l][1][i]*sn[m][2][i];
	    slpm = sn[l][1][i]*cs[m][2][i] - cs[l][1][i]*sn[m][2][i];
	    cstr3 += q[i]*(cs[k][0][i]*clpm - sn[k][0][i]*slpm);
	    sstr3 += q[i]*(sn[k][0][i]*clpm + cs[k][0][i]*slpm);
	    
	    clpm = cs[l][1][i]*cs[m][2][i] - sn[l][1][i]*sn[m][2][i];
	    slpm = -sn[l][1][i]*cs[m][2][i] - cs[l][1][i]*sn[m][2][i];
	    cstr4 += q[i]*(cs[k][0][i]*clpm - sn[k][0][i]*slpm);
	    sstr4 += q[i]*(sn[k][0][i]*clpm + cs[k][0][i]*slpm);
	  }
	  sfacrl[n] = cstr1;
	  sfacim[n++] = sstr1;
	  sfacrl[n] = cstr2;
	  sfacim[n++] = sstr2;
	  sfacrl[n] = cstr3;
	  sfacim[n++] = sstr3;
	  sfacrl[n] = cstr4;
	  sfacim[n++] = sstr4;
	}
      }
    }
  }
}

/* ----------------------------------------------------------------------
   pre-compute coefficients for each Ewald K-vector 
------------------------------------------------------------------------- */

void Ewald::coeffs()
{
  int k,l,m;
  double sqk,vterm;

  double unitkx = unitk[0];
  double unitky = unitk[1];
  double unitkz = unitk[2];
  double g_ewald_sq_inv = 1.0 / (g_ewald*g_ewald);
  double preu = 4.0*PI/volume;

  kcount = 0;

  // (k,0,0), (0,l,0), (0,0,m)

  for (m = 1; m <= kmax; m++) {
    sqk = (m*unitkx) * (m*unitkx);
    if (sqk <= gsqmx) {
      kxvecs[kcount] = m;
      kyvecs[kcount] = 0;
      kzvecs[kcount] = 0;
      ug[kcount] = preu*exp(-0.25*sqk*g_ewald_sq_inv)/sqk;
      eg[kcount][0] = 2.0*unitkx*m*ug[kcount];
      eg[kcount][1] = 0.0;
      eg[kcount][2] = 0.0;
      vterm = -2.0*(1.0/sqk + 0.25*g_ewald_sq_inv);
      vg[kcount][0] = 1.0 + vterm*(unitkx*m)*(unitkx*m);
      vg[kcount][1] = 1.0;
      vg[kcount][2] = 1.0;
      vg[kcount][3] = 0.0;
      vg[kcount][4] = 0.0;
      vg[kcount][5] = 0.0;
      kcount++;
    }
    sqk = (m*unitky) * (m*unitky);
    if (sqk <= gsqmx) {
      kxvecs[kcount] = 0;
      kyvecs[kcount] = m;
      kzvecs[kcount] = 0;
      ug[kcount] = preu*exp(-0.25*sqk*g_ewald_sq_inv)/sqk;
      eg[kcount][0] = 0.0;
      eg[kcount][1] = 2.0*unitky*m*ug[kcount];
      eg[kcount][2] = 0.0;
      vterm = -2.0*(1.0/sqk + 0.25*g_ewald_sq_inv);
      vg[kcount][0] = 1.0;
      vg[kcount][1] = 1.0 + vterm*(unitky*m)*(unitky*m);
      vg[kcount][2] = 1.0;
      vg[kcount][3] = 0.0;
      vg[kcount][4] = 0.0;
      vg[kcount][5] = 0.0;
      kcount++;
    }
    sqk = (m*unitkz) * (m*unitkz);
    if (sqk <= gsqmx) {
      kxvecs[kcount] = 0;
      kyvecs[kcount] = 0;
      kzvecs[kcount] = m;
      ug[kcount] = preu*exp(-0.25*sqk*g_ewald_sq_inv)/sqk;
      eg[kcount][0] = 0.0;
      eg[kcount][1] = 0.0;
      eg[kcount][2] = 2.0*unitkz*m*ug[kcount];
      vterm = -2.0*(1.0/sqk + 0.25*g_ewald_sq_inv);
      vg[kcount][0] = 1.0;
      vg[kcount][1] = 1.0;
      vg[kcount][2] = 1.0 + vterm*(unitkz*m)*(unitkz*m);
      vg[kcount][3] = 0.0;
      vg[kcount][4] = 0.0;
      vg[kcount][5] = 0.0;
      kcount++;
    }
  }

  // 1 = (k,l,0), 2 = (k,-l,0)

  for (k = 1; k <= kmax; k++) {
    for (l = 1; l <= kmax; l++) {
      sqk = (unitkx*k) * (unitkx*k) + (unitky*l) * (unitky*l);
      if (sqk <= gsqmx) {
	kxvecs[kcount] = k;
	kyvecs[kcount] = l;
	kzvecs[kcount] = 0;
	ug[kcount] = preu*exp(-0.25*sqk*g_ewald_sq_inv)/sqk;
	eg[kcount][0] = 2.0*unitkx*k*ug[kcount];
	eg[kcount][1] = 2.0*unitky*l*ug[kcount];
	eg[kcount][2] = 0.0;
	vterm = -2.0*(1.0/sqk + 0.25*g_ewald_sq_inv);
	vg[kcount][0] = 1.0 + vterm*(unitkx*k)*(unitkx*k);
	vg[kcount][1] = 1.0 + vterm*(unitky*l)*(unitky*l);
	vg[kcount][2] = 1.0;
	vg[kcount][3] = vterm*unitkx*k*unitky*l;
	vg[kcount][4] = 0.0;
	vg[kcount][5] = 0.0;
	kcount++;

	kxvecs[kcount] = k;
	kyvecs[kcount] = -l;
	kzvecs[kcount] = 0;
	ug[kcount] = preu*exp(-0.25*sqk*g_ewald_sq_inv)/sqk;
	eg[kcount][0] = 2.0*unitkx*k*ug[kcount];
	eg[kcount][1] = -2.0*unitky*l*ug[kcount];
	eg[kcount][2] = 0.0;
	vg[kcount][0] = 1.0 + vterm*(unitkx*k)*(unitkx*k);
	vg[kcount][1] = 1.0 + vterm*(unitky*l)*(unitky*l);
	vg[kcount][2] = 1.0;
	vg[kcount][3] = -vterm*unitkx*k*unitky*l;
	vg[kcount][4] = 0.0;
	vg[kcount][5] = 0.0;
	kcount++;;
      }
    }
  }

  // 1 = (0,l,m), 2 = (0,l,-m)

  for (l = 1; l <= kmax; l++) {
    for (m = 1; m <= kmax; m++) {
      sqk = (unitky*l) * (unitky*l) + (unitkz*m) * (unitkz*m);
      if (sqk <= gsqmx) {
	kxvecs[kcount] = 0;
	kyvecs[kcount] = l;
	kzvecs[kcount] = m;
	ug[kcount] = preu*exp(-0.25*sqk*g_ewald_sq_inv)/sqk;
	eg[kcount][0] =  0.0;
	eg[kcount][1] =  2.0*unitky*l*ug[kcount];
	eg[kcount][2] =  2.0*unitkz*m*ug[kcount];
	vterm = -2.0*(1.0/sqk + 0.25*g_ewald_sq_inv);
	vg[kcount][0] = 1.0;
	vg[kcount][1] = 1.0 + vterm*(unitky*l)*(unitky*l);
	vg[kcount][2] = 1.0 + vterm*(unitkz*m)*(unitkz*m);
	vg[kcount][3] = 0.0;
	vg[kcount][4] = 0.0;
	vg[kcount][5] = vterm*unitky*l*unitkz*m;
	kcount++;

	kxvecs[kcount] = 0;
	kyvecs[kcount] = l;
	kzvecs[kcount] = -m;
	ug[kcount] = preu*exp(-0.25*sqk*g_ewald_sq_inv)/sqk;
	eg[kcount][0] =  0.0;
	eg[kcount][1] =  2.0*unitky*l*ug[kcount];
	eg[kcount][2] = -2.0*unitkz*m*ug[kcount];
	vg[kcount][0] = 1.0;
	vg[kcount][1] = 1.0 + vterm*(unitky*l)*(unitky*l);
	vg[kcount][2] = 1.0 + vterm*(unitkz*m)*(unitkz*m);
	vg[kcount][3] = 0.0;
	vg[kcount][4] = 0.0;
	vg[kcount][5] = -vterm*unitky*l*unitkz*m;
	kcount++;
      }
    }
  }

  // 1 = (k,0,m), 2 = (k,0,-m)

  for (k = 1; k <= kmax; k++) {
    for (m = 1; m <= kmax; m++) {
      sqk = (unitkx*k) * (unitkx*k) + (unitkz*m) * (unitkz*m);
      if (sqk <= gsqmx) {
	kxvecs[kcount] = k;
	kyvecs[kcount] = 0;
	kzvecs[kcount] = m;
	ug[kcount] = preu*exp(-0.25*sqk*g_ewald_sq_inv)/sqk;
	eg[kcount][0] =  2.0*unitkx*k*ug[kcount];
	eg[kcount][1] =  0.0;
	eg[kcount][2] =  2.0*unitkz*m*ug[kcount];
	vterm = -2.0*(1.0/sqk + 0.25*g_ewald_sq_inv);
	vg[kcount][0] = 1.0 + vterm*(unitkx*k)*(unitkx*k);
	vg[kcount][1] = 1.0;
	vg[kcount][2] = 1.0 + vterm*(unitkz*m)*(unitkz*m);
	vg[kcount][3] = 0.0;
	vg[kcount][4] = vterm*unitkx*k*unitkz*m;
	vg[kcount][5] = 0.0;
	kcount++;

	kxvecs[kcount] = k;
	kyvecs[kcount] = 0;
	kzvecs[kcount] = -m;
	ug[kcount] = preu*exp(-0.25*sqk*g_ewald_sq_inv)/sqk;
	eg[kcount][0] =  2.0*unitkx*k*ug[kcount];
	eg[kcount][1] =  0.0;
	eg[kcount][2] = -2.0*unitkz*m*ug[kcount];
	vg[kcount][0] = 1.0 + vterm*(unitkx*k)*(unitkx*k);
	vg[kcount][1] = 1.0;
	vg[kcount][2] = 1.0 + vterm*(unitkz*m)*(unitkz*m);
	vg[kcount][3] = 0.0;
	vg[kcount][4] = -vterm*unitkx*k*unitkz*m;
	vg[kcount][5] = 0.0;
	kcount++;
      }
    }
  }

  // 1 = (k,l,m), 2 = (k,-l,m), 3 = (k,l,-m), 4 = (k,-l,-m)

  for (k = 1; k <= kmax; k++) {
    for (l = 1; l <= kmax; l++) {
      for (m = 1; m <= kmax; m++) {
	sqk = (unitkx*k) * (unitkx*k) + (unitky*l) * (unitky*l) + 
	  (unitkz*m) * (unitkz*m);
	if (sqk <= gsqmx) {
	  kxvecs[kcount] = k;
	  kyvecs[kcount] = l;
	  kzvecs[kcount] = m;
	  ug[kcount] = preu*exp(-0.25*sqk*g_ewald_sq_inv)/sqk;
	  eg[kcount][0] = 2.0*unitkx*k*ug[kcount];
	  eg[kcount][1] = 2.0*unitky*l*ug[kcount];
	  eg[kcount][2] = 2.0*unitkz*m*ug[kcount];
	  vterm = -2.0*(1.0/sqk + 0.25*g_ewald_sq_inv);
	  vg[kcount][0] = 1.0 + vterm*(unitkx*k)*(unitkx*k);
	  vg[kcount][1] = 1.0 + vterm*(unitky*l)*(unitky*l);
	  vg[kcount][2] = 1.0 + vterm*(unitkz*m)*(unitkz*m);
	  vg[kcount][3] = vterm*unitkx*k*unitky*l;
	  vg[kcount][4] = vterm*unitkx*k*unitkz*m;
	  vg[kcount][5] = vterm*unitky*l*unitkz*m;
	  kcount++;

	  kxvecs[kcount] = k;
	  kyvecs[kcount] = -l;
	  kzvecs[kcount] = m;
	  ug[kcount] = preu*exp(-0.25*sqk*g_ewald_sq_inv)/sqk;
	  eg[kcount][0] = 2.0*unitkx*k*ug[kcount];
	  eg[kcount][1] = -2.0*unitky*l*ug[kcount];
	  eg[kcount][2] = 2.0*unitkz*m*ug[kcount];
	  vg[kcount][0] = 1.0 + vterm*(unitkx*k)*(unitkx*k);
	  vg[kcount][1] = 1.0 + vterm*(unitky*l)*(unitky*l);
	  vg[kcount][2] = 1.0 + vterm*(unitkz*m)*(unitkz*m);
	  vg[kcount][3] = -vterm*unitkx*k*unitky*l;
	  vg[kcount][4] = vterm*unitkx*k*unitkz*m;
	  vg[kcount][5] = -vterm*unitky*l*unitkz*m;
	  kcount++;

	  kxvecs[kcount] = k;
	  kyvecs[kcount] = l;
	  kzvecs[kcount] = -m;
	  ug[kcount] = preu*exp(-0.25*sqk*g_ewald_sq_inv)/sqk;
	  eg[kcount][0] = 2.0*unitkx*k*ug[kcount];
	  eg[kcount][1] = 2.0*unitky*l*ug[kcount];
	  eg[kcount][2] = -2.0*unitkz*m*ug[kcount];
	  vg[kcount][0] = 1.0 + vterm*(unitkx*k)*(unitkx*k);
	  vg[kcount][1] = 1.0 + vterm*(unitky*l)*(unitky*l);
	  vg[kcount][2] = 1.0 + vterm*(unitkz*m)*(unitkz*m);
	  vg[kcount][3] = vterm*unitkx*k*unitky*l;
	  vg[kcount][4] = -vterm*unitkx*k*unitkz*m;
	  vg[kcount][5] = -vterm*unitky*l*unitkz*m;
	  kcount++;

	  kxvecs[kcount] = k;
	  kyvecs[kcount] = -l;
	  kzvecs[kcount] = -m;
	  ug[kcount] = preu*exp(-0.25*sqk*g_ewald_sq_inv)/sqk;
	  eg[kcount][0] = 2.0*unitkx*k*ug[kcount];
	  eg[kcount][1] = -2.0*unitky*l*ug[kcount];
	  eg[kcount][2] = -2.0*unitkz*m*ug[kcount];
	  vg[kcount][0] = 1.0 + vterm*(unitkx*k)*(unitkx*k);
	  vg[kcount][1] = 1.0 + vterm*(unitky*l)*(unitky*l);
	  vg[kcount][2] = 1.0 + vterm*(unitkz*m)*(unitkz*m);
	  vg[kcount][3] = -vterm*unitkx*k*unitky*l;
	  vg[kcount][4] = -vterm*unitkx*k*unitkz*m;
	  vg[kcount][5] = vterm*unitky*l*unitkz*m;
	  kcount++;;
	}
      }
    }
  }
}

/* ----------------------------------------------------------------------
   allocate memory that depends on # of K-vectors 
------------------------------------------------------------------------- */

void Ewald::allocate()
{
  kxvecs = new int[kmax3d];
  kyvecs = new int[kmax3d];
  kzvecs = new int[kmax3d];

  ug = new double[kmax3d];
  eg = memory->create_2d_double_array(kmax3d,3,"ewald:eg");
  vg = memory->create_2d_double_array(kmax3d,6,"ewald:vg");

  sfacrl = new double[kmax3d];
  sfacim = new double[kmax3d];
  sfacrl_all = new double[kmax3d];
  sfacim_all = new double[kmax3d];
}

/* ----------------------------------------------------------------------
   deallocate memory that depends on # of K-vectors 
------------------------------------------------------------------------- */

void Ewald::deallocate()
{
  delete [] kxvecs;
  delete [] kyvecs;
  delete [] kzvecs;
  
  delete [] ug;
  memory->destroy_2d_double_array(eg);
  memory->destroy_2d_double_array(vg);

  delete [] sfacrl;
  delete [] sfacim;
  delete [] sfacrl_all;
  delete [] sfacim_all;
}

/* ----------------------------------------------------------------------
   Slab-geometry correction term to dampen inter-slab interactions between
   periodically repeating slabs.  Yields good approximation to 2-D Ewald if 
   adequate empty space is left between repeating slabs (J. Chem. Phys. 
   111, 3155).  Slabs defined here to be parallel to the xy plane. 
------------------------------------------------------------------------- */

void Ewald::slabcorr(int eflag)
{
  // compute local contribution to global dipole moment
  
  double *q = atom->q;
  double **x = atom->x;
  int nlocal = atom->nlocal;

  double dipole = 0.0;
  for (int i = 0; i < nlocal; i++) dipole += q[i]*x[i][2];
  
  // sum local contributions to get global dipole moment

  double dipole_all;
  MPI_Allreduce(&dipole,&dipole_all,1,MPI_DOUBLE,MPI_SUM,world);

  // compute corrections
  
  double e_slabcorr = 2.0*PI*dipole_all*dipole_all/volume;
  
  if (eflag) energy += qqrd2e*e_slabcorr;

  // add on force corrections

  double ffact = -4.0*PI*dipole_all/volume; 
  double **f = atom->f;

  for (int i = 0; i < nlocal; i++) f[i][2] += qqrd2e*q[i]*ffact;
}

/* ----------------------------------------------------------------------
   memory usage of local arrays 
------------------------------------------------------------------------- */

double Ewald::memory_usage()
{
  double bytes = 3 * kmax3d * sizeof(int);
  bytes += (1 + 3 + 6) * kmax3d * sizeof(double);
  bytes += 4 * kmax3d * sizeof(double);
  bytes += nmax*3 * sizeof(double);
  bytes += 2 * (2*kmax+1)*3*nmax * sizeof(double);
  return bytes;
}
