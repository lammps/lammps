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

#include "mpi.h"
#include "string.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "ewald_n.h"
#include "math_vector.h"
#include "math_const.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "pair.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define SMALL 0.00001

#define KSPACE_ILLEGAL	"Illegal kspace_style ewald/n command"
#define KSPACE_ORDER	"Unsupported order in kspace_style ewald/n for"
#define KSPACE_MIX	"Unsupported mixing rule in kspace_style ewald/n for"

enum{GEOMETRIC,ARITHMETIC,SIXTHPOWER};         // same as in pair.cpp

//#define DEBUG

/* ---------------------------------------------------------------------- */

EwaldN::EwaldN(LAMMPS *lmp, int narg, char **arg) : KSpace(lmp, narg, arg)
{
  if (narg!=1) error->all(FLERR,KSPACE_ILLEGAL);
  accuracy_relative = fabs(atof(arg[0]));
  memset(function, 0, EWALD_NORDER*sizeof(int));
  kenergy = kvirial = NULL;
  cek_local = cek_global = NULL;
  ekr_local = NULL;
  hvec = NULL;
  kvec = NULL;
  B = NULL;
  first_output = 0;
  energy_self_peratom = NULL;
  virial_self_peratom = NULL;
  nmax = 0;
  q2 = 0;
  b2 = 0;
}

EwaldN::~EwaldN()
{
  deallocate();
  deallocate_peratom();
  delete [] ekr_local;
  delete [] B;
}


/* --------------------------------------------------------------------- */

void EwaldN::init()
{
  nkvec = nkvec_max = nevec = nevec_max = 0;
  nfunctions = nsums = sums = 0;
  nbox = -1;
  bytes = 0.0;

  if (!comm->me) {					// output message
    if (screen) fprintf(screen,"EwaldN initialization ...\n");
    if (logfile) fprintf(logfile,"EwaldN initialization ...\n");
  }

  if (domain->dimension == 2)				// check for errors
    error->all(FLERR,"Cannot use EwaldN with 2d simulation");
  if (slabflag == 0 && domain->nonperiodic > 0)
    error->all(FLERR,"Cannot use nonperiodic boundaries with EwaldN");
  if (slabflag == 1) {
    if (domain->xperiodic != 1 || domain->yperiodic != 1 || 
	domain->boundary[2][0] != 1 || domain->boundary[2][1] != 1)
      error->all(FLERR,"Incorrect boundaries with slab EwaldN");
  }

  scale = 1.0;
  //mumurd2e = force->mumurd2e;
  //dielectric = force->dielectric;
  mumurd2e = dielectric = 1.0;
  
  int tmp;
  Pair *pair = force->pair;
  int *ptr = pair ? (int *) pair->extract("ewald_order",tmp) : NULL;
  double *cutoff = pair ? (double *) pair->extract("cut_coul",tmp) : NULL;
  if (!(ptr||cutoff)) 
    error->all(FLERR,"KSpace style is incompatible with Pair style");
  int ewald_order = ptr ? *((int *) ptr) : 1<<1;
  int ewald_mix = ptr ? *((int *) pair->extract("ewald_mix",tmp)) : GEOMETRIC;
  memset(function, 0, EWALD_NFUNCS*sizeof(int));
  for (int i=0; i<=EWALD_NORDER; ++i)			// transcribe order
    if (ewald_order&(1<<i)) {				// from pair_style
      int n[] = EWALD_NSUMS, k = 0;
      char str[128];
      switch (i) {
	case 1:
	  k = 0; break;
	case 3:
	  k = 3; break;
	case 6:
	  if (ewald_mix==GEOMETRIC) { k = 1; break; }
	  else if (ewald_mix==ARITHMETIC) { k = 2; break; }
	  sprintf(str, "%s pair_style %s", KSPACE_MIX, force->pair_style);
	  error->all(FLERR,str);
	default:
	  sprintf(str, "%s pair_style %s", KSPACE_ORDER, force->pair_style);
	  error->all(FLERR,str);
      }
      nfunctions += function[k] = 1;
      nsums += n[k];
    }
    
    
  g_ewald = 0;
  pair->init();  // so B is defined
  init_coeffs();
  init_coeff_sums();  
    
  double qsum, qsqsum, bsbsum;
  qsum = qsqsum = bsbsum = 0.0;
  if (function[0]) {
    qsum = sum[0].x;
    qsqsum = sum[0].x2;
  }
  if (function[1]) {
    bsbsum = sum[1].x2;
  } 
  if (function[2]) {
    bsbsum = sum[2].x2;
  }    
    
    
  if (qsqsum == 0.0 && bsbsum == 0.0)
      error->all(FLERR,"Cannot use Ewald/n solver on system with no charge or LJ particles");
  if (fabs(qsum) > SMALL && comm->me == 0) {
      char str[128];
      sprintf(str,"System is not charge neutral, net charge = %g",qsum);
      error->warning(FLERR,str);
  }
    
  //set accuracy (force units) from accuracy_relative or accuracy_absolute
    
  if (accuracy_absolute >= 0.0) accuracy = accuracy_absolute;
  else accuracy = accuracy_relative * two_charge_force;
    
  // setup K-space resolution
    
  q2 = qsqsum * force->qqrd2e / force->dielectric;
  b2 = bsbsum; //Are these units right?
  bigint natoms = atom->natoms;
    
  if (function[0]) {       //Coulombic
    g_ewald = accuracy*sqrt(natoms*(*cutoff)*shape_det(domain->h)) / (2.0*q2);
    if (g_ewald >= 1.0)
        error->all(FLERR,"KSpace accuracy too large to estimate G vector");
    g_ewald = sqrt(-log(g_ewald)) / *cutoff;
  }
  else if (function[1] || function[2]) {    //Only LJ

    double *cutoffLJ = pair ? (double *) pair->extract("cut_LJ",tmp) : NULL;
    
    //Try Newton Solver
    
    //Use old method to get guess
    g_ewald = (1.35 - 0.15*log(accuracy))/ *cutoffLJ;
      
    double g_ewald_new = NewtonSolve(g_ewald,(*cutoffLJ),natoms,shape_det(domain->h),b2);
           
    if (g_ewald_new > 0.0) g_ewald = g_ewald_new;
    else error->warning(FLERR,"Ewald/n Newton solver failed, using old method to estimate g_ewald");  
    
    if (g_ewald >= 1.0)
        error->all(FLERR,"KSpace accuracy too large to estimate G vector");
  }
    
  if (!comm->me) {					// output results
      if (screen) fprintf(screen, "  G vector = %g\n", g_ewald);
      if (logfile) fprintf(logfile, "  G vector = %g\n", g_ewald);
  }  
    
    
  deallocate_peratom();
  peratom_allocate_flag = 0;
}


/* ----------------------------------------------------------------------
   adjust EwaldN coeffs, called initially and whenever volume has changed 
------------------------------------------------------------------------- */

void EwaldN::setup()
{
  volume = shape_det(domain->h)*slab_volfactor;		// cell volume
  memcpy(unit, domain->h_inv, sizeof(shape));		// wave vector units
  shape_scalar_mult(unit, 2.0*MY_PI);
  unit[2] /= slab_volfactor;
    
  //int nbox_old = nbox, nkvec_old = nkvec;
  if (accuracy>=1) {
    nbox = 0;
    error->all(FLERR,"KSpace accuracy too low"); 
  }
    
  bigint natoms = atom->natoms;
  double err;
  int kxmax = 1;
  int kymax = 1;
  int kzmax = 1;
  err = rms(kxmax,domain->h[0],natoms,q2,b2);
  while (err > accuracy) {
    kxmax++;
    err = rms(kxmax,domain->h[0],natoms,q2,b2);
  }
  err = rms(kymax,domain->h[1],natoms,q2,b2);
  while (err > accuracy) {
    kymax++;
    err = rms(kymax,domain->h[1],natoms,q2,b2);
  }
  err = rms(kzmax,domain->h[2]*slab_volfactor,natoms,q2,b2);
  while (err > accuracy) {
    kzmax++;
    err = rms(kzmax,domain->h[2]*slab_volfactor,natoms,q2,b2);
  } 
  nbox = MAX(kxmax,kymax);
  nbox = MAX(nbox,kzmax);
  double gsqxmx = unit[0]*unit[0]*kxmax*kxmax;
  double gsqymx = unit[1]*unit[1]*kymax*kymax;
  double gsqzmx = unit[2]*unit[2]*kzmax*kzmax;
  gsqmx = MAX(gsqxmx,gsqymx);
  gsqmx = MAX(gsqmx,gsqzmx);
  gsqmx *= 1.00001;
    
  reallocate();
  coefficients();					// compute coeffs
  init_coeffs();
  init_coeff_sums();
  init_self();

  if (!(first_output||comm->me)) {			// output on first
    first_output = 1;
    if (screen) fprintf(screen,
       	"  vectors: nbox = %d, nkvec = %d\n", nbox, nkvec);
    if (logfile) fprintf(logfile,
	"  vectors: nbox = %d, nkvec = %d\n", nbox, nkvec);
  }  
}

/* ----------------------------------------------------------------------
   compute RMS accuracy for a dimension
------------------------------------------------------------------------- */

double EwaldN::rms(int km, double prd, bigint natoms, double q2, double b2)
{  
  double value = 0.0;
    
  //Coulombic 
    
  double g2 = g_ewald*g_ewald;
    
  value += 2.0*q2*g_ewald/prd * 
    sqrt(1.0/(MY_PI*km*natoms)) * 
    exp(-MY_PI*MY_PI*km*km/(g2*prd*prd));
    
  //Lennard-Jones  
    
  double g7 = g2*g2*g2*g_ewald;
    
  value += 4.0*b2*g7/3.0 * 
    sqrt(1.0/(MY_PI*natoms)) * 
    (exp(-MY_PI*MY_PI*km*km/(g2*prd*prd)) *
    (MY_PI*km/(g_ewald*prd) + 1));
    
  return value;
}

void EwaldN::reallocate()				// allocate memory
{
  int ix, iy, iz;
  int nkvec_max = nkvec;
  vector h;

  nkvec = 0;						// determine size(kvec)
  int *kflag = new int[(nbox+1)*(2*nbox+1)*(2*nbox+1)];
  int *flag = kflag;

  for (ix=0; ix<=nbox; ++ix)
    for (iy=-nbox; iy<=nbox; ++iy)
      for (iz=-nbox; iz<=nbox; ++iz)
	if (!(ix||iy||iz)) *(flag++) = 0;
	else if ((!ix)&&(iy<0)) *(flag++) = 0;
	else if ((!(ix||iy))&&(iz<0)) *(flag++) = 0;	// use symmetry
	else {
	  h[0] = unit[0]*ix;
	  h[1] = unit[5]*ix+unit[1]*iy;
	  h[2] = unit[4]*ix+unit[3]*iy+unit[2]*iz;
	  if ((*(flag++) = h[0]*h[0]+h[1]*h[1]+h[2]*h[2]<=gsqmx)) ++nkvec;
	}
  
  if (nkvec>nkvec_max) {
    deallocate();					// free memory
    hvec = new hvector[nkvec];				// hvec
    bytes += (nkvec-nkvec_max)*sizeof(hvector);
    kvec = new kvector[nkvec];				// kvec
    bytes += (nkvec-nkvec_max)*sizeof(kvector);
    kenergy = new double[nkvec*nfunctions];		// kenergy
    bytes += (nkvec-nkvec_max)*nfunctions*sizeof(double);
    kvirial = new double[6*nkvec*nfunctions];		// kvirial
    bytes += 6*(nkvec-nkvec_max)*nfunctions*sizeof(double);
    cek_local = new complex[nkvec*nsums];		// cek_local
    bytes += (nkvec-nkvec_max)*nsums*sizeof(complex);
    cek_global = new complex[nkvec*nsums];		// cek_global
    bytes += (nkvec-nkvec_max)*nsums*sizeof(complex);
    nkvec_max = nkvec;
  }

  flag = kflag;						// create index and
  kvector *k = kvec;					// wave vectors
  hvector *hi = hvec;
  for (ix=0; ix<=nbox; ++ix)
    for (iy=-nbox; iy<=nbox; ++iy)
      for (iz=-nbox; iz<=nbox; ++iz)
	if (*(flag++)) {
	  hi->x = unit[0]*ix;
	  hi->y = unit[5]*ix+unit[1]*iy;
	  (hi++)->z = unit[4]*ix+unit[3]*iy+unit[2]*iz;
	  k->x = ix+nbox; k->y = iy+nbox; (k++)->z = iz+nbox; }

  delete [] kflag;
}


void EwaldN::reallocate_atoms()
{
  if (eflag_atom || vflag_atom)
    if (atom->nlocal > nmax) {
      deallocate_peratom();
      allocate_peratom();
      nmax = atom->nmax;
    }
    
  if ((nevec = atom->nmax*(2*nbox+1))<=nevec_max) return;
  delete [] ekr_local;
  ekr_local = new cvector[nevec];
  bytes += (nevec-nevec_max)*sizeof(cvector);
  nevec_max = nevec;
}


void EwaldN::allocate_peratom()				
{
  memory->create(energy_self_peratom,
      atom->nmax,EWALD_NFUNCS,"ewald/n:energy_self_peratom");
  memory->create(virial_self_peratom,
      atom->nmax,EWALD_NFUNCS,"ewald/n:virial_self_peratom");
}


void EwaldN::deallocate_peratom()			// free memory
{
  memory->destroy(energy_self_peratom);
  memory->destroy(virial_self_peratom);
}


void EwaldN::deallocate()				// free memory
{
  delete [] hvec;		hvec = NULL;
  delete [] kvec;		kvec = NULL;
  delete [] kenergy;		kenergy = NULL;
  delete [] kvirial;		kvirial = NULL;
  delete [] cek_local;		cek_local = NULL;
  delete [] cek_global;		cek_global = NULL;
}


void EwaldN::coefficients()				// set up pre-factors
{
  vector h;
  hvector *hi = hvec, *nh;
  double eta2 = 0.25/(g_ewald*g_ewald);
  double b1, b2, expb2, h1, h2, c1, c2;
  double *ke = kenergy, *kv = kvirial;
  int func0 = function[0], func12 = function[1]||function[2],
      func3 = function[3];

  for (nh = (hi = hvec)+nkvec; hi<nh; ++hi) {		// wave vectors
    memcpy(h, hi, sizeof(vector));
    expb2 = exp(-(b2 = (h2 = vec_dot(h, h))*eta2));
    if (func0) {					// qi*qj/r coeffs
      *(ke++) = c1 = expb2/h2;
      *(kv++) = c1-(c2 = 2.0*c1*(1.0+b2)/h2)*h[0]*h[0];
      *(kv++) = c1-c2*h[1]*h[1];			// lammps convention
      *(kv++) = c1-c2*h[2]*h[2];			// instead of voigt
      *(kv++) = -c2*h[1]*h[0];				
      *(kv++) = -c2*h[2]*h[0];
      *(kv++) = -c2*h[2]*h[1];
    }
    if (func12) {					// -Bij/r^6 coeffs
      b1 = sqrt(b2);					// minus sign folded
      h1 = sqrt(h2);					// into constants
      *(ke++) = c1 = -h1*h2*((c2=sqrt(MY_PI)*erfc(b1))+(0.5/b2-1.0)*expb2/b1);
      *(kv++) = c1-(c2 = 3.0*h1*(c2-expb2/b1))*h[0]*h[0];
      *(kv++) = c1-c2*h[1]*h[1];			// lammps convention
      *(kv++) = c1-c2*h[2]*h[2];			// instead of voigt
      *(kv++) = -c2*h[1]*h[0];
      *(kv++) = -c2*h[2]*h[0];
      *(kv++) = -c2*h[2]*h[1];
    }
    if (func3) {					// dipole coeffs
      *(ke++) = c1 = expb2/h2;
      *(kv++) = c1-(c2 = 2.0*c1*(1.0+b2)/h2)*h[0]*h[0];
      *(kv++) = c1-c2*h[1]*h[1];			// lammps convention
      *(kv++) = c1-c2*h[2]*h[2];			// instead of voigt
      *(kv++) = -c2*h[1]*h[0];				
      *(kv++) = -c2*h[2]*h[0];
      *(kv++) = -c2*h[2]*h[1];
    }
  }
}


void EwaldN::init_coeffs()				// local pair coeffs
{
  int tmp;
  int n = atom->ntypes;

  if (function[1]) {					// geometric 1/r^6
    double **b = (double **) force->pair->extract("B",tmp);
    delete [] B;
    B = new double[n+1];
    bytes += (n+1)*sizeof(double);
    for (int i=0; i<=n; ++i) B[i] = sqrt(fabs(b[i][i]));
  }
  if (function[2]) {					// arithmetic 1/r^6
    double **epsilon = (double **) force->pair->extract("epsilon",tmp);
    double **sigma = (double **) force->pair->extract("sigma",tmp);
    double eps_i, sigma_i, sigma_n, *bi = B = new double[7*n+7];
    double c[7] = {
      1.0, sqrt(6.0), sqrt(15.0), sqrt(20.0), sqrt(15.0), sqrt(6.0), 1.0};

    if (!(epsilon&&sigma))
      error->all(
	  FLERR,"epsilon or sigma reference not set by pair style in ewald/n");
    for (int i=0; i<=n; ++i) {
      eps_i = sqrt(epsilon[i][i]);
      sigma_i = sigma[i][i];
      sigma_n = 1.0;
      for (int j=0; j<7; ++j) {
	*(bi++) = sigma_n*eps_i*c[j]; sigma_n *= sigma_i;
      }
    }
  }
}


void EwaldN::init_coeff_sums()				// sums based on atoms
{
  if (sums) return;					// calculated only once
  sums = 1;

  Sum sum_local[EWALD_MAX_NSUMS];

  memset(sum_local, 0, EWALD_MAX_NSUMS*sizeof(Sum));
  if (function[0]) {					// 1/r
    double *q = atom->q, *qn = q+atom->nlocal;
    for (double *i=q; i<qn; ++i) {
      sum_local[0].x += i[0]; sum_local[0].x2 += i[0]*i[0]; }
  }
  if (function[1]) {					// geometric 1/r^6
    int *type = atom->type, *ntype = type+atom->nlocal;
    for (int *i=type; i<ntype; ++i) {
      sum_local[1].x += B[i[0]]; sum_local[1].x2 += B[i[0]]*B[i[0]]; }
  }
  if (function[2]) {					// arithmetic 1/r^6
    double *bi;
    int *type = atom->type, *ntype = type+atom->nlocal;
    for (int *i=type; i<ntype; ++i) {
      bi = B+7*i[0];
      sum_local[2].x2 += bi[0]*bi[6];
      for (int k=2; k<9; ++k) sum_local[k].x += *(bi++); 
    }
  }
  if (function[3]&&atom->mu) {				// dipole
    double *mu = atom->mu[0], *nmu = mu+4*atom->nlocal;
    for (double *i = mu; i < nmu; i += 4)
      sum_local[9].x2 += i[3]*i[3];
  }
  MPI_Allreduce(sum_local, sum, 2*EWALD_MAX_NSUMS, MPI_DOUBLE, MPI_SUM, world);
}


void EwaldN::init_self()
{
  double g1 = g_ewald, g2 = g1*g1, g3 = g1*g2;
  const double qscale = force->qqrd2e * scale;

  memset(energy_self, 0, EWALD_NFUNCS*sizeof(double));	// self energy
  memset(virial_self, 0, EWALD_NFUNCS*sizeof(double));

  if (function[0]) {					// 1/r
    virial_self[0] = -0.5*MY_PI*qscale/(g2*volume)*sum[0].x*sum[0].x;
    energy_self[0] = sum[0].x2*qscale*g1/sqrt(MY_PI)-virial_self[0];
  }
  if (function[1]) {					// geometric 1/r^6
    virial_self[1] = MY_PI*sqrt(MY_PI)*g3/(6.0*volume)*sum[1].x*sum[1].x;
    energy_self[1] = -sum[1].x2*g3*g3/12.0+virial_self[1];
  }
  if (function[2]) {					// arithmetic 1/r^6
    virial_self[2] = MY_PI*sqrt(MY_PI)*g3/(48.0*volume)*(sum[2].x*sum[8].x+
	sum[3].x*sum[7].x+sum[4].x*sum[6].x+0.5*sum[5].x*sum[5].x);
    energy_self[2] = -sum[2].x2*g3*g3/3.0+virial_self[2];
  }
  if (function[3]) {					// dipole
    virial_self[3] = 0;					// in surface
    energy_self[3] = sum[9].x2*mumurd2e*2.0*g3/3.0/sqrt(MY_PI)-virial_self[3];
  }
}


void EwaldN::init_self_peratom()
{
  if (!(vflag_atom || eflag_atom)) return; 

  double g1 = g_ewald, g2 = g1*g1, g3 = g1*g2;
  const double qscale = force->qqrd2e * scale;
  double *energy = energy_self_peratom[0];
  double *virial = virial_self_peratom[0];
  int nlocal = atom->nlocal;
    
  memset(energy, 0, EWALD_NFUNCS*nlocal*sizeof(double));
  memset(virial, 0, EWALD_NFUNCS*nlocal*sizeof(double));

  if (function[0]) {					// 1/r
    double *ei = energy;
    double *vi = virial;
    double ce = qscale*g1/sqrt(MY_PI);
    double cv = -0.5*MY_PI*qscale/(g2*volume);
    double *qi = atom->q, *qn = qi + nlocal;
    for (; qi < qn; qi++, vi += EWALD_NFUNCS, ei += EWALD_NFUNCS) {
      double q = *qi;
      *vi = cv*q*sum[0].x;
      *ei = ce*q*q-vi[0];
    }
  }
  if (function[1]) {					// geometric 1/r^6
    double *ei = energy+1;
    double *vi = virial+1;
    double ce = -g3*g3/12.0;
    double cv = MY_PI*sqrt(MY_PI)*g3/(6.0*volume);
    int *typei = atom->type, *typen = typei + atom->nlocal;
    for (; typei < typen; typei++, vi += EWALD_NFUNCS, ei += EWALD_NFUNCS) {
      double b = B[*typei];
      *vi = cv*b*sum[1].x;
      *ei = ce*b*b+vi[0];
    }
  }
  if (function[2]) {					// arithmetic 1/r^6
    double *bi;
    double *ei = energy+2;
    double *vi = virial+2;
    double ce = -g3*g3/3.0;
    double cv = 0.5*MY_PI*sqrt(MY_PI)*g3/(48.0*volume);
    int *typei = atom->type, *typen = typei + atom->nlocal;
    for (; typei < typen; typei++, vi += EWALD_NFUNCS, ei += EWALD_NFUNCS) {
      bi = B+7*typei[0]+7;       
      for (int k=2; k<9; ++k) *vi += cv*sum[k].x*(--bi)[0];

      /* PJV 20120225: 
         should this be this instead?  above implies an inverse dependence
	 seems to be the above way in original;  i recall having tested
	 arithmetic mixing in the conception phase, but an extra test would
	 be prudent (pattern repeats in multiple functions below)

      bi = B+7*typei[0];       
      for (int k=2; k<9; ++k) *vi += cv*sum[k].x*(bi++)[0];
      
      */

      *ei = ce*bi[0]*bi[6]+vi[0];
    }
  }
  if (function[3]&&atom->mu) {				// dipole
    double *ei = energy+3;
    double *vi = virial+3;
    double *imu = atom->mu[0], *nmu = imu+4*atom->nlocal;
    double ce = mumurd2e*2.0*g3/3.0/sqrt(MY_PI);
    for (; imu < nmu; imu += 4, vi += EWALD_NFUNCS, ei += EWALD_NFUNCS) {
      *vi = 0;						// in surface
      *ei = ce*imu[3]*imu[3]-vi[0];
    }
  }
}


/* ----------------------------------------------------------------------
   compute the EwaldN long-range force, energy, virial 
------------------------------------------------------------------------- */

void EwaldN::compute(int eflag, int vflag)
{
  if (!nbox) return;
    
  // set energy/virial flags
  // invoke allocate_peratom() if needed for first time
    
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = eflag_global = vflag_global = eflag_atom = vflag_atom = 0;  
    
  if (!peratom_allocate_flag && (eflag_atom || vflag_atom)) {
      allocate_peratom();
      peratom_allocate_flag = 1;
      nmax = atom->nmax;
  }
    
  reallocate_atoms();
  init_self_peratom();
  compute_ek();
  compute_force();
  compute_surface();
  compute_energy();
  compute_energy_peratom();
  compute_virial();
  compute_virial_peratom();
}


void EwaldN::compute_ek()
{
  cvector *ekr = ekr_local;
  int lbytes = (2*nbox+1)*sizeof(cvector);
  hvector *h = NULL;
  kvector *k, *nk = kvec+nkvec;
  cvector *z = new cvector[2*nbox+1];
  cvector z1, *zx, *zy, *zz, *zn = z+2*nbox;
  complex *cek, zxyz, zxy = COMPLEX_NULL, cx = COMPLEX_NULL;
  vector mui;
  double *x = atom->x[0], *xn = x+3*atom->nlocal, *q = atom->q, qi = 0.0;
  double bi = 0.0, ci[7];
  double *mu = atom->mu ? atom->mu[0] : NULL;
  int i, kx, ky, n = nkvec*nsums, *type = atom->type, tri = domain->triclinic;
  int func[EWALD_NFUNCS];

  memcpy(func, function, EWALD_NFUNCS*sizeof(int));
  memset(cek_local, 0, n*sizeof(complex));		// reset sums
  while (x<xn) {
    zx = (zy = (zz = z+nbox)+1)-2;
    C_SET(zz->x, 1, 0); C_SET(zz->y, 1, 0); C_SET(zz->z, 1, 0);	// z[0]
    if (tri) {						// triclinic z[1]
      C_ANGLE(z1.x, unit[0]*x[0]+unit[5]*x[1]+unit[4]*x[2]);
      C_ANGLE(z1.y, unit[1]*x[1]+unit[3]*x[2]);
      C_ANGLE(z1.z, x[2]*unit[2]); x += 3;
    }
    else {						// orthogonal z[1]
      C_ANGLE(z1.x, *(x++)*unit[0]);
      C_ANGLE(z1.y, *(x++)*unit[1]);
      C_ANGLE(z1.z, *(x++)*unit[2]);
    }
    for (; zz<zn; --zx, ++zy, ++zz) {			// set up z[k]=e^(ik.r)
      C_RMULT(zy->x, zz->x, z1.x);			// 3D k-vector
      C_RMULT(zy->y, zz->y, z1.y); C_CONJ(zx->y, zy->y);
      C_RMULT(zy->z, zz->z, z1.z); C_CONJ(zx->z, zy->z);
    }
    kx = ky = -1;
    cek = cek_local;
    if (func[0]) qi = *(q++);
    if (func[1]) bi = B[*type];
    if (func[2]) memcpy(ci, B+7*type[0], 7*sizeof(double));
    if (func[3]) { 
      memcpy(mui, mu, sizeof(vector));
      vec_scalar_mult(mui, mu[3]);
      mu += 4;
      h = hvec;
    }
    for (k=kvec; k<nk; ++k) {				// compute rho(k)
      if (ky!=k->y) {					// based on order in 
	if (kx!=k->x) cx = z[kx = k->x].x;		// reallocate
	C_RMULT(zxy, z[ky = k->y].y, cx);
      }
      C_RMULT(zxyz, z[k->z].z, zxy);
      if (func[0]) {
       	cek->re += zxyz.re*qi; (cek++)->im += zxyz.im*qi;
      }
      if (func[1]) {
       	cek->re += zxyz.re*bi; (cek++)->im += zxyz.im*bi;
      }
      if (func[2]) for (i=0; i<7; ++i) {
	cek->re += zxyz.re*ci[i]; (cek++)->im += zxyz.im*ci[i];
      }
      if (func[3]) {
	register double muk = mui[0]*h->x+mui[1]*h->y+mui[2]*h->z; ++h;
	cek->re += zxyz.re*muk; (cek++)->im += zxyz.im*muk;
      }
    }
    ekr = (cvector *) ((char *) memcpy(ekr, z, lbytes)+lbytes);
    ++type;
  }
  MPI_Allreduce(cek_local, cek_global, 2*n, MPI_DOUBLE, MPI_SUM, world);

  delete [] z;
}


void EwaldN::compute_force()
{
  kvector *k;
  hvector *h, *nh;
  cvector *z = ekr_local;
  vector sum[EWALD_MAX_NSUMS], mui = COMPLEX_NULL;
  complex *cek, zc, zx = COMPLEX_NULL, zxy = COMPLEX_NULL;
  double *f = atom->f[0], *fn = f+3*atom->nlocal, *q = atom->q, *t = NULL;
  double *mu = atom->mu ? atom->mu[0] : NULL;
  const double qscale = force->qqrd2e * scale;
  double *ke, c[EWALD_NFUNCS] = {
    8.0*MY_PI*qscale/volume, 2.0*MY_PI*sqrt(MY_PI)/(12.0*volume),
    2.0*MY_PI*sqrt(MY_PI)/(192.0*volume), 8.0*MY_PI*mumurd2e/volume};
  double kt = 4.0*pow(g_ewald, 3.0)/3.0/sqrt(MY_PI)/c[3];
  int i, kx, ky, lbytes = (2*nbox+1)*sizeof(cvector), *type = atom->type;
  int func[EWALD_NFUNCS];

  if (atom->torque) t = atom->torque[0];
  memcpy(func, function, EWALD_NFUNCS*sizeof(int));
  memset(sum, 0, EWALD_MAX_NSUMS*sizeof(vector));	// fj = -dE/dr =
  for (; f<fn; f+=3) {					//      -i*qj*fac*
    k = kvec;						//       Sum[conj(d)-d]
    kx = ky = -1;					// d = k*conj(ekj)*ek
    ke = kenergy;
    cek = cek_global;
    memset(sum, 0, EWALD_MAX_NSUMS*sizeof(vector));
    if (func[3]) { 
      register double di = mu[3] * c[3];
      mui[0] = di*(mu++)[0]; mui[1] = di*(mu++)[0]; mui[2] = di*(mu++)[0];
      mu++;
    }
    for (nh = (h = hvec)+nkvec; h<nh; ++h, ++k) {
      if (ky!=k->y) {					// based on order in
	if (kx!=k->x) zx = z[kx = k->x].x; 		// reallocate
	C_RMULT(zxy, z[ky = k->y].y, zx);
      }
      C_CRMULT(zc, z[k->z].z, zxy);
      if (func[0]) {					// 1/r
	register double im = *(ke++)*(zc.im*cek->re+cek->im*zc.re); ++cek;
	sum[0][0] += h->x*im; sum[0][1] += h->y*im; sum[0][2] += h->z*im;
      }
      if (func[1]) {					// geometric 1/r^6
	register double im = *(ke++)*(zc.im*cek->re+cek->im*zc.re); ++cek;
	sum[1][0] += h->x*im; sum[1][1] += h->y*im; sum[1][2] += h->z*im;
      }
      if (func[2]) {					// arithmetic 1/r^6
	register double im, c = *(ke++);
	for (i=2; i<9; ++i) {
	  im = c*(zc.im*cek->re+cek->im*zc.re); ++cek;
	  sum[i][0] += h->x*im; sum[i][1] += h->y*im; sum[i][2] += h->z*im;
	}
      }
      if (func[3]) {					// dipole
	register double im = *(ke++)*(zc.im*cek->re+
	    cek->im*zc.re)*(mui[0]*h->x+mui[1]*h->y+mui[2]*h->z); ++cek;
	sum[9][0] += h->x*im; sum[9][1] += h->y*im; sum[9][2] += h->z*im;
      }
    }
    if (func[0]) {					// 1/r
      register double qi = *(q++)*c[0];
      f[0] -= sum[0][0]*qi; f[1] -= sum[0][1]*qi; f[2] -= sum[0][2]*qi;
    }
    if (func[1]) {					// geometric 1/r^6
      register double bi = B[*type]*c[1];
      f[0] -= sum[1][0]*bi; f[1] -= sum[1][1]*bi; f[2] -= sum[1][2]*bi;
    }
    if (func[2]) {					// arithmetic 1/r^6
      register double *bi = B+7*type[0]+7;
      for (i=2; i<9; ++i) {
	register double c2 = (--bi)[0]*c[2];
	f[0] -= sum[i][0]*c2; f[1] -= sum[i][1]*c2; f[2] -= sum[i][2]*c2;
      }
    }
    if (func[3]) {					// dipole
      f[0] -= sum[9][0]; f[1] -= sum[9][1]; f[2] -= sum[9][2];
      *(t++) -= mui[1]*sum[0][2]+mui[2]*sum[0][1]-mui[0]*kt;	// torque
      *(t++) -= mui[2]*sum[0][0]+mui[0]*sum[0][2]-mui[1]*kt;
      *(t++) -= mui[0]*sum[0][1]+mui[1]*sum[0][0]-mui[2]*kt;
    }
    z = (cvector *) ((char *) z+lbytes);
    ++type;
  }
}


void EwaldN::compute_surface()
{
  if (!function[3]) return;
  if (!atom->mu) return;

  vector sum_local = VECTOR_NULL, sum_total;
  memset(sum_local, 0, sizeof(vector));
  double *i, *n, *mu = atom->mu[0];

  for (n = (i = mu) + 4*atom->nlocal; i < n; ++i) {
    register double di = i[3];
    sum_local[0] += di*(i++)[0];
    sum_local[1] += di*(i++)[0];
    sum_local[2] += di*(i++)[0];
  }
  MPI_Allreduce(sum_local, sum_total, 3, MPI_DOUBLE, MPI_SUM, world);

  energy_self[3] += virial_self[3];
  virial_self[3] =
    mumurd2e*(2.0*MY_PI*vec_dot(sum_total,sum_total)/(2.0*dielectric+1)/volume);
  energy_self[3] -= virial_self[3];

  if (!(vflag_atom || eflag_atom)) return;
 
  double *ei = energy_self_peratom[0]+3;
  double *vi = virial_self_peratom[0]+3;
  double cv = 2.0*mumurd2e*MY_PI/(2.0*dielectric+1)/volume;

  for (i = mu; i < n; i += 4, ei += EWALD_NFUNCS, vi += EWALD_NFUNCS) {
    *ei += *vi;
    *vi = cv*i[3]*(i[0]*sum_total[0]+i[1]*sum_total[1]+i[2]*sum_total[2]);
    *ei -= *vi;
  }
}


void EwaldN::compute_energy()
{
  energy = 0.0;
  if (!eflag_global) return;
  
  complex *cek = cek_global;
  double *ke = kenergy;
  const double qscale = force->qqrd2e * scale;
  double c[EWALD_NFUNCS] = {
    4.0*MY_PI*qscale/volume, 2.0*MY_PI*sqrt(MY_PI)/(24.0*volume),
    2.0*MY_PI*sqrt(MY_PI)/(192.0*volume), 4.0*MY_PI*mumurd2e/volume};
  double sum[EWALD_NFUNCS];
  int func[EWALD_NFUNCS];

  memcpy(func, function, EWALD_NFUNCS*sizeof(int));
  memset(sum, 0, EWALD_NFUNCS*sizeof(double));		// reset sums
  for (int k=0; k<nkvec; ++k) {				// sum over k vectors
    if (func[0]) {					// 1/r
      sum[0] += *(ke++)*(cek->re*cek->re+cek->im*cek->im); ++cek; }
    if (func[1]) {					// geometric 1/r^6
      sum[1] += *(ke++)*(cek->re*cek->re+cek->im*cek->im); ++cek; }
    if (func[2]) {					// arithmetic 1/r^6
      register double r =
	    (cek[0].re*cek[6].re+cek[0].im*cek[6].im)+
	    (cek[1].re*cek[5].re+cek[1].im*cek[5].im)+
	    (cek[2].re*cek[4].re+cek[2].im*cek[4].im)+
	0.5*(cek[3].re*cek[3].re+cek[3].im*cek[3].im); cek += 7;
      sum[2] += *(ke++)*r;
    }
    if (func[3]) {					// dipole
      sum[3] += *(ke++)*(cek->re*cek->re+cek->im*cek->im); ++cek; }
  }  
  for (int k=0; k<EWALD_NFUNCS; ++k) energy += c[k]*sum[k]-energy_self[k];
  if (slabflag) compute_slabcorr();
}


void EwaldN::compute_energy_peratom()
{
  if (!eflag_atom) return;
      
  kvector *k;
  hvector *h, *nh;
  cvector *z = ekr_local;
  vector  mui = VECTOR_NULL;
  double sum[EWALD_MAX_NSUMS];
  complex *cek, zc = COMPLEX_NULL, zx = COMPLEX_NULL, zxy = COMPLEX_NULL;
  double *q = atom->q;
  double *eatomj = eatom;
  double *mu = atom->mu ? atom->mu[0] : NULL;
  const double qscale = force->qqrd2e * scale;
  double *ke = kenergy;
  double c[EWALD_NFUNCS] = {
      4.0*MY_PI*qscale/volume, 2.0*MY_PI*sqrt(MY_PI)/(24.0*volume),
      2.0*MY_PI*sqrt(MY_PI)/(192.0*volume), 4.0*MY_PI*mumurd2e/volume};
  int i, kx, ky, lbytes = (2*nbox+1)*sizeof(cvector), *type = atom->type;
  int func[EWALD_NFUNCS];
  
  memcpy(func, function, EWALD_NFUNCS*sizeof(int));
  for (int j = 0; j < atom->nlocal; j++, ++eatomj) {
    k = kvec;						
    kx = ky = -1;					
    ke = kenergy;
    cek = cek_global;
    memset(sum, 0, EWALD_MAX_NSUMS*sizeof(double));
    if (func[3]) { 
      register double di = mu[3] * c[3];
      mui[0] = di*(mu++)[0]; mui[1] = di*(mu++)[0]; mui[2] = di*(mu++)[0];
      mu++;
    }
    for (nh = (h = hvec)+nkvec; h<nh; ++h, ++k) {
      if (ky!=k->y) {					// based on order in
	if (kx!=k->x) zx = z[kx = k->x].x; 		// reallocate
	C_RMULT(zxy, z[ky = k->y].y, zx);
      }
      C_CRMULT(zc, z[k->z].z, zxy);
      if (func[0]) {					// 1/r
	sum[0] += *(ke++)*(cek->re*zc.re - cek->im*zc.im); ++cek; }
      if (func[1]) {					// geometric 1/r^6
	sum[1] += *(ke++)*(cek->re*zc.re - cek->im*zc.im); ++cek; }
      if (func[2]) {					// arithmetic 1/r^6
	register double im, c = *(ke++);
	for (i=2; i<9; ++i) {
	  im = c*(cek->re*zc.re - cek->im*zc.im); ++cek;
	  sum[i] += im;
	}
      }
      if (func[3]) {					// dipole
	sum[9] += *(ke++)*(cek->re*zc.re +
		  cek->im*zc.im)*(mui[0]*h->x+mui[1]*h->y+mui[2]*h->z); ++cek;
      }
    }
    
    if (func[0]) {					// 1/r
      register double qj = *(q++)*c[0];
      *eatomj += sum[0]*qj - energy_self_peratom[j][0];
    }
    if (func[1]) {					// geometric 1/r^6
      register double bj = B[*type]*c[1];
      *eatomj += sum[1]*bj - energy_self_peratom[j][1];
    }
    if (func[2]) {					// arithmetic 1/r^6
      register double *bj = B+7*type[0]+7;
      for (i=2; i<9; ++i) {
	register double c2 = (--bj)[0]*c[2];
	*eatomj += 0.5*sum[i]*c2;
      }  
      *eatomj -= energy_self_peratom[j][2];           
    }
    if (func[3]) {					// dipole
      *eatomj += sum[9] - energy_self_peratom[j][3];
    }
    z = (cvector *) ((char *) z+lbytes);
    ++type;   
  }
}


#define swap(a, b) { register double t = a; a= b; b = t; }

void EwaldN::compute_virial()
{
  memset(virial, 0, sizeof(shape));
  if (!vflag_global) return;

  complex *cek = cek_global;
  double *kv = kvirial;
  const double qscale = force->qqrd2e * scale;
  double c[EWALD_NFUNCS] = {
    4.0*MY_PI*qscale/volume, 2.0*MY_PI*sqrt(MY_PI)/(24.0*volume),
    2.0*MY_PI*sqrt(MY_PI)/(192.0*volume), 4.0*MY_PI*mumurd2e/volume};
  shape sum[EWALD_NFUNCS];
  int func[EWALD_NFUNCS];

  memcpy(func, function, EWALD_NFUNCS*sizeof(int));
  memset(sum, 0, EWALD_NFUNCS*sizeof(shape));
  for (int k=0; k<nkvec; ++k) {				// sum over k vectors
    if (func[0]) { 					// 1/r
      register double r = cek->re*cek->re+cek->im*cek->im; ++cek;
      sum[0][0] += *(kv++)*r; sum[0][1] += *(kv++)*r; sum[0][2] += *(kv++)*r;
      sum[0][3] += *(kv++)*r; sum[0][4] += *(kv++)*r; sum[0][5] += *(kv++)*r;
    }
    if (func[1]) {					// geometric 1/r^6
      register double r = cek->re*cek->re+cek->im*cek->im; ++cek;
      sum[1][0] += *(kv++)*r; sum[1][1] += *(kv++)*r; sum[1][2] += *(kv++)*r;
      sum[1][3] += *(kv++)*r; sum[1][4] += *(kv++)*r; sum[1][5] += *(kv++)*r;
    }
    if (func[2]) {					// arithmetic 1/r^6
      register double r =
	    (cek[0].re*cek[6].re+cek[0].im*cek[6].im)+
	    (cek[1].re*cek[5].re+cek[1].im*cek[5].im)+
	    (cek[2].re*cek[4].re+cek[2].im*cek[4].im)+
	0.5*(cek[3].re*cek[3].re+cek[3].im*cek[3].im); cek += 7;
      sum[2][0] += *(kv++)*r; sum[2][1] += *(kv++)*r; sum[2][2] += *(kv++)*r;
      sum[2][3] += *(kv++)*r; sum[2][4] += *(kv++)*r; sum[2][5] += *(kv++)*r;
    }
    if (func[3]) {
      register double r = cek->re*cek->re+cek->im*cek->im; ++cek;
      sum[3][0] += *(kv++)*r; sum[3][1] += *(kv++)*r; sum[3][2] += *(kv++)*r;
      sum[3][3] += *(kv++)*r; sum[3][4] += *(kv++)*r; sum[3][5] += *(kv++)*r;
    }
  }
  for (int k=0; k<EWALD_NFUNCS; ++k)
    if (func[k]) {
      shape self = {virial_self[k], virial_self[k], virial_self[k], 0, 0, 0};
      shape_scalar_mult(sum[k], c[k]);
      shape_add(virial, sum[k]);
      shape_subtr(virial, self);
    }
}

void EwaldN::compute_virial_peratom()
{
  if (!vflag_atom) return;
      
  kvector *k;
  hvector *h, *nh;
  cvector *z = ekr_local;
  vector  mui = VECTOR_NULL;
  complex *cek, zc = COMPLEX_NULL, zx = COMPLEX_NULL, zxy = COMPLEX_NULL;
  double *kv;
  double *q = atom->q;
  double *vatomj = vatom[0];
  double *mu = atom->mu ? atom->mu[0] : NULL;
  const double qscale = force->qqrd2e * scale;
  double c[EWALD_NFUNCS] = {
    4.0*MY_PI*qscale/volume, 2.0*MY_PI*sqrt(MY_PI)/(24.0*volume),
    2.0*MY_PI*sqrt(MY_PI)/(192.0*volume), 4.0*MY_PI*mumurd2e/volume};
  shape sum[EWALD_MAX_NSUMS];
  int func[EWALD_NFUNCS];

  memcpy(func, function, EWALD_NFUNCS*sizeof(int));
  int i, kx, ky, lbytes = (2*nbox+1)*sizeof(cvector), *type = atom->type;
  for (int j = 0; j < atom->nlocal; j++, vatomj += 6) {					
    k = kvec;						
    kx = ky = -1;					
    kv = kvirial;
    cek = cek_global;
    memset(sum, 0, EWALD_MAX_NSUMS*sizeof(shape));
    if (func[3]) { 
      register double di = mu[3] * c[3];
      mui[0] = di*(mu++)[0]; mui[1] = di*(mu++)[1]; mui[2] = di*(mu++)[2];
      mu++;
    }
    for (nh = (h = hvec)+nkvec; h<nh; ++h, ++k) {
      if (ky!=k->y) {					// based on order in
	  if (kx!=k->x) zx = z[kx = k->x].x; 		// reallocate
	  C_RMULT(zxy, z[ky = k->y].y, zx);
      }
      C_CRMULT(zc, z[k->z].z, zxy);
      if (func[0]) {					// 1/r
	  register double r = cek->re*zc.re - cek->im*zc.im; ++cek;
	  sum[0][0] += *(kv++)*r;
	  sum[0][1] += *(kv++)*r;
	  sum[0][2] += *(kv++)*r;
	  sum[0][3] += *(kv++)*r;
	  sum[0][4] += *(kv++)*r;
	  sum[0][5] += *(kv++)*r;
      }
      if (func[1]) {					// geometric 1/r^6
	  register double r = cek->re*zc.re - cek->im*zc.im; ++cek;
	  sum[1][0] += *(kv++)*r;
	  sum[1][1] += *(kv++)*r;
	  sum[1][2] += *(kv++)*r;
	  sum[1][3] += *(kv++)*r;
	  sum[1][4] += *(kv++)*r;
	  sum[1][5] += *(kv++)*r;
      }
      if (func[2]) {					// arithmetic 1/r^6
	register double r;
	for (i=2; i<9; ++i) {
	  r = cek->re*zc.re - cek->im*zc.im; ++cek;
	  sum[i][0] += *(kv++)*r;
	  sum[i][1] += *(kv++)*r;
	  sum[i][2] += *(kv++)*r;
	  sum[i][3] += *(kv++)*r;
	  sum[i][4] += *(kv++)*r;
	  sum[i][5] += *(kv++)*r;
      kv -= 6;
	}
    kv += 6;
      }
      if (func[3]) {					// dipole
	 register double
	   r = (cek->re*zc.re - cek->im*zc.im)
	      *(mui[0]*h->x+mui[1]*h->y+mui[2]*h->z); ++cek;
	 sum[9][0] += *(kv++)*r;
	 sum[9][1] += *(kv++)*r;
	 sum[9][2] += *(kv++)*r;
	 sum[9][3] += *(kv++)*r;
	 sum[9][4] += *(kv++)*r;
	 sum[9][5] += *(kv++)*r;
      }
    }
		
    if (func[0]) {					// 1/r
      register double qi = *(q++)*c[0];
      for (int n = 0; n < 6; n++) vatomj[n] += sum[0][n]*qi;
    }
    if (func[1]) {					// geometric 1/r^6
      register double bi = B[*type]*c[1];
      for (int n = 0; n < 6; n++) vatomj[n] += sum[1][n]*bi;
    }
    if (func[2]) {					// arithmetic 1/r^6
      register double *bj = B+7*type[0]+7;
      for (i=2; i<9; ++i) {
	register double c2 = (--bj)[0]*c[2];
	for (int n = 0; n < 6; n++) vatomj[n] += 0.5*sum[i][n]*c2;
      }
    }
    if (func[3]) {					// dipole
      for (int n = 0; n < 6; n++) vatomj[n] += sum[9][n];
    }

    for (int k=0; k<EWALD_NFUNCS; ++k) {
      if (func[k]) {
	for (int n = 0; n < 3; n++) vatomj[n] -= virial_self_peratom[j][k]; 
      }
    }
    
    z = (cvector *) ((char *) z+lbytes);
    ++type;
  }
}


/* ----------------------------------------------------------------------
   Slab-geometry correction term to dampen inter-slab interactions between
   periodically repeating slabs.  Yields good approximation to 2-D EwaldN if 
   adequate empty space is left between repeating slabs (J. Chem. Phys. 
   111, 3155).  Slabs defined here to be parallel to the xy plane. 
------------------------------------------------------------------------- */

void EwaldN::compute_slabcorr()
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
  
  const double e_slabcorr = 2.0*MY_PI*dipole_all*dipole_all/volume;
  const double qscale = force->qqrd2e * scale;
  
  if (eflag_global) energy += qscale * e_slabcorr;
  
  // per-atom energy
  
  if (eflag_atom) {
    double efact = 2.0*MY_PI*dipole_all/volume; 
    for (int i = 0; i < nlocal; i++) eatom[i] += qscale * q[i]*x[i][2]*efact;
  }
  
  // add on force corrections
  
  double ffact = -4.0*MY_PI*dipole_all/volume; 
  double **f = atom->f;
  
  for (int i = 0; i < nlocal; i++) f[i][2] += qscale * q[i]*ffact;
}

/* ----------------------------------------------------------------------
  Newton solver used to find g_ewald for LJ systems
 ------------------------------------------------------------------------- */

double EwaldN::NewtonSolve(double x, double Rc, bigint natoms, double vol, double b2)
{
  double dx,tol;
  int maxit;

  maxit = 10000; //Maximum number of iterations
  tol = 0.00001; //Convergence tolerance
        
  //Begin algorithm
  
  for (int i = 0; i < maxit; i++) {
    dx = f(x,Rc,natoms,vol,b2) / derivf(x,Rc,natoms,vol,b2); 
    x = x - dx; //Update x
    if (fabs(dx) < tol) return x;
  }
  return -1;
}

/* ----------------------------------------------------------------------
 Calculate f(x)
 ------------------------------------------------------------------------- */

double EwaldN::f(double x, double Rc, bigint natoms, double vol, double b2)
{
  double a = Rc*x;
  double f = (4.0*MY_PI*b2*pow(x,4)/vol/sqrt(natoms)*erfc(a) *
    (6.0*pow(a,-5.0) + 6.0*pow(a,-3.0) + 3.0/a + a) - accuracy);
  return f;
}

/* ----------------------------------------------------------------------
 Calculate numerical derivative f'(x)
 ------------------------------------------------------------------------- */
            
double EwaldN::derivf(double x, double Rc, bigint natoms, double vol, double b2)
{  
  double h = 0.000001;  //Derivative step-size
  return (f(x + h,Rc,natoms,vol,b2) - f(x,Rc,natoms,vol,b2)) / h;
}       

