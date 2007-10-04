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

#include "mpi.h"
#include "string.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "ewald_n.h"
#include "math_vector.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "pair.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define KSPACE_ILLEGAL	"Illegal kspace_style ewald/n command"
#define KSPACE_ORDER	"Unsupported order in kspace_style ewald/n for"
#define KSPACE_MIX	"Unsupported mixing rule in kspace_style ewald/n for"

enum{GEOMETRIC,ARITHMETIC,SIXTHPOWER};         // same as in pair.cpp

//#define DEBUG

/* ---------------------------------------------------------------------- */

EwaldN::EwaldN(LAMMPS *lmp, int narg, char **arg) : KSpace(lmp, narg, arg)
{
  if (narg!=1) error->all(KSPACE_ILLEGAL);
  precision = fabs(atof(arg[0]));
  memset(function, 0, EWALD_NORDER*sizeof(int));
  kenergy = kvirial = NULL;
  cek_local = cek_global = NULL;
  ekr_local = NULL;
  hvec = NULL;
  kvec = NULL;
  B = NULL;
  first_output = 0;
}

EwaldN::~EwaldN()
{
  deallocate();
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
    error->all("Cannot use EwaldN with 2d simulation");
  if (slabflag == 0 && domain->nonperiodic > 0)
    error->all("Cannot use nonperiodic boundaries with EwaldN");
  if (slabflag == 1) {
    if (domain->xperiodic != 1 || domain->yperiodic != 1 || 
	domain->boundary[2][0] != 1 || domain->boundary[2][1] != 1)
      error->all("Incorrect boundaries with slab EwaldN");
  }

  qqrd2e = force->qqrd2e;				// check pair_style
  //mumurd2e = force->mumurd2e;
  //dielectric = force->dielectric;
  mumurd2e = dielectric = 1.0;
  
  Pair *pair = force->pair;
  int *ptr = pair ? (int *) pair->extract("ewald_order") : NULL;
  double *cutoff = pair ? (double *) pair->extract("cut_coul") : NULL;
  if (!(ptr||cutoff)) 
    error->all("KSpace style is incompatible with Pair style");
  int ewald_order = ptr ? *((int *) ptr) : 1<<1;
  int ewald_mix = ptr ? *((int *) pair->extract("ewald_mix")) : GEOMETRIC;
  memset(function, 0, EWALD_NFUNCS*sizeof(int));
  for (int i=0; i<=EWALD_NORDER; ++i)			// transcribe order
    if (ewald_order&(1<<i)) {				// from pair_style
      int n[] = EWALD_NSUMS, k;
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
	  error->all(str);
	default:
	  sprintf(str, "%s pair_style %s", KSPACE_ORDER, force->pair_style);
	  error->all(str);
      }
      nfunctions += function[k] = 1;
      nsums += n[k];
    }
  
  g_ewald = (1.35 - 0.15*log(precision))/ *cutoff;	// determine resolution
  g2_max = -4.0*g_ewald*g_ewald*log(precision);

  if (!comm->me) {					// output results
    if (screen) fprintf(screen, "  G vector = %g\n", g_ewald);
    if (logfile) fprintf(logfile, "  G vector = %g\n", g_ewald);
  }
}


/* ----------------------------------------------------------------------
   adjust EwaldN coeffs, called initially and whenever volume has changed 
------------------------------------------------------------------------- */

void EwaldN::setup()
{
  volume = shape_det(domain->h)*slab_volfactor;		// cell volume
  memcpy(unit, domain->h_inv, sizeof(shape));		// wave vector units
  shape_scalar_mult(unit, 2.0*M_PI);
  unit[2] /= slab_volfactor;

  //int nbox_old = nbox, nkvec_old = nkvec;
  if (precision>=1) nbox = 0;
  else {
    vector n = {1.0, 1.0, 1.0};				// based on cutoff
    vec_scalar_mult(n, g_ewald*sqrt(-log(precision))/M_PI);
    shape_vec_dot(n, n, domain->h);
    n[2] *= slab_volfactor;
    nbox = (int) n[0];
    if (nbox<(int) n[1]) nbox = (int) n[1];
    if (nbox<(int) n[2]) nbox = (int) n[2];
  }
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


void EwaldN::reallocate()				// allocate memory
{
  int ix, iy, iz;
  int nkvec_max = nkvec;
  vector h;

  nkvec = 0;						// determine size(kvec)
  int kflag[(nbox+1)*(2*nbox+1)*(2*nbox+1)], *flag = kflag;
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
	  if ((*(flag++) = h[0]*h[0]+h[1]*h[1]+h[2]*h[2]<=g2_max)) ++nkvec;
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

}


void EwaldN::reallocate_atoms()
{
  if ((nevec = atom->nmax*(2*nbox+1))<=nevec_max) return;
  delete [] ekr_local;
  ekr_local = new cvector[nevec];
  bytes += (nevec-nevec_max)*sizeof(cvector);
  nevec_max = nevec;
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
      *(ke++) = c1 = -h1*h2*((c2=sqrt(M_PI)*erfc(b1))+(0.5/b2-1.0)*expb2/b1);
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
  int n = atom->ntypes;

  if (function[1]) {					// geometric 1/r^6
    double **b = (double **) force->pair->extract("B");
    delete [] B;
    B = new double[n+1];
    bytes += (n+1)*sizeof(double);
    for (int i=0; i<=n; ++i) B[i] = sqrt(fabs(b[i][i]));
  }
  if (function[2]) {					// arithmetic 1/r^6
    double **epsilon = (double **) force->pair->extract("epsilon");
    double **sigma = (double **) force->pair->extract("sigma");
    if (!(epsilon&&sigma))
      error->all("epsilon or sigma reference not set by pair style in ewald/n");
    double eps_i, sigma_i, sigma_n, *bi = B = new double[7*n+7];
    double c[7] = {
      1.0, sqrt(6.0), sqrt(15.0), sqrt(20.0), sqrt(15.0), sqrt(6.0), 1.0};
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
  if (function[2]) {					// aritmetic 1/r^6
    double *bi;
    int *type = atom->type, *ntype = type+atom->nlocal;
    for (int *i=type; i<ntype; ++i) {
      bi = B+7*i[0];
      sum_local[2].x2 += bi[0]*bi[6];
      for (int k=2; k<9; ++k) sum_local[k].x += *(bi++);
    }
  }
  if (function[3]) {					// dipole
    int *type = atom->type, *ntype = type+atom->nlocal;
    double *dipole = atom->dipole;
    for (int *i=type; i<ntype; ++i)
      sum_local[9].x2 += dipole[i[0]]*dipole[i[0]];
  }
  MPI_Allreduce(sum_local, sum, 2*EWALD_MAX_NSUMS, MPI_DOUBLE, MPI_SUM, world);
}


void EwaldN::init_self()
{
  double g1 = g_ewald, g2 = g1*g1, g3 = g1*g2;

  memset(energy_self, 0, EWALD_NFUNCS*sizeof(double));	// self energy
  memset(virial_self, 0, EWALD_NFUNCS*sizeof(double));
  if (function[0]) {					// 1/r
    virial_self[0] = -0.5*M_PI*qqrd2e/(g2*volume)*sum[0].x*sum[0].x;
    energy_self[0] = sum[0].x2*qqrd2e*g1/sqrt(M_PI)-virial_self[0];
  }
  if (function[1]) {					// geometric 1/r^6
    virial_self[1] = M_PI*sqrt(M_PI)*g3/(6.0*volume)*sum[1].x*sum[1].x;
    energy_self[1] = -sum[1].x2*g3*g3/12.0+virial_self[1];
  }
  if (function[2]) {					// arithmetic 1/r^6
    virial_self[2] = M_PI*sqrt(M_PI)*g3/(48.0*volume)*(sum[2].x*sum[8].x+
	sum[3].x*sum[7].x+sum[4].x*sum[6].x+0.5*sum[5].x*sum[5].x);
    energy_self[2] = -sum[2].x2*g3*g3/3.0+virial_self[2];
  }
  if (function[3]) {					// dipole
    virial_self[3] = 0;					// in surface
    energy_self[3] = sum[9].x2*mumurd2e*2.0*g3/3.0/sqrt(M_PI)-virial_self[3];
  }
}

/* ----------------------------------------------------------------------
   compute the EwaldN long-range force, energy, virial 
------------------------------------------------------------------------- */

void EwaldN::compute(int eflag, int vflag)
{
  if (!nbox) return;
  reallocate_atoms();
  compute_ek();
  compute_force();
  compute_surface();
  compute_energy(eflag);
  compute_virial(vflag);
}


void EwaldN::compute_ek()
{
  cvector *ekr = ekr_local;
  int lbytes = (2*nbox+1)*sizeof(cvector);
  hvector *h;
  kvector *k, *nk = kvec+nkvec;
  cvector z1, z[2*nbox+1], *zx, *zy, *zz, *zn = z+2*nbox;
  complex *cek, zxyz, zxy, cx;
  vector mui;
  double *x = atom->x[0], *xn = x+3*atom->nlocal, *q = atom->q, qi, bi, ci[7];
  double *dipole = atom->dipole, *mu = atom->mu ? atom->mu[0] : NULL;
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
      memcpy(mui, mu, sizeof(vector)); mu += 3;
      vec_scalar_mult(mui, dipole[*type]);
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
}


void EwaldN::compute_force()
{
  kvector *k;
  hvector *h, *nh;
  cvector *z = ekr_local;
  vector sum[EWALD_MAX_NSUMS], mui;
  complex *cek, zc, zx, zxy;
  double *f = atom->f[0], *fn = f+3*atom->nlocal, *q = atom->q, *t = NULL;
  double *dipole = atom->dipole, *mu = atom->mu ? atom->mu[0] : NULL;
  double *ke, c[EWALD_NFUNCS] = {
    8.0*M_PI*qqrd2e/volume, 2.0*M_PI*sqrt(M_PI)/(12.0*volume),
    2.0*M_PI*sqrt(M_PI)/(192.0*volume), 8.0*M_PI*mumurd2e/volume};
  double kt = 4.0*pow(g_ewald, 3.0)/3.0/sqrt(M_PI)/c[3];
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
      register double di = dipole[*type]*c[3];
      mui[0] = di*(mu++)[0]; mui[1] = di*(mu++)[1]; mui[2] = di*(mu++)[2];
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

  vector sum_local = VECTOR_NULL, sum_total;
  double *mu = atom->mu ? atom->mu[0] : NULL, *dipole = atom->dipole;
  int *type = atom->type, *ntype = type+atom->nlocal;

  for (int *i=type; i<ntype; ++i) {
    register double di = dipole[i[0]];
    sum_local[0] += di*(mu++)[0];
    sum_local[1] += di*(mu++)[1];
    sum_local[2] += di*(mu++)[2];
  }
  MPI_Allreduce(sum_local, sum_total, 3, MPI_DOUBLE, MPI_SUM, world);

  energy_self[3] += virial_self[3];
  virial_self[3] =
    mumurd2e*(2.0*M_PI*vec_dot(sum_total,sum_total)/(2.0*dielectric+1)/volume);
  energy_self[3] -= virial_self[3];
}


void EwaldN::compute_energy(int eflag)
{
  energy = 0.0;
  if (!eflag) return;
  
  complex *cek = cek_global;
  double *ke = kenergy;
  double c[EWALD_NFUNCS] = {
    4.0*M_PI*qqrd2e/volume, 2.0*M_PI*sqrt(M_PI)/(24.0*volume),
    2.0*M_PI*sqrt(M_PI)/(192.0*volume), 4.0*M_PI*mumurd2e/volume};
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
  if (slabflag) compute_slabcorr(eflag);
}


#define swap(a, b) { register double t = a; a= b; b = t; }

void EwaldN::compute_virial(int vflag)
{
  memset(virial, 0, sizeof(shape));
  if (!vflag) return;

  complex *cek = cek_global;
  double *kv = kvirial;
  double c[EWALD_NFUNCS] = {
    4.0*M_PI*qqrd2e/volume, 2.0*M_PI*sqrt(M_PI)/(24.0*volume),
    2.0*M_PI*sqrt(M_PI)/(192.0*volume), 4.0*M_PI*mumurd2e/volume};
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

/* ----------------------------------------------------------------------
   Slab-geometry correction term to dampen inter-slab interactions between
   periodically repeating slabs.  Yields good approximation to 2-D EwaldN if 
   adequate empty space is left between repeating slabs (J. Chem. Phys. 
   111, 3155).  Slabs defined here to be parallel to the xy plane. 
------------------------------------------------------------------------- */

void EwaldN::compute_slabcorr(int eflag)
{
  // compute local contribution to global dipole moment
  
  double *q = atom->q;
  double *x = atom->x[0]-1, *xn = x+3*atom->nlocal-1;
  double dipole = 0.0, dipole_all;
  
  while ((x+=3)<xn) dipole += *x * *(q++);
  MPI_Allreduce(&dipole, &dipole_all, 1, MPI_DOUBLE, MPI_SUM, world);
  
  double ffact = -4.0*M_PI*qqrd2e*dipole_all/volume;	// force correction
  double *f = atom->f[0]-1, *fn = f+3*atom->nlocal-3;
  
  q = atom->q;
  while ((f+=3)<fn) *f += ffact* *(q++);

  if (eflag) 						// energy correction
    energy += qqrd2e*(2.0*M_PI*dipole_all*dipole_all/volume);
}

