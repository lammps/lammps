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
   Contributing authors: Pieter in 't Veld (SNL), Stan Moore (SNL)
------------------------------------------------------------------------- */

#include "ewald_disp.h"
#include <mpi.h>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "math_vector.h"
#include "math_const.h"
#include "math_special.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "pair.h"
#include "domain.h"
#include "memory.h"
#include "error.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathSpecial;

#define SMALL 0.00001

//#define DEBUG

/* ---------------------------------------------------------------------- */

EwaldDisp::EwaldDisp(LAMMPS *lmp) : KSpace(lmp),
  kenergy(NULL), kvirial(NULL), energy_self_peratom(NULL), virial_self_peratom(NULL),
  ekr_local(NULL), hvec(NULL), kvec(NULL), B(NULL), cek_local(NULL), cek_global(NULL)
{
  ewaldflag = dispersionflag = dipoleflag = 1;

  memset(function, 0, EWALD_NFUNCS*sizeof(int));
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
  M2 = 0;
}

void EwaldDisp::settings(int narg, char **arg)
{
  if (narg!=1) error->all(FLERR,"Illegal kspace_style ewald/n command");
  accuracy_relative = fabs(force->numeric(FLERR,arg[0]));
}


/* ---------------------------------------------------------------------- */

EwaldDisp::~EwaldDisp()
{
  deallocate();
  deallocate_peratom();
  delete [] ekr_local;
  delete [] B;
}

/* --------------------------------------------------------------------- */

void EwaldDisp::init()
{
  nkvec = nkvec_max = nevec = nevec_max = 0;
  nfunctions = nsums = sums = 0;
  nbox = -1;
  bytes = 0.0;

  if (!comm->me) {
    if (screen) fprintf(screen,"EwaldDisp initialization ...\n");
    if (logfile) fprintf(logfile,"EwaldDisp initialization ...\n");
  }

  triclinic_check();
  if (domain->dimension == 2)
    error->all(FLERR,"Cannot use EwaldDisp with 2d simulation");
  if (slabflag == 0 && domain->nonperiodic > 0)
    error->all(FLERR,"Cannot use non-periodic boundaries with EwaldDisp");
  if (slabflag == 1) {
    if (domain->xperiodic != 1 || domain->yperiodic != 1 ||
        domain->boundary[2][0] != 1 || domain->boundary[2][1] != 1)
      error->all(FLERR,"Incorrect boundaries with slab EwaldDisp");
  }

  scale = 1.0;
  mumurd2e = force->qqrd2e;
  dielectric = force->dielectric;

  int tmp;
  Pair *pair = force->pair;
  int *ptr = pair ? (int *) pair->extract("ewald_order",tmp) : NULL;
  double *cutoff = pair ? (double *) pair->extract("cut_coul",tmp) : NULL;
  if (!(ptr||cutoff))
    error->all(FLERR,"KSpace style is incompatible with Pair style");
  int ewald_order = ptr ? *((int *) ptr) : 1<<1;
  int ewald_mix = ptr ? *((int *) pair->extract("ewald_mix",tmp)) : Pair::GEOMETRIC;
  memset(function, 0, EWALD_NFUNCS*sizeof(int));
  for (int i=0; i<=EWALD_NORDER; ++i)                        // transcribe order
    if (ewald_order&(1<<i)) {                                // from pair_style
      int n[] = EWALD_NSUMS, k = 0;
      switch (i) {
        case 1:
          k = 0; break;
        case 3:
          k = 3; break;
        case 6:
          if (ewald_mix==Pair::GEOMETRIC) { k = 1; break; }
          else if (ewald_mix==Pair::ARITHMETIC) { k = 2; break; }
          error->all(FLERR,
                     "Unsupported mixing rule in kspace_style ewald/disp");
        default:
          error->all(FLERR,"Unsupported order in kspace_style ewald/disp");
      }
      nfunctions += function[k] = 1;
      nsums += n[k];
    }

  if (!gewaldflag) g_ewald = g_ewald_6 = 1.0;
  pair->init();  // so B is defined
  init_coeffs();
  init_coeff_sums();
  if (function[0]) qsum_qsq();
  else qsqsum = qsum = 0.0;
  natoms_original = atom->natoms;
  if (!gewaldflag) g_ewald = g_ewald_6 = 0.0;

  // turn off coulombic if no charge

  if (function[0] && qsqsum == 0.0) {
    function[0] = 0;
    nfunctions -= 1;
    nsums -= 1;
  }

  double bsbsum = 0.0;
  M2 = 0.0;
  if (function[1]) bsbsum = sum[1].x2;
  if (function[2]) bsbsum = sum[2].x2;

  if (function[3]) M2 = sum[9].x2;

  if (function[3] && strcmp(update->unit_style,"electron") == 0)
    error->all(FLERR,"Cannot (yet) use 'electron' units with dipoles");

  if (qsqsum == 0.0 && bsbsum == 0.0 && M2 == 0.0)
      error->all(FLERR,"Cannot use Ewald/disp solver "
                 "on system with no charge, dipole, or LJ particles");
  if (fabs(qsum) > SMALL && comm->me == 0) {
      char str[128];
      sprintf(str,"System is not charge neutral, net charge = %g",qsum);
      error->warning(FLERR,str);
  }

  if (!function[1] && !function[2]) dispersionflag = 0;
  if (!function[3]) dipoleflag = 0;

  pair_check();

  // set accuracy (force units) from accuracy_relative or accuracy_absolute

  if (accuracy_absolute >= 0.0) accuracy = accuracy_absolute;
  else accuracy = accuracy_relative * two_charge_force;

  // setup K-space resolution

  q2 = qsqsum * force->qqrd2e;
  M2 *= mumurd2e;
  b2 = bsbsum; //Are these units right?
  bigint natoms = atom->natoms;

  if (!gewaldflag) {
    if (function[0]) {
      if (accuracy <= 0.0)
        error->all(FLERR,"KSpace accuracy must be > 0");
      if (q2 == 0.0)
        error->all(FLERR,"Must use 'kspace_modify gewald' for uncharged system");
      g_ewald = accuracy*sqrt(natoms*(*cutoff)*shape_det(domain->h)) / (2.0*q2);
      if (g_ewald >= 1.0) g_ewald = (1.35 - 0.15*log(accuracy))/(*cutoff);
      else g_ewald = sqrt(-log(g_ewald)) / (*cutoff);
    } else if (function[3]) {
      //Try Newton Solver
      //Use old method to get guess
      g_ewald = (1.35 - 0.15*log(accuracy))/ *cutoff;
      double g_ewald_new =
        NewtonSolve(g_ewald,(*cutoff),natoms,shape_det(domain->h),M2);
      if (g_ewald_new > 0.0) g_ewald = g_ewald_new;
      else error->warning(FLERR,"Ewald/disp Newton solver failed, "
                          "using old method to estimate g_ewald");
    } else if (function[1] || function[2]) {
      //Try Newton Solver
      //Use old method to get guess
      g_ewald = (1.35 - 0.15*log(accuracy))/ *cutoff;

      double g_ewald_new =
        NewtonSolve(g_ewald,(*cutoff),natoms,shape_det(domain->h),b2);
      if (g_ewald_new > 0.0) g_ewald = g_ewald_new;
      else error->warning(FLERR,"Ewald/disp Newton solver failed, "
                          "using old method to estimate g_ewald");
    }
  }

  if (!comm->me) {
      if (screen) fprintf(screen, "  G vector = %g,   accuracy = %g\n", g_ewald,accuracy);
      if (logfile) fprintf(logfile, "  G vector = %g   accuracy = %g\n", g_ewald,accuracy);
  }

  g_ewald_6 = g_ewald;
  deallocate_peratom();
  peratom_allocate_flag = 0;
}

/* ----------------------------------------------------------------------
   adjust EwaldDisp coeffs, called initially and whenever volume has changed
------------------------------------------------------------------------- */

void EwaldDisp::setup()
{
  volume = shape_det(domain->h)*slab_volfactor;
  memcpy(unit, domain->h_inv, sizeof(shape));
  shape_scalar_mult(unit, 2.0*MY_PI);
  unit[2] /= slab_volfactor;

  // int nbox_old = nbox, nkvec_old = nkvec;

  if (accuracy >= 1) {
    nbox = 0;
    error->all(FLERR,"KSpace accuracy too low");
  }

  bigint natoms = atom->natoms;
  double err;
  int kxmax = 1;
  int kymax = 1;
  int kzmax = 1;
  err = rms(kxmax,domain->h[0],natoms,q2,b2,M2);
  while (err > accuracy) {
    kxmax++;
    err = rms(kxmax,domain->h[0],natoms,q2,b2,M2);
  }
  err = rms(kymax,domain->h[1],natoms,q2,b2,M2);
  while (err > accuracy) {
    kymax++;
    err = rms(kymax,domain->h[1],natoms,q2,b2,M2);
  }
  err = rms(kzmax,domain->h[2]*slab_volfactor,natoms,q2,b2,M2);
  while (err > accuracy) {
    kzmax++;
    err = rms(kzmax,domain->h[2]*slab_volfactor,natoms,q2,b2,M2);
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
  coefficients();
  init_coeffs();
  init_coeff_sums();
  init_self();

  if (!(first_output||comm->me)) {
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

double EwaldDisp::rms(int km, double prd, bigint natoms,
                      double q2, double b2, double M2)
{
  double value = 0.0;
  if (natoms == 0) natoms = 1; // avoid division by zero

  // Coulombic

  double g2 = g_ewald*g_ewald;

  value += 2.0*q2*g_ewald/prd *
    sqrt(1.0/(MY_PI*km*natoms)) *
    exp(-MY_PI*MY_PI*km*km/(g2*prd*prd));

  // Lennard-Jones

  double g7 = g2*g2*g2*g_ewald;

  value += 4.0*b2*g7/3.0 *
    sqrt(1.0/(MY_PI*natoms)) *
    (exp(-MY_PI*MY_PI*km*km/(g2*prd*prd)) *
    (MY_PI*km/(g_ewald*prd) + 1));

  // dipole

  value += 8.0*MY_PI*M2/volume*g_ewald *
    sqrt(2.0*MY_PI*km*km*km/(15.0*natoms)) *
    exp(-pow(MY_PI*km/(g_ewald*prd),2.0));

  return value;
}

void EwaldDisp::reallocate()
{
  int ix, iy, iz;
  int nkvec_max = nkvec;
  vector h;

  nkvec = 0;
  int *kflag = new int[(nbox+1)*(2*nbox+1)*(2*nbox+1)];
  int *flag = kflag;

  for (ix=0; ix<=nbox; ++ix)
    for (iy=-nbox; iy<=nbox; ++iy)
      for (iz=-nbox; iz<=nbox; ++iz)
        if (!(ix||iy||iz)) *(flag++) = 0;
        else if ((!ix)&&(iy<0)) *(flag++) = 0;
        else if ((!(ix||iy))&&(iz<0)) *(flag++) = 0;        // use symmetry
        else {
          h[0] = unit[0]*ix;
          h[1] = unit[5]*ix+unit[1]*iy;
          h[2] = unit[4]*ix+unit[3]*iy+unit[2]*iz;
          if ((*(flag++) = h[0]*h[0]+h[1]*h[1]+h[2]*h[2]<=gsqmx)) ++nkvec;
        }

  if (nkvec>nkvec_max) {
    deallocate();                                        // free memory
    hvec = new hvector[nkvec];                                // hvec
    bytes += (nkvec-nkvec_max)*sizeof(hvector);
    kvec = new kvector[nkvec];                                // kvec
    bytes += (nkvec-nkvec_max)*sizeof(kvector);
    kenergy = new double[nkvec*nfunctions];                // kenergy
    bytes += (nkvec-nkvec_max)*nfunctions*sizeof(double);
    kvirial = new double[6*nkvec*nfunctions];                // kvirial
    bytes += 6*(nkvec-nkvec_max)*nfunctions*sizeof(double);
    cek_local = new complex[nkvec*nsums];                // cek_local
    bytes += (nkvec-nkvec_max)*nsums*sizeof(complex);
    cek_global = new complex[nkvec*nsums];                // cek_global
    bytes += (nkvec-nkvec_max)*nsums*sizeof(complex);
    nkvec_max = nkvec;
  }

  flag = kflag;                                           // create index and
  kvector *k = kvec;                                      // wave vectors
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

/* ---------------------------------------------------------------------- */

void EwaldDisp::reallocate_atoms()
{
  if (eflag_atom || vflag_atom)
    if (atom->nmax > nmax) {
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

/* ---------------------------------------------------------------------- */

void EwaldDisp::allocate_peratom()
{
  memory->create(energy_self_peratom,
      atom->nmax,EWALD_NFUNCS,"ewald/n:energy_self_peratom");
  memory->create(virial_self_peratom,
      atom->nmax,EWALD_NFUNCS,"ewald/n:virial_self_peratom");
}

/* ---------------------------------------------------------------------- */

void EwaldDisp::deallocate_peratom()                        // free memory
{
  if (energy_self_peratom) {
    memory->destroy(energy_self_peratom);
    energy_self_peratom = NULL;
  }

  if (virial_self_peratom) {
    memory->destroy(virial_self_peratom);
    virial_self_peratom = NULL;
  }
}

/* ---------------------------------------------------------------------- */

void EwaldDisp::deallocate()                                // free memory
{
  delete [] hvec;                hvec = NULL;
  delete [] kvec;                kvec = NULL;
  delete [] kenergy;                kenergy = NULL;
  delete [] kvirial;                kvirial = NULL;
  delete [] cek_local;                cek_local = NULL;
  delete [] cek_global;                cek_global = NULL;
}

/* ---------------------------------------------------------------------- */

void EwaldDisp::coefficients()
{
  vector h;
  hvector *hi = hvec, *nh;
  double eta2 = 0.25/(g_ewald*g_ewald);
  double b1, b2, expb2, h1, h2, c1, c2;
  double *ke = kenergy, *kv = kvirial;
  int func0 = function[0], func12 = function[1]||function[2],
      func3 = function[3];

  for (nh = (hi = hvec)+nkvec; hi<nh; ++hi) {                // wave vectors
    memcpy(h, hi, sizeof(vector));
    expb2 = exp(-(b2 = (h2 = vec_dot(h, h))*eta2));
    if (func0) {                                        // qi*qj/r coeffs
      *(ke++) = c1 = expb2/h2;
      *(kv++) = c1-(c2 = 2.0*c1*(1.0+b2)/h2)*h[0]*h[0];
      *(kv++) = c1-c2*h[1]*h[1];                        // lammps convention
      *(kv++) = c1-c2*h[2]*h[2];                        // instead of voigt
      *(kv++) = -c2*h[1]*h[0];
      *(kv++) = -c2*h[2]*h[0];
      *(kv++) = -c2*h[2]*h[1];
    }
    if (func12) {                                        // -Bij/r^6 coeffs
      b1 = sqrt(b2);                                        // minus sign folded
      h1 = sqrt(h2);                                        // into constants
      *(ke++) = c1 = -h1*h2*((c2=MY_PIS*erfc(b1))+(0.5/b2-1.0)*expb2/b1);
      *(kv++) = c1-(c2 = 3.0*h1*(c2-expb2/b1))*h[0]*h[0];
      *(kv++) = c1-c2*h[1]*h[1];                        // lammps convention
      *(kv++) = c1-c2*h[2]*h[2];                        // instead of voigt
      *(kv++) = -c2*h[1]*h[0];
      *(kv++) = -c2*h[2]*h[0];
      *(kv++) = -c2*h[2]*h[1];
    }
    if (func3) {                                        // dipole coeffs
      *(ke++) = c1 = expb2/h2;
      *(kv++) = c1-(c2 = 2.0*c1*(1.0+b2)/h2)*h[0]*h[0];
      *(kv++) = c1-c2*h[1]*h[1];                        // lammps convention
      *(kv++) = c1-c2*h[2]*h[2];                        // instead of voigt
      *(kv++) = -c2*h[1]*h[0];
      *(kv++) = -c2*h[2]*h[0];
      *(kv++) = -c2*h[2]*h[1];
    }
  }
}

/* ---------------------------------------------------------------------- */

void EwaldDisp::init_coeffs()
{
  int tmp;
  int n = atom->ntypes;

  if (function[1]) {                                        // geometric 1/r^6
    double **b = (double **) force->pair->extract("B",tmp);
    delete [] B;
    B = new double[n+1];
    B[0] = 0.0;
    bytes += (n+1)*sizeof(double);
    for (int i=1; i<=n; ++i) B[i] = sqrt(fabs(b[i][i]));
  }
  if (function[2]) {                                        // arithmetic 1/r^6
    double **epsilon = (double **) force->pair->extract("epsilon",tmp);
    double **sigma = (double **) force->pair->extract("sigma",tmp);
    double eps_i, sigma_i, sigma_n, *bi = B = new double[7*n+7];
    double c[7] = {
      1.0, sqrt(6.0), sqrt(15.0), sqrt(20.0), sqrt(15.0), sqrt(6.0), 1.0};

    if (!(epsilon&&sigma))
      error->all(
          FLERR,"Epsilon or sigma reference not set by pair style in ewald/n");
    for (int j=0; j<7; ++j)
      *(bi++) = 0.0;
    for (int i=1; i<=n; ++i) {
      eps_i = sqrt(epsilon[i][i]);
      sigma_i = sigma[i][i];
      sigma_n = 1.0;
      for (int j=0; j<7; ++j) {
        *(bi++) = sigma_n*eps_i*c[j]; sigma_n *= sigma_i;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void EwaldDisp::init_coeff_sums()
{
  if (sums) return;                            // calculated only once
  sums = 1;

  Sum sum_local[EWALD_MAX_NSUMS];

  memset(sum_local, 0, EWALD_MAX_NSUMS*sizeof(Sum));
  memset(sum,       0, EWALD_MAX_NSUMS*sizeof(Sum));

  // now perform qsum and qsq via parent qsum_qsq()

  sum_local[0].x = 0.0;
  sum_local[0].x2 = 0.0;

  //if (function[0]) {                                        // 1/r
  //  double *q = atom->q, *qn = q+atom->nlocal;
  //  for (double *i=q; i<qn; ++i) {
  //    sum_local[0].x += i[0]; sum_local[0].x2 += i[0]*i[0]; }
  //}

  if (function[1]) {                                        // geometric 1/r^6
    int *type = atom->type, *ntype = type+atom->nlocal;
    for (int *i=type; i<ntype; ++i) {
      sum_local[1].x += B[i[0]]; sum_local[1].x2 += B[i[0]]*B[i[0]]; }
  }
  if (function[2]) {                                        // arithmetic 1/r^6
    double *bi;
    int *type = atom->type, *ntype = type+atom->nlocal;
    for (int *i=type; i<ntype; ++i) {
      bi = B+7*i[0];
      sum_local[2].x2 += bi[0]*bi[6];
      for (int k=2; k<9; ++k) sum_local[k].x += *(bi++);
    }
  }
  if (function[3]&&atom->mu) {                                // dipole
    double *mu = atom->mu[0], *nmu = mu+4*atom->nlocal;
    for (double *i = mu; i < nmu; i += 4)
      sum_local[9].x2 += i[3]*i[3];
  }
  MPI_Allreduce(sum_local, sum, 2*EWALD_MAX_NSUMS, MPI_DOUBLE, MPI_SUM, world);
}

/* ---------------------------------------------------------------------- */

void EwaldDisp::init_self()
{
  double g1 = g_ewald, g2 = g1*g1, g3 = g1*g2;
  const double qscale = force->qqrd2e * scale;

  memset(energy_self, 0, EWALD_NFUNCS*sizeof(double));        // self energy
  memset(virial_self, 0, EWALD_NFUNCS*sizeof(double));

  if (function[0]) {                                        // 1/r
    virial_self[0] = -0.5*MY_PI*qscale/(g2*volume)*qsum*qsum;
    energy_self[0] = qsqsum*qscale*g1/MY_PIS-virial_self[0];
  }
  if (function[1]) {                                        // geometric 1/r^6
    virial_self[1] = MY_PI*MY_PIS*g3/(6.0*volume)*sum[1].x*sum[1].x;
    energy_self[1] = -sum[1].x2*g3*g3/12.0+virial_self[1];
  }
  if (function[2]) {                                        // arithmetic 1/r^6
    virial_self[2] = MY_PI*MY_PIS*g3/(48.0*volume)*(sum[2].x*sum[8].x+
        sum[3].x*sum[7].x+sum[4].x*sum[6].x+0.5*sum[5].x*sum[5].x);
    energy_self[2] = -sum[2].x2*g3*g3/3.0+virial_self[2];
  }
  if (function[3]) {                                        // dipole
    virial_self[3] = 0;                                        // in surface
    energy_self[3] = sum[9].x2*mumurd2e*2.0*g3/3.0/MY_PIS-virial_self[3];
  }
}

/* ---------------------------------------------------------------------- */

void EwaldDisp::init_self_peratom()
{
  if (!(vflag_atom || eflag_atom)) return;

  double g1 = g_ewald, g2 = g1*g1, g3 = g1*g2;
  const double qscale = force->qqrd2e * scale;
  double *energy = energy_self_peratom[0];
  double *virial = virial_self_peratom[0];
  int nlocal = atom->nlocal;

  memset(energy, 0, EWALD_NFUNCS*nlocal*sizeof(double));
  memset(virial, 0, EWALD_NFUNCS*nlocal*sizeof(double));

  if (function[0]) {                                        // 1/r
    double *ei = energy;
    double *vi = virial;
    double ce = qscale*g1/MY_PIS;
    double cv = -0.5*MY_PI*qscale/(g2*volume);
    double *qi = atom->q, *qn = qi + nlocal;
    for (; qi < qn; qi++, vi += EWALD_NFUNCS, ei += EWALD_NFUNCS) {
      double q = *qi;
      *vi = cv*q*qsum;
      *ei = ce*q*q-vi[0];
    }
  }
  if (function[1]) {                                        // geometric 1/r^6
    double *ei = energy+1;
    double *vi = virial+1;
    double ce = -g3*g3/12.0;
    double cv = MY_PI*MY_PIS*g3/(6.0*volume);
    int *typei = atom->type, *typen = typei + atom->nlocal;
    for (; typei < typen; typei++, vi += EWALD_NFUNCS, ei += EWALD_NFUNCS) {
      double b = B[*typei];
      *vi = cv*b*sum[1].x;
      *ei = ce*b*b+vi[0];
    }
  }
  if (function[2]) {                                        // arithmetic 1/r^6
    double *bi;
    double *ei = energy+2;
    double *vi = virial+2;
    double ce = -g3*g3/3.0;
    double cv = 0.5*MY_PI*MY_PIS*g3/(48.0*volume);
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
  if (function[3]&&atom->mu) {                                // dipole
    double *ei = energy+3;
    double *vi = virial+3;
    double *imu = atom->mu[0], *nmu = imu+4*atom->nlocal;
    double ce = mumurd2e*2.0*g3/3.0/MY_PIS;
    for (; imu < nmu; imu += 4, vi += EWALD_NFUNCS, ei += EWALD_NFUNCS) {
      *vi = 0;                                                // in surface
      *ei = ce*imu[3]*imu[3]-vi[0];
    }
  }
}

/* ----------------------------------------------------------------------
   compute the EwaldDisp long-range force, energy, virial
------------------------------------------------------------------------- */

void EwaldDisp::compute(int eflag, int vflag)
{
  if (!nbox) return;

  // set energy/virial flags
  // invoke allocate_peratom() if needed for first time

  ev_init(eflag,vflag);

  if (!peratom_allocate_flag && (eflag_atom || vflag_atom)) {
      allocate_peratom();
      peratom_allocate_flag = 1;
      nmax = atom->nmax;
  }

  reallocate_atoms();
  init_self_peratom();
  compute_ek();
  compute_force();
  //compute_surface(); // assume conducting metal (tinfoil) boundary conditions

  // update qsum and qsqsum, if atom count has changed and energy needed

  if ((eflag_global || eflag_atom) && atom->natoms != natoms_original) {
    if (function[0]) qsum_qsq();
    natoms_original = atom->natoms;
  }

  compute_energy();
  compute_energy_peratom();
  compute_virial();
  compute_virial_dipole();
  compute_virial_peratom();

  if (slabflag) compute_slabcorr();
}


void EwaldDisp::compute_ek()
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
  memset(cek_local, 0, n*sizeof(complex));                // reset sums
  while (x<xn) {
    zx = (zy = (zz = z+nbox)+1)-2;
    C_SET(zz->x, 1, 0); C_SET(zz->y, 1, 0); C_SET(zz->z, 1, 0);        // z[0]
    if (tri) {                                                // triclinic z[1]
      C_ANGLE(z1.x, unit[0]*x[0]+unit[5]*x[1]+unit[4]*x[2]);
      C_ANGLE(z1.y, unit[1]*x[1]+unit[3]*x[2]);
      C_ANGLE(z1.z, x[2]*unit[2]); x += 3;
    }
    else {                                                // orthogonal z[1]
      C_ANGLE(z1.x, *(x++)*unit[0]);
      C_ANGLE(z1.y, *(x++)*unit[1]);
      C_ANGLE(z1.z, *(x++)*unit[2]);
    }
    for (; zz<zn; --zx, ++zy, ++zz) {                  // set up z[k]=e^(ik.r)
      C_RMULT(zy->x, zz->x, z1.x);                        // 3D k-vector
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
      mu += 4;
      h = hvec;
    }
    for (k=kvec; k<nk; ++k) {                                // compute rho(k)
      if (ky!=k->y) {                                   // based on order in
        if (kx!=k->x) cx = z[kx = k->x].x;                // reallocate
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
        double muk = mui[0]*h->x+mui[1]*h->y+mui[2]*h->z; ++h;
        cek->re += zxyz.re*muk; (cek++)->im += zxyz.im*muk;
      }
    }
    ekr = (cvector *) ((char *) memcpy(ekr, z, lbytes)+lbytes);
    ++type;
  }
  MPI_Allreduce(cek_local, cek_global, 2*n, MPI_DOUBLE, MPI_SUM, world);

  delete [] z;
}

/* ---------------------------------------------------------------------- */

void EwaldDisp::compute_force()
{
  kvector *k;
  hvector *h, *nh;
  cvector *z = ekr_local;
  vector sum[EWALD_MAX_NSUMS], mui = COMPLEX_NULL;
  complex *cek, zc, zx = COMPLEX_NULL, zxy = COMPLEX_NULL;
  complex *cek_coul;
  double *f = atom->f[0], *fn = f+3*atom->nlocal, *q = atom->q, *t = NULL;
  double *mu = atom->mu ? atom->mu[0] : NULL;
  const double qscale = force->qqrd2e * scale;
  double *ke, c[EWALD_NFUNCS] = {
    8.0*MY_PI*qscale/volume, 2.0*MY_PI*MY_PIS/(12.0*volume),
    2.0*MY_PI*MY_PIS/(192.0*volume), 8.0*MY_PI*mumurd2e/volume};
  int i, kx, ky, lbytes = (2*nbox+1)*sizeof(cvector), *type = atom->type;
  int func[EWALD_NFUNCS];

  if (atom->torque) t = atom->torque[0];
  memcpy(func, function, EWALD_NFUNCS*sizeof(int));
  memset(sum, 0, EWALD_MAX_NSUMS*sizeof(vector));        // fj = -dE/dr =
  for (; f<fn; f+=3) {                                    //      -i*qj*fac*
    k = kvec;                                         //       Sum[conj(d)-d]
    kx = ky = -1;                                        // d = k*conj(ekj)*ek
    ke = kenergy;
    cek = cek_global;
    memset(sum, 0, EWALD_MAX_NSUMS*sizeof(vector));
    if (func[3]) {
      double di = c[3];
      mui[0] = di*(mu++)[0]; mui[1] = di*(mu++)[0]; mui[2] = di*(mu++)[0];
      mu++;
    }
    for (nh = (h = hvec)+nkvec; h<nh; ++h, ++k) {
      if (ky!=k->y) {                                   // based on order in
        if (kx!=k->x) zx = z[kx = k->x].x;                 // reallocate
        C_RMULT(zxy, z[ky = k->y].y, zx);
      }
      C_CRMULT(zc, z[k->z].z, zxy);
      if (func[0]) {                                        // 1/r
        double im = *(ke++)*(zc.im*cek->re+cek->im*zc.re);
        if (func[3]) cek_coul = cek;
        ++cek;
        sum[0][0] += h->x*im; sum[0][1] += h->y*im; sum[0][2] += h->z*im;
      }
      if (func[1]) {                                        // geometric 1/r^6
        double im = *(ke++)*(zc.im*cek->re+cek->im*zc.re); ++cek;
        sum[1][0] += h->x*im; sum[1][1] += h->y*im; sum[1][2] += h->z*im;
      }
      if (func[2]) {                                        // arithmetic 1/r^6
        double im, c = *(ke++);
        for (i=2; i<9; ++i) {
          im = c*(zc.im*cek->re+cek->im*zc.re); ++cek;
          sum[i][0] += h->x*im; sum[i][1] += h->y*im; sum[i][2] += h->z*im;
        }
      }
      if (func[3]) {                                        // dipole
        double im = *(ke)*(zc.im*cek->re+
            cek->im*zc.re)*(mui[0]*h->x+mui[1]*h->y+mui[2]*h->z);
        double im2 = *(ke)*(zc.re*cek->re-
            cek->im*zc.im);
        sum[9][0] += h->x*im; sum[9][1] += h->y*im; sum[9][2] += h->z*im;
        t[0] += -mui[1]*h->z*im2 + mui[2]*h->y*im2;        // torque
        t[1] += -mui[2]*h->x*im2 + mui[0]*h->z*im2;
        t[2] += -mui[0]*h->y*im2 + mui[1]*h->x*im2;
        if (func[0]) {                                      // charge-dipole
          double qi = *(q)*c[0];
          im = - *(ke)*(zc.re*cek_coul->re -
              cek_coul->im*zc.im)*(mui[0]*h->x+mui[1]*h->y+mui[2]*h->z);
          im += *(ke)*(zc.re*cek->re - cek->im*zc.im)*qi;
          sum[9][0] += h->x*im; sum[9][1] += h->y*im; sum[9][2] += h->z*im;

          im2 =  *(ke)*(zc.re*cek_coul->im + cek_coul->re*zc.im);
          im2 += -*(ke)*(zc.re*cek->im - cek->im*zc.re);
          t[0] += -mui[1]*h->z*im2 + mui[2]*h->y*im2;        // torque
          t[1] += -mui[2]*h->x*im2 + mui[0]*h->z*im2;
          t[2] += -mui[0]*h->y*im2 + mui[1]*h->x*im2;
        }
        ++cek;
        ke++;
      }
    }
    if (func[0]) {                                        // 1/r
      double qi = *(q++)*c[0];
      f[0] -= sum[0][0]*qi; f[1] -= sum[0][1]*qi; f[2] -= sum[0][2]*qi;
    }
    if (func[1]) {                                        // geometric 1/r^6
      double bi = B[*type]*c[1];
      f[0] -= sum[1][0]*bi; f[1] -= sum[1][1]*bi; f[2] -= sum[1][2]*bi;
    }
    if (func[2]) {                                        // arithmetic 1/r^6
      double *bi = B+7*type[0]+7;
      for (i=2; i<9; ++i) {
        double c2 = (--bi)[0]*c[2];
        f[0] -= sum[i][0]*c2; f[1] -= sum[i][1]*c2; f[2] -= sum[i][2]*c2;
      }
    }
    if (func[3]) {                                        // dipole
      f[0] -= sum[9][0]; f[1] -= sum[9][1]; f[2] -= sum[9][2];
    }
    z = (cvector *) ((char *) z+lbytes);
    ++type;
    t += 3;
  }
}

/* ---------------------------------------------------------------------- */

void EwaldDisp::compute_surface()
{
  // assume conducting metal (tinfoil) boundary conditions, so this function is
  // not called because dielectric at the boundary --> infinity, which makes all
  // the terms here zero.

  if (!function[3]) return;
  if (!atom->mu) return;

  vector sum_local = VECTOR_NULL, sum_total;
  memset(sum_local, 0, sizeof(vector));
  double *i, *n, *mu = atom->mu[0];

  for (n = (i = mu) + 4*atom->nlocal; i < n; ++i) {
    sum_local[0] += (i++)[0];
    sum_local[1] += (i++)[0];
    sum_local[2] += (i++)[0];
  }
  MPI_Allreduce(sum_local, sum_total, 3, MPI_DOUBLE, MPI_SUM, world);

  virial_self[3] =
    mumurd2e*(2.0*MY_PI*vec_dot(sum_total,sum_total)/(2.0*dielectric+1)/volume);
  energy_self[3] -= virial_self[3];

  if (!(vflag_atom || eflag_atom)) return;

  double *ei = energy_self_peratom[0]+3;
  double *vi = virial_self_peratom[0]+3;
  double cv = 2.0*mumurd2e*MY_PI/(2.0*dielectric+1)/volume;

  for (i = mu; i < n; i += 4, ei += EWALD_NFUNCS, vi += EWALD_NFUNCS) {
    *vi = cv*(i[0]*sum_total[0]+i[1]*sum_total[1]+i[2]*sum_total[2]);
    *ei -= *vi;
  }
}

/* ---------------------------------------------------------------------- */

void EwaldDisp::compute_energy()
{
  energy = 0.0;
  if (!eflag_global) return;

  complex *cek = cek_global;
  complex *cek_coul;
  double *ke = kenergy;
  const double qscale = force->qqrd2e * scale;
  double c[EWALD_NFUNCS] = {
    4.0*MY_PI*qscale/volume, 2.0*MY_PI*MY_PIS/(24.0*volume),
    2.0*MY_PI*MY_PIS/(192.0*volume), 4.0*MY_PI*mumurd2e/volume};
  double sum[EWALD_NFUNCS];
  int func[EWALD_NFUNCS];

  memcpy(func, function, EWALD_NFUNCS*sizeof(int));
  memset(sum, 0, EWALD_NFUNCS*sizeof(double));                // reset sums
  for (int k=0; k<nkvec; ++k) {                       // sum over k vectors
    if (func[0]) {                                        // 1/r
      sum[0] += *(ke++)*(cek->re*cek->re+cek->im*cek->im);
      if (func[3]) cek_coul = cek;
      ++cek;
    }
    if (func[1]) {                                        // geometric 1/r^6
      sum[1] += *(ke++)*(cek->re*cek->re+cek->im*cek->im); ++cek; }
    if (func[2]) {                                        // arithmetic 1/r^6
      double r =
            (cek[0].re*cek[6].re+cek[0].im*cek[6].im)+
            (cek[1].re*cek[5].re+cek[1].im*cek[5].im)+
            (cek[2].re*cek[4].re+cek[2].im*cek[4].im)+
        0.5*(cek[3].re*cek[3].re+cek[3].im*cek[3].im); cek += 7;
      sum[2] += *(ke++)*r;
    }
    if (func[3]) {                                        // dipole
      sum[3] += *(ke)*(cek->re*cek->re+cek->im*cek->im);
      if (func[0]) {                                      // charge-dipole
        sum[3] += *(ke)*2.0*(cek->re*cek_coul->im - cek->im*cek_coul->re);
      }
      ke++;
      ++cek;
    }
  }
  for (int k=0; k<EWALD_NFUNCS; ++k) energy += c[k]*sum[k]-energy_self[k];
}

/* ---------------------------------------------------------------------- */

void EwaldDisp::compute_energy_peratom()
{
  if (!eflag_atom) return;

  kvector *k;
  hvector *h, *nh;
  cvector *z = ekr_local;
  vector  mui = VECTOR_NULL;
  double sum[EWALD_MAX_NSUMS];
  complex *cek, zc = COMPLEX_NULL, zx = COMPLEX_NULL, zxy = COMPLEX_NULL;
  complex *cek_coul;
  double *q = atom->q;
  double *eatomj = eatom;
  double *mu = atom->mu ? atom->mu[0] : NULL;
  const double qscale = force->qqrd2e * scale;
  double *ke = kenergy;
  double c[EWALD_NFUNCS] = {
      4.0*MY_PI*qscale/volume, 2.0*MY_PI*MY_PIS/(24.0*volume),
      2.0*MY_PI*MY_PIS/(192.0*volume), 4.0*MY_PI*mumurd2e/volume};
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
      double di = c[3];
      mui[0] = di*(mu++)[0]; mui[1] = di*(mu++)[0]; mui[2] = di*(mu++)[0];
      mu++;
    }
    for (nh = (h = hvec)+nkvec; h<nh; ++h, ++k) {
      if (ky!=k->y) {                              // based on order in
        if (kx!=k->x) zx = z[kx = k->x].x;                 // reallocate
        C_RMULT(zxy, z[ky = k->y].y, zx);
      }
      C_CRMULT(zc, z[k->z].z, zxy);
      if (func[0]) {                                        // 1/r
        sum[0] += *(ke++)*(cek->re*zc.re - cek->im*zc.im);
        if (func[3]) cek_coul = cek;
        ++cek;
      }
      if (func[1]) {                                        // geometric 1/r^6
        sum[1] += *(ke++)*(cek->re*zc.re - cek->im*zc.im); ++cek; }
      if (func[2]) {                                        // arithmetic 1/r^6
        double im, c = *(ke++);
        for (i=2; i<9; ++i) {
          im = c*(cek->re*zc.re - cek->im*zc.im); ++cek;
          sum[i] += im;
        }
      }
      if (func[3]) {                                        // dipole
        double muk = (mui[0]*h->x+mui[1]*h->y+mui[2]*h->z);
        sum[9] += *(ke)*(cek->re*zc.re - cek->im*zc.im)*muk;
        if (func[0]) {                                      // charge-dipole
          double qj = *(q)*c[0];
          sum[9] += *(ke)*(cek_coul->im*zc.re + cek_coul->re*zc.im)*muk;
          sum[9] -= *(ke)*(cek->re*zc.im + cek->im*zc.re)*qj;
        }
        ++cek;
        ke++;
      }
    }

    if (func[0]) {                                        // 1/r
      double qj = *(q++)*c[0];
      *eatomj += sum[0]*qj - energy_self_peratom[j][0];
    }
    if (func[1]) {                                        // geometric 1/r^6
      double bj = B[*type]*c[1];
      *eatomj += sum[1]*bj - energy_self_peratom[j][1];
    }
    if (func[2]) {                                        // arithmetic 1/r^6
      double *bj = B+7*type[0]+7;
      for (i=2; i<9; ++i) {
        double c2 = (--bj)[0]*c[2];
        *eatomj += 0.5*sum[i]*c2;
      }
      *eatomj -= energy_self_peratom[j][2];
    }
    if (func[3]) {                                        // dipole
      *eatomj += sum[9] - energy_self_peratom[j][3];
    }
    z = (cvector *) ((char *) z+lbytes);
    ++type;
  }
}

/* ---------------------------------------------------------------------- */

#define swap(a, b) { register double t = a; a= b; b = t; }

void EwaldDisp::compute_virial()
{
  memset(virial, 0, sizeof(shape));
  if (!vflag_global) return;

  complex *cek = cek_global;
  complex *cek_coul;
  double *kv = kvirial;
  const double qscale = force->qqrd2e * scale;
  double c[EWALD_NFUNCS] = {
    4.0*MY_PI*qscale/volume, 2.0*MY_PI*MY_PIS/(24.0*volume),
    2.0*MY_PI*MY_PIS/(192.0*volume), 4.0*MY_PI*mumurd2e/volume};
  shape sum[EWALD_NFUNCS];
  int func[EWALD_NFUNCS];

  memcpy(func, function, EWALD_NFUNCS*sizeof(int));
  memset(sum, 0, EWALD_NFUNCS*sizeof(shape));
  for (int k=0; k<nkvec; ++k) {                      // sum over k vectors
    if (func[0]) {                                         // 1/r
      double r = cek->re*cek->re+cek->im*cek->im;
      if (func[3]) cek_coul = cek;
      ++cek;
      sum[0][0] += *(kv++)*r; sum[0][1] += *(kv++)*r; sum[0][2] += *(kv++)*r;
      sum[0][3] += *(kv++)*r; sum[0][4] += *(kv++)*r; sum[0][5] += *(kv++)*r;
    }
    if (func[1]) {                                        // geometric 1/r^6
      double r = cek->re*cek->re+cek->im*cek->im; ++cek;
      sum[1][0] += *(kv++)*r; sum[1][1] += *(kv++)*r; sum[1][2] += *(kv++)*r;
      sum[1][3] += *(kv++)*r; sum[1][4] += *(kv++)*r; sum[1][5] += *(kv++)*r;
    }
    if (func[2]) {                                        // arithmetic 1/r^6
      double r =
            (cek[0].re*cek[6].re+cek[0].im*cek[6].im)+
            (cek[1].re*cek[5].re+cek[1].im*cek[5].im)+
            (cek[2].re*cek[4].re+cek[2].im*cek[4].im)+
        0.5*(cek[3].re*cek[3].re+cek[3].im*cek[3].im); cek += 7;
      sum[2][0] += *(kv++)*r; sum[2][1] += *(kv++)*r; sum[2][2] += *(kv++)*r;
      sum[2][3] += *(kv++)*r; sum[2][4] += *(kv++)*r; sum[2][5] += *(kv++)*r;
    }
    if (func[3]) {
      double r = cek->re*cek->re+cek->im*cek->im;
      sum[3][0] += *(kv++)*r; sum[3][1] += *(kv++)*r; sum[3][2] += *(kv++)*r;
      sum[3][3] += *(kv++)*r; sum[3][4] += *(kv++)*r; sum[3][5] += *(kv++)*r;
      if (func[0]) {                                      // charge-dipole
        kv -= 6;
        double r = 2.0*(cek->re*cek_coul->im - cek->im*cek_coul->re);
        sum[3][0] += *(kv++)*r; sum[3][1] += *(kv++)*r; sum[3][2] += *(kv++)*r;
        sum[3][3] += *(kv++)*r; sum[3][4] += *(kv++)*r; sum[3][5] += *(kv++)*r;
      }
      ++cek;
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

/* ---------------------------------------------------------------------- */

void EwaldDisp::compute_virial_dipole()
{
  if (!function[3]) return;
  if (!vflag_atom && !vflag_global) return;
  kvector *k;
  hvector *h, *nh;
  cvector *z = ekr_local;
  vector mui = COMPLEX_NULL;
  double sum[6];
  double sum_total[6];
  complex *cek, zc, zx = COMPLEX_NULL, zxy = COMPLEX_NULL;
  complex *cek_coul;
  double *mu = atom->mu ? atom->mu[0] : NULL;
  double *vatomj = NULL;
  if (vflag_atom && vatom) vatomj = vatom[0];
  const double qscale = force->qqrd2e * scale;
  double *ke, c[EWALD_NFUNCS] = {
    8.0*MY_PI*qscale/volume, 2.0*MY_PI*MY_PIS/(12.0*volume),
    2.0*MY_PI*MY_PIS/(192.0*volume), 8.0*MY_PI*mumurd2e/volume};
  int i, kx, ky, lbytes = (2*nbox+1)*sizeof(cvector), *type = atom->type;
  int func[EWALD_NFUNCS];

  memcpy(func, function, EWALD_NFUNCS*sizeof(int));
  memset(&sum[0], 0, 6*sizeof(double));
  memset(&sum_total[0], 0, 6*sizeof(double));
  for (int j = 0; j < atom->nlocal; j++) {
    k = kvec;
    kx = ky = -1;
    ke = kenergy;
    cek = cek_global;
    memset(&sum[0], 0, 6*sizeof(double));
    if (func[3]) {
      double di = c[3];
      mui[0] = di*(mu++)[0]; mui[1] = di*(mu++)[0]; mui[2] = di*(mu++)[0];
      mu++;
    }
    for (nh = (h = hvec)+nkvec; h<nh; ++h, ++k) {
      if (ky!=k->y) {                                   // based on order in
        if (kx!=k->x) zx = z[kx = k->x].x;                 // reallocate
        C_RMULT(zxy, z[ky = k->y].y, zx);
      }
      C_CRMULT(zc, z[k->z].z, zxy);
      double im = 0.0;
      if (func[0]) {                                        // 1/r
        ke++;
        if (func[3]) cek_coul = cek;
        ++cek;
      }
      if (func[1]) {                                        // geometric 1/r^6
        ke++;
        ++cek;
      }
      if (func[2]) {                                        // arithmetic 1/r^6
        ke++;
        for (i=2; i<9; ++i) {
          ++cek;
        }
      }
      if (func[3]) {                                        // dipole
        im = *(ke)*(zc.re*cek->re - cek->im*zc.im);
        if (func[0]) {                                      // charge-dipole
          im += *(ke)*(zc.im*cek_coul->re + cek_coul->im*zc.re);
        }
        sum[0] -= mui[0]*h->x*im;
        sum[1] -= mui[1]*h->y*im;
        sum[2] -= mui[2]*h->z*im;
        sum[3] -= mui[0]*h->y*im;
        sum[4] -= mui[0]*h->z*im;
        sum[5] -= mui[1]*h->z*im;
        ++cek;
        ke++;
      }
    }

    if (vflag_global)
      for (int n = 0; n < 6; n++)
        sum_total[n] -= sum[n];

    if (vflag_atom)
      for (int n = 0; n < 6; n++)
        vatomj[n] -= sum[n];

    z = (cvector *) ((char *) z+lbytes);
    ++type;
    if (vflag_atom) vatomj += 6;
  }

  if (vflag_global) {
    MPI_Allreduce(&sum_total[0],&sum[0],6,MPI_DOUBLE,MPI_SUM,world);
    for (int n = 0; n < 6; n++)
      virial[n] += sum[n];
  }
}

/* ---------------------------------------------------------------------- */

void EwaldDisp::compute_virial_peratom()
{
  if (!vflag_atom) return;

  kvector *k;
  hvector *h, *nh;
  cvector *z = ekr_local;
  vector  mui = VECTOR_NULL;
  complex *cek, zc = COMPLEX_NULL, zx = COMPLEX_NULL, zxy = COMPLEX_NULL;
  complex *cek_coul;
  double *kv;
  double *q = atom->q;
  double *vatomj = vatom ? vatom[0] : NULL;
  double *mu = atom->mu ? atom->mu[0] : NULL;
  const double qscale = force->qqrd2e * scale;
  double c[EWALD_NFUNCS] = {
    4.0*MY_PI*qscale/volume, 2.0*MY_PI*MY_PIS/(24.0*volume),
    2.0*MY_PI*MY_PIS/(192.0*volume), 4.0*MY_PI*mumurd2e/volume};
  shape sum[EWALD_MAX_NSUMS];
  int func[EWALD_NFUNCS];

  memcpy(func, function, EWALD_NFUNCS*sizeof(int));
  int i, kx, ky, lbytes = (2*nbox+1)*sizeof(cvector), *type = atom->type;
  for (int j = 0; j < atom->nlocal; j++) {
    k = kvec;
    kx = ky = -1;
    kv = kvirial;
    cek = cek_global;
    memset(sum, 0, EWALD_MAX_NSUMS*sizeof(shape));
    if (func[3]) {
      double di = c[3];
      mui[0] = di*(mu++)[0]; mui[1] = di*(mu++)[0]; mui[2] = di*(mu++)[0];
      mu++;
    }
    for (nh = (h = hvec)+nkvec; h<nh; ++h, ++k) {
      if (ky!=k->y) {                                // based on order in
          if (kx!=k->x) zx = z[kx = k->x].x;                 // reallocate
          C_RMULT(zxy, z[ky = k->y].y, zx);
      }
      C_CRMULT(zc, z[k->z].z, zxy);
      if (func[0]) {                                        // 1/r
          if (func[3]) cek_coul = cek;
          double r = cek->re*zc.re - cek->im*zc.im; ++cek;
          sum[0][0] += *(kv++)*r;
          sum[0][1] += *(kv++)*r;
          sum[0][2] += *(kv++)*r;
          sum[0][3] += *(kv++)*r;
          sum[0][4] += *(kv++)*r;
          sum[0][5] += *(kv++)*r;
      }
      if (func[1]) {                                        // geometric 1/r^6
          double r = cek->re*zc.re - cek->im*zc.im; ++cek;
          sum[1][0] += *(kv++)*r;
          sum[1][1] += *(kv++)*r;
          sum[1][2] += *(kv++)*r;
          sum[1][3] += *(kv++)*r;
          sum[1][4] += *(kv++)*r;
          sum[1][5] += *(kv++)*r;
      }
      if (func[2]) {                                        // arithmetic 1/r^6
        double r;
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
      if (func[3]) {                                        // dipole
         double muk = (mui[0]*h->x+mui[1]*h->y+mui[2]*h->z);
         double
           r = (cek->re*zc.re - cek->im*zc.im)*muk;
         sum[9][0] += *(kv++)*r;
         sum[9][1] += *(kv++)*r;
         sum[9][2] += *(kv++)*r;
         sum[9][3] += *(kv++)*r;
         sum[9][4] += *(kv++)*r;
         sum[9][5] += *(kv++)*r;
         if (func[0]) {                                      // charge-dipole
           kv -= 6;
           double qj = *(q)*c[0];
           r = (cek_coul->im*zc.re + cek_coul->re*zc.im)*muk;
           r += -(cek->re*zc.im + cek->im*zc.re)*qj;
           sum[9][0] += *(kv++)*r; sum[9][1] += *(kv++)*r; sum[9][2] += *(kv++)*r;
           sum[9][3] += *(kv++)*r; sum[9][4] += *(kv++)*r; sum[9][5] += *(kv++)*r;
         }
         ++cek;
      }
    }

    if (func[0]) {                                        // 1/r
      double qi = *(q++)*c[0];
      for (int n = 0; n < 6; n++) vatomj[n] += sum[0][n]*qi;
    }
    if (func[1]) {                                        // geometric 1/r^6
      double bi = B[*type]*c[1];
      for (int n = 0; n < 6; n++) vatomj[n] += sum[1][n]*bi;
    }
    if (func[2]) {                                        // arithmetic 1/r^6
      double *bj = B+7*type[0]+7;
      for (i=2; i<9; ++i) {
        double c2 = (--bj)[0]*c[2];
        for (int n = 0; n < 6; n++) vatomj[n] += 0.5*sum[i][n]*c2;
      }
    }
    if (func[3]) {                                        // dipole
      for (int n = 0; n < 6; n++) vatomj[n] += sum[9][n];
    }

    for (int k=0; k<EWALD_NFUNCS; ++k) {
      if (func[k]) {
        for (int n = 0; n < 3; n++) vatomj[n] -= virial_self_peratom[j][k];
      }
    }

    z = (cvector *) ((char *) z+lbytes);
    ++type;
    vatomj += 6;
  }
}

/* ----------------------------------------------------------------------
   Slab-geometry correction term to dampen inter-slab interactions between
   periodically repeating slabs.  Yields good approximation to 2D Ewald if
   adequate empty space is left between repeating slabs (J. Chem. Phys.
   111, 3155).  Slabs defined here to be parallel to the xy plane. Also
   extended to non-neutral systems (J. Chem. Phys. 131, 094107).
------------------------------------------------------------------------- */

void EwaldDisp::compute_slabcorr()
{
  // compute local contribution to global dipole moment

  double *q = atom->q;
  double **x = atom->x;
  double zprd = domain->zprd;
  int nlocal = atom->nlocal;

  double dipole = 0.0;
  for (int i = 0; i < nlocal; i++) dipole += q[i]*x[i][2];

  if (function[3] && atom->mu) {
    double **mu = atom->mu;
    for (int i = 0; i < nlocal; i++) dipole += mu[i][2];
  }

  // sum local contributions to get global dipole moment

  double dipole_all;
  MPI_Allreduce(&dipole,&dipole_all,1,MPI_DOUBLE,MPI_SUM,world);

  // need to make non-neutral systems and/or
  //  per-atom energy translationally invariant

  double dipole_r2 = 0.0;
  if (eflag_atom || fabs(qsum) > SMALL) {

    if (function[3] && atom->mu)
      error->all(FLERR,"Cannot (yet) use kspace slab correction with "
        "long-range dipoles and non-neutral systems or per-atom energy");

    for (int i = 0; i < nlocal; i++)
      dipole_r2 += q[i]*x[i][2]*x[i][2];

    // sum local contributions

    double tmp;
    MPI_Allreduce(&dipole_r2,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    dipole_r2 = tmp;
  }

  // compute corrections

  const double e_slabcorr = MY_2PI*(dipole_all*dipole_all -
    qsum*dipole_r2 - qsum*qsum*zprd*zprd/12.0)/volume;
  const double qscale = force->qqrd2e * scale;

  if (eflag_global) energy += qscale * e_slabcorr;

  // per-atom energy

  if (eflag_atom) {
    double efact = qscale * MY_2PI/volume;
    for (int i = 0; i < nlocal; i++)
      eatom[i] += efact * q[i]*(x[i][2]*dipole_all - 0.5*(dipole_r2 +
        qsum*x[i][2]*x[i][2]) - qsum*zprd*zprd/12.0);
  }

  // add on force corrections

  double ffact = qscale * (-4.0*MY_PI/volume);
  double **f = atom->f;

  for (int i = 0; i < nlocal; i++)
    f[i][2] += ffact * q[i]*(dipole_all - qsum*x[i][2]);

  // add on torque corrections

  if (function[3] && atom->mu && atom->torque) {
    double **mu = atom->mu;
    double **torque = atom->torque;
    for (int i = 0; i < nlocal; i++) {
      torque[i][0] += ffact * dipole_all * mu[i][1];
      torque[i][1] += -ffact * dipole_all * mu[i][0];
    }
  }
}

/* ----------------------------------------------------------------------
   Newton solver used to find g_ewald for LJ systems
------------------------------------------------------------------------- */

double EwaldDisp::NewtonSolve(double x, double Rc,
                              bigint natoms, double vol, double b2)
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
    if (x < 0 || x != x) // solver failed
      return -1;
  }
  return -1;
}

/* ----------------------------------------------------------------------
 Calculate f(x)
 ------------------------------------------------------------------------- */

double EwaldDisp::f(double x, double Rc, bigint natoms, double vol, double b2)
{
  double a = Rc*x;
  double f = 0.0;

  if (function[3]) { // dipole
    double rg2 = a*a;
    double rg4 = rg2*rg2;
    double rg6 = rg4*rg2;
    double Cc = 4.0*rg4 + 6.0*rg2 + 3.0;
    double Dc = 8.0*rg6 + 20.0*rg4 + 30.0*rg2 + 15.0;
    f = (b2/(sqrt(vol*powint(x,4)*powint(Rc,9)*natoms)) *
      sqrt(13.0/6.0*Cc*Cc + 2.0/15.0*Dc*Dc - 13.0/15.0*Cc*Dc) *
      exp(-rg2)) - accuracy;
  } else if (function[1] || function[2]) { // LJ
    f = (4.0*MY_PI*b2*powint(x,4)/vol/sqrt((double)natoms)*erfc(a) *
      (6.0*powint(a,-5) + 6.0*powint(a,-3) + 3.0/a + a) - accuracy);
  }

  return f;
}

/* ----------------------------------------------------------------------
 Calculate numerical derivative f'(x)
 ------------------------------------------------------------------------- */

double EwaldDisp::derivf(double x, double Rc,
                         bigint natoms, double vol, double b2)
{
  double h = 0.000001;  //Derivative step-size
  return (f(x + h,Rc,natoms,vol,b2) - f(x,Rc,natoms,vol,b2)) / h;
}
