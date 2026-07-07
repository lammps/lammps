// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

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

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "math_extra.h"
#include "math_special.h"
#include "memory.h"
#include "pair.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace EwaldConst;
using namespace MathConst;
using namespace MathSpecial;
using namespace MathExtra;

static constexpr double SMALL = 0.00001;

struct LAMMPS_NS::complex { double re, im; };
struct LAMMPS_NS::cvector { complex x, y, z; };
struct LAMMPS_NS::hvector { double x, y, z; };
struct LAMMPS_NS::kvector { long x, y, z; };

// complex arithmetic helpers operating on the local complex struct.
// the destination may alias an input, so the multiplies copy first.

static inline void c_rmult(complex &d, const complex &x, const complex &y) {  // d = x*y
  complex t = x;
  d.re = t.re*y.re - t.im*y.im;
  d.im = t.re*y.im + t.im*y.re;
}

static inline void c_crmult(complex &d, const complex &x, const complex &y) { // d = x*conj(y)
  complex t = x;
  d.re = t.re*y.re - t.im*y.im;
  d.im = -t.re*y.im - t.im*y.re;
}

static inline void c_set(complex &d, double re, double im) {
  d.re = re;
  d.im = im;
}

static inline void c_conj(complex &d, const complex &x) {                     // d = conj(x)
  d.re = x.re;
  d.im = -x.im;
}

static inline void c_angle(complex &d, double angle) {                        // d = e^(i*angle)
  d.re = cos(angle);
  d.im = sin(angle);
}

static inline void shape_add(double *dest, const double *src) {                // h_a+h_b
  dest[0] += src[0]; dest[1] += src[1]; dest[2] += src[2];
  dest[3] += src[3]; dest[4] += src[4]; dest[5] += src[5]; }

static inline void shape_subtr(double *dest, const double *src) {                // h_a-h_b
  dest[0] -= src[0]; dest[1] -= src[1]; dest[2] -= src[2];
  dest[3] -= src[3]; dest[4] -= src[4]; dest[5] -= src[5]; }

static inline double shape_det(double *s) {
  return s[0]*s[1]*s[2]; }

static inline void shape_scalar_mult(double *dest, double f) {                // f*h
  dest[0] *= f; dest[1] *= f; dest[2] *= f;
  dest[3] *= f; dest[4] *= f; dest[5] *= f; }

/* ---------------------------------------------------------------------- */

EwaldDisp::EwaldDisp(LAMMPS *lmp) : KSpace(lmp),
  kenergy(nullptr), kvirial(nullptr), energy_self_peratom(nullptr), virial_self_peratom(nullptr),
  ekr(nullptr), hvec(nullptr), kvec(nullptr), B(nullptr), sfac(nullptr), sfac_all(nullptr)
{
  ewaldflag = dispersionflag = dipoleflag = 1;

  memset(termflag, 0, EWALD_NTERMS*sizeof(int));
  first_output = 0;
  nmax = 0;
  q2 = 0;
  b2 = 0;
  M2 = 0;
}

/* ---------------------------------------------------------------------- */

void EwaldDisp::settings(int narg, char **arg)
{
  if (narg!=1) error->all(FLERR,"Illegal kspace_style {} command", force->kspace_style);
  accuracy_relative = fabs(utils::numeric(FLERR,arg[0],false,lmp));
  if (accuracy_relative > 1.0)
    error->all(FLERR, "Invalid relative accuracy {:g} for kspace_style {}",
               accuracy_relative, force->kspace_style);
}

/* ---------------------------------------------------------------------- */

EwaldDisp::~EwaldDisp()
{
  deallocate();
  deallocate_peratom();
  memory->destroy(ekr);
  memory->destroy(B);
}

/* --------------------------------------------------------------------- */

void EwaldDisp::init()
{
  kcount = kcount_max = nevec = nevec_max = 0;
  nterms = nsums = coeff_sums_flag = 0;
  kmax = -1;

  if (!comm->me) utils::logmesg(lmp,"EwaldDisp initialization ...\n");

  triclinic_check();
  if (domain->dimension == 2)
    error->all(FLERR,"Cannot use EwaldDisp with 2d simulation");
  if (slabflag == 0 && domain->nonperiodic > 0)
    error->all(FLERR,"Cannot use non-periodic boundaries with EwaldDisp");
  if (slabflag == 1) {
    if (domain->xperiodic != 1 || domain->yperiodic != 1 ||
        domain->boundary[2][0] != 1 || domain->boundary[2][1] != 1)
      error->all(FLERR,"Incorrect boundaries with slab EwaldDisp");
    if (domain->triclinic && (domain->yz != 0.0 || domain->xz != 0.0))
      error->all(FLERR,"Triclinic slab (EW3DC) correction requires xz = yz = 0 "
                 "(the slab normal must be the z axis); xy tilt is allowed");
  }

  scale = 1.0;
  mumurd2e = force->qqrd2e;
  dielectric = force->dielectric;

  int tmp;
  Pair *pair = force->pair;
  int *ptr = pair ? (int *) pair->extract("ewald_order",tmp) : nullptr;
  double *cutoff = pair ? (double *) pair->extract("cut_coul",tmp) : nullptr;
  if (!(ptr||cutoff))
    error->all(FLERR,"KSpace style is incompatible with Pair style");
  int ewald_order = ptr ? *ptr : EWALD_COUL;
  int ewald_mix = ptr ? *((int *) pair->extract("ewald_mix",tmp)) : Pair::GEOMETRIC;

  // transcribe the interaction orders handled by the pair style into termflag

  static constexpr int nsums_term[EWALD_NTERMS] = {1, 1, 7, 1};  // sum channels per term

  memset(termflag, 0, EWALD_NTERMS*sizeof(int));
  for (int i = 0; i <= EWALD_MAXORDER; ++i)
    if (ewald_order & (1 << i)) {
      int t = -1;
      switch (i) {
        case 1:
          t = TERM_COUL;
          break;
        case 3:
          t = TERM_DIPOLE;
          break;
        case 6:
          // honor kspace_modify mix/disp (base-class mixflag): 1 = force geometric,
          // 2 = none, 0 = follow the pair style's rule.  ewald/disp has no eigenvalue
          // splitting (the "none" path of pppm/disp), so mix/disp none is rejected.
          if (mixflag == 2)
            error->all(FLERR,"kspace_modify mix/disp none is not supported by "
                             "kspace_style ewald/disp");
          if (mixflag == 1 || ewald_mix == Pair::GEOMETRIC) t = TERM_DISP_GEOM;
          else if (ewald_mix == Pair::ARITHMETIC) t = TERM_DISP_ARITH;
          else error->all(FLERR,"Unsupported mixing rule in kspace_style ewald/disp");
          break;
        default:
          error->all(FLERR,"Unsupported order in kspace_style ewald/disp");
      }
      termflag[t] = 1;
      nterms++;
      nsums += nsums_term[t];
    }

  if (!gewaldflag) g_ewald = 1.0;
  if (!gewaldflag_6) g_ewald_6 = 1.0;
  pair->init();  // so B is defined
  init_coeffs();
  init_coeff_sums();
  if (termflag[TERM_COUL]) qsum_qsq();
  else qsqsum = qsum = 0.0;
  natoms_original = atom->natoms;
  if (!gewaldflag) g_ewald = 0.0;
  if (!gewaldflag_6) g_ewald_6 = 0.0;

  // turn off coulombic if no charge

  if (termflag[TERM_COUL] && qsqsum == 0.0) {
    termflag[TERM_COUL] = 0;
    nterms -= 1;
    nsums -= 1;
  }

  double bsbsum = 0.0;
  M2 = 0.0;
  if (termflag[TERM_DISP_GEOM]) bsbsum = coeff_sum[1].x2;
  if (termflag[TERM_DISP_ARITH]) bsbsum = coeff_sum[2].x2;

  if (termflag[TERM_DIPOLE]) M2 = coeff_sum[9].x2;

  if (termflag[TERM_DIPOLE] && strcmp(update->unit_style,"electron") == 0)
    error->all(FLERR,"Cannot (yet) use 'electron' units with dipoles");

  if (qsqsum == 0.0 && bsbsum == 0.0 && M2 == 0.0)
      error->all(FLERR,"Cannot use Ewald/disp solver on system without "
                 "charged, dipole, or LJ particles");
  if (fabs(qsum) > SMALL && comm->me == 0)
    error->warning(FLERR,"System is not charge neutral, net charge = {:.8g}" + utils::errorurl(29),qsum);

  if (!termflag[TERM_DISP_GEOM] && !termflag[TERM_DISP_ARITH]) dispersionflag = 0;
  if (!termflag[TERM_DIPOLE]) dipoleflag = 0;

  // compute two charge force

  two_charge();

  // extract short-range Coulombic cutoff from pair style

  pair_check();

  // set accuracy (force units) from accuracy_relative or accuracy_absolute

  if (accuracy_absolute >= 0.0) accuracy = accuracy_absolute;
  else accuracy = accuracy_relative * two_charge_force;

  // setup K-space resolution

  q2 = qsqsum * force->qqrd2e;
  M2 *= mumurd2e;
  b2 = bsbsum;
  bigint natoms = atom->natoms;

  if (!gewaldflag) {
    if (termflag[TERM_COUL]) {
      if (accuracy <= 0.0)
        error->all(FLERR,"KSpace accuracy must be > 0");
      if (q2 == 0.0)
        error->all(FLERR,"Must use 'kspace_modify gewald' for uncharged system");
      g_ewald = accuracy*sqrt(natoms*(*cutoff)*shape_det(domain->h)) / (2.0*q2);
      if (g_ewald >= 1.0) g_ewald = (1.35 - 0.15*log(accuracy))/(*cutoff);
      else g_ewald = sqrt(-log(g_ewald)) / (*cutoff);
    } else if (termflag[TERM_DIPOLE]) {
      // use the old method as initial guess for the Newton solver
      g_ewald = (1.35 - 0.15*log(accuracy))/ *cutoff;
      double g_ewald_new =
        NewtonSolve(g_ewald,(*cutoff),natoms,shape_det(domain->h),M2);
      if (g_ewald_new > 0.0) g_ewald = g_ewald_new;
      else error->warning(FLERR,"Ewald/disp Newton solver failed, "
                          "using old method to estimate g_ewald");
    } else if (termflag[TERM_DISP_GEOM] || termflag[TERM_DISP_ARITH]) {
      // use the old method as initial guess for the Newton solver
      g_ewald = (1.35 - 0.15*log(accuracy))/ *cutoff;

      double g_ewald_new =
        NewtonSolve(g_ewald,(*cutoff),natoms,shape_det(domain->h),b2);
      if (g_ewald_new > 0.0) g_ewald = g_ewald_new;
      else error->warning(FLERR,"Ewald/disp Newton solver failed, "
                          "using old method to estimate g_ewald");
    }
  }

  if (comm->me == 0)
    utils::logmesg(lmp,"  G vector = {:.8g},   accuracy = {:.8g}\n",
                   g_ewald,accuracy);

  // apply coulomb g_ewald to dispersion unless it is explicitly set

  if (!gewaldflag_6) g_ewald_6 = g_ewald;
  deallocate_peratom();
  peratom_allocate_flag = 0;
}

/* ----------------------------------------------------------------------
   adjust EwaldDisp coeffs, called initially and whenever volume has changed
------------------------------------------------------------------------- */

void EwaldDisp::setup()
{
  volume = shape_det(domain->h)*slab_volfactor;
  memcpy(unitk, domain->h_inv, 6*sizeof(double));
  shape_scalar_mult(unitk, 2.0*MY_PI);
  unitk[2] /= slab_volfactor;

  if (accuracy >= 1) {
    kmax = 0;
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
  kmax = MAX(kxmax,kymax);
  kmax = MAX(kmax,kzmax);
  double gsqxmx = unitk[0]*unitk[0]*kxmax*kxmax;
  double gsqymx = unitk[1]*unitk[1]*kymax*kymax;
  double gsqzmx = unitk[2]*unitk[2]*kzmax*kzmax;
  gsqmx = MAX(gsqxmx,gsqymx);
  gsqmx = MAX(gsqmx,gsqzmx);
  gsqmx *= 1.00001;

  reallocate();
  coeffs();
  init_coeffs();
  init_coeff_sums();
  init_self();

  if (!(first_output||comm->me)) {
    first_output = 1;
    utils::logmesg(lmp,"  vectors: kmax = {}, kcount = {}\n", kmax, kcount);
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

  // Lennard-Jones reciprocal-space RMS force error (uses the dispersion
  // splitting parameter g_ewald_6).  Kolafa-Perram / Deserno-Holm analysis
  // extended to r^-6: with the large-h asymptotic of the dispersion influence
  // coefficient B(h) ~ 24 g^5 exp(-h^2/(4g^2))/h^2 and h_max = 2 pi km/prd,
  //   dF_recip ~ (b2 g^6 / sqrt(N V)) sqrt(prd/(2 pi km)) exp(-pi^2 km^2/(g^2 prd^2)).

  double g2_6 = g_ewald_6*g_ewald_6;
  double g6 = g2_6*g2_6*g2_6;

  value += b2*g6/sqrt((double)natoms*volume) *
    sqrt(prd/(2.0*MY_PI*km)) *
    exp(-MY_PI*MY_PI*km*km/(g2_6*prd*prd));

  // dipole

  value += 8.0*MY_PI*M2/volume*g_ewald *
    sqrt(2.0*MY_PI*km*km*km/(15.0*natoms)) *
    exp(-pow(MY_PI*km/(g_ewald*prd),2.0));

  return value;
}

/* ----------------------------------------------------------------------
   build the list of k-vectors within the cutoff and the per-k-vector
   arrays sized by it
------------------------------------------------------------------------- */

void EwaldDisp::reallocate()
{
  int *kflag;
  memory->create(kflag, (kmax+1)*(2*kmax+1)*(2*kmax+1), "ewald/disp:kflag");

  // flag which k-vectors are within the cutoff, exploiting inversion
  // symmetry: keep the half-space ix > 0, plus iy > 0 in the ix == 0
  // plane, plus iz > 0 on the ix == iy == 0 axis

  kcount = 0;
  int n = 0;
  for (int ix = 0; ix <= kmax; ix++) {
    for (int iy = -kmax; iy <= kmax; iy++) {
      for (int iz = -kmax; iz <= kmax; iz++) {
        if ((ix == 0) && (iy == 0) && (iz == 0)) kflag[n] = 0;
        else if ((ix == 0) && (iy < 0)) kflag[n] = 0;
        else if ((ix == 0) && (iy == 0) && (iz < 0)) kflag[n] = 0;
        else {
          double h[3];
          h[0] = unitk[0]*ix;
          h[1] = unitk[5]*ix + unitk[1]*iy;
          h[2] = unitk[4]*ix + unitk[3]*iy + unitk[2]*iz;
          kflag[n] = (h[0]*h[0] + h[1]*h[1] + h[2]*h[2] <= gsqmx) ? 1 : 0;
          if (kflag[n]) kcount++;
        }
        n++;
      }
    }
  }

  if (kcount > kcount_max) {
    deallocate();
    memory->create(hvec, kcount, "ewald/disp:hvec");
    memory->create(kvec, kcount, "ewald/disp:kvec");
    memory->create(kenergy, kcount*nterms, "ewald/disp:kenergy");
    memory->create(kvirial, 6*kcount*nterms, "ewald/disp:kvirial");
    memory->create(sfac, kcount*nsums, "ewald/disp:sfac");
    memory->create(sfac_all, kcount*nsums, "ewald/disp:sfac_all");
    kcount_max = kcount;
  }

  // store the flagged k-vectors: cartesian components in hvec, integer
  // indices shifted by kmax (i.e. always >= 0) in kvec.  the resulting
  // list is ordered x-major, then y, then z; eik_dot_r() and the
  // compute_xxx() loops rely on that order.

  n = 0;
  int k = 0;
  for (int ix = 0; ix <= kmax; ix++) {
    for (int iy = -kmax; iy <= kmax; iy++) {
      for (int iz = -kmax; iz <= kmax; iz++) {
        if (kflag[n]) {
          hvec[k].x = unitk[0]*ix;
          hvec[k].y = unitk[5]*ix + unitk[1]*iy;
          hvec[k].z = unitk[4]*ix + unitk[3]*iy + unitk[2]*iz;
          kvec[k].x = ix + kmax;
          kvec[k].y = iy + kmax;
          kvec[k].z = iz + kmax;
          k++;
        }
        n++;
      }
    }
  }

  memory->destroy(kflag);
}

/* ---------------------------------------------------------------------- */

void EwaldDisp::reallocate_atoms()
{
  if ((eflag_atom || vflag_atom) && atom->nmax > nmax) {
    deallocate_peratom();
    allocate_peratom();
    nmax = atom->nmax;
  }

  nevec = atom->nmax*(2*kmax+1);
  if (nevec <= nevec_max) return;
  memory->destroy(ekr);
  memory->create(ekr, nevec, "ewald/disp:ekr");
  nevec_max = nevec;
}

/* ---------------------------------------------------------------------- */

void EwaldDisp::allocate_peratom()
{
  memory->create(energy_self_peratom,
      atom->nmax,EWALD_NTERMS,"ewald/disp:energy_self_peratom");
  memory->create(virial_self_peratom,
      atom->nmax,EWALD_NTERMS,"ewald/disp:virial_self_peratom");
}

/* ---------------------------------------------------------------------- */

void EwaldDisp::deallocate_peratom()
{
  memory->destroy(energy_self_peratom);
  memory->destroy(virial_self_peratom);
}

/* ---------------------------------------------------------------------- */

void EwaldDisp::deallocate()
{
  memory->destroy(hvec);
  memory->destroy(kvec);
  memory->destroy(kenergy);
  memory->destroy(kvirial);
  memory->destroy(sfac);
  memory->destroy(sfac_all);
}

/* ----------------------------------------------------------------------
   energy and virial prefactors for each k-vector and enabled term
------------------------------------------------------------------------- */

void EwaldDisp::coeffs()
{
  const double eta2 = 0.25/(g_ewald*g_ewald);
  const double eta2_6 = 0.25/(g_ewald_6*g_ewald_6);    // dispersion uses its own g_ewald_6
  int ne = 0, nv = 0;                          // running indices into kenergy/kvirial

  for (int k = 0; k < kcount; k++) {
    const double h[3] = {hvec[k].x, hvec[k].y, hvec[k].z};
    const double h2 = dot3(h, h);
    const double bk2 = h2*eta2;
    const double expb2 = exp(-bk2);

    if (termflag[TERM_COUL]) {                          // qi*qj/r coeffs
      const double c1 = expb2/h2;
      const double c2 = 2.0*c1*(1.0 + bk2)/h2;
      kenergy[ne++] = c1;
      kvirial[nv++] = c1 - c2*h[0]*h[0];                // lammps convention
      kvirial[nv++] = c1 - c2*h[1]*h[1];                // instead of voigt
      kvirial[nv++] = c1 - c2*h[2]*h[2];
      kvirial[nv++] = -c2*h[1]*h[0];
      kvirial[nv++] = -c2*h[2]*h[0];
      kvirial[nv++] = -c2*h[2]*h[1];
    }
    if (termflag[TERM_DISP_GEOM] || termflag[TERM_DISP_ARITH]) {   // -Bij/r^6 coeffs
      const double bk2_6 = h2*eta2_6;                   // dispersion split width (g_ewald_6)
      const double expb2_6 = exp(-bk2_6);
      const double b1 = sqrt(bk2_6);                    // minus sign folded
      const double h1 = sqrt(h2);                       // into constants
      const double ce = MY_PIS*erfc(b1);
      const double c1 = -h1*h2*(ce + (0.5/bk2_6 - 1.0)*expb2_6/b1);
      const double c2 = 3.0*h1*(ce - expb2_6/b1);
      kenergy[ne++] = c1;
      kvirial[nv++] = c1 - c2*h[0]*h[0];                // lammps convention
      kvirial[nv++] = c1 - c2*h[1]*h[1];                // instead of voigt
      kvirial[nv++] = c1 - c2*h[2]*h[2];
      kvirial[nv++] = -c2*h[1]*h[0];
      kvirial[nv++] = -c2*h[2]*h[0];
      kvirial[nv++] = -c2*h[2]*h[1];
    }
    if (termflag[TERM_DIPOLE]) {                        // dipole coeffs
      const double c1 = expb2/h2;
      const double c2 = 2.0*c1*(1.0 + bk2)/h2;
      kenergy[ne++] = c1;
      kvirial[nv++] = c1 - c2*h[0]*h[0];                // lammps convention
      kvirial[nv++] = c1 - c2*h[1]*h[1];                // instead of voigt
      kvirial[nv++] = c1 - c2*h[2]*h[2];
      kvirial[nv++] = -c2*h[1]*h[0];
      kvirial[nv++] = -c2*h[2]*h[0];
      kvirial[nv++] = -c2*h[2]*h[1];
    }
  }
}

/* ----------------------------------------------------------------------
   per-type dispersion coefficients, extracted from the pair style
------------------------------------------------------------------------- */

void EwaldDisp::init_coeffs()
{
  int tmp;
  const int n = atom->ntypes;

  if (termflag[TERM_DISP_GEOM]) {                       // geometric 1/r^6
    auto **b = (double **) force->pair->extract("B",tmp);
    memory->destroy(B);
    memory->create(B, n+1, "ewald/disp:B");
    B[0] = 0.0;
    for (int i = 1; i <= n; ++i) B[i] = sqrt(fabs(b[i][i]));
  }
  if (termflag[TERM_DISP_ARITH]) {                      // arithmetic 1/r^6
    auto **epsilon = (double **) force->pair->extract("epsilon",tmp);
    auto **sigma = (double **) force->pair->extract("sigma",tmp);
    if (!(epsilon && sigma))
      error->all(FLERR,"Epsilon or sigma reference not set by pair style in ewald/disp");
    memory->destroy(B);
    memory->create(B, 7*n+7, "ewald/disp:B");

    // the seven per-type coefficients of the binomial expansion of
    // (0.5*(sigma_i+sigma_j))^6: sqrt(eps_i)*c[j]*sigma_i^j, j = 0..6

    const double c[7] = {
      1.0, sqrt(6.0), sqrt(15.0), sqrt(20.0), sqrt(15.0), sqrt(6.0), 1.0};

    for (int j = 0; j < 7; ++j) B[j] = 0.0;
    for (int i = 1; i <= n; ++i) {
      const double eps_i = sqrt(epsilon[i][i]);
      const double sigma_i = sigma[i][i];
      double sigma_n = 1.0;
      for (int j = 0; j < 7; ++j) {
        B[7*i+j] = sigma_n*eps_i*c[j];
        sigma_n *= sigma_i;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   global sums of the per-type dispersion coefficients and dipole moments
------------------------------------------------------------------------- */

void EwaldDisp::init_coeff_sums()
{
  if (coeff_sums_flag) return;                          // calculated only once
  coeff_sums_flag = 1;

  Sum sum_local[EWALD_MAX_NSUMS];

  memset(sum_local, 0, EWALD_MAX_NSUMS*sizeof(Sum));
  memset(coeff_sum, 0, EWALD_MAX_NSUMS*sizeof(Sum));

  // sum_local[0] stays zero: the charge sums are computed by qsum_qsq()

  const int nlocal = atom->nlocal;
  const int *type = atom->type;

  if (termflag[TERM_DISP_GEOM]) {                       // geometric 1/r^6
    for (int i = 0; i < nlocal; i++) {
      sum_local[1].x += B[type[i]];
      sum_local[1].x2 += B[type[i]]*B[type[i]];
    }
  }
  if (termflag[TERM_DISP_ARITH]) {                      // arithmetic 1/r^6
    for (int i = 0; i < nlocal; i++) {
      const double *bi = B + 7*type[i];
      sum_local[2].x2 += bi[0]*bi[6];
      for (int k = 2; k < 9; k++) sum_local[k].x += bi[k-2];
    }
  }
  if (termflag[TERM_DIPOLE] && atom->mu) {              // dipole
    double **mu = atom->mu;
    for (int i = 0; i < nlocal; i++) sum_local[9].x2 += mu[i][3]*mu[i][3];
  }
  MPI_Allreduce(sum_local, coeff_sum, 2*EWALD_MAX_NSUMS, MPI_DOUBLE, MPI_SUM, world);
}

/* ---------------------------------------------------------------------- */

void EwaldDisp::init_self()
{
  double g1 = g_ewald, g2 = g1*g1, g3 = g1*g2;
  double g3_6 = g_ewald_6*g_ewald_6*g_ewald_6;    // dispersion uses its own g_ewald_6
  const double qscale = force->qqrd2e * scale;

  memset(energy_self, 0, EWALD_NTERMS*sizeof(double));  // self energy
  memset(virial_self, 0, EWALD_NTERMS*sizeof(double));

  if (termflag[TERM_COUL]) {                            // 1/r
    virial_self[TERM_COUL] = -0.5*MY_PI*qscale/(g2*volume)*qsum*qsum;
    energy_self[TERM_COUL] = qsqsum*qscale*g1/MY_PIS-virial_self[TERM_COUL];
  }
  if (termflag[TERM_DISP_GEOM]) {                       // geometric 1/r^6
    virial_self[TERM_DISP_GEOM] =
      MY_PI*MY_PIS*g3_6/(6.0*volume)*coeff_sum[1].x*coeff_sum[1].x;
    energy_self[TERM_DISP_GEOM] =
      -coeff_sum[1].x2*g3_6*g3_6/12.0+virial_self[TERM_DISP_GEOM];
  }
  if (termflag[TERM_DISP_ARITH]) {                      // arithmetic 1/r^6
    virial_self[TERM_DISP_ARITH] =
      MY_PI*MY_PIS*g3_6/(48.0*volume)*(coeff_sum[2].x*coeff_sum[8].x+
        coeff_sum[3].x*coeff_sum[7].x+coeff_sum[4].x*coeff_sum[6].x+
        0.5*coeff_sum[5].x*coeff_sum[5].x);
    energy_self[TERM_DISP_ARITH] =
      -coeff_sum[2].x2*g3_6*g3_6/3.0+virial_self[TERM_DISP_ARITH];
  }
  if (termflag[TERM_DIPOLE]) {                          // dipole
    virial_self[TERM_DIPOLE] = 0;
    energy_self[TERM_DIPOLE] =
      coeff_sum[9].x2*mumurd2e*2.0*g3/3.0/MY_PIS-virial_self[TERM_DIPOLE];
  }
}

/* ---------------------------------------------------------------------- */

void EwaldDisp::init_self_peratom()
{
  if (!(vflag_atom || eflag_atom)) return;

  double g1 = g_ewald, g2 = g1*g1, g3 = g1*g2;
  double g3_6 = g_ewald_6*g_ewald_6*g_ewald_6;    // dispersion uses its own g_ewald_6
  const double qscale = force->qqrd2e * scale;
  const int nlocal = atom->nlocal;

  memset(energy_self_peratom[0], 0, EWALD_NTERMS*nlocal*sizeof(double));
  memset(virial_self_peratom[0], 0, EWALD_NTERMS*nlocal*sizeof(double));

  if (termflag[TERM_COUL]) {                            // 1/r
    const double *q = atom->q;
    const double ce = qscale*g1/MY_PIS;
    const double cv = -0.5*MY_PI*qscale/(g2*volume);
    for (int i = 0; i < nlocal; i++) {
      virial_self_peratom[i][TERM_COUL] = cv*q[i]*qsum;
      energy_self_peratom[i][TERM_COUL] =
        ce*q[i]*q[i]-virial_self_peratom[i][TERM_COUL];
    }
  }
  if (termflag[TERM_DISP_GEOM]) {                       // geometric 1/r^6
    const int *type = atom->type;
    const double ce = -g3_6*g3_6/12.0;
    const double cv = MY_PI*MY_PIS*g3_6/(6.0*volume);
    for (int i = 0; i < nlocal; i++) {
      const double b = B[type[i]];
      virial_self_peratom[i][TERM_DISP_GEOM] = cv*b*coeff_sum[1].x;
      energy_self_peratom[i][TERM_DISP_GEOM] =
        ce*b*b+virial_self_peratom[i][TERM_DISP_GEOM];
    }
  }
  if (termflag[TERM_DISP_ARITH]) {                      // arithmetic 1/r^6
    const int *type = atom->type;
    const double ce = -g3_6*g3_6/3.0;
    const double cv = 0.5*MY_PI*MY_PIS*g3_6/(48.0*volume);
    for (int i = 0; i < nlocal; i++) {
      const double *bi = B + 7*type[i];
      for (int k = 2; k < 9; k++)
        virial_self_peratom[i][TERM_DISP_ARITH] += cv*coeff_sum[k].x*bi[8-k];
      energy_self_peratom[i][TERM_DISP_ARITH] =
        ce*bi[0]*bi[6]+virial_self_peratom[i][TERM_DISP_ARITH];
    }
  }
  if (termflag[TERM_DIPOLE] && atom->mu) {              // dipole
    double **mu = atom->mu;
    const double ce = mumurd2e*2.0*g3/3.0/MY_PIS;
    for (int i = 0; i < nlocal; i++) {
      virial_self_peratom[i][TERM_DIPOLE] = 0.0;
      energy_self_peratom[i][TERM_DIPOLE] = ce*mu[i][3]*mu[i][3];
    }
  }
}

/* ----------------------------------------------------------------------
   compute the EwaldDisp long-range force, energy, virial
------------------------------------------------------------------------- */

void EwaldDisp::compute(int eflag, int vflag)
{
  if (!kmax) return;

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

  // tinfoil (conducting metal) boundary conditions are assumed, so there
  // is no surface dipole correction term

  eik_dot_r();
  compute_force();

  // update qsum and qsqsum, if atom count has changed and energy needed

  if ((eflag_global || eflag_atom) && atom->natoms != natoms_original) {
    if (termflag[TERM_COUL]) qsum_qsq();
    natoms_original = atom->natoms;
  }

  compute_energy();
  compute_energy_peratom();
  compute_virial();
  compute_virial_dipole();
  compute_virial_peratom();

  if (slabflag) slabcorr();
}

/* ----------------------------------------------------------------------
   compute the per-atom powers e^(i*k*r) and the structure factors
------------------------------------------------------------------------- */

void EwaldDisp::eik_dot_r()
{
  const int nz = 2*kmax+1;             // per-atom e^(ikr) table entries
  const int nlocal = atom->nlocal;
  const int tri = domain->triclinic;
  double **x = atom->x;
  const double *q = atom->q;
  const int *type = atom->type;
  double **mu = atom->mu;

  memset(sfac, 0, (sizeof(complex)*kcount*nsums));        // reset sums

  for (int i = 0; i < nlocal; i++) {

    // z[kmax+j] with j = -kmax..kmax holds e^(i*j*k0.r) per dimension for
    // this atom, built by the recursion z[j] = z[j-1]*z[1]; negative j
    // entries are the conjugates (only needed for y and z since the
    // k-vector list has x >= 0)

    cvector *z = &ekr[i*nz];
    cvector z1;

    c_set(z[kmax].x, 1, 0);                             // j = 0
    c_set(z[kmax].y, 1, 0);
    c_set(z[kmax].z, 1, 0);
    if (tri) {                                          // j = 1, triclinic
      c_angle(z1.x, unitk[0]*x[i][0]+unitk[5]*x[i][1]+unitk[4]*x[i][2]);
      c_angle(z1.y, unitk[1]*x[i][1]+unitk[3]*x[i][2]);
      c_angle(z1.z, x[i][2]*unitk[2]);
    } else {                                            // j = 1, orthogonal
      c_angle(z1.x, x[i][0]*unitk[0]);
      c_angle(z1.y, x[i][1]*unitk[1]);
      c_angle(z1.z, x[i][2]*unitk[2]);
    }
    for (int j = 1; j <= kmax; j++) {
      c_rmult(z[kmax+j].x, z[kmax+j-1].x, z1.x);
      c_rmult(z[kmax+j].y, z[kmax+j-1].y, z1.y);
      c_conj(z[kmax-j].y, z[kmax+j].y);
      c_rmult(z[kmax+j].z, z[kmax+j-1].z, z1.z);
      c_conj(z[kmax-j].z, z[kmax+j].z);
    }

    double qi = 0.0, bi = 0.0;
    const double *ci = nullptr, *mui = nullptr;
    if (termflag[TERM_COUL]) qi = q[i];
    if (termflag[TERM_DISP_GEOM]) bi = B[type[i]];
    if (termflag[TERM_DISP_ARITH]) ci = &B[7*type[i]];
    if (termflag[TERM_DIPOLE]) mui = mu[i];

    int n = 0;                         // running channel index into sfac
    int kxold = -1, kyold = -1;
    complex cx = {0, 0}, zxy = {0, 0}, zxyz;

    for (int k = 0; k < kcount; k++) {                  // compute rho(k)
      const int kx = kvec[k].x, ky = kvec[k].y, kz = kvec[k].z;

      // kvec is ordered x-major, then y, then z (see reallocate()), so
      // e^(i*kx*x) and e^(i*(kx*x+ky*y)) can be reused as long as kx
      // and ky do not change

      if (ky != kyold) {
        if (kx != kxold) { cx = z[kx].x; kxold = kx; }
        c_rmult(zxy, z[ky].y, cx);
        kyold = ky;
      }
      c_rmult(zxyz, z[kz].z, zxy);

      if (termflag[TERM_COUL]) {
        sfac[n].re += zxyz.re*qi;
        sfac[n].im += zxyz.im*qi;
        n++;
      }
      if (termflag[TERM_DISP_GEOM]) {
        sfac[n].re += zxyz.re*bi;
        sfac[n].im += zxyz.im*bi;
        n++;
      }
      if (termflag[TERM_DISP_ARITH]) {
        for (int m = 0; m < 7; m++) {
          sfac[n].re += zxyz.re*ci[m];
          sfac[n].im += zxyz.im*ci[m];
          n++;
        }
      }
      if (termflag[TERM_DIPOLE]) {
        const double muk = mui[0]*hvec[k].x+mui[1]*hvec[k].y+mui[2]*hvec[k].z;
        sfac[n].re += zxyz.re*muk;
        sfac[n].im += zxyz.im*muk;
        n++;
      }
    }
  }
  MPI_Allreduce(sfac, sfac_all, 2*kcount*nsums, MPI_DOUBLE, MPI_SUM, world);
}

/* ----------------------------------------------------------------------
   convert the structure factors into forces (and torques for dipoles):
   fj = -dE/dr = -i*qj*fac*Sum[conj(d)-d] with d = k*conj(ekj)*ek
------------------------------------------------------------------------- */

void EwaldDisp::compute_force()
{
  const int nz = 2*kmax+1;
  const int nlocal = atom->nlocal;
  double **f = atom->f;
  double **tq = atom->torque;
  const double *q = atom->q;
  const int *type = atom->type;
  double **mu = atom->mu;
  const double qscale = force->qqrd2e * scale;
  const double c[EWALD_NTERMS] = {
    8.0*MY_PI*qscale/volume, 2.0*MY_PI*MY_PIS/(12.0*volume),
    2.0*MY_PI*MY_PIS/(192.0*volume), 8.0*MY_PI*mumurd2e/volume};
  double eksum[EWALD_MAX_NSUMS][3];
  double mui[3] = {0.0,0.0,0.0};

  for (int i = 0; i < nlocal; i++) {
    const cvector *z = &ekr[i*nz];
    int n = 0, ne = 0;                 // running indices into sfac_all/kenergy
    int ncoul = 0;                     // coulomb channel index for charge-dipole
    int kxold = -1, kyold = -1;
    complex zx = {0, 0}, zxy = {0, 0}, zc;

    memset(eksum, 0, EWALD_MAX_NSUMS*3*sizeof(double));
    if (termflag[TERM_DIPOLE]) {
      mui[0] = c[TERM_DIPOLE]*mu[i][0];
      mui[1] = c[TERM_DIPOLE]*mu[i][1];
      mui[2] = c[TERM_DIPOLE]*mu[i][2];
    }

    for (int k = 0; k < kcount; k++) {
      const int kx = kvec[k].x, ky = kvec[k].y, kz = kvec[k].z;
      const double hx = hvec[k].x, hy = hvec[k].y, hz = hvec[k].z;

      // see eik_dot_r() for the reuse of zx and zxy

      if (ky != kyold) {
        if (kx != kxold) { zx = z[kx].x; kxold = kx; }
        c_rmult(zxy, z[ky].y, zx);
        kyold = ky;
      }
      c_crmult(zc, z[kz].z, zxy);

      if (termflag[TERM_COUL]) {                        // 1/r
        const double im = kenergy[ne++]*(zc.im*sfac_all[n].re+sfac_all[n].im*zc.re);
        if (termflag[TERM_DIPOLE]) ncoul = n;
        n++;
        eksum[0][0] += hx*im; eksum[0][1] += hy*im; eksum[0][2] += hz*im;
      }
      if (termflag[TERM_DISP_GEOM]) {                   // geometric 1/r^6
        const double im = kenergy[ne++]*(zc.im*sfac_all[n].re+sfac_all[n].im*zc.re);
        n++;
        eksum[1][0] += hx*im; eksum[1][1] += hy*im; eksum[1][2] += hz*im;
      }
      if (termflag[TERM_DISP_ARITH]) {                  // arithmetic 1/r^6
        const double ck = kenergy[ne++];
        for (int m = 2; m < 9; m++) {
          const double im = ck*(zc.im*sfac_all[n].re+sfac_all[n].im*zc.re);
          n++;
          eksum[m][0] += hx*im; eksum[m][1] += hy*im; eksum[m][2] += hz*im;
        }
      }
      if (termflag[TERM_DIPOLE]) {                      // dipole
        const double muk = mui[0]*hx+mui[1]*hy+mui[2]*hz;
        double im = kenergy[ne]*(zc.im*sfac_all[n].re+sfac_all[n].im*zc.re)*muk;
        double im2 = kenergy[ne]*(zc.re*sfac_all[n].re-sfac_all[n].im*zc.im);
        eksum[9][0] += hx*im; eksum[9][1] += hy*im; eksum[9][2] += hz*im;
        tq[i][0] += -mui[1]*hz*im2 + mui[2]*hy*im2;     // torque
        tq[i][1] += -mui[2]*hx*im2 + mui[0]*hz*im2;
        tq[i][2] += -mui[0]*hy*im2 + mui[1]*hx*im2;
        if (termflag[TERM_COUL]) {                      // charge-dipole
          const double qi = q[i]*c[TERM_COUL];
          im = -kenergy[ne]*(zc.re*sfac_all[ncoul].re-sfac_all[ncoul].im*zc.im)*muk;
          im += kenergy[ne]*(zc.re*sfac_all[n].re-sfac_all[n].im*zc.im)*qi;
          eksum[9][0] += hx*im; eksum[9][1] += hy*im; eksum[9][2] += hz*im;

          im2 = kenergy[ne]*(zc.re*sfac_all[ncoul].im+sfac_all[ncoul].re*zc.im);
          im2 += -kenergy[ne]*(zc.re*sfac_all[n].im-sfac_all[n].im*zc.re);
          tq[i][0] += -mui[1]*hz*im2 + mui[2]*hy*im2;   // torque
          tq[i][1] += -mui[2]*hx*im2 + mui[0]*hz*im2;
          tq[i][2] += -mui[0]*hy*im2 + mui[1]*hx*im2;
        }
        n++;
        ne++;
      }
    }

    if (termflag[TERM_COUL]) {                          // 1/r
      const double qi = q[i]*c[TERM_COUL];
      f[i][0] -= eksum[0][0]*qi;
      f[i][1] -= eksum[0][1]*qi;
      f[i][2] -= eksum[0][2]*qi;
    }
    if (termflag[TERM_DISP_GEOM]) {                     // geometric 1/r^6
      const double bi = B[type[i]]*c[TERM_DISP_GEOM];
      f[i][0] -= eksum[1][0]*bi;
      f[i][1] -= eksum[1][1]*bi;
      f[i][2] -= eksum[1][2]*bi;
    }
    if (termflag[TERM_DISP_ARITH]) {                    // arithmetic 1/r^6
      const double *bi = &B[7*type[i]];
      for (int m = 2; m < 9; m++) {
        const double c2 = bi[8-m]*c[TERM_DISP_ARITH];
        f[i][0] -= eksum[m][0]*c2;
        f[i][1] -= eksum[m][1]*c2;
        f[i][2] -= eksum[m][2]*c2;
      }
    }
    if (termflag[TERM_DIPOLE]) {                        // dipole
      f[i][0] -= eksum[9][0];
      f[i][1] -= eksum[9][1];
      f[i][2] -= eksum[9][2];
    }
  }
}

/* ---------------------------------------------------------------------- */

void EwaldDisp::compute_energy()
{
  energy = 0.0;
  if (!eflag_global) return;

  const double qscale = force->qqrd2e * scale;
  const double c[EWALD_NTERMS] = {
    4.0*MY_PI*qscale/volume, 2.0*MY_PI*MY_PIS/(24.0*volume),
    2.0*MY_PI*MY_PIS/(192.0*volume), 4.0*MY_PI*mumurd2e/volume};
  double esum[EWALD_NTERMS];
  int n = 0, ne = 0, ncoul = 0;

  memset(esum, 0, EWALD_NTERMS*sizeof(double));         // reset sums
  for (int k = 0; k < kcount; k++) {                    // sum over k vectors
    if (termflag[TERM_COUL]) {                          // 1/r
      esum[TERM_COUL] +=
        kenergy[ne++]*(sfac_all[n].re*sfac_all[n].re+sfac_all[n].im*sfac_all[n].im);
      if (termflag[TERM_DIPOLE]) ncoul = n;
      n++;
    }
    if (termflag[TERM_DISP_GEOM]) {                     // geometric 1/r^6
      esum[TERM_DISP_GEOM] +=
        kenergy[ne++]*(sfac_all[n].re*sfac_all[n].re+sfac_all[n].im*sfac_all[n].im);
      n++;
    }
    if (termflag[TERM_DISP_ARITH]) {                    // arithmetic 1/r^6
      const double r =
            (sfac_all[n+0].re*sfac_all[n+6].re+sfac_all[n+0].im*sfac_all[n+6].im)+
            (sfac_all[n+1].re*sfac_all[n+5].re+sfac_all[n+1].im*sfac_all[n+5].im)+
            (sfac_all[n+2].re*sfac_all[n+4].re+sfac_all[n+2].im*sfac_all[n+4].im)+
        0.5*(sfac_all[n+3].re*sfac_all[n+3].re+sfac_all[n+3].im*sfac_all[n+3].im);
      n += 7;
      esum[TERM_DISP_ARITH] += kenergy[ne++]*r;
    }
    if (termflag[TERM_DIPOLE]) {                        // dipole
      esum[TERM_DIPOLE] +=
        kenergy[ne]*(sfac_all[n].re*sfac_all[n].re+sfac_all[n].im*sfac_all[n].im);
      if (termflag[TERM_COUL]) {                        // charge-dipole
        esum[TERM_DIPOLE] += kenergy[ne]*2.0*
          (sfac_all[n].re*sfac_all[ncoul].im-sfac_all[n].im*sfac_all[ncoul].re);
      }
      ne++;
      n++;
    }
  }
  for (int t = 0; t < EWALD_NTERMS; t++) energy += c[t]*esum[t]-energy_self[t];
}

/* ---------------------------------------------------------------------- */

void EwaldDisp::compute_energy_peratom()
{
  if (!eflag_atom) return;

  const int nz = 2*kmax+1;
  const int nlocal = atom->nlocal;
  const double *q = atom->q;
  const int *type = atom->type;
  double **mu = atom->mu;
  const double qscale = force->qqrd2e * scale;
  const double c[EWALD_NTERMS] = {
    4.0*MY_PI*qscale/volume, 2.0*MY_PI*MY_PIS/(24.0*volume),
    2.0*MY_PI*MY_PIS/(192.0*volume), 4.0*MY_PI*mumurd2e/volume};
  double esum[EWALD_MAX_NSUMS];
  double mui[3] = {0.0,0.0,0.0};

  for (int i = 0; i < nlocal; i++) {
    const cvector *z = &ekr[i*nz];
    int n = 0, ne = 0, ncoul = 0;
    int kxold = -1, kyold = -1;
    complex zx = {0, 0}, zxy = {0, 0}, zc = {0, 0};

    memset(esum, 0, EWALD_MAX_NSUMS*sizeof(double));
    if (termflag[TERM_DIPOLE]) {
      mui[0] = c[TERM_DIPOLE]*mu[i][0];
      mui[1] = c[TERM_DIPOLE]*mu[i][1];
      mui[2] = c[TERM_DIPOLE]*mu[i][2];
    }

    for (int k = 0; k < kcount; k++) {
      const int kx = kvec[k].x, ky = kvec[k].y, kz = kvec[k].z;

      // see eik_dot_r() for the reuse of zx and zxy

      if (ky != kyold) {
        if (kx != kxold) { zx = z[kx].x; kxold = kx; }
        c_rmult(zxy, z[ky].y, zx);
        kyold = ky;
      }
      c_crmult(zc, z[kz].z, zxy);

      if (termflag[TERM_COUL]) {                        // 1/r
        esum[0] += kenergy[ne++]*(sfac_all[n].re*zc.re-sfac_all[n].im*zc.im);
        if (termflag[TERM_DIPOLE]) ncoul = n;
        n++;
      }
      if (termflag[TERM_DISP_GEOM]) {                   // geometric 1/r^6
        esum[1] += kenergy[ne++]*(sfac_all[n].re*zc.re-sfac_all[n].im*zc.im);
        n++;
      }
      if (termflag[TERM_DISP_ARITH]) {                  // arithmetic 1/r^6
        const double ck = kenergy[ne++];
        for (int m = 2; m < 9; m++) {
          esum[m] += ck*(sfac_all[n].re*zc.re-sfac_all[n].im*zc.im);
          n++;
        }
      }
      if (termflag[TERM_DIPOLE]) {                      // dipole
        const double muk = mui[0]*hvec[k].x+mui[1]*hvec[k].y+mui[2]*hvec[k].z;
        esum[9] += kenergy[ne]*(sfac_all[n].re*zc.re-sfac_all[n].im*zc.im)*muk;
        if (termflag[TERM_COUL]) {                      // charge-dipole
          const double qj = q[i]*c[TERM_COUL];
          esum[9] += kenergy[ne]*
            (sfac_all[ncoul].im*zc.re+sfac_all[ncoul].re*zc.im)*muk;
          esum[9] -= kenergy[ne]*(sfac_all[n].re*zc.im+sfac_all[n].im*zc.re)*qj;
        }
        n++;
        ne++;
      }
    }

    if (termflag[TERM_COUL]) {                          // 1/r
      const double qj = q[i]*c[TERM_COUL];
      eatom[i] += esum[0]*qj - energy_self_peratom[i][TERM_COUL];
    }
    if (termflag[TERM_DISP_GEOM]) {                     // geometric 1/r^6
      const double bj = B[type[i]]*c[TERM_DISP_GEOM];
      eatom[i] += esum[1]*bj - energy_self_peratom[i][TERM_DISP_GEOM];
    }
    if (termflag[TERM_DISP_ARITH]) {                    // arithmetic 1/r^6
      const double *bj = &B[7*type[i]];
      for (int m = 2; m < 9; m++) {
        const double c2 = bj[8-m]*c[TERM_DISP_ARITH];
        eatom[i] += 0.5*esum[m]*c2;
      }
      eatom[i] -= energy_self_peratom[i][TERM_DISP_ARITH];
    }
    if (termflag[TERM_DIPOLE]) {                        // dipole
      eatom[i] += esum[9] - energy_self_peratom[i][TERM_DIPOLE];
    }
  }
}

/* ---------------------------------------------------------------------- */

void EwaldDisp::compute_virial()
{
  memset(virial, 0, 6*sizeof(double));
  if (!vflag_global) return;

  const double qscale = force->qqrd2e * scale;
  const double c[EWALD_NTERMS] = {
    4.0*MY_PI*qscale/volume, 2.0*MY_PI*MY_PIS/(24.0*volume),
    2.0*MY_PI*MY_PIS/(192.0*volume), 4.0*MY_PI*mumurd2e/volume};
  double vsum[EWALD_NTERMS][6];
  int n = 0, nv = 0, ncoul = 0;

  memset(vsum, 0, EWALD_NTERMS*6*sizeof(double));
  for (int k = 0; k < kcount; k++) {                    // sum over k vectors
    if (termflag[TERM_COUL]) {                          // 1/r
      const double r = sfac_all[n].re*sfac_all[n].re+sfac_all[n].im*sfac_all[n].im;
      if (termflag[TERM_DIPOLE]) ncoul = n;
      n++;
      for (int m = 0; m < 6; m++) vsum[TERM_COUL][m] += kvirial[nv+m]*r;
      nv += 6;
    }
    if (termflag[TERM_DISP_GEOM]) {                     // geometric 1/r^6
      const double r = sfac_all[n].re*sfac_all[n].re+sfac_all[n].im*sfac_all[n].im;
      n++;
      for (int m = 0; m < 6; m++) vsum[TERM_DISP_GEOM][m] += kvirial[nv+m]*r;
      nv += 6;
    }
    if (termflag[TERM_DISP_ARITH]) {                    // arithmetic 1/r^6
      const double r =
            (sfac_all[n+0].re*sfac_all[n+6].re+sfac_all[n+0].im*sfac_all[n+6].im)+
            (sfac_all[n+1].re*sfac_all[n+5].re+sfac_all[n+1].im*sfac_all[n+5].im)+
            (sfac_all[n+2].re*sfac_all[n+4].re+sfac_all[n+2].im*sfac_all[n+4].im)+
        0.5*(sfac_all[n+3].re*sfac_all[n+3].re+sfac_all[n+3].im*sfac_all[n+3].im);
      n += 7;
      for (int m = 0; m < 6; m++) vsum[TERM_DISP_ARITH][m] += kvirial[nv+m]*r;
      nv += 6;
    }
    if (termflag[TERM_DIPOLE]) {                        // dipole
      const double r = sfac_all[n].re*sfac_all[n].re+sfac_all[n].im*sfac_all[n].im;
      for (int m = 0; m < 6; m++) vsum[TERM_DIPOLE][m] += kvirial[nv+m]*r;
      if (termflag[TERM_COUL]) {                        // charge-dipole
        const double r2 =
          2.0*(sfac_all[n].re*sfac_all[ncoul].im-sfac_all[n].im*sfac_all[ncoul].re);
        for (int m = 0; m < 6; m++) vsum[TERM_DIPOLE][m] += kvirial[nv+m]*r2;
      }
      nv += 6;
      n++;
    }
  }

  for (int t = 0; t < EWALD_NTERMS; t++)
    if (termflag[t]) {
      double self[6] = {virial_self[t], virial_self[t], virial_self[t], 0, 0, 0};
      shape_scalar_mult(vsum[t], c[t]);
      shape_add(virial, vsum[t]);
      shape_subtr(virial, self);
    }
}

/* ---------------------------------------------------------------------- */

void EwaldDisp::compute_virial_dipole()
{
  if (!termflag[TERM_DIPOLE]) return;
  if (!vflag_atom && !vflag_global) return;

  const int nz = 2*kmax+1;
  const int nlocal = atom->nlocal;
  double **mu = atom->mu;
  const double qscale = force->qqrd2e * scale;
  const double c[EWALD_NTERMS] = {
    8.0*MY_PI*qscale/volume, 2.0*MY_PI*MY_PIS/(12.0*volume),
    2.0*MY_PI*MY_PIS/(192.0*volume), 8.0*MY_PI*mumurd2e/volume};
  double vsum[6], vsum_total[6];
  double mui[3];

  // offsets of the dipole channel within the per-k-vector blocks of
  // sfac_all (nsums channels) and kenergy (nterms entries); the coulomb
  // channel, when enabled, is always the first one

  const int ndip = termflag[TERM_COUL]+termflag[TERM_DISP_GEOM]+7*termflag[TERM_DISP_ARITH];
  const int nedip = termflag[TERM_COUL]+termflag[TERM_DISP_GEOM]+termflag[TERM_DISP_ARITH];

  memset(vsum_total, 0, 6*sizeof(double));
  for (int i = 0; i < nlocal; i++) {
    const cvector *z = &ekr[i*nz];
    int kxold = -1, kyold = -1;
    complex zx = {0, 0}, zxy = {0, 0}, zc;

    memset(vsum, 0, 6*sizeof(double));
    mui[0] = c[TERM_DIPOLE]*mu[i][0];
    mui[1] = c[TERM_DIPOLE]*mu[i][1];
    mui[2] = c[TERM_DIPOLE]*mu[i][2];

    for (int k = 0; k < kcount; k++) {
      const int kx = kvec[k].x, ky = kvec[k].y, kz = kvec[k].z;

      // see eik_dot_r() for the reuse of zx and zxy

      if (ky != kyold) {
        if (kx != kxold) { zx = z[kx].x; kxold = kx; }
        c_rmult(zxy, z[ky].y, zx);
        kyold = ky;
      }
      c_crmult(zc, z[kz].z, zxy);

      const complex &sd = sfac_all[k*nsums+ndip];
      const double ke = kenergy[k*nterms+nedip];
      double im = ke*(zc.re*sd.re - sd.im*zc.im);
      if (termflag[TERM_COUL]) {                        // charge-dipole
        const complex &sc = sfac_all[k*nsums];
        im += ke*(zc.im*sc.re + sc.im*zc.re);
      }
      vsum[0] -= mui[0]*hvec[k].x*im;
      vsum[1] -= mui[1]*hvec[k].y*im;
      vsum[2] -= mui[2]*hvec[k].z*im;
      vsum[3] -= mui[0]*hvec[k].y*im;
      vsum[4] -= mui[0]*hvec[k].z*im;
      vsum[5] -= mui[1]*hvec[k].z*im;
    }

    if (vflag_global)
      for (int m = 0; m < 6; m++)
        vsum_total[m] -= vsum[m];

    if (vflag_atom)
      for (int m = 0; m < 6; m++)
        vatom[i][m] -= vsum[m];
  }

  if (vflag_global) {
    MPI_Allreduce(vsum_total,vsum,6,MPI_DOUBLE,MPI_SUM,world);
    for (int m = 0; m < 6; m++)
      virial[m] += vsum[m];
  }
}

/* ---------------------------------------------------------------------- */

void EwaldDisp::compute_virial_peratom()
{
  if (!vflag_atom) return;

  const int nz = 2*kmax+1;
  const int nlocal = atom->nlocal;
  const double *q = atom->q;
  const int *type = atom->type;
  double **mu = atom->mu;
  const double qscale = force->qqrd2e * scale;
  const double c[EWALD_NTERMS] = {
    4.0*MY_PI*qscale/volume, 2.0*MY_PI*MY_PIS/(24.0*volume),
    2.0*MY_PI*MY_PIS/(192.0*volume), 4.0*MY_PI*mumurd2e/volume};
  double vsum[EWALD_MAX_NSUMS][6];
  double mui[3] = {0.0,0.0,0.0};

  for (int i = 0; i < nlocal; i++) {
    const cvector *z = &ekr[i*nz];
    int n = 0, nv = 0, ncoul = 0;
    int kxold = -1, kyold = -1;
    complex zx = {0, 0}, zxy = {0, 0}, zc = {0, 0};

    memset(vsum, 0, EWALD_MAX_NSUMS*6*sizeof(double));
    if (termflag[TERM_DIPOLE]) {
      mui[0] = c[TERM_DIPOLE]*mu[i][0];
      mui[1] = c[TERM_DIPOLE]*mu[i][1];
      mui[2] = c[TERM_DIPOLE]*mu[i][2];
    }

    for (int k = 0; k < kcount; k++) {
      const int kx = kvec[k].x, ky = kvec[k].y, kz = kvec[k].z;

      // see eik_dot_r() for the reuse of zx and zxy

      if (ky != kyold) {
        if (kx != kxold) { zx = z[kx].x; kxold = kx; }
        c_rmult(zxy, z[ky].y, zx);
        kyold = ky;
      }
      c_crmult(zc, z[kz].z, zxy);

      if (termflag[TERM_COUL]) {                        // 1/r
        if (termflag[TERM_DIPOLE]) ncoul = n;
        const double r = sfac_all[n].re*zc.re-sfac_all[n].im*zc.im;
        n++;
        for (int m = 0; m < 6; m++) vsum[0][m] += kvirial[nv+m]*r;
        nv += 6;
      }
      if (termflag[TERM_DISP_GEOM]) {                   // geometric 1/r^6
        const double r = sfac_all[n].re*zc.re-sfac_all[n].im*zc.im;
        n++;
        for (int m = 0; m < 6; m++) vsum[1][m] += kvirial[nv+m]*r;
        nv += 6;
      }
      if (termflag[TERM_DISP_ARITH]) {                  // arithmetic 1/r^6
        // the seven sub-channels share the same six kvirial entries
        for (int s = 2; s < 9; s++) {
          const double r = sfac_all[n].re*zc.re-sfac_all[n].im*zc.im;
          n++;
          for (int m = 0; m < 6; m++) vsum[s][m] += kvirial[nv+m]*r;
        }
        nv += 6;
      }
      if (termflag[TERM_DIPOLE]) {                      // dipole
        const double muk = mui[0]*hvec[k].x+mui[1]*hvec[k].y+mui[2]*hvec[k].z;
        double r = (sfac_all[n].re*zc.re-sfac_all[n].im*zc.im)*muk;
        for (int m = 0; m < 6; m++) vsum[9][m] += kvirial[nv+m]*r;
        if (termflag[TERM_COUL]) {                      // charge-dipole
          const double qj = q[i]*c[TERM_COUL];
          r = (sfac_all[ncoul].im*zc.re+sfac_all[ncoul].re*zc.im)*muk;
          r += -(sfac_all[n].re*zc.im+sfac_all[n].im*zc.re)*qj;
          for (int m = 0; m < 6; m++) vsum[9][m] += kvirial[nv+m]*r;
        }
        nv += 6;
        n++;
      }
    }

    if (termflag[TERM_COUL]) {                          // 1/r
      const double qi = q[i]*c[TERM_COUL];
      for (int m = 0; m < 6; m++) vatom[i][m] += vsum[0][m]*qi;
    }
    if (termflag[TERM_DISP_GEOM]) {                     // geometric 1/r^6
      const double bi = B[type[i]]*c[TERM_DISP_GEOM];
      for (int m = 0; m < 6; m++) vatom[i][m] += vsum[1][m]*bi;
    }
    if (termflag[TERM_DISP_ARITH]) {                    // arithmetic 1/r^6
      const double *bj = &B[7*type[i]];
      for (int s = 2; s < 9; s++) {
        const double c2 = bj[8-s]*c[TERM_DISP_ARITH];
        for (int m = 0; m < 6; m++) vatom[i][m] += 0.5*vsum[s][m]*c2;
      }
    }
    if (termflag[TERM_DIPOLE]) {                        // dipole
      for (int m = 0; m < 6; m++) vatom[i][m] += vsum[9][m];
    }

    for (int t = 0; t < EWALD_NTERMS; t++)
      if (termflag[t])
        for (int m = 0; m < 3; m++) vatom[i][m] -= virial_self_peratom[i][t];
  }
}

/* ----------------------------------------------------------------------
   Slab-geometry correction term to dampen inter-slab interactions between
   periodically repeating slabs.  Yields good approximation to 2D Ewald if
   adequate empty space is left between repeating slabs (J. Chem. Phys.
   111, 3155).  Slabs defined here to be parallel to the xy plane. Also
   extended to non-neutral systems (J. Chem. Phys. 131, 094107).
------------------------------------------------------------------------- */

void EwaldDisp::slabcorr()
{
  // compute local contribution to global dipole moment

  double *q = atom->q;
  double **x = atom->x;
  double zprd_slab = domain->zprd*slab_volfactor;
  int nlocal = atom->nlocal;

  double dipole = 0.0;
  for (int i = 0; i < nlocal; i++) dipole += q[i]*x[i][2];

  if (termflag[TERM_DIPOLE] && atom->mu) {
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

    if (termflag[TERM_DIPOLE] && atom->mu)
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
    qsum*dipole_r2 - qsum*qsum*zprd_slab*zprd_slab/12.0)/volume;
  const double qscale = force->qqrd2e * scale;

  if (eflag_global) energy += qscale * e_slabcorr;

  // per-atom energy

  if (eflag_atom) {
    double efact = qscale * MY_2PI/volume;
    for (int i = 0; i < nlocal; i++)
      eatom[i] += efact * q[i]*(x[i][2]*dipole_all - 0.5*(dipole_r2 +
        qsum*x[i][2]*x[i][2]) - qsum*zprd_slab*zprd_slab/12.0);
  }

  // add on force corrections

  double ffact = qscale * (-4.0*MY_PI/volume);
  double **f = atom->f;

  for (int i = 0; i < nlocal; i++)
    f[i][2] += ffact * q[i]*(dipole_all - qsum*x[i][2]);

  // add on torque corrections

  if (termflag[TERM_DIPOLE] && atom->mu && atom->torque) {
    double **mu = atom->mu;
    double **torque = atom->torque;
    for (int i = 0; i < nlocal; i++) {
      torque[i][0] += ffact * dipole_all * mu[i][1];
      torque[i][1] += -ffact * dipole_all * mu[i][0];
    }
  }
}

/* ----------------------------------------------------------------------
   estimate the memory used by the allocated arrays
------------------------------------------------------------------------- */

double EwaldDisp::memory_usage()
{
  double bytes = 0.0;
  bytes += (double)kcount_max*(sizeof(hvector)+sizeof(kvector));     // hvec, kvec
  bytes += (double)kcount_max*nterms*7.0*sizeof(double);             // kenergy, kvirial
  bytes += (double)kcount_max*nsums*2.0*sizeof(complex);             // sfac, sfac_all
  bytes += (double)nevec_max*sizeof(cvector);                        // ekr
  if (termflag[TERM_DISP_GEOM]) bytes += (double)(atom->ntypes+1)*sizeof(double);     // B
  if (termflag[TERM_DISP_ARITH]) bytes += 7.0*(atom->ntypes+1)*sizeof(double);        // B
  if (peratom_allocate_flag)
    bytes += 2.0*nmax*EWALD_NTERMS*sizeof(double);    // energy/virial_self_peratom
  return bytes;
}

/* ----------------------------------------------------------------------
   Newton solver used to find g_ewald for LJ systems
------------------------------------------------------------------------- */

double EwaldDisp::NewtonSolve(double x, double Rc,
                              bigint natoms, double vol, double b2)
{
  const int maxit = 10000;       // maximum number of iterations
  const double tol = 0.00001;    // convergence tolerance
  double dx;

  for (int i = 0; i < maxit; i++) {
    double dfx = derivf(x,Rc,natoms,vol,b2);
    if (dfx == 0.0 || dfx != dfx) // flat/invalid derivative
      return -1;
    dx = f(x,Rc,natoms,vol,b2) / dfx;
    // Damped (backtracking) Newton: the error estimate is monotonic in g with a
    // unique positive root, but a raw Newton step can overshoot to g <= 0 at
    // large cutoffs (where f is flat near the initial guess) and abort.  Halve
    // the step until it keeps x positive so the solver converges robustly.
    while (x - dx <= 0.0) dx *= 0.5;
    x = x - dx; //Update x
    if (fabs(dx) < tol) return x;
    if (x != x) // solver failed
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

  if (termflag[TERM_DIPOLE]) {                          // dipole
    double rg2 = a*a;
    double rg4 = rg2*rg2;
    double rg6 = rg4*rg2;
    double Cc = 4.0*rg4 + 6.0*rg2 + 3.0;
    double Dc = 8.0*rg6 + 20.0*rg4 + 30.0*rg2 + 15.0;
    f = (b2/(sqrt(vol*powint(x,4)*powint(Rc,9)*natoms)) *
      sqrt(13.0/6.0*Cc*Cc + 2.0/15.0*Dc*Dc - 13.0/15.0*Cc*Dc) *
      exp(-rg2)) - accuracy;
  } else if (termflag[TERM_DISP_GEOM] || termflag[TERM_DISP_ARITH]) {  // LJ
    // Kolafa-Perram real-space RMS force error for r^-6, extended by
    // Isele-Holder et al. (J. Chem. Phys. 137, 174107 (2012), Eq. 20).  The
    // published Eq. 20 has a dimensional typo in the denominator (sqrt(N) V Rc);
    // the correct form is sqrt(N V Rc), which makes the estimate dimensionally
    // consistent and match measured bulk errors to ~Kolafa-Perram accuracy.
    double a2 = a*a;
    f = (b2*sqrt(MY_PI)*powint(x,5)/sqrt(vol*Rc*(double)natoms) *
      (6.0/(a2*a2*a2) + 6.0/(a2*a2) + 3.0/a2 + 1.0) * exp(-a2) - accuracy);
  }

  return f;
}

/* ----------------------------------------------------------------------
 Calculate numerical derivative f'(x)
 ------------------------------------------------------------------------- */

double EwaldDisp::derivf(double x, double Rc,
                         bigint natoms, double vol, double b2)
{
  double h = 0.000001;           // derivative step-size
  return (f(x + h,Rc,natoms,vol,b2) - f(x,Rc,natoms,vol,b2)) / h;
}
