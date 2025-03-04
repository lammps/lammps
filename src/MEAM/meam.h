/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_MEAM_H
#define LMP_MEAM_H

#include "exceptions.h"

#include <cmath>
#include <string>

constexpr int MAXELT = 8;
// This should really be true for sanity reasons (see MEAM::alloyparams),
// but breaks too many potential files - including our tests. So disable for now.
constexpr bool STRICT_IJ_ORDER = false;

namespace LAMMPS_NS {
class Memory;

namespace MEAM_NS {

typedef enum { FCC, BCC, HCP, DIM, DIA, DIA3, B1, C11, L12, B2, CH4, LIN, ZIG, TRI, SC } lattice_t;

class MEAM {
 public:
  MEAM(Memory *mem);
  virtual ~MEAM();

  int copymode;
  int msmeamflag;

 protected:
  Memory *memory;

  // cutforce = force cutoff
  // cutforcesq = force cutoff squared

  double cutforce, cutforcesq;

  // Ec_meam = cohesive energy
  // re_meam = nearest-neighbor distance
  // B_meam = bulk modulus
  // ielt_meam = atomic number of element
  // A_meam = adjustable parameter
  // alpha_meam = sqrt(9*Omega*B/Ec)
  // rho0_meam = density scaling parameter
  // delta_meam = heat of formation for alloys
  // beta[0-3]_meam = electron density constants
  // t[0-3]_meam = coefficients on densities in Gamma computation
  // rho_ref_meam = background density for reference structure
  // ibar_meam(i) = selection parameter for Gamma function for elt i,
  // lattce_meam(i,j) = lattce configuration for elt i or alloy (i,j)
  // neltypes = maximum number of element type defined
  // eltind = index number of pair (similar to Voigt notation; ij = ji)
  // phir = pair potential function array
  // phirar[1-6] = spline coeffs
  // attrac_meam = attraction parameter in Rose energy
  // repuls_meam = repulsion parameter in Rose energy
  // nn2_meam = 1 if second nearest neighbors are to be computed, else 0
  // zbl_meam = 1 if zbl potential for small r to be use, else 0
  // emb_lin_neg = 1 if linear embedding function for rhob to be used, else 0
  // bkgd_dyn = 1 if reference densities follows Dynamo, else 0
  // Cmin_meam, Cmax_meam = min and max values in screening cutoff
  // rc_meam = cutoff distance for meam
  // delr_meam = cutoff region for meam
  // ebound_meam = factor giving maximum boundary of sceen fcn ellipse
  // augt1 = flag for whether t1 coefficient should be augmented
  // ialloy = flag for newer alloy formulation (as in dynamo code)
  // mix_ref_t = flag to recover "old" way of computing t in reference config
  // erose_form = selection parameter for form of E_rose function
  // gsmooth_factor = factor determining length of G smoothing region
  // vind[23]D = Voight notation index maps for 2 and 3D
  // v2D,v3D = array of factors to apply for Voight notation

  // MS-MEAM parameters

  // msmeamflag = flag to activate MS-MEAM (public; above)
  // betam[1-3]_meam = MS-MEAM electron density constants
  // tm[1-3]_meam = MS-MEAM coefficients on densities in Gamma computation

  // nr,dr = pair function discretization parameters
  // nrar,rdrar = spline coeff array parameters

  // theta = angle between three atoms in line, zigzag, and trimer reference structures
  // stheta_meam = sin(theta/2) in radian used in line, zigzag, and trimer reference structures
  // ctheta_meam = cos(theta/2) in radian used in line, zigzag, and trimer reference structures

  double Ec_meam[MAXELT][MAXELT], re_meam[MAXELT][MAXELT];
  double A_meam[MAXELT], alpha_meam[MAXELT][MAXELT], rho0_meam[MAXELT];
  double delta_meam[MAXELT][MAXELT];
  double beta0_meam[MAXELT], beta1_meam[MAXELT];
  double beta2_meam[MAXELT], beta3_meam[MAXELT];
  double t0_meam[MAXELT], t1_meam[MAXELT];
  double t2_meam[MAXELT], t3_meam[MAXELT];
  double rho_ref_meam[MAXELT];
  int ibar_meam[MAXELT], ielt_meam[MAXELT];
  lattice_t lattce_meam[MAXELT][MAXELT];
  int nn2_meam[MAXELT][MAXELT];
  int zbl_meam[MAXELT][MAXELT];
  int eltind[MAXELT][MAXELT];
  int neltypes;

  double **phir;

  double **phirar, **phirar1, **phirar2, **phirar3, **phirar4, **phirar5, **phirar6;

  double attrac_meam[MAXELT][MAXELT], repuls_meam[MAXELT][MAXELT];

  double Cmin_meam[MAXELT][MAXELT][MAXELT];
  double Cmax_meam[MAXELT][MAXELT][MAXELT];
  double rc_meam, delr_meam, ebound_meam[MAXELT][MAXELT];
  int augt1, ialloy, mix_ref_t, erose_form;
  int emb_lin_neg, bkgd_dyn;
  double gsmooth_factor;

  int vind2D[3][3], vind3D[3][3][3];    // x-y-z to Voigt-like index
  int v2D[6], v3D[10];                  // multiplicity of Voigt index (i.e. [1] -> xy+yx = 2

  int nr, nrar;
  double dr, rdrar;

  // MS-MEAM parameters

  double t1m_meam[MAXELT], t2m_meam[MAXELT], t3m_meam[MAXELT];
  double beta1m_meam[MAXELT], beta2m_meam[MAXELT], beta3m_meam[MAXELT];

 public:
  int nmax;
  double *rho, *rho0, *rho1, *rho2, *rho3, *frhop;
  double *gamma, *dgamma1, *dgamma2, *dgamma3, *arho2b;
  double **arho1, **arho2, **arho3, **arho3b, **t_ave, **tsq_ave;

  // MS-MEAM arrays

  double **arho1m, **arho2m, *arho2mb, **arho3m, **arho3mb;

  int maxneigh;
  double *scrfcn, *dscrfcn, *fcpair;

  //angle for trimer, zigzag, line reference structures
  double stheta_meam[MAXELT][MAXELT];
  double ctheta_meam[MAXELT][MAXELT];

 protected:
  // meam_funcs.cpp
  double G_gam(const double gamma, const int ibar, int &errorflag) const;
  double dG_gam(const double gamma, const int ibar, double &dG) const;
  bool rhobar12(const double r, const int a, const int b, double &rhobar1, double &rhobar2) const;
  double embedding(const double A, const double Ec, const double rhobar, double &dF) const;
  double invert_eam(const double r, const int a, const int b, const double Eu, const double F1, const double F2) const;
  double phi_meam(double, int, int) const;
  double phi_2nn_series(const double scrn, const int Z1, const int Z2, const int a, const int b,
                        const double r, const double arat) const;

 protected:
  void getscreen(int i, double *scrfcn, double *dscrfcn, double *fcpair, double **x, int numneigh,
                 int *firstneigh, int numneigh_full, int *firstneigh_full, int ntype, int *type,
                 int *fmap);
  void calc_rho1(int i, int ntype, int *type, int *fmap, double **x, int numneigh, int *firstneigh,
                 double *scrfcn, double *fcpair);

  void alloyparams();
  void compute_pair_meam();
  void compute_reference_density();
  void get_tavref(double *, double *, double *, double *, double *, double *, double, double,
                  double, double, double, double, double, double, double, int, int, lattice_t) const;
  double get_sijk(double, int, int, int) const;
  void get_densref(double, int, int, double *, double *, double *, double *, double *, double *,
                   double *, double *, double *, double *, double *, double *, double *, double *) const; // last 6 args for msmeam
  void interpolate_meam(int);

 public:
  void setup_library(int nelt, lattice_t *lat, int *ielement, double *atwt, double *alpha,
                     double *b0, double *b1, double *b2, double *b3, double *alat, double *esub,
                     double *asub, double *t0, double *t1, double *t2, double *t3,
                     double *rozero, int *ibar);
  void setup_library_ms(int nelt, double *b1m, double *b2m, double *b3m,
                        double *t1m, double *t2m, double *t3m);
  void setup_param(int which, double value, int nindex, int *index /*index(3)*/,
                   int *errorflag);
  void setup_finish(double *cutmax);

  void density_precompute();
  void density_setup(int atom_nmax, int nall, int n_neigh);
  void density_local(int i, int ntype, int *type, int *fmap, double **x, int numneigh,
                     int *firstneigh, int numneigh_full, int *firstneigh_full, int fnoffset);

  void eval_energy(int nlocal, int eflag_either, int eflag_global, int eflag_atom,
                   double *eng_vdwl, double *eatom, int ntype, int *type, int *fmap,
                   double **scale, int &errorflag);
  void eval_force(int i, int eflag_global, int eflag_atom, int vflag_global, int vflag_atom,
                  double *eng_vdwl, double *eatom, int ntype, int *type, int *fmap, double **scale,
                  double **x, int numneigh, int *firstneigh, int numneigh_full,
                  int *firstneigh_full, int fnoffset, double **f, double **vatom, double *virial);
};

class MEAMException : public LAMMPSException {
 public:
  MEAMException(const std::string &msg) :
      LAMMPSException(msg)
  {
  }
};

//-----------------------------------------------------------------------------
// Reference lattice definition
// Any value not set will be 0/false/nullptr, make sure this is a sensible default

typedef struct {
  // Name of the lattice in potential files
  const char * name;
  // lattice is valid for single element references
  bool single;
  // factor from alat to re
  double re;
  // Number of nearest neighbors (=coordination number)
  int Zij;
  // Number of second-nearest neighbors
  int Zij2;
  // number of atoms that screen the 2NN bond
  int Nscr2;
  // distance ratio R1/R2 (a2nn in dynamo)
  double ratio_2nn;
  // ratio_2nn depends on sin(theta)
  bool ratio_2nn_angular;
  // shape function getter, assume s is initialized to 0.0
  void (*shpfcn)(const double sthe, const double cthe, double (&s)[3]);
  // i-j cohesive energy from formation energy, use averaging if not specified
  double (*ecoh)(const double Eii, const lattice_t ilat, const double Ejj, const double delta);
  // special: reference structure includes third nearest neighbors
  bool nn3;
} reference_lattice_t;

constexpr int MAXLAT = 15;
extern const reference_lattice_t lattice_defs[MAXLAT];

//-----------------------------------------------------------------------------
// Functions we need for compat

static
inline bool iszero(const double f)
{
  return fabs(f) < 1e-20;
}

static
inline bool isone(const double f)
{
  return fabs(f - 1.0) < 1e-20;
}

//-----------------------------------------------------------------------------
// Helper functions

static
inline double fdiv_zero(const double n, const double d)
{
  if (iszero(d)) return 0.0;
  return n / d;
}


//-----------------------------------------------------------------------------
// Pure Functions where inlining is known to give performance benefit: give the
// compiler visibility at call site, but leave final inline decision to optimizer

// Cutoff function
//
static
double fcut(const double xi)
{
  double a;
  if (xi >= 1.0)
    return 1.0;
  else if (xi <= 0.0)
    return 0.0;
  else {
    // ( 1.d0 - (1.d0 - xi)**4 )**2, but with better codegen
    a = 1.0 - xi;
    a *= a;
    a *= a;
    a = 1.0 - a;
    return a * a;
  }
}

// Cutoff function and derivative
//
static
double dfcut(const double xi, double &dfc)
{
  double a, a3, a4, a1m4;
  if (xi >= 1.0) {
    dfc = 0.0;
    return 1.0;
  } else if (xi <= 0.0) {
    dfc = 0.0;
    return 0.0;
  } else {
    a = 1.0 - xi;
    a3 = a * a * a;
    a4 = a * a3;
    a1m4 = 1.0 - a4;

    dfc = 8 * a1m4 * a3;
    return a1m4 * a1m4;
  }
}

// Screening ellipse, excluding multiple screening
//
static
double Csijk(const double C, const double Cmin, const double Cmax)
{
  double x;
  x = (C - Cmin) / (Cmax - Cmin);
  return fcut(x);
}

// Derivative of Cikj w.r.t. rij
//     Inputs: rij,rij2,rik2,rjk2
//
static
double dCfunc(const double rij2, const double rik2, const double rjk2)
{
  double rij4, a, asq, b, denom;

  rij4 = rij2 * rij2;
  a = rik2 - rjk2;
  b = rik2 + rjk2;
  asq = a * a;
  denom = rij4 - asq;
  denom = denom * denom;
  return -4 * (-2 * rij2 * asq + rij4 * b + asq * b) / denom;
}

// Derivative of Cikj w.r.t. rik and rjk
//     Inputs: rij,rij2,rik2,rjk2
//
static
void dCfunc2(const double rij2, const double rik2, const double rjk2, double &dCikj1,
             double &dCikj2)
{
  double rij4, rik4, rjk4, a, denom;

  rij4 = rij2 * rij2;
  rik4 = rik2 * rik2;
  rjk4 = rjk2 * rjk2;
  a = rik2 - rjk2;
  denom = rij4 - a * a;
  denom = denom * denom;
  dCikj1 = 4 * rij2 * (rij4 + rik4 + 2 * rik2 * rjk2 - 3 * rjk4 - 2 * rij2 * a) / denom;
  dCikj2 = 4 * rij2 * (rij4 - 3 * rik4 + 2 * rik2 * rjk2 + rjk4 + 2 * rij2 * a) / denom;
}

//-----------------------------------------------------------------------------
// Pure Functions where no benefit from inlining is expected
extern bool str_to_lat(const std::string & str, bool single, lattice_t& lat);
extern double zbl(const double r, const int z1, const int z2);
extern double erose(const double r, const double re, const double alpha, const double Ec,
                    const double repuls, const double attrac, const int form);

extern void get_shpfcn(const lattice_t latt, const double sthe, const double cthe,
                       double (&s)[3]);

extern int get_Zij(const lattice_t latt);
extern int get_Zij2(const lattice_t latt, const double cmin, const double cmax, const double sthe,
                    double &arat, double &S);
extern int get_Zij2_b2nn(const lattice_t latt, const double cmin, const double cmax, double &S);


}    // namespace MEAM_NS
}    // namespace LAMMPS_NS
#endif
