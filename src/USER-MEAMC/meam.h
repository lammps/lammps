#ifndef LMP_MEAM_H
#define LMP_MEAM_H

#include "memory.h"
#include <cmath>

#define maxelt 5

namespace LAMMPS_NS {

typedef enum { FCC, BCC, HCP, DIM, DIA, B1, C11, L12, B2 } lattice_t;

class MEAM
{
public:
  MEAM(Memory* mem);
  ~MEAM();

private:
  Memory* memory;

  // cutforce = force cutoff
  // cutforcesq = force cutoff squared

  double cutforce, cutforcesq;

  // Ec_meam = cohesive energy
  // re_meam = nearest-neighbor distance
  // Omega_meam = atomic volume
  // B_meam = bulk modulus
  // Z_meam = number of first neighbors for reference structure
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

  // nr,dr = pair function discretization parameters
  // nrar,rdrar = spline coeff array parameters

  double Ec_meam[maxelt][maxelt], re_meam[maxelt][maxelt];
  double Omega_meam[maxelt], Z_meam[maxelt];
  double A_meam[maxelt], alpha_meam[maxelt][maxelt], rho0_meam[maxelt];
  double delta_meam[maxelt][maxelt];
  double beta0_meam[maxelt], beta1_meam[maxelt];
  double beta2_meam[maxelt], beta3_meam[maxelt];
  double t0_meam[maxelt], t1_meam[maxelt];
  double t2_meam[maxelt], t3_meam[maxelt];
  double rho_ref_meam[maxelt];
  int ibar_meam[maxelt], ielt_meam[maxelt];
  lattice_t lattce_meam[maxelt][maxelt];
  int nn2_meam[maxelt][maxelt];
  int zbl_meam[maxelt][maxelt];
  int eltind[maxelt][maxelt];
  int neltypes;

  double** phir;

  double **phirar, **phirar1, **phirar2, **phirar3, **phirar4, **phirar5, **phirar6;

  double attrac_meam[maxelt][maxelt], repuls_meam[maxelt][maxelt];

  double Cmin_meam[maxelt][maxelt][maxelt];
  double Cmax_meam[maxelt][maxelt][maxelt];
  double rc_meam, delr_meam, ebound_meam[maxelt][maxelt];
  int augt1, ialloy, mix_ref_t, erose_form;
  int emb_lin_neg, bkgd_dyn;
  double gsmooth_factor;

  int vind2D[3][3], vind3D[3][3][3];                  // x-y-z to Voigt-like index
  int v2D[6], v3D[10];                                // multiplicity of Voigt index (i.e. [1] -> xy+yx = 2

  int nr, nrar;
  double dr, rdrar;

public:
  int nmax;
  double *rho, *rho0, *rho1, *rho2, *rho3, *frhop;
  double *gamma, *dgamma1, *dgamma2, *dgamma3, *arho2b;
  double **arho1, **arho2, **arho3, **arho3b, **t_ave, **tsq_ave;

  int maxneigh;
  double *scrfcn, *dscrfcn, *fcpair;

protected:
  // meam_funcs.cpp

  //-----------------------------------------------------------------------------
  // Cutoff function
  //
  static double fcut(const double xi) {
    double a;
    if (xi >= 1.0)
      return 1.0;
    else if (xi <= 0.0)
      return 0.0;
    else {
      // ( 1.d0 - (1.d0 - xi)**4 )**2, but with better codegen
      a = 1.0 - xi;
      a *= a; a *= a;
      a = 1.0 - a;
      return a * a;
    }
  }

  //-----------------------------------------------------------------------------
  // Cutoff function and derivative
  //
  static double dfcut(const double xi, double& dfc) {
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
      a1m4 = 1.0-a4;

      dfc = 8 * a1m4 * a3;
      return a1m4*a1m4;
    }
  }

  //-----------------------------------------------------------------------------
  // Derivative of Cikj w.r.t. rij
  //     Inputs: rij,rij2,rik2,rjk2
  //
  static double dCfunc(const double rij2, const double rik2, const double rjk2) {
    double rij4, a, asq, b,denom;

    rij4 = rij2 * rij2;
    a = rik2 - rjk2;
    b = rik2 + rjk2;
    asq = a*a;
    denom = rij4 - asq;
    denom = denom * denom;
    return -4 * (-2 * rij2 * asq + rij4 * b + asq * b) / denom;
  }

  //-----------------------------------------------------------------------------
  // Derivative of Cikj w.r.t. rik and rjk
  //     Inputs: rij,rij2,rik2,rjk2
  //
  static void dCfunc2(const double rij2, const double rik2, const double rjk2,
               double& dCikj1, double& dCikj2) {
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

  double G_gam(const double gamma, const int ibar, int &errorflag) const;
  double dG_gam(const double gamma, const int ibar, double &dG) const;
  static double zbl(const double r, const int z1, const int z2);
  double embedding(const double A, const double Ec, const double rhobar, double& dF) const;
  static double erose(const double r, const double re, const double alpha, const double Ec, const double repuls, const double attrac, const int form);

  static void get_shpfcn(const lattice_t latt, double (&s)[3]);
  static int get_Zij(const lattice_t latt);
  static int get_Zij2(const lattice_t latt, const double cmin, const double cmax, double &a, double &S);
protected:
  void meam_checkindex(int, int, int, int*, int*);
  void getscreen(int i, double* scrfcn, double* dscrfcn, double* fcpair, double** x, int numneigh,
                 int* firstneigh, int numneigh_full, int* firstneigh_full, int ntype, int* type, int* fmap);
  void calc_rho1(int i, int ntype, int* type, int* fmap, double** x, int numneigh, int* firstneigh,
                 double* scrfcn, double* fcpair);

  void alloyparams();
  void compute_pair_meam();
  double phi_meam(double, int, int);
  void compute_reference_density();
  void get_tavref(double*, double*, double*, double*, double*, double*, double, double, double, double,
                  double, double, double, int, int, lattice_t);
  void get_sijk(double, int, int, int, double*);
  void get_densref(double, int, int, double*, double*, double*, double*, double*, double*, double*, double*);
  void interpolate_meam(int);

public:
  void meam_setup_global(int nelt, lattice_t* lat, double* z, int* ielement, double* atwt, double* alpha,
                         double* b0, double* b1, double* b2, double* b3, double* alat, double* esub,
                         double* asub, double* t0, double* t1, double* t2, double* t3, double* rozero,
                         int* ibar);
  void meam_setup_param(int which, double value, int nindex, int* index /*index(3)*/, int* errorflag);
  void meam_setup_done(double* cutmax);
  void meam_dens_setup(int atom_nmax, int nall, int n_neigh);
  void meam_dens_init(int i, int ntype, int* type, int* fmap, double** x, int numneigh, int* firstneigh,
                      int numneigh_full, int* firstneigh_full, int fnoffset);
  void meam_dens_final(int nlocal, int eflag_either, int eflag_global, int eflag_atom, double* eng_vdwl,
                       double* eatom, int ntype, int* type, int* fmap, int& errorflag);
  void meam_force(int i, int eflag_either, int eflag_global, int eflag_atom, int vflag_atom, double* eng_vdwl,
                  double* eatom, int ntype, int* type, int* fmap, double** x, int numneigh, int* firstneigh,
                  int numneigh_full, int* firstneigh_full, int fnoffset, double** f, double** vatom);
};

// Functions we need for compat

static inline bool iszero(const double f) {
  return fabs(f) < 1e-20;
}

template <typename TYPE, size_t maxi, size_t maxj>
static inline void setall2d(TYPE (&arr)[maxi][maxj], const TYPE v) {
  for (size_t i = 0; i < maxi; i++)
    for (size_t j = 0; j < maxj; j++)
      arr[i][j] = v;
}

template <typename TYPE, size_t maxi, size_t maxj, size_t maxk>
static inline void setall3d(TYPE (&arr)[maxi][maxj][maxk], const TYPE v) {
  for (size_t i = 0; i < maxi; i++)
    for (size_t j = 0; j < maxj; j++)
      for (size_t k = 0; k < maxk; k++)
        arr[i][j][k] = v;
}

// Helper functions

static inline double fdiv_zero(const double n, const double d) {
  if (iszero(d))
    return 0.0;
  return n / d;
}

}
#endif
