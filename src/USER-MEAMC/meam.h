#ifndef LMP_MEAM_H
#define LMP_MEAM_H

#include "memory.h"
#include <stdlib.h>

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
  int vind2D[3][3], vind3D[3][3][3];
  int v2D[6], v3D[10];

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
  void meam_checkindex(int, int, int, int*, int*);
  void G_gam(double, int, double*, int*);
  void dG_gam(double, int, double*, double*);
  void getscreen(int i, double* scrfcn, double* dscrfcn, double* fcpair, double** x, int numneigh,
                 int* firstneigh, int numneigh_full, int* firstneigh_full, int ntype, int* type, int* fmap);
  void screen(int i, int j, double** x, double rijsq, double* sij, int numneigh_full, int* firstneigh_full,
              int ntype, int* type, int* fmap);
  void calc_rho1(int i, int ntype, int* type, int* fmap, double** x, int numneigh, int* firstneigh,
                 double* scrfcn, double* fcpair);
  void dsij(int i, int j, int k, int jn, int numneigh, double rij2, double* dsij1, double* dsij2, int ntype,
            int* type, int* fmap, double** x, double* scrfcn, double* fcpair);
  void fcut(double, double*);
  void dfcut(double, double*, double*);
  void dCfunc(double, double, double, double*);
  void dCfunc2(double, double, double, double*, double*);

  void alloyparams();
  void compute_pair_meam();
  double phi_meam(double, int, int);
  void compute_reference_density();
  void get_shpfcn(double*, lattice_t);
  void get_tavref(double*, double*, double*, double*, double*, double*, double, double, double, double,
                  double, double, double, int, int, lattice_t);
  void get_Zij(int*, lattice_t);
  void get_Zij2(int*, double*, double*, lattice_t, double, double);
  void get_sijk(double, int, int, int, double*);
  void get_densref(double, int, int, double*, double*, double*, double*, double*, double*, double*, double*);
  double zbl(double, int, int);
  double erose(double, double, double, double, double, double, int);
  void interpolate_meam(int);
  double compute_phi(double, int, int);

public:
  void meam_setup_global(int nelt, lattice_t* lat, double* z, int* ielement, double* atwt, double* alpha,
                         double* b0, double* b1, double* b2, double* b3, double* alat, double* esub,
                         double* asub, double* t0, double* t1, double* t2, double* t3, double* rozero,
                         int* ibar);
  void meam_setup_param(int which, double value, int nindex, int* index /*index(3)*/, int* errorflag);
  void meam_setup_done(double* cutmax);
  void meam_dens_setup(int atom_nmax, int nall, int n_neigh);
  void meam_dens_init(int i, int ntype, int* type, int* fmap, double** x, int numneigh, int* firstneigh,
                      int numneigh_full, int* firstneigh_full, int fnoffset, int* errorflag);
  void meam_dens_final(int nlocal, int eflag_either, int eflag_global, int eflag_atom, double* eng_vdwl,
                       double* eatom, int ntype, int* type, int* fmap, int* errorflag);
  void meam_force(int i, int eflag_either, int eflag_global, int eflag_atom, int vflag_atom, double* eng_vdwl,
                  double* eatom, int ntype, int* type, int* fmap, double** x, int numneigh, int* firstneigh,
                  int numneigh_full, int* firstneigh_full, int fnoffset, double** f, double** vatom,
                  int* errorflag);
};

// Functions we need for compat

#define iszero(f) (fabs(f) < 1e-20)

#define setall2d(arr, v)                                                                                     \
  {                                                                                                          \
    for (int __i = 0; __i < maxelt; __i++)                                                                   \
      for (int __j = 0; __j < maxelt; __j++)                                                                 \
        arr[__i][__j] = v;                                                                                   \
  }
#define setall3d(arr, v)                                                                                     \
  {                                                                                                          \
    for (int __i = 0; __i < maxelt; __i++)                                                                   \
      for (int __j = 0; __j < maxelt; __j++)                                                                 \
        for (int __k = 0; __k < maxelt; __k++)                                                               \
          arr[__i][__j][__k] = v;                                                                            \
  }
};
#endif
