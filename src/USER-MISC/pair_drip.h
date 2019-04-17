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

#ifdef PAIR_CLASS

PairStyle(drip,PairDRIP)

#else

#ifndef LMP_PAIR_DRIP_H
#define LMP_PAIR_DRIP_H

#include "pair.h"
#include "my_page.h"
#include <cmath>

namespace LAMMPS_NS {

#define DIM 3
typedef double V3[3];


class PairDRIP : public Pair {
 public:
  PairDRIP(class LAMMPS *);
  virtual ~PairDRIP();

  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);
  void init_style();
//  int pack_forward_comm(int, int *, double *, int, int *);
//  void unpack_forward_comm(int, int, double *);

 protected:
  double cutmax;                   // max cutoff for all species
  int me;
  int maxlocal;                    // size of numneigh, firstneigh arrays
  int pgsize;                      // size of neighbor page
  int oneatom;                     // max # of neighbors for one atom
  MyPage<int> *ipage;              // neighbor list pages


  struct Param {
    int ielement,jelement;
    double C0,C2,C4,C,delta,lambda,A,z0,B,eta,rhocut,rcut;
    double rhocutsq, rcutsq;
    double delta2inv,z06;
  };
  Param *params;         // parameter set for I-J interactions
  char **elements;       // names of unique elements
  int **elem2param;      // mapping from element pairs to parameters
  int *map;              // mapping from atom types to elements
  int nelements;         // # of unique elements
  int nparams;           // # of stored parameter sets
  int maxparam;          // max # of parameter sets
  int nmax;              // max # of atoms
  int ** nearest3neigh;  // nearest 3 neighbors of atoms

  void read_file( char * );
  void allocate();

  // DRIP specific functions
  double calc_attractive(int const i, int const j, Param& p,
      double const rsq, double const * rvec, double * const fi, double * const fj);

 double calc_repulsive(int const i, int const j,
     Param& p, double const rsq, double const * rvec,
     int const nbi1, int const nbi2, int const nbi3, double const * ni,
     V3 const * dni_dri, V3 const * dni_drnb1, V3 const * dni_drnb2,
     V3 const * dni_drnb3, double * const fi, double * const fj);


  void find_nearest3neigh();


 void calc_normal(int const i, int& k1, int& k2, int& k3,
     double * const normal, V3 *const dn_dri, V3 *const dn_drk1,
     V3 *const dn_drk2, V3 *const dn_drk3);



void get_drhosqij( double const* rij, double const* ni,
    V3 const* dni_dri, V3 const* dni_drn1,
    V3 const* dni_drn2, V3 const* dni_drn3,
    double* const drhosq_dri, double* const drhosq_drj,
    double* const drhosq_drn1, double* const drhosq_drn2,
    double* const drhosq_drn3);


  double td(double C0, double C2, double C4, double delta,
      double const* const rvec, double r,
      const double* const n,
      double& rho_sq, double& dtd);

  double dihedral(
      const int i, const int j, Param& p, double const rhosq, double& d_drhosq,
      double* const d_dri, double* const d_drj,
      double* const d_drk1, double* const d_drk2, double* const d_drk3,
      double* const d_drl1, double* const d_drl2, double* const d_drl3);

  double deriv_cos_omega( double const* rk, double const* ri,
      double const* rj, double const* rl, double* const dcos_drk,
      double* const dcos_dri, double* const dcos_drj, double* const dcos_drl);

  double tap(double r, double cutoff, double& dtap);

  double tap_rho(double rhosq, double cut_rhosq, double& drhosq);

void deriv_cross( double const* rk, double const* rl, double const* rm,
    double* const cross, V3 *const dcross_drk,
    V3 *const dcross_drl, V3 *const dcross_drm);



  inline double dot(double const* x, double const* y) {
    return x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
  }


  inline void mat_dot_vec(V3 const* X, double const* y, double* const z)
  {
    for (int k = 0; k < 3; k++) {
      z[k] = X[k][0] * y[0] + X[k][1] * y[1] + X[k][2] * y[2];
    }
  }


};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

*/

