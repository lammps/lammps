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
   Contributing author: Mingjian Wen (University of Minnesota)
   e-mail: wenxx151@umn.edu

   This implements the DRIP model as described in
   M. Wen, S. Carr, S. Fang, E. Kaxiras, and E. B. Tadmor,
   Phys. Rev. B, 98, 235404 (2018).
------------------------------------------------------------------------- */


#ifdef PAIR_CLASS

PairStyle(drip, PairDRIP)

#else

#ifndef LMP_PAIR_DRIP_H
#define LMP_PAIR_DRIP_H

#include "pair.h"
#include "my_page.h"
#include <cmath>

namespace LAMMPS_NS {

#define DIM 3
typedef double   V3[3];


class PairDRIP : public Pair {
public:
  PairDRIP(class LAMMPS *);
  virtual ~PairDRIP();

  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);
  void init_style();

protected:
  struct Param
  {
    int    ielement, jelement;
    double C0, C2, C4, C, delta, lambda, A, z0, B, eta, rhocut, rcut;
    double rhocutsq, rcutsq;
  };
  Param *params;         // parameter set for I-J interactions
  int **nearest3neigh;   // nearest 3 neighbors of atoms
  char **elements;       // names of unique elements
  int **elem2param;      // mapping from element pairs to parameters
  int *map;              // mapping from atom types to elements
  int nelements;         // # of unique elements
  double cutmax;         // max cutoff for all species

  void read_file(char *);
  void allocate();

  // DRIP specific functions
  double calc_attractive(int const, int const, Param&, double const,
      double const *, double *const, double *const);

  double calc_repulsive(int const, int const, Param&, double const,
      double const *, double const *, V3 const *, V3 const *, V3 const *,
      V3 const *, double *const, double *const);

  void find_nearest3neigh();

  void calc_normal(int const, double *const, V3 *const, V3 *const, V3 *const,
      V3 *const);

  void get_drhosqij(double const *, double const *, V3 const *, V3 const *,
      V3 const *, V3 const *, double *const, double *const, double *const,
      double *const, double *const);

  double td(double, double, double, double, double const *const, double,
      const double *const, double&, double&);

  double dihedral(const int, const int, Param&, double const, double&,
      double *const, double *const, double *const, double *const, double *const,
      double *const, double *const, double *const);

  double deriv_cos_omega(double const *, double const *, double const *,
      double const *, double *const, double *const, double *const,
      double *const);

  double tap(double, double, double&);

  double tap_rho(double, double, double&);

  void deriv_cross(double const *, double const *, double const *,
      double *const, V3 *const, V3 *const, V3 *const);

  // inline functions
  inline double dot(double const *x, double const *y)
  {
    return x[0]*y[0]+x[1]*y[1]+x[2]*y[2];
  }

  inline void mat_dot_vec(V3 const *X, double const *y, double *const z)
  {
    for (int k = 0; k < 3; k++) {
      z[k] = X[k][0]*y[0]+X[k][1]*y[1]+X[k][2]*y[2];
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

E: No enough neighbors to construct normal

Cannot find three neighbors within cutoff of the target atom.
Check the configuration.

*/
