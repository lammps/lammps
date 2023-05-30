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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(resquared,PairRESquared);
// clang-format on
#else

#ifndef LMP_PAIR_RESQUARED_H
#define LMP_PAIR_RESQUARED_H

#include "pair.h"

namespace LAMMPS_NS {

class PairRESquared : public Pair {
 public:
  PairRESquared(LAMMPS *lmp);
  ~PairRESquared() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;

 protected:
  enum { SPHERE_SPHERE, SPHERE_ELLIPSE, ELLIPSE_SPHERE, ELLIPSE_ELLIPSE };

  double cut_global;
  double **cut;

  double **shape1;              // per-type radii in x, y and z
  double **shape2;              // per-type radii in x, y and z SQUARED
  double *lshape;               // product of the radii
  double **well;                // well depth scaling along each axis
  double **epsilon, **sigma;    // epsilon and sigma values for atom-type pairs

  int **form;
  double **lj1, **lj2, **lj3, **lj4;
  double **offset;
  int *setwell;
  class AtomVecEllipsoid *avec;

  // per-particle temporaries for RE-squared calculation

  struct RE2Vars {
    // per particle precomputations for energy, force, torque

    double A[3][3];        // Rotation matrix (lab->body)
    double aTe[3][3];      // A'*E
    double gamma[3][3];    // A'*S^2*A

    // per particle precomputations for torque

    double sa[3][3];          // S^2*A;
    double lA[3][3][3];       // -A*rotation generator (x,y, or z)
    double lAtwo[3][3][3];    // A'*S^2*lA
    double lAsa[3][3][3];     // lAtwo+lA'*sa
  };

  void allocate();

  void precompute_i(const int i, RE2Vars &ws);
  double det_prime(const double m[3][3], const double m2[3][3]);
  double resquared_analytic(const int i, const int j, const RE2Vars &wi, const RE2Vars &wj,
                            const double *r, const double rsq, double *fforce, double *ttor,
                            double *rtor);
  double resquared_lj(const int i, const int j, const RE2Vars &wi, const double *r,
                      const double rsq, double *fforce, double *ttor, bool calc_torque);

  double cr60;        // 60^1/3
  double b_alpha;     // 45/56
  double solv_f_a;    // 3.0/(4.0*PI*-36)
  double solv_f_r;    // 3.0/(4.0*PI*2025)
};

}    // namespace LAMMPS_NS
#endif
#endif
