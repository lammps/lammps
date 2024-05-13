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
PairStyle(gayberne,PairGayBerne);
// clang-format on
#else

#ifndef LMP_PAIR_GAYBERNE_H
#define LMP_PAIR_GAYBERNE_H

#include "pair.h"

namespace LAMMPS_NS {

class PairGayBerne : public Pair {
 public:
  PairGayBerne(LAMMPS *lmp);
  ~PairGayBerne() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  void write_data(FILE *) override;
  void write_data_all(FILE *) override;

 protected:
  enum { SPHERE_SPHERE, SPHERE_ELLIPSE, ELLIPSE_SPHERE, ELLIPSE_ELLIPSE };

  double cut_global;
  double **cut;

  double gamma, upsilon, mu;    // Gay-Berne parameters
  double **shape1;              // per-type radii in x, y and z
  double **shape2;              // per-type radii in x, y and z SQUARED
  double *lshape;               // precalculation based on the shape
  double **well;                // well depth scaling along each axis ^ -1.0/mu
  double **epsilon, **sigma;    // epsilon and sigma values for atom-type pairs

  int **form;
  double **lj1, **lj2, **lj3, **lj4;
  double **offset;
  int *setwell;
  class AtomVecEllipsoid *avec;

  void allocate();
  double gayberne_analytic(const int i, const int j, double a1[3][3], double a2[3][3],
                           double b1[3][3], double b2[3][3], double g1[3][3], double g2[3][3],
                           double *r12, const double rsq, double *fforce, double *ttor,
                           double *rtor);
  double gayberne_lj(const int i, const int j, double a1[3][3], double b1[3][3], double g1[3][3],
                     double *r12, const double rsq, double *fforce, double *ttor);
  void compute_eta_torque(double m[3][3], double m2[3][3], double *s, double ans[3][3]);
};

}    // namespace LAMMPS_NS
#endif
#endif
