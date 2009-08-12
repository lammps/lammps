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

#ifndef PAIR_GAYBERNE_H
#define PAIR_GAYBERNE_H

#include "pair.h"

namespace LAMMPS_NS {

class PairGayBerne : public Pair {
 public:
  PairGayBerne(LAMMPS *lmp);
  virtual ~PairGayBerne();
  virtual void compute(int, int);
  virtual void settings(int, char **);
  void coeff(int, char **);
  virtual void init_style();
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);

 protected:
  double cut_global;
  double **cut;

  double gamma,upsilon,mu;   // Gay-Berne parameters
  double **shape;            // radii in x, y and z SQUARED
  double *lshape;            // precalculation based on the shape
  double **well;             // well depth scaling along each axis ^ -1.0/mu
  double **epsilon,**sigma;  // epsilon and sigma values for atom-type pairs

  int **form;
  double **lj1,**lj2,**lj3,**lj4;
  double **offset;
  int *setwell;

  void allocate();
  double gayberne_analytic(const int i, const int j, double a1[3][3],
                           double a2[3][3], double b1[3][3], double b2[3][3],
                           double g1[3][3], double g2[3][3], double *r12, 
                           const double rsq, double *fforce, double *ttor, 
                           double *rtor);
  double gayberne_lj(const int i, const int j, double a1[3][3],
                     double b1[3][3],double g1[3][3],double *r12, 
		     const double rsq, double *fforce, double *ttor);
  void compute_eta_torque(double m[3][3], double m2[3][3],
                          double *s, double ans[3][3]);
};

}
#endif
