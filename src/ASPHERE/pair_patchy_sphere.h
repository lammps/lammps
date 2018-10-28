/* -*- c++ -*- ----------------------------------------------------------
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

PairStyle(patchysphere,PairPatchySphere)

#else

#ifndef LMP_PAIR_PATCHY_SPHERE_H
#define LMP_PAIR_PATCHY_SPHERE_H

#include "pair.h"

namespace LAMMPS_NS {

class PairPatchySphere : public Pair {
 public:
  PairPatchySphere(LAMMPS *);
  virtual ~PairPatchySphere();
  virtual void compute(int, int);
  virtual void settings(int, char **);
  void coeff(int, char **);
  virtual void init_style();
  double init_one(int, int);
	//void write_restart(FILE *);
	//void read_restart(FILE *);
	//void write_restart_settings(FILE *);
	//void read_restart_settings(FILE *);
	//void write_data(FILE *);
	//void write_data_all(FILE *);

 protected:
  // enumerate the "encounter" possibilities.
  enum{SPHERE_SPHERE,SPHERE_PATCHY,PATCHY_SPHERE,PATCHY_PATCHY};

  double cut_global;
  double **cut;
  double **inner_cut;
  double **offset;
  double **epsilon, **sigma;


  // Whole bunch of stuff related to bond vectors...
  double **epsilon_b; // Patch binding strength.
  int n_bond_pairs;   // Number of bond pairs (combinations of alpha, beta, ...)
  int **bond_pairs;   // Enumerate bond pairs.

  double **phi_m, **theta_m; // Theta and phi scale factors.

  int n_bonds;            // Number of bond vectors.
  double **bond_vectors;  // Bond vectors.

  class AtomVecEllipsoid *avec;

  bool debug_output;


  void allocate();
  double bond_vec_analytic(const int i, const int j, double a1[3][3],
                           double a2[3][3], double *r12, const double rsq,
                           double *fforce, double *ttor, double *rtor);


  double bond_vec_theta_part(const int i, const int j,
                             const int itype, const int jtype,
                             const double *r12, const double *b_alp,
                             const double *b_bet, double *fforce, double *ttor,
                             double *rtor);
  double bond_vec_phi_part(const int i, const int j,
                           const int itype, const int jtype,
                           const double *r12, const double rsq,
                           const double *b_alp, const double *b_bet,
                           const double *b_gam, const double *b_eps,
                           double *fforce, double *ttor, double *rtor);

};


} // namespace

#endif // LMP_PAIR_PATCHY_SPHERE_H
#endif // PAIR_CLASS
