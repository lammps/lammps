/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(ilp/graphene/hbn/opt,PairILPGrapheneHBNOpt);
// clang-format on
#else

#ifndef LMP_PAIR_ILP_GRAPHENE_HBN_OPT_H
#define LMP_PAIR_ILP_GRAPHENE_HBN_OPT_H

#include "pair.h"

namespace LAMMPS_NS {

class PairILPGrapheneHBNOpt : public Pair {
 public:
  PairILPGrapheneHBNOpt(class LAMMPS *);
  virtual ~PairILPGrapheneHBNOpt();

  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);
  void init_style();
  void calc_single_normal(int i, int *ILP_neigh, int nneigh, double *normal, double (*dnormdri)[3],
                          double (*dnormdrk)[3][3]);
  void update_internal_list();
  double single(int, int, int, int, double, double, double, double &);
  template <int, int, int> void eval();
  static constexpr int NPARAMS_PER_LINE = 13;

 protected:
  int me;
  int maxlocal;    // size of numneigh, firstneigh arrays
  int tap_flag;    // flag to turn on/off taper function

  struct Param {
    double z0, alpha, epsilon, C, delta, d, sR, reff, C6, S;
    double delta2inv, seff, lambda, rcut, seffinv;
    int ielement, jelement;
  };
  Param *params;    // parameter set for I-J interactions
  int nmax;         // max # of atoms

  double cut_global;
  double cut_normal;
  double **cutILPsq;    // mapping the cutoff from element pairs to parameters
  double **offset;
  int *layered_neigh;
  int **first_layered_neigh;
  int *num_intra, *num_inter, *num_vdw;
  int inum_max, jnum_max;
  void read_file(char *);
  void allocate();
};

}    // namespace LAMMPS_NS

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
