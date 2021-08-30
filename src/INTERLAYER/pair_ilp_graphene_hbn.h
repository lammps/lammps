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
PairStyle(ilp/graphene/hbn,PairILPGrapheneHBN);
// clang-format on
#else

#ifndef LMP_PAIR_ILP_GRAPHENE_HBN_H
#define LMP_PAIR_ILP_GRAPHENE_HBN_H

#include "pair.h"

namespace LAMMPS_NS {

class PairILPGrapheneHBN : public Pair {
 public:
  PairILPGrapheneHBN(class LAMMPS *);
  virtual ~PairILPGrapheneHBN();

  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);
  void init_style();
  void ILP_neigh();
  void calc_normal();
  void calc_FRep(int, int);
  void calc_FvdW(int, int);
  double single(int, int, int, int, double, double, double, double &);

  static constexpr int NPARAMS_PER_LINE = 13;

 protected:
  int me;
  int maxlocal;            // size of numneigh, firstneigh arrays
  int pgsize;              // size of neighbor page
  int oneatom;             // max # of neighbors for one atom
  MyPage<int> *ipage;      // neighbor list pages
  int *ILP_numneigh;       // # of pair neighbors for each atom
  int **ILP_firstneigh;    // ptr to 1st neighbor of each atom
  int tap_flag;            // flag to turn on/off taper function

  struct Param {
    double z0, alpha, epsilon, C, delta, d, sR, reff, C6, S;
    double delta2inv, seff, lambda, rcut;
    int ielement, jelement;
  };
  Param *params;    // parameter set for I-J interactions
  int nmax;         // max # of atoms

  double cut_global;
  double cut_normal;
  double **cutILPsq;    // mapping the cutoff from element pairs to parameters
  double **offset;
  double **normal;
  double ***dnormdri;
  double ****dnormal;

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
