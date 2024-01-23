/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
-------------------------------------------------------------------------
   pair_LocalDensity written by:
   Tanmoy Sanyal and M. Scott Shell from UC Santa Barbara
   David Rosenberger: TU Darmstadt
-------------------------------------------------------------------------*/

#ifdef PAIR_CLASS
// clang-format off
PairStyle(local/density,PairLocalDensity);
// clang-format on
#else

#ifndef LMP_PAIR_LOCAL_DENSITY_H
#define LMP_PAIR_LOCAL_DENSITY_H

#include "pair.h"

namespace LAMMPS_NS {

class PairLocalDensity : public Pair {
 public:
  PairLocalDensity(class LAMMPS *);
  ~PairLocalDensity() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  double single(int, int, int, int, double, double, double, double &) override;

  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  double memory_usage() override;

 protected:
  //------------------------------------------------------------------------
  //This information is read from the tabulated input file

  int nLD, nrho;                            // number of LD types
  int **a, **b;                             // central and neigh atom filters
  double *uppercut, *lowercut;              // upper and lower cutoffs
  double *uppercutsq, *lowercutsq;          // square of above cutoffs
  double *c0, *c2, *c4, *c6;                // coeffs for indicator function
  double *rho_min, *rho_max, *delta_rho;    // min, max & grid-size for LDs
  double **rho, **frho;                     // LD and LD function tables

  //------------------------------------------------------------------------

  double ***frho_spline;    // splined LD potentials
  double cutmax;            // max cutoff for all elements
  double cutforcesq;        // square of global upper cutoff

  int nmax;             // max size of per-atom arrays
  double **localrho;    // per-atom LD
  double **fp;          // per-atom LD potential function derivative

  void allocate();

  // read tabulated input file
  void parse_file(char *);

  // convert array to spline
  void array2spline();

  // cubic spline interpolation
  void interpolate_cbspl(int, double, double *, double **);
};

}    // namespace LAMMPS_NS

#endif
#endif
