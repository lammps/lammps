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
/* ----------------------------------------------------------------------
   Contributing author: David Immel (d.immel@fz-juelich.de, FZJ, Germany)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(lambda/zone/apip,PairLambdaZoneAPIP);
// clang-format on
#else

#ifndef LMP_PAIR_LAMBDA_ZONE_APIP_H
#define LMP_PAIR_LAMBDA_ZONE_APIP_H

#include "pair.h"

namespace LAMMPS_NS {

class PairLambdaZoneAPIP : public Pair {
  friend class FixLambdaAPIP;

 public:
  PairLambdaZoneAPIP(class LAMMPS *);
  ~PairLambdaZoneAPIP() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void *extract(const char *, int &) override;

 protected:
  // pro forma pair style variables
  double **cut;

  // lambda calculation variables
  double timer;            ///< accumulated compute time
  double time_per_atom;    ///< compute time for one atom
  int n_calculations;      ///< number of accumulated computations

  double cut_global;    ///< used cutoff
  double cut_lo;    ///< distance at which the cutoff function of the transition zone decays from 1
  double cut_hi;    ///< distance at which the cutoff function of the transition zone is 0
  double cut_width;           ///< cut_hi - cut_lo
  double cut_hi_sq;           ///< cut_hi_sq * cut_hi_sq
  double lambda_non_group;    ///< lambda for atoms that are not in the group of the fix
  int groupbit;               ///< group for which lambda is calculated

  // variables for calculation
  double *lambda_ta;    ///< time averaged lambda input
  int nmax_ta;          ///< size of lambda_ta

  virtual void allocate();
  void calculate_time_per_atom();
  void calculate_lambda();
  double switching_function_poly_distance(double);
};

}    // namespace LAMMPS_NS

#endif
#endif
