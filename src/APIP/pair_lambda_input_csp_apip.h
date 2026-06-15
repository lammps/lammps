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
PairStyle(lambda/input/csp/apip,PairLambdaInputCSPAPIP);
// clang-format on
#else

#ifndef LMP_PAIR_LAMBDA_INPUT_CSP_APIP_H
#define LMP_PAIR_LAMBDA_INPUT_CSP_APIP_H

#include "pair_lambda_input_apip.h"

namespace LAMMPS_NS {

class PairLambdaInputCSPAPIP : virtual public PairLambdaInputAPIP {

 public:
  PairLambdaInputCSPAPIP(class LAMMPS *);
  ~PairLambdaInputCSPAPIP() override;
  void settings(int, char **) override;

 protected:
  // csp variables
  double cut_csp_sq;    ///< squared cutoff
  int maxneigh;         ///< number of atoms for which distsq and nearest are allocated
  int nnn;              ///< number of nearest neighbours that are used in the csp calculation
  int nnn_buffer;       ///< number of additional (to nnn) stored nearest neighbours
  double *distsq;       ///< distance sq to each neighbor
  int *nearest;         ///< atom indices of neighbors

  int calculate_lambda_input() override;
  void select(int, int, double *);
  void select2(int, int, double *, int *);
};

}    // namespace LAMMPS_NS

#endif
#endif
