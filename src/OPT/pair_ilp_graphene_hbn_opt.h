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

#include "pair_ilp_graphene_hbn.h"

namespace LAMMPS_NS {

class PairILPGrapheneHBNOpt : public PairILPGrapheneHBN {
 public:
  PairILPGrapheneHBNOpt(class LAMMPS *);
  ~PairILPGrapheneHBNOpt() override;

  void compute(int, int) override;
  void init_style() override;

 protected:
  void computeILP(int, int);
  inline void cross_deriv(double *pv, double (*dpv)[3][3], double (*vet)[3], int j, int k, int l);
  inline void calc_dnormal(double (*dnormal)[3][3], double (*dn1)[3][3], double *n1, double nn,
                           double nn2);
  inline void ev_tally_buffer(int i, int j, double *ei, double *vi, double *ej, double *vj,
                              int nlocal, int newton_pair, double evdwl, double ecoul, double fx,
                              double fy, double fz, double delx, double dely, double delz);
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
