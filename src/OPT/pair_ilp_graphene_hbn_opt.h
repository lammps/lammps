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

class PairILPGrapheneHBNOpt : virtual public PairILPGrapheneHBN {
 public:
  PairILPGrapheneHBNOpt(class LAMMPS *);
  ~PairILPGrapheneHBNOpt() override;

  void compute(int, int) override;
  void init_style() override;

 protected:
  void update_internal_list();
  template <int MAX_NNEIGH>
  void calc_normal(int i, int *ILP_neigh, int nneigh, double *normal, double (*dnormdri)[3],
                   double (*dnormdrk)[3][3]);
  template <int MAX_NNEIGH, int EFLAG, int VFLAG_EITHER, int TAP_FLAG, int VARIANT = ILP_GrhBN>
  void eval();
  int *layered_neigh;
  int **first_layered_neigh;
  int *special_type;
  int *num_intra, *num_inter, *num_vdw;
  int inum_max, jnum_max;
};

}    // namespace LAMMPS_NS

#endif
#endif
