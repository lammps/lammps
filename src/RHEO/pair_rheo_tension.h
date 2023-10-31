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
PairStyle(rheo/tension,PairRHEOTension)
// clang-format on
#else

#ifndef LMP_PAIR_RHEO_TENSION_H
#define LMP_PAIR_RHEO_TENSION_H

#include "pair.h"

namespace LAMMPS_NS {

class PairRHEOTension : public Pair {
 public:
  PairRHEOTension(class LAMMPS *);
  ~PairRHEOTension() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void setup() override;
  void init_style() override;
  double init_one(int, int) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;

 protected:
  int nmax_store;
  double **nt, *ct;
  double **alpha;
  double h, hsq, hinv, hinv3;

  void allocate();

  class ComputeRHEOKernel *compute_kernel;
  class FixRHEO *fix_rheo;
};

}    // namespace LAMMPS_NS

#endif
#endif
