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

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(fabric,ComputeFabric);
// clang-format on
#else

#ifndef LMP_COMPUTE_FABRIC_H
#define LMP_COMPUTE_FABRIC_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeFabric : public Compute {
 public:
  ComputeFabric(class LAMMPS *, int, char **);
  ~ComputeFabric();
  void init();
  void init_list(int, class NeighList *);
  void compute_vector();
  double compute_scalar();

 private:
  int ntensors, pstyle, cutstyle;
  double nc;
  int *tensor_style;
  int **type_filter;
  class NeighList *list;

  int cn_flag, br_flag, fn_flag, ft_flag;
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute fabric radius style requires atom attribute radius

Self-explanatory.

E: No pair style is defined for compute fabric

Self-explanatory.

E: Pair style does not support compute fabric normal or tangential force

Pair style must be single enabled to calculate the normal or tangential force tensors

E: Pair style does not calculate tangential forces for compute fabric

The tangential force tensor can only be calculated for granular pair styles with tangential forces

E: Compute fabric does not support pair styles that extend beyond contact

Granular pair styles that extend beyond contact such as JKR or DMT are not supported

*/
