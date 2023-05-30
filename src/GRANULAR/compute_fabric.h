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
  ~ComputeFabric() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void compute_vector() override;
  double compute_scalar() override;

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
