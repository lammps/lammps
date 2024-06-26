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
ComputeStyle(pod/global,ComputePODGlobal);
// clang-format on
#else

#ifndef LMP_COMPUTE_POD_GLOBAL_H
#define LMP_COMPUTE_POD_GLOBAL_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputePODGlobal : public Compute {
 public:
  ComputePODGlobal(class LAMMPS *, int, char **);
  ~ComputePODGlobal() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void compute_array() override;
  double memory_usage() override;
  void lammpsNeighborList(double **x, int **firstneigh, tagint *atomid, int *atomtype,
                          int *numneigh, double rcutsq, int i);
  void map_element2type(int narg, char **arg, int nelements);

 private:
  class NeighList *list;
  class EAPOD *podptr;
  double **pod;
  double cutmax;
  int nij;
  int nijmax;

  double *tmpmem;    // temporary memory
  double *rij;       // (xj - xi) for all pairs (I, J)
  char **elements;
  int *map;
  int *ai;    // IDs of atoms I for all pairs (I, J)
  int *aj;    // IDs of atoms J for all pairs (I, J)
  int *ti;    // types of atoms I for all pairs (I, J)
  int *tj;    // types of atoms J  for all pairs (I, J)
};

}    // namespace LAMMPS_NS

#endif
#endif
