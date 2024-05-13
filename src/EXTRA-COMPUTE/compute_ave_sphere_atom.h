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
ComputeStyle(ave/sphere/atom,ComputeAveSphereAtom);
// clang-format on
#else

#ifndef LMP_COMPUTE_AVE_SPHERE_ATOM_H
#define LMP_COMPUTE_AVE_SPHERE_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeAveSphereAtom : public Compute {
 public:
  ComputeAveSphereAtom(class LAMMPS *, int, char **);
  ~ComputeAveSphereAtom() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void compute_peratom() override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  double memory_usage() override;

 protected:
  int nmax;
  double cutoff, cutsq, volume;
  class NeighList *list;

  double **result;
};

}    // namespace LAMMPS_NS

#endif
#endif
