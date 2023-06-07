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

/*------------------------------------------------------------------------
  Contributing Authors : Romain Vermorel (LFCR), Laurent Joly (ULyon)
  --------------------------------------------------------------------------*/

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(stress/mop/profile,ComputeStressMopProfile);
// clang-format on
#else

#ifndef LMP_COMPUTE_STRESS_MOP_PROFILE_H
#define LMP_COMPUTE_STRESS_MOP_PROFILE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeStressMopProfile : public Compute {
 public:
  ComputeStressMopProfile(class LAMMPS *, int, char **);
  ~ComputeStressMopProfile() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void compute_array() override;

 private:
  void compute_pairs();
  void compute_bonds();
  void setup_bins();

  int nvalues, dir;
  int *which;

  int bondflag;

  int originflag;
  double origin, delta, offset, invdelta;
  int nbins;
  double **coord, **coordp;
  double **values_local, **values_global;
  double **bond_local, **bond_global;
  double **local_contribution;

  double dt, nktv2p, ftm2v;
  double area;
  class NeighList *list;
};

}    // namespace LAMMPS_NS

#endif
#endif
