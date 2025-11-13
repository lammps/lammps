/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(pace,ComputePACE);
// clang-format on
#else

#ifndef LMP_COMPUTE_PACE_H
#define LMP_COMPUTE_PACE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputePACE : public Compute {
 public:
  ComputePACE(class LAMMPS *, int, char **);
  ~ComputePACE() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void compute_array() override;
  double memory_usage() override;

 private:
  int natoms, lastcol, ncoeff, ndims_force, ndims_virial;
  int bikflag, bik_rows, dgradflag, dgrad_rows;
  double cutmax, **pace, **paceall;
  class NeighList *list;

  Compute *c_pe, *c_virial;
  std::string id_virial;

  struct ACECimpl *acecimpl;
  std::vector<int> number_of_functions, type_offsets;
};

}    // namespace LAMMPS_NS

#endif
#endif
