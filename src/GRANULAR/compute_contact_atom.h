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
ComputeStyle(contact/atom,ComputeContactAtom);
// clang-format on
#else

#ifndef LMP_COMPUTE_CONTACT_ATOM_H
#define LMP_COMPUTE_CONTACT_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeContactAtom : public Compute {
 public:
  ComputeContactAtom(class LAMMPS *, int, char **);
  ~ComputeContactAtom() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void compute_peratom() override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  double memory_usage() override;

 private:
  int nmax;

  char *group2;
  int jgroup, jgroupbit;

  class NeighList *list;
  double *contact;
};

}    // namespace LAMMPS_NS

#endif
#endif
