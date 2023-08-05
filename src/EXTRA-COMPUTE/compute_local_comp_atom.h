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
ComputeStyle(local/comp/atom,ComputeLocalCompAtom);
// clang-format on
#else

#ifndef LMP_COMPUTE_LOCAL_COMP_ATOM_H
#define LMP_COMPUTE_LOCAL_COMP_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeLocalCompAtom : public Compute {
 public:
  ComputeLocalCompAtom(class LAMMPS *, int, char **);
  ~ComputeLocalCompAtom() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void compute_peratom() override;
  double memory_usage() override;

 protected:
  int nmax, ntypes;
  double cutoff;            // global cutoff distance
  double cutsq;             // cutoff**2
  class NeighList *list;    // neighbor list
  double **result;          // peratom array of local compositions
};

}    // namespace LAMMPS_NS

#endif
#endif
