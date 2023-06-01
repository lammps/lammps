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
ComputeStyle(local_composition/atom,ComputeLocalCompositionAtom);
// clang-format on
#else

#ifndef LMP_COMPUTE_LOCAL_COMPOSITION_ATOM_H
#define LMP_COMPUTE_LOCAL_COMPOSITION_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeLocalCompositionAtom : public Compute {
 public:
  ComputeLocalCompositionAtom(class LAMMPS *, int, char **);
  ~ComputeLocalCompositionAtom() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void compute_peratom() override;
  double memory_usage() override;

 protected:
  int nmax;
  double cutoff;            // global cutoff distance
  double cutsq;             // cutoff**2
  double volume;            // local volume
  double nelements;         // number of elements
  int *map;                 // map types to [0,nelements)
  class NeighList *list;    // neighbor list

  double **result;          // peratom array of local compositions
};

}    // namespace LAMMPS_NS

#endif
#endif
