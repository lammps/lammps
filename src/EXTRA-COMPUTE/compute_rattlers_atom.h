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
ComputeStyle(rattlers/atom,ComputeRattlersAtom);
// clang-format on
#else

#ifndef LMP_COMPUTE_RATTLERS_ATOM_H
#define LMP_COMPUTE_RATTLERS_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeRattlersAtom : public Compute {
 public:
  ComputeRattlersAtom(class LAMMPS *, int, char **);
  ~ComputeRattlersAtom() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void compute_peratom() override;
  double compute_scalar() override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;

 private:
  int pstyle, cutstyle;
  int ncontacts_rattler, max_tries, nmax, invoked_peratom;
  int *ncontacts;
  double *rattler;
  class NeighList *list;

};

}    // namespace LAMMPS_NS

#endif
#endif
