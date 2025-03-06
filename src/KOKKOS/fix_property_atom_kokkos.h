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

#ifdef FIX_CLASS
// clang-format off
FixStyle(property/atom/kk,FixPropertyAtomKokkos);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_PROPERTY_ATOM_KOKKOS_H
#define LMP_FIX_PROPERTY_ATOM_KOKKOS_H

#include "fix_property_atom.h"
#include "atom_vec_kokkos.h"

namespace LAMMPS_NS {

class FixPropertyAtomKokkos : public FixPropertyAtom {
 public:
  FixPropertyAtomKokkos(class LAMMPS *, int, char **);
  void post_constructor() override;
  ~FixPropertyAtomKokkos() override;
  void grow_arrays(int) override;

  void sync(ExecutionSpace space, unsigned int mask);
  void modified(ExecutionSpace space, unsigned int mask);
  void sync_overlapping_device(ExecutionSpace space, unsigned int mask);

 private:
  int dvector_flag;
};

}

#endif
#endif

