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
FixStyle(colvars/kk,FixColvarsKokkos);
FixStyle(colvars/kk/device,FixColvarsKokkos);
FixStyle(colvars/kk/host,FixColvarsKokkos);
// clang-format on
#else

#ifndef LMP_FIX_COLVARS_KOKKOS_H
#define LMP_FIX_COLVARS_KOKKOS_H

#include "fix_colvars.h"

namespace LAMMPS_NS {

class FixColvarsKokkos : public FixColvars {

 public:
  FixColvarsKokkos(class LAMMPS *, int, char **);

  void post_force(int) override;
  void end_of_step() override;

};

}    // namespace LAMMPS_NS

#endif
#endif
