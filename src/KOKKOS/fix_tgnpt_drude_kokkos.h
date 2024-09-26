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
FixStyle(tgnpt/drude/kk,FixTGNPTDrudeKokkos<LMPDeviceType>);
FixStyle(tgnpt/drude/kk/device,FixTGNPTDrudeKokkos<LMPDeviceType>);
FixStyle(tgnpt/drude/kk/host,FixTGNPTDrudeKokkos<LMPHostType>);
// clang-format on
#else

#ifndef LMP_FIX_TGNPT_DRUDE_KOKKOS_H
#define LMP_FIX_TGNPT_DRUDE_KOKKOS_H

#include "fix_tgnh_drude_kokkos.h"

namespace LAMMPS_NS {

class FixTGNPTDrudeKokkos : public FixTGNHDrudeKokkos {
 public:
  FixTGNPTDrudeKokkos(class LAMMPS *, int, char **);
};

}    // namespace LAMMPS_NS

#endif
#endif
