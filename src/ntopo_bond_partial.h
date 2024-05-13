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

#ifdef NTOPO_CLASS
// clang-format off
NTopoStyle(NTOPO_BOND_PARTIAL,NTopoBondPartial);
// clang-format on
#else

#ifndef LMP_TOPO_BOND_PARTIAL_H
#define LMP_TOPO_BOND_PARTIAL_H

#include "ntopo.h"

namespace LAMMPS_NS {

class NTopoBondPartial : public NTopo {
 public:
  NTopoBondPartial(class LAMMPS *);
  void build() override;
};

}    // namespace LAMMPS_NS

#endif
#endif
