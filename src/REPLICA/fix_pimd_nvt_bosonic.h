/* ----------------------------------------------------------------------
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
FixStyle(pimd/nvt/bosonic,FixPIMDBNVT);
// clang-format on
#else

#ifndef FIX_PIMDB_NVT_H
#define FIX_PIMDB_NVT_H

#include "fix_pimd_nvt.h"

namespace LAMMPS_NS {

class FixPIMDBNVT : public FixPIMDNVT {
 public:
  FixPIMDBNVT(class LAMMPS *, int, char **);
  ~FixPIMDBNVT();
  double compute_vector(int) override;

 protected:
  void prepare_coordinates() override;
  void spring_force() override;
  void pre_spring_force_estimators() override;

 private:
  class BosonicExchange *bosonic_exchange;
  double prim;
};

}    // namespace LAMMPS_NS

#endif
#endif
