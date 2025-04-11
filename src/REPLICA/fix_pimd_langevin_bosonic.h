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
FixStyle(pimd/langevin/bosonic, FixPIMDBLangevin);
// clang-format on
#else

#ifndef FIX_PIMDB_LANGEVIN_H
#define FIX_PIMDB_LANGEVIN_H

#include "fix_pimd_langevin.h"

namespace LAMMPS_NS {

class FixPIMDBLangevin : public FixPIMDLangevin {
 public:
  FixPIMDBLangevin(class LAMMPS *, int, char **);
  ~FixPIMDBLangevin();

  double compute_vector(int) override;
  void compute_spring_energy() override;
  void compute_t_prim() override;

  char **filtered_args;
  int filtered_narg;

 protected:
  void prepare_coordinates() override;
  void spring_force() override;

 private:
  const int nbosons;
  bool synch_energies;
  class BosonicExchange *bosonic_exchange;
  double **f_tag_order;
  char **filter_args(
      int, char **);    // for hold memory of filtered arguments when calling the parent constructor
};

}    // namespace LAMMPS_NS

#endif
#endif
