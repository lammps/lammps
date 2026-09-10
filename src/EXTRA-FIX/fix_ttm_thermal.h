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
FixStyle(ttm/thermal,FixTTMThermal);
// clang-format on
#else

#ifndef LMP_FIX_TTM_THERMAL_H
#define LMP_FIX_TTM_THERMAL_H

#include "fix_ttm.h"

namespace LAMMPS_NS {

class FixTTMThermal : public FixTTM {
 public:
  FixTTMThermal(class LAMMPS *, int, char **);
  ~FixTTMThermal() override;

  void post_constructor() override;
  void init() override;

  void post_force(int) override;
  void end_of_step() override;

  double compute_vector(int) override;
  double memory_usage() override;

 protected:
  double inductive_power;

  std::string e_property_file;

  double ***gamma_p_grid;
  double ***inductive_response_grid;
  double ***c_e_grid;
  double ***k_e_grid;

  void allocate_grid() override;
  void deallocate_grid() override;
  virtual void read_electron_properties(const std::string &);
};

}    // namespace LAMMPS_NS

#endif
#endif
