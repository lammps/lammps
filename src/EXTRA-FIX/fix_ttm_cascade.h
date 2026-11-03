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
FixStyle(ttm/cascade,FixTTMCascade);
// clang-format on
#else

#ifndef LMP_FIX_TTM_CASCADE_H
#define LMP_FIX_TTM_CASCADE_H

#include "fix_ttm_grid.h"

namespace LAMMPS_NS {

class FixTTMCascade : public FixTTMGrid {
 public:
  FixTTMCascade(class LAMMPS *, int, char **);
  ~FixTTMCascade() override;
  void post_force(int) override;
  void end_of_step() override;
  double compute_vector(int) override;

 protected:
  bool cutoff_active, offset_active, cetable_active, ketable_active;
  double time_offset;
  double ***thermal_conductivity_grid;

  // tabular specific heat data
  std::vector<double> temp_ce_values;
  std::vector<double> ce_values;
  std::vector<double> dtemp_ce_values;
  std::vector<double> dce_values;
  std::vector<double> ce_integral_values;

  // tabular thermal conductivity data
  std::vector<double> temp_ke_values;
  std::vector<double> ke_values;
  std::vector<double> dtemp_ke_values;
  std::vector<double> dke_values;

  void allocate_grid() override;
  void deallocate_grid() override;

  void tableinterpreader(const std::string &filename, const std::string &keyword);
  double linearinterpolation(double temp, const std::string &keyword);
  double integrated_ce(double temp);

  // heat flux gradient
  double heat_flux_gradient(int ix, int iy, int iz, double dxinv, double dyinv, double dzinv);
};

}    // namespace LAMMPS_NS

#endif
#endif
