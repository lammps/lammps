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

#ifndef LMP_PACKAGE_REGISTRY_H
#define LMP_PACKAGE_REGISTRY_H

#include <map>
#include <string>

namespace LAMMPS_NS {

// Process-global map of every package style keyword to the name of the package
// it belongs to (independent of whether that package is enabled in this build).
// Used by LAMMPS::match_style() / utils::check_packages_for_style() to produce
// "style X is part of package Y which is not enabled" diagnostics, and by the
// hybrid styles to recognize known style names.  This data is fixed at build
// time (no factory functions, no plugins), so a single instance is populated
// once per process by the generated register_package_styles().

struct PackageStyleLists {
  std::map<std::string, std::string> angle_styles;
  std::map<std::string, std::string> atom_styles;
  std::map<std::string, std::string> body_styles;
  std::map<std::string, std::string> bond_styles;
  std::map<std::string, std::string> command_styles;
  std::map<std::string, std::string> compute_styles;
  std::map<std::string, std::string> dihedral_styles;
  std::map<std::string, std::string> dump_styles;
  std::map<std::string, std::string> fix_styles;
  std::map<std::string, std::string> improper_styles;
  std::map<std::string, std::string> integrate_styles;
  std::map<std::string, std::string> kspace_styles;
  std::map<std::string, std::string> minimize_styles;
  std::map<std::string, std::string> pair_styles;
  std::map<std::string, std::string> reader_styles;
  std::map<std::string, std::string> region_styles;
};

// the single process-global instance (Meyers singleton)
inline PackageStyleLists &package_styles()
{
  static PackageStyleLists lists;
  return lists;
}

// fill package_styles() from the generated style markers; called once per
// process from register_builtin_styles()
void register_package_styles();

}    // namespace LAMMPS_NS

#endif
