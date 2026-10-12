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

#include "creator_registry.h"

#include <mutex>

namespace LAMMPS_NS {

// forward declarations of the per-category registration functions defined in
// the generated style_*.cpp files.  Each is added here as its category is
// migrated to the global registry.

void register_pair_styles();
void register_bond_styles();
void register_angle_styles();
void register_dihedral_styles();
void register_improper_styles();
void register_kspace_styles();
void register_fix_styles();
void register_compute_styles();
void register_integrate_styles();
void register_minimize_styles();
void register_region_styles();
void register_dump_styles();
void register_command_styles();
void register_atom_styles();
void register_body_styles();
void register_reader_styles();

// fills the global package-style -> package-name map (package_registry.h)
void register_package_styles();

/* ----------------------------------------------------------------------
   register all built-in styles exactly once per process.  The std::call_once
   guard makes this safe to call from every LAMMPS constructor (including
   concurrent instances) while running the registration only the first time.
------------------------------------------------------------------------- */

void register_builtin_styles()
{
  static std::once_flag flag;
  std::call_once(flag, []() {
    // per-category register_*_styles() calls are added here during migration
    register_pair_styles();
    register_bond_styles();
    register_angle_styles();
    register_dihedral_styles();
    register_improper_styles();
    register_kspace_styles();
    register_fix_styles();
    register_compute_styles();
    register_integrate_styles();
    register_minimize_styles();
    register_region_styles();
    register_dump_styles();
    register_command_styles();
    register_atom_styles();
    register_body_styles();
    register_reader_styles();

    // style-keyword -> package-name map used for diagnostics
    register_package_styles();
  });
}

}    // namespace LAMMPS_NS
