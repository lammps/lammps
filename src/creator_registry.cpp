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

    // style-keyword -> package-name map used for diagnostics
    register_package_styles();
  });
}

}    // namespace LAMMPS_NS
