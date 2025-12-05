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

#ifndef LMP_JSON_FWD_H
#define LMP_JSON_FWD_H

// Forward declarations for header-only JSON class
// For use in headers

#include "nlohmann/json_fwd.hpp" // IWYU pragma: export

namespace LAMMPS_NS {
using json = ::nlohmann_lmp::basic_json<nlohmann_lmp::ordered_map>;
}
#endif
