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

/* ----------------------------------------------------------------------
   Contributing authors: Lars Pastewka (University of Freiburg), Guillaume Fraux (EPFL)
------------------------------------------------------------------------- */

#ifndef LMP_NETCDF_UNITS_H
#define LMP_NETCDF_UNITS_H

#if defined(LMP_HAS_NETCDF) || defined(LMP_HAS_PNETCDF)

#include <string>

namespace LAMMPS_NS {
class Error;

namespace NetCDFUnits {
  // type of quantity for per-atom values (used to get the unit)
  enum Quantity {
    UNKNOWN = 0,
    TIME,
    DISTANCE,
    VELOCITY,
    FORCE,
    DIPOLE_MOMENT,
  };

  // for compatibility with older NetCDF versions
  static constexpr int LMP_MAX_VAR_DIMS = 1024;

  // get the name of the unit for the given `quantity` in the given LAMMPS
  // `unit_style` any error will be reported through `error`
  std::string get_unit_for(const char *unit_style, int quantity, Error *error);
}    // namespace NetCDFUnits
}    // namespace LAMMPS_NS

#endif
#endif
