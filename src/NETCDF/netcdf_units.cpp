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

/* ----------------------------------------------------------------------
   Contributing author: Lars Pastewka (University of Freiburg)
------------------------------------------------------------------------- */

#if defined(LMP_HAS_NETCDF) || defined(LMP_HAS_PNETCDF)

#include "netcdf_units.h"

#include "error.h"

using namespace LAMMPS_NS;

std::string NetCDFUnits::get_unit_for(const char *unit_style, int quantity, Error *error)
{
  if (!strcmp(unit_style, "lj")) {
    if (quantity == Quantity::UNKNOWN) {
      return "";
    } else {
      return "lj";
    }
  } else if (!strcmp(unit_style, "real")) {
    switch (quantity) {
      case Quantity::UNKNOWN:
        return "";
      case Quantity::TIME:
        return "femtosecond";
      case Quantity::DISTANCE:
        return "angstrom";
      case Quantity::VELOCITY:
        return "angstrom/femtosecond";
      case Quantity::FORCE:
        return "(Kcal/mol)/angstrom)";
      case Quantity::DIPOLE_MOMENT:
        return "e * angstrom";
    }
  } else if (!strcmp(unit_style, "metal")) {
    switch (quantity) {
      case Quantity::UNKNOWN:
        return "";
      case Quantity::TIME:
        return "picosecond";
      case Quantity::DISTANCE:
        return "angstrom";
      case Quantity::VELOCITY:
        return "angstrom/picosecond";
      case Quantity::FORCE:
        return "eV/angstrom";
      case Quantity::DIPOLE_MOMENT:
        return "e * angstrom";
    }
  } else if (!strcmp(unit_style, "si")) {
    switch (quantity) {
      case Quantity::UNKNOWN:
        return "";
      case Quantity::TIME:
        return "second";
      case Quantity::DISTANCE:
        return "meter";
      case Quantity::VELOCITY:
        return "meter/second";
      case Quantity::FORCE:
        return "Newton";
      case Quantity::DIPOLE_MOMENT:
        return "Coulomb * meter";
    }
  } else if (!strcmp(unit_style, "cgs")) {
    switch (quantity) {
      case Quantity::UNKNOWN:
        return "";
      case Quantity::TIME:
        return "second";
      case Quantity::DISTANCE:
        return "centimeter";
      case Quantity::VELOCITY:
        return "centimeter/second";
      case Quantity::FORCE:
        return "dynes";
      case Quantity::DIPOLE_MOMENT:
        return "statcoul * cm";
    }
  } else if (!strcmp(unit_style, "electron")) {
    switch (quantity) {
      case Quantity::UNKNOWN:
        return "";
      case Quantity::TIME:
        return "femtoseconds";
      case Quantity::DISTANCE:
        return "Bohr";
      case Quantity::VELOCITY:
        return "Bohr/atomic time units";
      case Quantity::FORCE:
        return "Hartree/Bohr";
      case Quantity::DIPOLE_MOMENT:
        return "Debye";
    }
  } else if (!strcmp(unit_style, "micro")) {
    switch (quantity) {
      case Quantity::UNKNOWN:
        return "";
      case Quantity::TIME:
        return "microseconds";
      case Quantity::DISTANCE:
        return "micrometers";
      case Quantity::VELOCITY:
        return "micrometers/microsecond";
      case Quantity::FORCE:
        return "picogram * micrometer/microsecond^2";
      case Quantity::DIPOLE_MOMENT:
        return "picocoulomb * micrometer";
    }
  } else if (!strcmp(unit_style, "nano")) {
    switch (quantity) {
      case Quantity::UNKNOWN:
        return "";
      case Quantity::TIME:
        return "nanoseconds";
      case Quantity::DISTANCE:
        return "nanometers";
      case Quantity::VELOCITY:
        return "nanometers/nanosecond";
      case Quantity::FORCE:
        return "attogram * nanometer/nanosecond^2";
      case Quantity::DIPOLE_MOMENT:
        return "e * nanometer";
    }
  }

  error->all(FLERR, "Unsupported unit style: {}", unit_style);
  return "";
}

#endif
