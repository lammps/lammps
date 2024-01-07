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

#ifndef LMP_ELEMENTS_H
#define LMP_ELEMENTS_H

/*! \file elements.h */

#include "lmptype.h"
#include "error.h"

#include <string>

namespace LAMMPS_NS {
namespace elements {

  std::string symbol(int atomic_number, Error *error);
  std::string name(int atomic_number, Error *error);
  std::string cpkHexColor(int atomic_number, Error *error);
  
  double atomic_mass(int atomic_number, Error *error); // units u
  double vdw_radius(int atomic_number, Error *error); // units pm
  double covalent_radius(int atomic_number, Error *error); // units pm
  
  int atomic_number_with_symbol(const std::string &symbol, Error *error);
  int atomic_number_with_closest_mass(double mass, Error *error);

} // namespace elements
} // namespace LAMMPS_NS
#endif

