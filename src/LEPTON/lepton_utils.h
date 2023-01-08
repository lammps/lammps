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

#include <string>

// forward declarations

namespace LAMMPS_NS {
class LAMMPS;
}

// utility functions and classes

namespace LeptonUtils {

/// remove whitespace and quotes from expression string
std::string condense(const std::string &);
/// substitute LAMMPS variable references with their value
std::string substitute(const std::string &, LAMMPS_NS::LAMMPS *);

}    // namespace LeptonUtils
