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
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "LMP_Lepton.h"

#include <cctype>

/// remove whitespace and quotes from expression string
std::string LMP_Lepton::condense(const std::string & in)
{
  std::string out;
  for (const auto &c : in)
    if (!isspace(c) && (c != '"') && (c != '\'')) out.push_back(c);
  return out;
}


