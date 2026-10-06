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

#ifndef LMP_PIMD_PARTITION_H
#define LMP_PIMD_PARTITION_H

#include "lmptype.h"

#include <vector>

namespace LAMMPS_NS {
class LAMMPS;

namespace PIMDUtils {
  void check_atom_consistency(LAMMPS *, const char *, int);
  void collect_atom_vectors(LAMMPS *, const char *, const std::vector<tagint> &, double **,
                            std::vector<double> &);
}    // namespace PIMDUtils
}    // namespace LAMMPS_NS

#endif
