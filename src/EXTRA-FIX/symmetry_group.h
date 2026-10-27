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

#ifndef LMP_SYMMETRY_GROUP_H
#define LMP_SYMMETRY_GROUP_H

#include "json_fwd.h"
#include "pointers.h"

#include <string>
#include <vector>

namespace LAMMPS_NS {

// One affine operator in fractional coordinates of the simulation box: s' = R*s + t.
struct SymmOp {
  double R[3][3];
  double t[3];
  double Rinv[3][3];
};

// The operators of a space group, read from a JSON format symmetry file.
// The identity is always the first operator.
class SymmetryGroup : protected Pointers {
 public:
  SymmetryGroup(class LAMMPS *);

  void read(const char *filename);

  [[nodiscard]] int n_ops() const { return ops.size(); }
  [[nodiscard]] const std::string &name() const { return group_name; }
  [[nodiscard]] const SymmOp &op(int i) const { return ops[i]; }

  // convert an operator in coordinate triplet notation, e.g. "-y,x-y,z+1/3",
  // throws std::invalid_argument for malformed input
  static SymmOp parse_xyz(const std::string &);

 private:
  std::string group_name;
  std::vector<SymmOp> ops;
  int ncell[3];    // unit cells in the box from the "supercell" entry

  void populate(const json &data);
  void apply_supercell();
  void compute_inverses();
  void validate_group();
};

}    // namespace LAMMPS_NS

#endif
