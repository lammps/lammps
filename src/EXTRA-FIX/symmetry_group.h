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

namespace LAMMPS_NS {

// One affine operator in fractional coordinates: s' = R*s + t.
struct SymmOp {
  double R[3][3];
  double t[3];
  double Rinv[3][3];
};

// One orbit as listed in the symmetry file: one atom ID per operator, the
// image of the first atom under that operator.  For an atom on a special
// position, the IDs repeat.
struct Orbit {
  std::vector<tagint> image_tag;    // length n_ops
};

class SymmetryGroup : protected Pointers {
 public:
  SymmetryGroup(class LAMMPS *);
  ~SymmetryGroup() override = default;

  void read(const char *filename);

  [[nodiscard]] int n_ops() const { return ops.size(); }
  [[nodiscard]] int n_orbits() const { return orbits.size(); }
  [[nodiscard]] const std::string &name() const { return group_name; }

  [[nodiscard]] const SymmOp &op(int i) const { return ops[i]; }
  [[nodiscard]] const Orbit &orbit(int i) const { return orbits[i]; }

 private:
  std::string group_name;
  std::vector<SymmOp> ops;
  std::vector<Orbit> orbits;

  void populate(const json &data);
  void compute_inverses();
  void validate_group();
};

}    // namespace LAMMPS_NS

#endif
