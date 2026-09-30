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

// Lattice families recognized in the symmetry-data file. Used by FixSymmetry
// to validate the simulation box at init.
enum LatticeFamily {
  LATTICE_TRICLINIC,
  LATTICE_MONOCLINIC,
  LATTICE_ORTHORHOMBIC,
  LATTICE_TETRAGONAL,
  LATTICE_HEXAGONAL,
  LATTICE_TRIGONAL,
  LATTICE_CUBIC
};

// One affine operator in fractional coordinates: s' = R*s + t.
struct SymmOp {
  double R[3][3];
  double t[3];
  double Rinv[3][3];
};

// One orbit: the asymmetric representative and an entry per non-identity
// operator giving the tag of the image atom. For an atom on a special
// (Wyckoff) position, the orbit is shorter than the group order; the
// corresponding image_tag slots simply repeat asym_tag, and the op
// indices that map asym onto itself are listed in `stabilizer`.
struct Orbit {
  tagint asym_tag;
  std::vector<tagint> image_tag;    // length n_ops; image_tag[0] == asym_tag

  // Site-symmetry: list of op indices k > 0 such that ops[k] fixes the
  // asym atom (i.e. R_k * s_asym + t_k == s_asym, mod lattice). Empty for
  // general positions. Precomputed projector data B and Binv are valid
  // only when this is non-empty.
  std::vector<int> stabilizer;
  double B[3][3];       // sum_k (R_k - I)^T (R_k - I) for k in stabilizer
  double Binv[3][3];    // Moore-Penrose pseudo-inverse of B via jacobi3 SVD trick
};

class SymmetryGroup : protected Pointers {
 public:
  SymmetryGroup(class LAMMPS *);
  ~SymmetryGroup() override = default;

  void read(const char *filename);

  [[nodiscard]] int n_ops() const { return ops.size(); }
  [[nodiscard]] int n_orbits() const { return orbits.size(); }
  [[nodiscard]] LatticeFamily lattice() const { return lattice_family; }
  [[nodiscard]] const std::string &name() const { return group_name; }

  [[nodiscard]] const SymmOp &op(int i) const { return ops[i]; }
  [[nodiscard]] const Orbit &orbit(int i) const { return orbits[i]; }

 private:
  std::string group_name;
  LatticeFamily lattice_family;
  std::vector<SymmOp> ops;
  std::vector<Orbit> orbits;

  void populate(const json &data);
  void compute_inverses();
  void compute_wyckoff_projectors();
  void validate_group();
};

}    // namespace LAMMPS_NS

#endif
