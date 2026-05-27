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

#ifdef FIX_CLASS
// clang-format off
FixStyle(symmetry,FixSymmetry);
// clang-format on
#else

#ifndef LMP_FIX_SYMMETRY_H
#define LMP_FIX_SYMMETRY_H

#include "fix.h"

#include <unordered_map>
#include <utility>
#include <vector>

namespace LAMMPS_NS {

class Respa;
class SymmetryGroup;

class FixSymmetry : public Fix {
 public:
  FixSymmetry(class LAMMPS *, int, char **);
  ~FixSymmetry() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  void end_of_step() override;

 private:
  SymmetryGroup *grp;
  double tol;
  char *symfile;

  // tag -> (orbit_id, op_id_within_orbit); op_id == 0 is the asymmetric representative
  std::unordered_map<tagint, std::pair<int, int>> tag_info;

  // scratch buffers for per-orbit MPI_Allreduce, sized 3 * n_orbits
  double *fsum_local, *fsum_global;
  double *vsum_local, *vsum_global;
  double *xasym_local, *xasym_global;
  int n_orbits;

  // Per-orbit Wyckoff constant c = sum_k (R_k - I)^T (t_k - L_k), where the
  // L_k are the lattice translations selected by the user's initial asym
  // position. Used to drive the Lagrange-multiplier projection in
  // end_of_step. Zero (and unused) for orbits with empty stabilizer.
  std::vector<double> wyckoff_c;

  void validate_box();
  void build_orbit_map();
  void wyckoff_project_locals();
};

}    // namespace LAMMPS_NS

#endif
#endif
