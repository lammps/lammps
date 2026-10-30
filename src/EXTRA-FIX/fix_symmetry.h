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

#include <string>
#include <unordered_map>
#include <vector>

namespace LAMMPS_NS {

class SymmetryGroup;

class FixSymmetry : public Fix {
 public:
  FixSymmetry(class LAMMPS *, int, char **);
  ~FixSymmetry() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void post_neighbor() override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  void min_setup(int) override;
  void min_post_neighbor() override;
  void min_post_force(int) override;
  bigint dof(int) override;

 private:
  SymmetryGroup *grp;
  double tol;
  char *symfile;

  // Symmetry orbits: sets of atoms related by the operators of the group.
  // All orbit data is identical on all MPI ranks.  The first member of an
  // orbit is its representative.  Every member is the image of the
  // representative under its canonical operator (the identity for the
  // representative).  If the representative sits on a special position,
  // the operators that leave it in place (its site symmetry) are listed in
  // "stabilizer" and the orbit has fewer members than the group has operators.
  struct OrbitInfo {
    int nmembers;                   // number of distinct atoms in the orbit
    int fixed_dim;                  // number of free coordinates of the representative
    bool special;                   // true if the representative is on a special position
    double proj[3][3];              // projector onto the free coordinates (fractional)
    std::vector<int> stabilizer;    // non-identity operators keeping the representative
  };
  std::vector<OrbitInfo> orbits;
  int n_orbits;
  int n_special;          // number of orbits on special positions
  int n_trans;            // number of translations of the whole system the group allows
  std::string summary;    // last summary printed to the screen and log file

  // atom ID -> global member index; member index -> orbit and canonical operator
  std::unordered_map<tagint, int> member_of_tag;
  std::vector<int> member_orbit;
  std::vector<int> member_op;

  // per-atom cache of the member index (or -1), rebuilt after reneighboring
  int *atom2member;
  int nmax_cache;

  // per-atom lattice translation selecting the periodic copy of an image
  // atom, valid until the next reneighboring
  std::vector<double> image_shift;
  bool shift_valid;

  // Cartesian forms of the operators and their inverses: M = h R h^-1
  std::vector<double> mcart, mcart_inv;
  double h_cart[6];

  // communication buffers: 9 doubles per orbit for the sums of forces and
  // velocities and the position of the representative
  double *buf_local, *buf_global;

  void symmetrize(bool);
  void detect_orbits();
  void build_cache();
  void update_cartesian_ops();
  void check_metric(const double H[3][3]);
  void check_fix_order();
  void print_summary();
  void warn_conflicting_fixes();
};

}    // namespace LAMMPS_NS

#endif
#endif
