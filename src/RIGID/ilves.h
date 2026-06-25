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

/* ----------------------------------------------------------------------
   ILVES constraint solver: abstract base class.

   Adapted from GROMACS 2021 ILVES (LGPL-2.1), src/gromacs/mdlib/ilves.{cpp,h}.
   The reusable algorithm (right-hand-side / left-hand-side assembly, position
   update, Lagrange-multiplier accumulation, bond partitioning and matrix
   weights) is preserved; the interface is adapted to LAMMPS (positions as
   double**, minimum-image via Domain, no SIMD/MPI/FEP in this single-rank
   port).  See ilves_graph.h for full attribution.
------------------------------------------------------------------------- */

#ifndef LMP_ILVES_H
#define LMP_ILVES_H

#include "ilves_compat.h"
#include "ilves_molecule.h"
#include "ilves_schur_solver.h"

#include <array>
#include <memory>
#include <utility>
#include <vector>

namespace LAMMPS_NS {

class LAMMPS;
class Domain;

namespace ILVES {

// per-dimension arrays of bond vectors (x[b]-x[a]) over all constraints
using VecReal = std::vector<real, AlignedAllocator<real>>;
using BondVecs = std::array<VecReal, DIM>;

class Ilves {
 public:
  /**
   * Build the ILVES solver for the given constraint list.
   *
   * @param lmp LAMMPS instance (used for the minimum-image convention).
   * @param nbonds Number of constraints.
   * @param catom1 catom1[k]/catom2[k] are the atom indices of constraint k.
   * @param catom2 See catom1.
   * @param cdist cdist[k] is the target length of constraint k.
   * @param invmass Per-atom inverse mass (1/m, no unit conversion), indexed by
   * atom index.  The solver keeps the pointer; it must stay alive and valid.
   * @param nthreads Number of OpenMP threads (1 for the serial port).
   * @param upper_tri True for the symmetric (Cholesky/LDLT) variant.
   */
  Ilves(LAMMPS *lmp, int nbonds, const int *catom1, const int *catom2, const real *cdist,
        const real *invmass, int nthreads, bool upper_tri);

  virtual ~Ilves() = default;

  /**
   * Iterate xprime (predicted positions) until the constraints are satisfied to
   * the relative tolerance tol or maxiter Newton iterations are reached.  x are
   * the reference positions (start of the step), used for the constant bond
   * vectors x_ab.  Returns {converged, number-of-iterations}.  The accumulated
   * Lagrange multipliers are left in current_lagr for add_constraint_forces().
   */
  virtual std::pair<bool, int> solve(double **x, double **xprime, real tol, int maxiter) = 0;

  /**
   * After a solve, add the constraint forces to f:
   *   f[atom1[k]] += current_lagr[k] * inv_dtfsq * x_ab[k]
   *   f[atom2[k]] -= current_lagr[k] * inv_dtfsq * x_ab[k]
   * with inv_dtfsq = 1 / (dt^2 * ftm2v), matching fix shake.
   */
  void add_constraint_forces(double **f, real inv_dtfsq) const;

  int num_constraints() const { return mol->bonds.num; }
  const Molecule *molecule() const { return mol.get(); }

 protected:
  LAMMPS *lmp;
  Domain *domain;

  int nthreads;

  std::unique_ptr<Molecule> mol;
  std::unique_ptr<SchurLinearSolver> schur_solver;

  // weights of the entries of the lhs (one vector per partition)
  std::vector<VecReal> part_lhs_weights;
  // current approximation of the Lagrange multipliers (one vector per partition)
  std::vector<VecReal> current_lagr;

  // x_ab[d][k]      = (x[b]-x[a])[d]      using reference positions x
  // xprime_ab[d][k] = (xprime[b]-xprime[a])[d] using predicted positions xprime
  BondVecs x_ab;
  BondVecs xprime_ab;

  // minimum-image difference out = xa - xb (closest periodic image)
  void min_image_sub(const double *xa, const double *xb, double *out) const;

  real make_rhs_scalar(double **x, double **xprime, bool compute_x_ab, real *rhs, int gstart,
                       int gend, int lstart);
  real make_rhs(int partition, double **x, double **xprime, bool compute_x_ab);

  void make_lhs_scalar(int partition, const BondVecs &xab1, const BondVecs &xab2, int lrowstart);
  void make_lhs(int partition, const BondVecs &xab1, const BondVecs &xab2);

  void update_current_lagr(int partition, bool first_time);
  void update_positions(int partition, double **xprime) const;

 private:
  bool disjoint_mol(int submol_max_size) const;
  void make_weights();
};

}    // namespace ILVES
}    // namespace LAMMPS_NS

#endif
