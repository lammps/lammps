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
   increment, Lagrange-multiplier accumulation, bond partitioning and matrix
   weights) is preserved; the interface is adapted to LAMMPS (positions as
   double**, no SIMD/FEP).

   The Newton iteration is global (convergence is the all-reduced maximum
   relative violation), so the fix drives the loop and the MPI communication
   uniformly across all ranks, calling these per-rank primitives:
     prepare()   -> assemble g(x) and the Jacobian + factorization; returns the
                    local max violation
     step(dx)    -> solve the linear system; accumulate per-atom position
                    increments into dx (home + ghost)
     recompute() -> accumulate the multipliers; reassemble g(x); return the
                    local max violation
   Between step() and recompute() the fix reverse-sums dx to the owning ranks,
   applies it, and forward-communicates the predicted positions to the ghosts.
   See ilves_graph.h for full attribution.
------------------------------------------------------------------------- */

#ifndef LMP_ILVES_H
#define LMP_ILVES_H

#include "ilves_molecule.h"
#include "ilves_schur_solver.h"

#include <array>
#include <memory>
#include <vector>

namespace LAMMPS_NS {

class LAMMPS;

namespace ILVES {

// dimension constants (matching the GROMACS DIM / XX / YY / ZZ spellings)
enum { XX = 0, YY = 1, ZZ = 2, DIM = 3 };

// per-dimension arrays of bond vectors (x[b]-x[a]) over all constraints
using VecDouble = std::vector<double>;
using BondVecs = std::array<VecDouble, DIM>;

class Ilves {
 public:
  Ilves(LAMMPS *lmp, int nbonds, const int *catom1, const int *catom2, const double *cdist,
        const double *invmass);

  virtual ~Ilves() = default;

  // assemble g(x) using reference positions x and predicted positions xprime;
  // returns the local max relative (squared) bond-length violation
  virtual double prepare(double **x, double **xprime) = 0;

  // one Newton step: solve the linear system, accumulate position increments
  // (for both atoms of every constraint, home or ghost) into dx
  virtual void step(double **dx) = 0;

  // accumulate the multipliers, then reassemble g(x); first_iter selects the
  // initial multiplier handling.  returns the local max relative violation
  double recompute(double **x, double **xprime, bool first_iter);

  // add the constraint contribution to the 6-component global virial:
  // sum over owned constraints of -lambda*inv_dtfsq * r (x) r
  void add_global_virial(double *v6, double inv_dtfsq) const;

  int num_constraints() const { return mol->bonds.num; }
  const Molecule *molecule() const { return mol.get(); }

 protected:
  LAMMPS *lmp;

  std::unique_ptr<Molecule> mol;
  std::unique_ptr<SchurLinearSolver> schur_solver;

  // weights of the entries of the lhs (one vector per partition)
  std::vector<VecDouble> part_lhs_weights;
  // current approximation of the Lagrange multipliers (one vector per partition)
  std::vector<VecDouble> current_lagr;

  // x_ab[d][k]      = (x[b]-x[a])[d]      using reference positions x
  // xprime_ab[d][k] = (xprime[b]-xprime[a])[d] using predicted positions xprime
  BondVecs x_ab;
  BondVecs xprime_ab;

  double make_rhs_scalar(double **x, double **xprime, bool compute_x_ab, double *rhs, int gstart,
                       int gend, int lstart);
  double make_rhs(int partition, double **x, double **xprime, bool compute_x_ab);

  void make_lhs_scalar(int partition, const BondVecs &xab1, const BondVecs &xab2, int lrowstart);
  void make_lhs(int partition, const BondVecs &xab1, const BondVecs &xab2);

  void update_current_lagr(int partition, bool first_time);
  // accumulate this iteration's position increments into dx (home + ghost)
  void accumulate_increment(int partition, double **dx) const;

 private:
  void make_weights();
};

}    // namespace ILVES
}    // namespace LAMMPS_NS

#endif
