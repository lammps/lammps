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

#ifdef INTEGRATE_CLASS
// clang-format off
IntegrateStyle(verlet/kk,VerletKokkos);
IntegrateStyle(verlet/kk/device,VerletKokkos);
IntegrateStyle(verlet/kk/host,VerletKokkos);
// clang-format on
#else

// clang-format off
#ifndef LMP_VERLET_KOKKOS_H
#define LMP_VERLET_KOKKOS_H

#include "verlet.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

class VerletKokkos : public Verlet {
 public:
  VerletKokkos(class LAMMPS *, int, char **);

  void init() override;
  void setup(int) override;
  void setup_minimal(int) override;
  void run(int) override;
  void force_clear() override;

  // A fix with a pre_force() kernel over all owned and ghost atoms (fix
  // OXDNA/LRF/kk) can take over zeroing the forces and torques, which saves
  // the two kernel launches of force_clear() on every step.  The fix requests
  // this in its init(); setup() checks whether it is safe for this run.
  void request_force_clear_by_fix(class Fix *fix) { force_clear_fix = fix; }
  int force_clear_by_fix(const class Fix *fix) const
    { return force_clear_fused && (fix == force_clear_fix); }

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    f(i,0) += f_merge_copy(i,0);
    f(i,1) += f_merge_copy(i,1);
    f(i,2) += f_merge_copy(i,2);
  }

 protected:
  DAT::t_kkacc_1d_3 f_merge_copy,f;
  int fuse_force_clear,fuse_integrate;
  class Fix *force_clear_fix;
  int force_clear_fused;

  void check_force_clear_fix();

  void fuse_check(int, int);
  int overlap_possible();
  int host_force_styles(uint64_t * = nullptr);
};
}

#endif
#endif

