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
FixStyle(rheo/oxidation,FixRHEOOxidation)
// clang-format on
#else

#ifndef LMP_FIX_RHEO_OXIDATION_H
#define LMP_FIX_RHEO_OXIDATION_H

#include "fix.h"

#include <vector>

namespace LAMMPS_NS {

class FixRHEOOxidation : public Fix {
 public:
  FixRHEOOxidation(class LAMMPS *, int, char **);
  ~FixRHEOOxidation() override;
  int setmask() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void setup_pre_force(int) override;
  void post_integrate() override;
  void pre_force(int) override;
  void post_force(int) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int *nbond;
  double rsurf, cut;

 private:
  int btype, index_nb;
  double cutsq;

  class NeighList *list;
  class ComputeRHEOSurface *compute_surface;
  class FixRHEO *fix_rheo;
};

}    // namespace LAMMPS_NS

#endif
#endif
