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
FixStyle(align/neighbor,FixAlignNeighbor);
// clang-format on
#else

#ifndef LMP_FIX_ALIGN_NEIGHBOR_H
#define LMP_FIX_ALIGN_NEIGHBOR_H

#include "fix.h"

namespace LAMMPS_NS {

class FixAlignNeighbor : public Fix {
 public:
  FixAlignNeighbor(class LAMMPS *, int, char **);
  int setmask() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void setup(int) override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;

 private:
  int mode;                // per-atom property that stores the orientation
  int symmetry;            // polar or nematic alignment
  double magnitude;        // prefactor of the alignment torque
  double cutoff, cutsq;    // alignment cutoff and its square
  int ilevel_respa;        // rRESPA level at which the torque is applied
  class NeighList *list;

  void post_force_dipole();
};

}    // namespace LAMMPS_NS
#endif
#endif
