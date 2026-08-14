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
FixStyle(graphics/periodic,FixGraphicsPeriodic);
// clang-format on
#else

#ifndef LMP_FIX_GRAPHICS_PERIODIC_H
#define LMP_FIX_GRAPHICS_PERIODIC_H

#include "fix.h"

namespace LAMMPS_NS {
class Compute;
class ComputeChunkAtom;

class FixGraphicsPeriodic : public Fix {
 public:
  FixGraphicsPeriodic(class LAMMPS *, int, char **);
  ~FixGraphicsPeriodic() override;
  int setmask() override;
  void setup(int) override;
  void end_of_step() override;

  int image(int *&, double **&) override;

 protected:
  bool atomflag, bondflag;
  int pxlo, pxhi, pylo, pyhi, pzlo, pzhi;
  double radius;

  int numobjs;
  int *imgobjs;
  double **imgparms;
};
}    // namespace LAMMPS_NS
#endif
#endif
