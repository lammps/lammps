/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(nve/sphere/bpm,FixNVESphereBPM)
// clang-format on
#else

#ifndef LMP_FIX_NVE_SPHERE_BPM_H
#define LMP_FIX_NVE_SPHERE_BPM_H

#include "fix_nve.h"

namespace LAMMPS_NS {

class FixNVESphereBPM : public FixNVE {
 public:
  FixNVESphereBPM(class LAMMPS *, int, char **);

  void init() override;
  void initial_integrate(int) override;
  void final_integrate() override;

 protected:
  double inertia, inv_inertia;
  int extra;
  int dlm;
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix nve/sphere/bpm disc requires 2d simulation

UNDOCUMENTED

E: Fix nve/sphere/bpm requires atom style sphere/bpm

Self-explanatory.

E: Fix nve/sphere/bpm update dipole requires atom attribute mu

An atom style with this attribute is needed.

E: Fix nve/sphere/bpm requires extended particles

This fix can only be used for particles of a finite size.


*/
