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
IntegrateStyle(verlet,Verlet);
// clang-format on
#else

#ifndef LMP_VERLET_H
#define LMP_VERLET_H

#include "integrate.h"

namespace LAMMPS_NS {

class Verlet : public Integrate {
 public:
  Verlet(class LAMMPS *, int, char **);
  void init() override;
  void setup(int flag) override;
  void setup_minimal(int) override;
  void run(int) override;
  void force_clear() override;
  void cleanup() override;

 protected:
  int triclinic;    // 0 if domain is orthog, 1 if triclinic
  int torqueflag, extraflag;
};

}    // namespace LAMMPS_NS

#endif
#endif
