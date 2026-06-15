/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_KSPACE_ZERO2_H
#define LMP_KSPACE_ZERO2_H

#include "kspace.h"

namespace LAMMPS_NS {

class KSpaceZero2 : public KSpace {
 public:
  KSpaceZero2(class LAMMPS *);

  void init() override;
  void setup() override;
  void settings(int, char **) override;

  void compute(int, int) override;
};
}    // namespace LAMMPS_NS
#endif
