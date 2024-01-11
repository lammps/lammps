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

#ifdef IMPROPER_CLASS
// clang-format off
ImproperStyle(ring,ImproperRing);
// clang-format on
#else

#ifndef LMP_IMPROPER_RING_H
#define LMP_IMPROPER_RING_H

#include "improper.h"

namespace LAMMPS_NS {

class ImproperRing : public Improper {
 public:
  ImproperRing(class LAMMPS *);
  ~ImproperRing() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_data(FILE *) override;

 protected:
  double *k, *chi;

  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif
