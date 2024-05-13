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
ImproperStyle(fourier,ImproperFourier);
// clang-format on
#else

#ifndef LMP_IMPROPER_FOURIER_H
#define LMP_IMPROPER_FOURIER_H

#include "improper.h"

namespace LAMMPS_NS {

class ImproperFourier : public Improper {
 public:
  ImproperFourier(class LAMMPS *);
  ~ImproperFourier() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_data(FILE *) override;

 protected:
  double *k, *C0, *C1, *C2;
  int *all;
  void addone(const int &i1, const int &i2, const int &i3, const int &i4, const int &type,
              const int &evflag, const int &eflag, const double &vb1x, const double &vb1y,
              const double &vb1z, const double &vb2x, const double &vb2y, const double &vb2z,
              const double &vb3x, const double &vb3y, const double &vb3z);
  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif
