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

#ifndef LMP_RANMARS_H
#define LMP_RANMARS_H

#include "pointers.h"

namespace LAMMPS_NS {

class RanMars : protected Pointers {
 public:
  // number of values stored by get_state(): a layout marker, u[1..97], i97, j97,
  //   c, cd, cm, and the flag and value of the cached second gaussian number
  static constexpr int STATE_SIZE = 1 + 97 + 2 + 3 + 2;

  RanMars(class LAMMPS *, int);
  ~RanMars() override;

  // uniform random number in (0,1), defined here so that it can be inlined.
  // the two wrap-arounds are written without branches: they are taken at random
  // and would be mispredicted about half the time.  adding 0.0 when they are
  // not taken leaves the values unchanged, so the sequence is the same.
  double uniform()
  {
    double uni = u[i97] - u[j97];
    uni += static_cast<double>(uni < 0.0);
    u[i97] = uni;
    i97--;
    if (i97 == 0) i97 = 97;
    j97--;
    if (j97 == 0) j97 = 97;
    c -= cd;
    c += cm * static_cast<double>(c < 0.0);
    uni -= c;
    uni += static_cast<double>(uni < 0.0);
    return uni;
  }

  double gaussian();
  double gaussian(double mu, double sigma);
  double rayleigh(double sigma);
  double besselexp(double theta, double alpha, double cp);
  void select_subset(bigint, int, int *, int *);
  void get_state(double *);
  void set_state(double *);
  static int state_size(const double *);

 private:
  // number of values in states stored before the cached gaussian was included
  static constexpr int LEGACY_STATE_SIZE = 98 + 2 + 3;

  int save;
  double second;
  double *u;
  int i97, j97;
  double c, cd, cm;
};

}    // namespace LAMMPS_NS

#endif
