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

#ifdef REGION_CLASS
// clang-format off
RegionStyle(intersect,RegIntersect);
// clang-format on
#else

#ifndef LMP_REGION_INTERSECT_H
#define LMP_REGION_INTERSECT_H

#include "region.h"

namespace LAMMPS_NS {

class RegIntersect : public Region {
 public:
  RegIntersect(class LAMMPS *, int, char **);
  ~RegIntersect();
  void init();
  int inside(double, double, double);
  int surface_interior(double *, double);
  int surface_exterior(double *, double);
  void shape_update();
  void pretransform();
  void set_velocity();
  void length_restart_string(int &);
  void write_restart(FILE *);
  int restart(char *, int &);
  void reset_vel();

 private:
  char **idsub;
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Region intersect region ID does not exist

Self-explanatory.

E: Region union region ID does not exist

One or more of the region IDs specified by the region union command
does not exist.

*/
