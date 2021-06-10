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
RegionStyle(block,RegBlock);
// clang-format on
#else

#ifndef LMP_REGION_BLOCK_H
#define LMP_REGION_BLOCK_H

#include "region.h"

namespace LAMMPS_NS {

class RegBlock : public Region {
  friend class FixPour;

 public:
  RegBlock(class LAMMPS *, int, char **);
  ~RegBlock();
  int inside(double, double, double);
  int surface_interior(double *, double);
  int surface_exterior(double *, double);

 protected:
  double xlo, xhi, ylo, yhi, zlo, zhi;
  double corners[6][4][3];
  double face[6][3];

  double find_closest_point(int, double *, double &, double &, double &);
  int inside_face(double *, int);
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Cannot use region INF or EDGE when box does not exist

Regions that extend to the box boundaries can only be used after the
create_box command has been used.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
