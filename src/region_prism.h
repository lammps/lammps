/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef REGION_CLASS

RegionStyle(prism,RegPrism)

#else

#ifndef LMP_REGION_PRISM_H
#define LMP_REGION_PRISM_H

#include "region.h"

namespace LAMMPS_NS {

class RegPrism : public Region {
  friend class CreateBox;

 public:
  RegPrism(class LAMMPS *, int, char **);
  ~RegPrism();
  int inside(double, double, double);
  int surface_interior(double *, double);
  int surface_exterior(double *, double);

 private:
  double xlo,xhi,ylo,yhi,zlo,zhi;
  double xy,xz,yz;
  double h[3][3],hinv[3][3];
  int dimension;
  double a[3],b[3],c[3];       // edge vectors of region
  double clo[3],chi[3];        // opposite corners of prism
  double face[6][3];           // unit normals of 6 prism faces
  double corners[8][3];        // 8 corner pts of prism
  int tri[12][3];              // 3 corner pts of 12 triangles (2 per face)

  void find_nearest(double *, double &, double &, double &);
  int inside_tri(double *, double *, double *, double *, double *);
  void point_on_line_segment(double *, double *, double *, double *);
  double closest(double *, double *, double *, double);

  void subtract(double *, double *, double *);
  void cross(double *, double *, double *);
  double dotproduct(double *, double *);
  void normalize(double *);
};

}

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
