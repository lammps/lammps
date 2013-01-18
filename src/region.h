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

#ifndef LMP_REGION_H
#define LMP_REGION_H

#include "pointers.h"

namespace LAMMPS_NS {

class Region : protected Pointers {
 public:
  char *id,*style;
  int interior;                     // 1 for interior, 0 for exterior
  int scaleflag;                    // 1 for lattice, 0 for box
  double xscale,yscale,zscale;      // scale factors for box/lattice units
  double extent_xlo,extent_xhi;     // bounding box on region
  double extent_ylo,extent_yhi;
  double extent_zlo,extent_zhi;
  int bboxflag;                     // 1 if bounding box is computable
  int varshape;                     // 1 if region shape changes over time

  // contact = particle near region surface

  struct Contact {
    double r;                 // distance between particle & surf, r > 0.0
    double delx,dely,delz;    // vector from surface pt to particle
  };
  Contact *contact;           // list of contacts
  int cmax;                   // max # of contacts possible with region

  Region(class LAMMPS *, int, char **);
  virtual ~Region();
  virtual void init();
  virtual int dynamic_check();

  // called by other classes to check point versus region

  int match(double, double, double);
  int surface(double, double, double, double);

  // implemented by each region, not called by other classes

  virtual int inside(double, double, double) = 0;
  virtual int surface_interior(double *, double) = 0;
  virtual int surface_exterior(double *, double) = 0;
  virtual void shape_update() {}

 protected:
  void add_contact(int, double *, double, double, double);
  void options(int, char **);

 private:
  int dynamic;        // 1 if region position/orientation changes over time
  int moveflag,rotateflag;   // 1 if position/orientation changes

  double point[3],axis[3],runit[3];
  char *xstr,*ystr,*zstr,*tstr;
  int xvar,yvar,zvar,tvar;
  double dx,dy,dz,theta;
  bigint lastshape,lastdynamic;

  void forward_transform(double &, double &, double &);
  void inverse_transform(double &, double &, double &);
  void rotate(double &, double &, double &, double);
};

}

#endif

/* ERROR/WARNING messages:

E: Variable name for region does not exist

Self-explanatory.

E: Variable for region is invalid style

Only equal-style variables can be used.

E: Variable for region is not equal style

Self-explanatory.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Region union or intersect cannot be dynamic

The sub-regions can be dynamic, but not the combined region.

E: Use of region with undefined lattice

If units = lattice (the default) for the region command, then a
lattice must first be defined via the lattice command.

E: Region cannot have 0 length rotation vector

Self-explanatory.

*/
