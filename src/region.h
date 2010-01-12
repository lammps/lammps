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
  int dynamic;                      // 1 if region changes over time

  // contact = particle near region surface

  struct Contact {
    double r;                 // distance between particle & surf, r > 0.0
    double delx,dely,delz;    // vector from surface pt to particle
  };
  Contact *contact;           // list of contacts
  int cmax;                   // max # of contacts possible with region
 
  Region(class LAMMPS *, int, char **);
  virtual ~Region();
  void init();
  virtual int dynamic_check();
  int match(double, double, double);
  int surface(double, double, double, double);

  virtual int inside(double, double, double) = 0;
  virtual int surface_interior(double *, double) = 0;
  virtual int surface_exterior(double *, double) = 0;

 protected:
  void options(int, char **);
  void add_contact(int, double *, double, double, double);

 private:
  int time_origin;
  double dt,period,omega_rotate;
  double vx,vy,vz;
  double ax,ay,az;
  double point[3],axis[3],runit[3];

  void rotate(double &, double &, double &, double);
};

}

#endif
