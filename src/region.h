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

#ifndef LMP_REGION_H
#define LMP_REGION_H

#include "pointers.h"    // IWYU pragma: export

namespace LAMMPS_NS {

class Region : protected Pointers {
 public:
  char *id, *style;
  Region **reglist;
  int interior;                     // 1 for interior, 0 for exterior
  int scaleflag;                    // 1 for lattice, 0 for box
  double xscale, yscale, zscale;    // scale factors for box/lattice units
  double extent_xlo, extent_xhi;    // bounding box on region
  double extent_ylo, extent_yhi;
  double extent_zlo, extent_zhi;
  int bboxflag;                // 1 if bounding box is computable
  int varshape;                // 1 if region shape changes over time
  int dynamic;                 // 1 if position/orient changes over time
  int moveflag, rotateflag;    // 1 if position/orientation changes
  int openflag;                // 1 if any face is open
  int open_faces[6];           // flags for which faces are open

  int copymode;    // 1 if copy of original class

  // contact = particle near region surface (for soft interactions)
  // touch = particle touching region surface (for granular interactions)

  struct Contact {
    double r;                   // distance between particle & surf, r > 0.0
    double delx, dely, delz;    // vector from surface pt to particle
    double radius;              // curvature of region at contact point
    int iwall;                  // unique id of wall for storing shear history
    int varflag;                // 1 if wall can be variable-controlled
  };
  Contact *contact;    // list of contacts
  int cmax;            // max # of contacts possible with region
  int tmax;            // max # of touching contacts possible

  // motion attributes of region
  // public so can be accessed by other classes

  double dx, dy, dz, theta;    // current displacement and orientation
  double v[3];                 // translational velocity
  double rpoint[3];            // current origin of rotation axis
  double omega[3];             // angular velocity
  double rprev;                // speed of time-dependent radius, if applicable
  double xcenter[3];    // translated/rotated center of cylinder/sphere (only used if varshape)
  double prev[5];       // stores displacement (X3), angle and if
                        //  necessary, region variable size (e.g. radius)
                        //  at previous time step
  int vel_timestep;     // store timestep at which set_velocity was called
                        //   prevents multiple fix/wall/gran/region calls
  int nregion;          // For union and intersect
  int size_restart;

  Region(class LAMMPS *, int, char **);
  ~Region() override;
  virtual void init();
  int dynamic_check();

  // called by other classes to check point versus region

  void prematch();
  int match(double, double, double);
  int surface(double, double, double, double);

  virtual void set_velocity();
  void velocity_contact(double *, double *, int);
  virtual void write_restart(FILE *);
  virtual int restart(char *, int &);
  virtual void length_restart_string(int &);
  virtual void reset_vel();

  // implemented by each region, not called by other classes

  virtual int inside(double, double, double) = 0;
  virtual int surface_interior(double *, double) = 0;
  virtual int surface_exterior(double *, double) = 0;
  virtual void shape_update() {}
  virtual void pretransform();
  virtual void set_velocity_shape() {}
  virtual void velocity_contact_shape(double *, double *) {}

 protected:
  void add_contact(int, double *, double, double, double);
  void options(int, char **);
  void point_on_line_segment(double *, double *, double *, double *);
  void forward_transform(double &, double &, double &);
  double point[3], runit[3];

 private:
  char *xstr, *ystr, *zstr, *tstr;
  int xvar, yvar, zvar, tvar;
  double axis[3];

  void inverse_transform(double &, double &, double &);
  void rotate(double &, double &, double &, double);
};

}    // namespace LAMMPS_NS

#endif
