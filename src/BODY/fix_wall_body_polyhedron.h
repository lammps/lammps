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

#ifdef FIX_CLASS
// clang-format off
FixStyle(wall/body/polyhedron,FixWallBodyPolyhedron);
// clang-format on
#else

#ifndef LMP_FIX_WALL_BODY_POLYHERON_H
#define LMP_FIX_WALL_BODY_POLYHERON_H

#include "fix.h"

namespace LAMMPS_NS {

class FixWallBodyPolyhedron : public Fix {
 public:
  FixWallBodyPolyhedron(class LAMMPS *, int, char **);
  ~FixWallBodyPolyhedron() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void post_force(int) override;
  void reset_dt() override;
  double memory_usage() override;

  int image(int *&, double **&) override;


 protected:
  int wallstyle, wiggle, axis;
  double kn, c_n, c_t;
  double lo, hi;
  double amplitude, period, omega;
  double dt;
  int time_origin;

  class AtomVecBody *avec;
  class BodyRoundedPolyhedron *bptr;

  double **discrete;    // list of all sub-particles for all bodies
  int ndiscrete;        // number of discretes in list
  int dmax;             // allocated size of discrete list
  int *dnum;            // number of discretes per line, 0 if uninit
  int *dfirst;          // index of first discrete per each line
  int nmax;             // allocated size of dnum,dfirst vectors

  double **edge;    // list of all edge for all bodies
  int nedge;        // number of edge in list
  int edmax;        // allocated size of edge list
  int *ednum;       // number of edges per line, 0 if uninit
  int *edfirst;     // index of first edge per each line
  int ednummax;     // allocated size of ednum,edfirst vectors

  double **face;    // list of all edge for all bodies
  int nface;        // number of faces in list
  int facmax;       // allocated size of face list
  int *facnum;      // number of faces per line, 0 if uninit
  int *facfirst;    // index of first face per each line
  int facnummax;    // allocated size of facnum,facfirst vectors

  double *enclosing_radius;    // enclosing radii for all bodies
  double *rounded_radius;      // rounded radii for all bodies

  // dump image data

  int numwalls;
  int *imgobjs;
  double **imgparms;

  // store particle interactions

  void body2space(int);



  void wall_force(int i, const double *xp, const double *n, double sd, const double *vwall,
                  double **x, double **v, double **angmom, double **f, double **torque);
  void sum_torque(double *xm, double *x, double fx, double fy, double fz, double *torque);
  void total_velocity(const double *p, double *xcm, double *vcm, double *angmom, double *inertia,
                      double *quat, double *vi);
  void distance(const double *x2, const double *x1, double &r);
};

}    // namespace LAMMPS_NS

#endif
#endif
