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
FixStyle(wall/body/polygon,FixWallBodyPolygon);
// clang-format on
#else

#ifndef LMP_FIX_WALL_BODY_POLYGON_H
#define LMP_FIX_WALL_BODY_POLYGON_H

#include "fix.h"

namespace LAMMPS_NS {

class FixWallBodyPolygon : public Fix {
 public:
  FixWallBodyPolygon(class LAMMPS *, int, char **);
  ~FixWallBodyPolygon() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void post_force(int) override;
  void reset_dt() override;

  struct Contact {
    int ibody, jbody;     // body (i.e. atom) indices (not tags)
    int vertex;           // vertex of the first polygon
    int edge;             // edge of the second polygon
    double xv[3];         // coordinates of the vertex
    double xe[3];         // coordinates of the projection of the vertex on the edge
    double separation;    // separation at contact
  };

 protected:
  int wallstyle, pairstyle, wiggle, axis;
  double kn;     // normal repulsion strength
  double c_n;    // normal damping coefficient
  double c_t;    // tangential damping coefficient
  double lo, hi, cylradius;
  double amplitude, period, omega;
  double dt;
  int time_origin;

  class AtomVecBody *avec;
  class BodyRoundedPolygon *bptr;

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

  double *enclosing_radius;    // enclosing radii for all bodies
  double *rounded_radius;      // rounded radii for all bodies

  void body2space(int);

  int vertex_against_wall(int ibody, double wall_pos, double **x, double **f, double **torque,
                          int side, Contact *contact_list, int &num_contacts, double *facc);

  int compute_distance_to_wall(double *x0, double rradi, double wall_pos, int side, double &d,
                               double hi[3], int &contact);
  double contact_separation(const Contact &c1, const Contact &c2);
  void contact_forces(Contact &contact, double j_a, double **x, double **v, double **angmom,
                      double **f, double **torque, double *vwall, double *facc);
  void sum_torque(double *xm, double *x, double fx, double fy, double fz, double *torque);
  void total_velocity(double *p, double *xcm, double *vcm, double *angmom, double *inertia,
                      double *quat, double *vi);
  void distance(const double *x2, const double *x1, double &r);
};

}    // namespace LAMMPS_NS

#endif
#endif
