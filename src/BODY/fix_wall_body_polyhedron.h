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
  virtual ~FixWallBodyPolyhedron();
  int setmask();
  void init();
  void setup(int);
  virtual void post_force(int);
  void reset_dt();

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
  double kn, c_n, c_t;
  double lo, hi, cylradius;
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

  void body2space(int);

  int edge_against_wall(int ibody, double wall_pos, int side, double *vwall, double **x, double **f,
                        double **torque, Contact *contact_list, int &num_contacts, double *facc);
  int sphere_against_wall(int i, double wall_pos, int side, double *vwall, double **x, double **v,
                          double **f, double **angmom, double **torque);

  int compute_distance_to_wall(int ibody, int edge_index, double *xmi, double rounded_radius_i,
                               double wall_pos, int side, double *vwall, int &contact);
  double contact_separation(const Contact &c1, const Contact &c2);
  void contact_forces(int ibody, double j_a, double *xi, double *xj, double delx, double dely,
                      double delz, double fx, double fy, double fz, double **x, double **v,
                      double **angmom, double **f, double **torque, double *vwall);

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

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix wall/body/polyhedron requires atom style body rounded/polyhedron

Self-explanatory.

E: Cannot use wall in periodic dimension

Self-explanatory.

E: Cannot wiggle and shear fix wall/body/polygon

Cannot specify both options at the same time.

E: Invalid wiggle direction for fix wall/body/polygon

Self-explanatory.

E: Fix wall/body/polygon is incompatible with Pair style

Must use a body pair style to define the parameters needed for
this fix.

*/
