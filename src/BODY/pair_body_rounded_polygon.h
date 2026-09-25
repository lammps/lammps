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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(body/rounded/polygon,PairBodyRoundedPolygon);
// clang-format on
#else

#ifndef LMP_PAIR_BODY_ROUNDED_POLYGON_H
#define LMP_PAIR_BODY_ROUNDED_POLYGON_H

#include "pair.h"

#include <vector>

namespace LAMMPS_NS {

class PairBodyRoundedPolygon : public Pair {
 public:
  PairBodyRoundedPolygon(class LAMMPS *);
  ~PairBodyRoundedPolygon() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  double memory_usage() override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;

  struct Contact {
    int ibody, jbody;     // body (i.e. atom) indices (not tags)
    int vertex;           // vertex of the first polygon
    int edge;             // edge of the second polygon
    double xv[3];         // coordinates of the vertex
    double xe[3];         // coordinates of the projection of the vertex on the edge
    double separation;    // separation at contact
    double fe[3];         // elastic force on the vertex, the edge gets -fe
    double fc[3];         // unscaled cohesive force on the vertex, the edge gets -fc
  };

  // scratch space for the interaction of a pair of bodies, one per thread

  struct Scratch {
    std::vector<Contact> contacts;    // vertex-edge contacts
    std::vector<int> vertex_done;     // flags for the vertices already interacted with
  };

 protected:
  double **k_n;        // normal repulsion strength
  double **k_na;       // normal attraction strength
  double c_n;          // normal damping coefficient
  double c_t;          // tangential damping coefficient
  double mu;           // normal friction coefficient during gross sliding
  double delta_ua;     // contact line (area for 3D models) modification factor
  double cut_inner;    // cutoff for interaction between vertex-edge surfaces

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
  double *maxrad;              // per-type maximum radius (enclosing + rounded)
  double w_ja;                 // work done by the force added by the j_a scaling
  double w_diss;               // work done by damping and friction
  double **fnc;                // per-body force and torque not deriving from the energy
  int nmax_fnc;                // allocated size of fnc
  char *id_fix_store;          // ID of fix STORE/ATOM with fnc of the previous step
  class FixStoreAtom *fix_store;

  void work_nonconservative();

  Scratch scratch;    // scratch space of the serial compute()

  void allocate();
  void body2space(int);

  // interaction between two bodies
  void pair_interaction(int i, int j, double delx, double dely, double delz, double rsq,
                        double **x, double **v, double **angmom, double **f, double **torque,
                        double **fnc, Scratch &s, double &evdwl, double *facc);
  // sphere-sphere interaction
  void sphere_against_sphere(int i, int j, double delx, double dely, double delz, double rsq,
                             double k_n, double k_na, double **v, double **f, double **fnc,
                             double &evdwl, double *facc);
  // vertex-edge interaction
  int vertex_against_edge(int i, int j, double k_n, double k_na, double **x, double **f,
                          double **torque, tagint *tag, Scratch &s, double &evdwl,
                          double *facc);
  // compute distance between a point and an edge from another body
  int compute_distance_to_vertex(int ibody, int edge_index, double *xmi, double rounded_radius,
                                 double *x0, double x0_rounded_radius, double cut_inner, double &d,
                                 double hi[3], double &t, int &contact);
  // compute contact forces if contact points are detected
  void contact_forces(Contact &contact, double j_a, int friction, double **x, double **v,
                      double **angmom, double **f, double **torque, double **fnc, double &evdwl,
                      double *facc);
  // normal force and energy at a given surface separation
  double normal_force(double R, double k_n, double k_na, double &fe, double &fc);

  // compute the separation between two contacts
  double contact_separation(const Contact &c1, const Contact &c2);

  // accumulate torque to a body given a force at a given point
  void sum_torque(double *xm, double *x, double fx, double fy, double fz, double *torque);
  // helper functions
  int opposite_sides(double *x1, double *x2, double *a, double *b);
  void total_velocity(double *p, double *xcm, double *vcm, double *angmom, double *inertia,
                      double *quat, double *vi);
  inline void distance(const double *x2, const double *x1, double &r);
};

}    // namespace LAMMPS_NS

#endif
#endif
