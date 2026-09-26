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
PairStyle(body/rounded/polyhedron,PairBodyRoundedPolyhedron);
// clang-format on
#else

#ifndef LMP_PAIR_BODY_ROUNDED_POLYHEDRON_H
#define LMP_PAIR_BODY_ROUNDED_POLYHEDRON_H

#include "pair.h"

#include <vector>

namespace LAMMPS_NS {

class PairBodyRoundedPolyhedron : public Pair {
 public:
  PairBodyRoundedPolyhedron(class LAMMPS *);
  ~PairBodyRoundedPolyhedron() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  double memory_usage() override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  void reset_dt() override;

  virtual void kernel_force(double R, int itype, int jtype, double &energy, double &fpair);
  // elastic and cohesive parts of the normal force, returns the energy
  double normal_force(double R, int itype, int jtype, double &fe, double &fc);

  struct Contact {
    int ibody, jbody;     // body (i.e. atom) indices (not tags)
    int type;             // 0 = VERTEX-FACE; 1 = EDGE-EDGE
    double fx, fy, fz;    // unscaled cohesive forces at contact
    double xi[3];         // coordinates of the contact point on ibody
    double xj[3];         // coordinates of the contact point on jbody
    double separation;    // contact surface separation
    double r;             // distance used to normalize xi - xj into the force direction,
                          // negative when jbody has crossed a face or an edge of ibody
    int unique;
  };

  // scratch space for the interaction of a pair of bodies, one per thread

  struct Scratch {
    std::vector<Contact> contacts;    // contacts between the two bodies
    std::vector<int> vertex_done;     // flags for the vertices already interacted with
    std::vector<double> reach;        // extent of the vertices towards the other body
    int ibody, jbody;                 // the two bodies of the pair
    double reach_min_i, reach_min_j;  // minimum reach of an interacting feature
    double *shear;                    // tangential displacement of the pair, or nullptr
    int shear_i;                      // body the tangential displacement refers to
    int touched;                      // 1 if the tangential displacement was updated
  };

 protected:
  double **k_n;        // normal repulsion strength
  double **k_na;       // normal attraction strength
  double c_n;          // normal damping coefficient
  double c_t;          // tangential damping coefficient
  double mu;           // normal friction coefficient during gross sliding
  double A_ua;         // characteristic contact area
  double cut_inner;    // cutoff for interaction between vertex-edge surfaces
  double **k_t;        // tangential stiffness of the contact history
  int history;         // 1 if the tangential displacement of the contacts is stored
  double dt;           // time step, for the update of the tangential displacement
  char *id_history;    // ID of fix NEIGH_HISTORY with the tangential displacements
  class FixNeighHistory *fix_history;

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
  void sphere_against_sphere(int ibody, int jbody, int itype, int jtype, double delx, double dely,
                             double delz, double rsq, double **x, double **v, double **angmom,
                             double **f, double **torque, double **fnc, Scratch &s,
                             double &evdwl, double *facc);
  // sphere-edge interaction
  void sphere_against_edge(int ibody, int jbody, int itype, int jtype, double **x, double **v,
                           double **f, double **torque, double **angmom, double **fnc,
                           Scratch &s, double &evdwl, double *facc);
  // sphere-face interaction
  void sphere_against_face(int ibody, int jbody, int itype, int jtype, double **x, double **v,
                           double **f, double **torque, double **angmom, double **fnc,
                           Scratch &s, double &evdwl, double *facc);
  // whether two edges, or the edges at two vertices, interact as edges
  int vertex_near(const Scratch &s, int ibody, int ni) const;
  int edge_near(const Scratch &s, int ibody, int ne) const;
  int face_near(const Scratch &s, int ibody, int nf) const;
  int edges_interact(int ibody, int ei, int jbody, int ej);
  int vertex_edges_interact(int ibody, int ni, int jbody, int ej, int nj);
  // vertex-edge and vertex-vertex interactions
  void vertex_against_edge(int ibody, int jbody, int itype, int jtype, double **x, double **v,
                           double **f, double **torque, double **angmom, double **fnc,
                           Scratch &s, double &evdwl, double *facc);
  void vertex_against_vertex(int ibody, int jbody, int itype, int jtype, double **x, double **v,
                             double **f, double **torque, double **angmom, double **fnc,
                             Scratch &s, double &evdwl, double *facc);
  // edge-edge interactions
  int edge_against_edge(int ibody, int jbody, int itype, int jtype, double **x, double **v,
                        double **f, double **torque, double **angmom, double **fnc, Scratch &s,
                        double &evdwl, double *facc);
  // edge-face interactions
  int edge_against_face(int ibody, int jbody, int itype, int jtype, double **x, double **v,
                        double **f, double **torque, double **angmom, double **fnc, Scratch &s,
                        double &evdwl, double *facc);

  // a face vs. a single edge
  int interaction_face_to_edge(int ibody, int face_index, double *xmi, double rounded_radius_i,
                               int jbody, int edge_index, double *xmj, double rounded_radius_j,
                               int itype, int jtype, double cut_inner, double **v, double **f,
                               double **torque, double **angmom, double **fnc, Scratch &s,
                               double &energy, double *facc);
  // an edge vs. an edge from another body
  int interaction_edge_to_edge(int ibody, int edge_index_i, double *xmi, double rounded_radius_i,
                               int jbody, int edge_index_j, double *xmj, double rounded_radius_j,
                               int itype, int jtype, double cut_inner, double **v, double **f,
                               double **torque, double **angmom, double **fnc, Scratch &s,
                               double &energy, double *facc);

  // compute contact forces if contact points are detected
  void contact_forces(int ibody, int jbody, double *xi, double *xj, double delx, double dely,
                      double delz, double r, double **x, double **v, double **angmom,
                      double **f, double **torque, double **fnc, double *facc);
  // contact point between the rounded surfaces of two bodies
  void contact_point(const double *pi, const double *pj, const double *n, double rradi,
                     double rradj, double *pc);
  // damping and friction forces at a contact point
  void damping_friction(int ibody, int jbody, const double *pc, const double *n, double fne,
                        int damping, int friction, double **x, double **v, double **angmom,
                        double **f, double **torque, double **fnc, int iref, double *facc,
                        Scratch *hs = nullptr);
  // friction force at the contact with the largest overlap
  void friction_force(Contact &contact, int itype, int jtype, double **x, double **v,
                      double **angmom, double **f, double **torque, double **fnc, int iref,
                      double *facc, Scratch &s);
  // friction force from a tangential spring with the contact history
  void tangential_spring(int ibody, int jbody, const double *n, const double *vt, double fne,
                         Scratch &s, double *fs);
  // create or delete the placeholder of fix NEIGH_HISTORY
  void history_dummy_fix(int flag);

  // compute force and torque between two bodies given a pair of interacting points
  void pair_force_and_torque(int ibody, int jbody, double *pi, double *pj, double r,
                             double contact_dist, int itype, int jtype, double **x, double **v,
                             double **f, double **torque, double **angmom, double **fnc,
                             int jflag, double &energy, double *facc);

  // rescale the cohesive forces if a contact area is detected
  void rescale_cohesive_forces(double **x, double **f, double **torque, double **fnc,
                               std::vector<Contact> &contacts, int itype, int jtype, int iref,
                               double *facc);

  // compute the separation between two contacts
  double contact_separation(const Contact &c1, const Contact &c2);

  // detect the unique contact points (as there may be double counts)
  void find_unique_contacts(std::vector<Contact> &contacts);

  // accumulate torque to a body given a force at a given point
  void sum_torque(double *xm, double *x, double fx, double fy, double fz, double *torque);

  // find the intersection point (if any) between an edge and a face
  int edge_face_intersect(double *x1, double *x2, double *x3, double *a, double *b, double *hi1,
                          double *hi2, double &d1, double &d2, int &inside_a, int &inside_b);
  // find the face of a body with the largest signed distance to a point
  double nearest_face(int ibody, double *xmi, const double *q, double *n);
  // helper functions
  int opposite_sides(double *n, double *x0, double *a, double *b);
  void project_pt_plane(const double *q, const double *p, const double *n, double *q_proj,
                        double &d);
  void project_pt_plane(const double *q, const double *x1, const double *x2, const double *x3,
                        double *q_proj, double &d, int &inside);
  void project_pt_line(const double *q, const double *xi1, const double *xi2, double *h, double &d,
                       double &t);
  void inside_polygon(int ibody, int face_index, double *xmi, const double *q1, const double *q2,
                      int &inside1, int &inside2);

  void distance_bt_edges(const double *x1, const double *x2, const double *x3, const double *x4,
                         double *h1, double *h2, double &t1, double &t2, double &r);
  void total_velocity(double *p, double *xcm, double *vcm, double *angmom, double *inertia,
                      double *quat, double *vi);
  void sanity_check();
};

}    // namespace LAMMPS_NS

#endif
#endif
