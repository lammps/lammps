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
FixStyle(surface/global,FixSurfaceGlobal)
// clang-format on
#else

#ifndef LMP_FIX_SURFACE_GLOBAL_H
#define LMP_FIX_SURFACE_GLOBAL_H

#include "fix_surface.h"

#include <map>
#include <unordered_set>
#include <vector>

namespace LAMMPS_NS {

namespace Granular_NS {
  class GranularModel;
}

class FixSurfaceGlobal : public FixSurface {
 public:
  // neighbor lists for spheres with surfs and shear history
  // accessed by fix shear/history

  class NeighList *list;
  class NeighList *listhistory;

  FixSurfaceGlobal(class LAMMPS *, int, char **);
  ~FixSurfaceGlobal() override;
  int setmask() override;
  void post_constructor() override;

  void init() override;
  void setup_pre_neighbor() override;
  void initial_integrate(int) override;
  void pre_neighbor() override;
  void post_force(int) override;

  int modify_param(int, char **x) override;
  double memory_usage() override;

  void *extract(const char *, int &) override;
  int image(int *&, double **&) override;

 private:
  int neigh_style;

  int dimension, firsttime, use_history;
  double dt, skin;
  double flatthresh;
  double Twall;
  int tvar;
  char *tstr;

  std::string id_fix;

  // per-surf properties

  int maxsurftype;
  int maxsurfmol;
  double **xsurf, **vsurf, **omegasurf, *radsurf;

  // granular models

  struct ModelTypes {
    int plo, phi;
    int slo, shi;
  };

  ModelTypes *modeltypes;
  class Granular_NS::GranularModel **models;          // list of command-line models
  class Granular_NS::GranularModel ***types2model;    // model for each particle/surf type

  int nmodel, maxmodel;
  int history, size_history, heat_flag;

  // neighbor params

  double triggersq;
  int last_setup_bins;
  class NBinManual *nb;
  class NStencilManual *ns;

  // settings for motion applied to specific surf types

  struct Motion {
    int active;
    int mstyle;
    int vxflag, vyflag, vzflag;
    int axflag, ayflag, azflag;
    int xvar, yvar, zvar;
    int vxvar, vyvar, vzvar;
    double vx, vy, vz;
    double ax, ay, az;
    double dx, dy, dz;
    double period;
    double point[3], axis[3], unit[3];
    double omega;
    char *xvarstr, *yvarstr, *zvarstr;
    char *vxvarstr, *vyvarstr, *vzvarstr;
    double time_origin;
  };

  struct Motion *motions;    // list of defined motions, can be flagged inactive
  int nmotion, maxmotion;    // # of defined motions versus allocated size
  int anymove;               // 1 if any surf motion is enabled
  int anymove_variable;      // 1 if any surf motion is style VARIABLE

  int *mol2motion;    // assingment of surf mols (1 to Nmol) to motions
                       // -1 = non-moving surf mol

  double **points_original, **xsurf_original;
  double **points_lastneigh;
  int *pointmove;

  // storage of granular history info

  class FixNeighHistory *fix_history;
  double *zeroes;

  // rigid body masses for use in granular interactions

  class Fix *fix_rigid;    // ptr to rigid body fix, NULL if none
  double *mass_rigid;      // rigid mass for owned+ghost atoms
  int nmax;                // allocated size of mass_rigid

  // data structs for extracting surfs from molecule or STL files

  Point *points;                 // global list of unique points
  Line *lines;                   // global list of lines
  Tri *tris;                     // global list of tris
  int npoints, nlines, ntris;    // count of each
  int nedges;                    // count of unique tri edges
  int maxpoints;                 // allocated length of points
  int nsurf;                     // count of lines or tris for 2d/3d
  int nsurf_ghost;               // count of ghost lines or tris for 2d/3d
  int *idsurf;                   // id of corresponding local lines or tris

  // ragged 2d arrays for 2d connectivity

  int **neigh_p1;     // indices of other lines connected to endpt 1
  int **pwhich_p1;    // which point (0/1) on other line is endpt 1
  int **nside_p1;     // consistency of other line normal
                      //   SAME_SIDE or OPPOSITE_SIDE
  int **aflag_p1;     // is this line + other line CONCAVE or CONVEX surf
                      //   surf = on normal side of this line
  int **fflag_p1;     // is this line + other line FLAT or NONFLAT

  int **neigh_p2;     // ditto for connections to endpt 2
  int **pwhich_p2;    // ditto for endpt 2
  int **nside_p2;     // ditto for endpt 2
  int **aflag_p2;     // ditto for endpt 2
  int **fflag_p2;     // ditto for endpt 2

  // ragged 2d arrays for 3d edge connectivity

  int **neigh_e1;     // indices of other tris connected to edge 1
  int **ewhich_e1;    // which edge (0/1/2) on other tri is edge 1
  int **nside_e1;     // consistency of other tri normal
                      //   SAME_SIDE or OPPOSITE_SIDE
  int **aflag_e1;     // is this tri + other tri a FLAT,CONCAVE,CONVEX surf
                      //   surf = on normal side of this tri
  int **fflag_e1;     // is this tri + other tri FLAT or NONFLAT
  int **neigh_e2;     // ditto for connections to edge 2
  int **ewhich_e2;    // ditto for edge 2
  int **nside_e2;     // ditto for edge 2
  int **aflag_e2;     // ditto for edge 2
  int **fflag_e2;     // ditto for edge 2
  int **neigh_e3;     // ditto for connections to edge 3
  int **ewhich_e3;    // ditto for edge 3
  int **nside_e3;     // ditto for edge 3
  int **aflag_e3;     // ditto for edge 3
  int **fflag_e3;     // ditto for edge 3

  // ragged 2d arrays for 3d corner connectivity

  int **neigh_c1;     // indices of other tris connected to cpt 1
  int **cwhich_c1;    // which corner point (0/1/2) on other tri is cpt 1
  int **nside_c1;     // consistency of other tri normal
                      //   SAME_SIDE or OPPOSITE_SIDE, only meaningful for FLAT
  int **fflag_c1;     // is this tri + other tri FLAT or NONFLAT
  int **neigh_c2;     // indices of other tris connected to cpt 21
  int **cwhich_c2;    // which corner point (0/1/2) on other tri is cpt 2
  int **nside_c2;     // ditto for corner 2
  int **fflag_c2;     // ditto for corner 2
  int **neigh_c3;     // indices of tris connected to cpt 3
  int **cwhich_c3;    // which corner point (0/1/2) on other tri is cpt 3
  int **nside_c3;     // ditto for corner 3
  int **fflag_c3;     // ditto for corner 3

  // per-surface 2d/3d connectivity

  Connect2d *connect2d;    // 2d connection info
  Connect3d *connect3d;    // 3d connection info

  std::vector<ContactSurf> contact_surfs;
  std::map<int, int> contacts_map;

  // data for DumpImage

  int *imflag;
  double **imdata;
  int imax;

  // private methods

  void check2d();
  void check3d();
  void connectivity2d_complete();
  void connectivity3d_complete();
  void check_molecules();
  void stats2d();
  void stats3d();

  void surface_attributes();

  // contact processing and force calculation

  void prewalk_connections2d();
  void prewalk_connections3d();
  void walk_connections2d(std::vector<int> &, std::unordered_set<int> &);
  void walk_connections3d(std::vector<int> &, std::unordered_set<int> &);
  double calculate_2d_forces(std::vector<int> &);
  double calculate_3d_forces(std::vector<int> &);
  void calculate_3d_edge_force(int, double *, double *, double *, double *);
  double dist_away_from_edge(double *, double *, double *, double *, double *);

  // surface movement

  int modify_param_move(Motion *, int, char **);
  void move_linear(int, int);
  void move_wiggle(int, int);
  void move_rotate(int, int);
  void move_transrotate(int, int);
  void move_rotate_point(int, double *, double *, double, double);
  void move_variable(int, int);
};

}    // namespace LAMMPS_NS

#endif
#endif
