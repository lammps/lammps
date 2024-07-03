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

#include <stdio.h>
#include "fix.h"

namespace LAMMPS_NS {

class FixSurfaceGlobal : public Fix {
 public:
  // neighbor lists for spheres with surfs and shear history
  // accessed by fix shear/history

  class NeighList *list;
  class NeighList *listhistory;

  FixSurfaceGlobal(class LAMMPS *, int, char **);
  ~FixSurfaceGlobal();
  int setmask() override;
  void post_constructor() override;
  void init() override;
  void setup_pre_neighbor() override;
  void initial_integrate(int) override;
  void pre_neighbor() override;
  void post_force(int) override;

  int modify_param(int, char **x) override;
  void reset_dt() override;
  double memory_usage() override;

  void *extract(const char *, int &) override;
  int image(int *&, double **&) override;

 private:
  int dimension,firsttime,pairstyle,history;
  int shearupdate;
  double kn,kt,gamman,gammat,xmu;
  double dt;
  double skin;

  double **xsurf,**vsurf,**omegasurf,*radsurf;;

  // motion

  int mstyle;
  double rperiod,omega_rotate,time_origin,triggersq;
  double xscale,yscale,zscale;
  double rpoint[3],raxis[3],runit[3];
  double **points_original,**xsurf_original;
  double **points_lastneigh;

  // storage of granular history info

  class FixNeighHistory *fix_history;
  double *zeroes;

  // rigid body masses for use in granular interactions

  class Fix *fix_rigid;    // ptr to rigid body fix, NULL if none
  double *mass_rigid;      // rigid mass for owned+ghost atoms
  int nmax;                // allocated size of mass_rigid

  // data structs for extracting surfs from molecule files

  struct Point {
    double x[3];
  };

  struct Line {
    int mol,type;           // molID and type of the line
    int p1,p2;              // indices of points in line segment
                            // rhand rule: Z x (p2-p1) = outward normal
  };

  struct Tri {
    int mol,type;           // modID and type of the triangle
    int p1,p2,p3;           // indices of points in triangle
                            // rhand rule: (p2-p1) x (p3-p1) = outward normal
    double norm[3];         // unit normal to tri plane
  };

  Point *points;              // global list of points
  Line *lines;                // global list of lines
  Tri *tris;                  // global list of tris
  int npoints,nlines,ntris;   // count of each
  int nsurf;                  // count of lines or tris for 2d/3d

  int **clist;                // ragged 2d array for global corner pt lists

  // 2d/3d connectivity

  struct Connect2d {      // line connectivity
    tagint neigh_p1;      // ID of line connected to pt1, 0 if none
    tagint neigh_p2;      // ditto for pt2
    int flags;            // flags for end pt coupling
                          // NOTE: document what is stored in flags
  };

  struct Connect3d {      // tri connectivity
    tagint neigh_e1;      // ID of tri connected to 1-2 edge, 0 if none
    tagint neigh_e2;      // ditto for 2-3 edge
    tagint neigh_e3;      // ditto for 3-1 edge
    tagint *neigh_c1;     // IDs of tris connected to 1st corner pt
    tagint *neigh_c2;     // ditto for 2nd corner pt
    tagint *neigh_c3;     // ditto for 3rd corner pt
    int nc1,nc2,nc3;      // # of tris connected to each corner pt
    int flags;            // flags for edge and corner pt coupling
                          // NOTE: document what is stored in flags
  };

  Connect2d *connect2d;             // 2d connection info
  Connect3d *connect3d;             // 3d connection info

  // data for DumpImage

  int *imflag;
  double **imdata;
  int imax;

  // private methods

  void hooke(int, int, double, double, double, double *, double *, double);
  void hooke_history(int, int, double, double, double, double, double, double,
                     double *, double *, double, double *);

  int overlap_sphere_line(int, int, double *, double *, double &);
  int endpt_neigh_check(int, int, int);
  int overlap_sphere_tri(int, int, double *, double *, double &);
  int nearest_point_line(double *, double *, double *, double *);
  int edge_neigh_check(int, int, int);
  int corner_neigh_check(int, int, int);

  void extract_from_molecules(char *);
  void connectivity2d_global();
  void connectivity3d_global();
  void set_attributes();
  void move_init();
  void move_clear();
};

}    // namespace LAMMPS_NS

#endif
#endif
