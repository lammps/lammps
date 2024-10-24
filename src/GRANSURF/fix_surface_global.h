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

namespace Granular_NS {
  class GranularModel;
}

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
  int dimension,firsttime,use_history;
  double dt, skin;

  // for granular model choices

  class Granular_NS::GranularModel *model;
  int history, size_history, heat_flag;

  double Twall;
  int tvar;
  char *tstr;

  double **xsurf,**vsurf,**omegasurf,*radsurf;

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
    //double norm[3];         // unit normal to line = Z x (p2-p1)
                              // not currently set or used
  };

  struct Tri {
    int mol,type;           // modID and type of the triangle
    int p1,p2,p3;           // indices of points in triangle
    double norm[3];         // unit normal to tri plane = (p2-p1) x (p3-p1)
  };

  Point *points;              // global list of unique points
  Line *lines;                // global list of lines
  Tri *tris;                  // global list of tris
  int npoints,nlines,ntris;   // count of each
  int nsurf;                  // count of lines or tris for 2d/3d

  int **plist;                // ragged 2d array for global line end pt lists
  int **elist;                // ragged 2d array for global tri edge lists
  int **clist;                // ragged 2d array for global tri corner pt lists

  // 2d/3d connectivity

  struct Connect2d {      // line connectivity
    int np1,np2;          // # of lines connected to pts 1,2 (including self)
    int *neigh_p1;        // indices of all lines connected to pt1 (if np1 > 1)
    int *neigh_p2;        // ditto for pt2
    int flags;            // future flags for end pt coupling
  };

  struct Connect3d {      // tri connectivity
    int ne1,ne2,ne3;      // # of tris connected to edges 1,2,3 (including self)
    int nc1,nc2,nc3;      // # of tris connected to corner pts 1,2,3 (including self)
    int *neigh_e1;        // indices of all tris connected to 1-2 edge (if ne1 > 1)
    int *neigh_e2;        // ditto for 2-3 edge
    int *neigh_e3;        // ditto for 3-1 edge
    int *neigh_c1;        // indices of all tris connected to corner pt 1 (if nc1 > 1)
    int *neigh_c2;        // ditto for corner pt 2
    int *neigh_c3;        // ditto for corner pt 3
    int flags;            // future flags for edge and corner pt coupling
  };

  Connect2d *connect2d;             // 2d connection info
  Connect3d *connect3d;             // 3d connection info

  // data for DumpImage

  int *imflag;
  double **imdata;
  int imax;

  // private methods

  int overlap_sphere_line(int, int, double *, double *, double &);
  int endpt_neigh_check(int, int, int);
  int overlap_sphere_tri(int, int, double *, double *, double &);
  int nearest_point_line(double *, double *, double *, double *);
  int edge_neigh_check(int, int, int);
  int corner_neigh_check(int, int, int);

  void extract_from_molecules(char *);
  void extract_from_stlfile(char *);
  void connectivity2d_global();
  void connectivity3d_global();
  void surface_attributes();
  void move_init();
  void move_clear();
};

}    // namespace LAMMPS_NS

#endif
#endif
