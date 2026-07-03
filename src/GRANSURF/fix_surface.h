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

#ifndef LMP_FIX_SURFACE_H
#define LMP_FIX_SURFACE_H

#include "fix.h"
#include <map>
#include <tuple>

namespace LAMMPS_NS {

class FixSurface : public Fix {
 public:
  // connectivity info between lines or tris

  struct Connect2d {    // line connectivity

    // counts, not including self
    int np1, np2;    // # of lines connected to endpts 1/2

    int external_pt[2];    // whether p1 and p2 are external

    // pairs of endpoint connections
    int *neigh_p1;     // indices (or IDs) of lines connected to endpt 1
    int *neigh_p2;     // ditto for connections to endpt 2
    int *pwhich_p1;    // which point (0,1) on other line is endpt 1
    int *pwhich_p2;    // ditto for endpt 2
    int *nside_p1;     // consistency of other line normal
    int *nside_p2;     // ditto for endpt 2
                       //   SAME_SIDE = 2 normals are on same side of surf
                       //   OPPOSITE_SIDE = opposite sides of surf
    int *aflag_p1;     // is this line + other line a CONCAVE or CONVEX surf
    int *aflag_p2;     // ditto for endpt 2
                       //   surf = on normal side of this line
                       //   aflag = CONCAVE, CONVEX
    int *fflag_p1;     // is this line + other line are FLAT or  NONFLAT
    int *fflag_p2;     // ditto for endpt 2
  };

  struct Connect3d {    // tri connectivity

    // counts, not including self
    // also not including edge-connected tris for corner pts
    int ne1, ne2, ne3;    // # of tris connected to edges 1,2,3
    int nc1, nc2, nc3;    // # of tris connected to corner pts 1,2,3

    int external_cor[3];     // whether c1, c2, and c3 are external
    int external_edge[3];    // whether e1, e2, and e3 are external

    // pairs of edge connections
    int *neigh_e1;     // indices (or IDs) of tris connected to edge 1
    int *neigh_e2;     // ditto for connections to edge 2
    int *neigh_e3;     // ditto for connections to edge 3
    int *ewhich_e1;    // which edge (0,1,2) on other tri shares edge 1
    int *ewhich_e2;    // ditto for edge 2
    int *ewhich_e3;    // ditto for edge 3
    int *nside_e1;     // consistency of other tri normal
    int *nside_e2;     // ditto for edge 2
    int *nside_e3;     // ditto for edge 3
                       //   SAME_SIDE = 2 normals are on same side of surf
                       //   OPPOSITE_SIDE = opposite sides of surf
    int *aflag_e1;     // is this tri + other tri a FLAT,CONCAVE,CONVEX surf
    int *aflag_e2;     // ditto for edge 2
    int *aflag_e3;     // ditto for edge 3
                       //   surf = on normal side of this tri
                       //   aflag = FLAT, CONCAVE, CONVEX
    int *fflag_e1;     // is this line + other line are FLAT or  NONFLAT
    int *fflag_e2;     // ditto for edge 2
    int *fflag_e3;     // ditto for edge 3

    // pairs of corner pt connections
    int *neigh_c1;     // indices (or IDs) of tris connected to corner pt 1
    int *neigh_c2;     // ditto for connections to corner pt 2
    int *neigh_c3;     // ditto for connections to corner pt 3
    int *cwhich_c1;    // which corner pt (0,1,2) on other tri shares corner pt 1
    int *cwhich_c2;    // ditto for corner pt 2
    int *cwhich_c3;    // ditto for corner pt 3
    int *nside_c1;     // consistency of normals, only meaningful for FLAT
    int *nside_c2;     // ditto for corner 2
    int *nside_c3;     // ditto for corner 3
    int *fflag_c1;     // are tris FLAT or NONFLAT
    int *fflag_c2;     // ditto for edge 2
    int *fflag_c3;     // ditto for edge 3
  };

  // struct for storing contact data

  struct ContactSurf {
    std::vector<int> cindex;
    std::vector<int> caflag;
    std::vector<int> cwhich;
    int index, neigh_index, type, flag, nside, external, priority;
    int convex_index, ck1, ck2, caflag1, caflag2;
    double overlap, rsq_com, weight_contribution, rmag;
    double dr[3], surf_norm[3], dr_force[3];
  };

  // comparators for sorting lists of particle-surface contacts; the presort
  // variant compares absolute normal-displacement alignment because the
  // surface normal signs are only assigned by the connection pre-walk

  static bool contact_presort(const ContactSurf &, const ContactSurf &);
  static bool contact_sort(const ContactSurf &, const ContactSurf &);

  // shared classification values for line/tri connection flags

  enum { NONFLAT = 0, FLAT = 1 };
  enum { CONCAVE = 0, CONVEX = 1 };
  enum { SAME_SIDE = 0, OPPOSITE_SIDE = 1 };

  FixSurface(class LAMMPS *, int, char **);
  ~FixSurface() override = default;

 protected:
  // shared computation of connection flags between two lines (2d) or
  // two tris (3d); the callers perform the endpoint/edge matching

  void point_connection2d(const double *inorm, const double *jnorm, int iwhich, int jwhich,
                          double flatthresh, int &fflag, int &nside, int &aflag);
  void edge_connection3d(const double *inorm, const double *jnorm, const double *iedge,
                         int jpfirst, int jpsecond, double flatthresh, int &fflag, int &ewhich,
                         int &nside, int &aflag);
  void corner_connection3d(const double *inorm, const double *jnorm, int cwhich,
                           double flatthresh, int &fflag, int &nside);

  // surfs read from molecule or STL files

  struct Point {
    double x[3];
  };

  struct Line {
    int mol, type;     // molID and type of the line
    int p1, p2;        // indices of points in line segment
    double norm[3];    // unit normal to line = Z x (p2-p1)
  };

  struct Tri {
    int mol, type;     // modID and type of the triangle
    int p1, p2, p3;    // indices of points in triangle
    double norm[3];    // unit normal to tri plane = (p2-p1) x (p3-p1)
  };

  // methods common to both global and local surfs

  void extract_from_molecule(char *, std::map<std::tuple<double, double, double, int>, int> *,
                             int &, int &, Point *&, int &, Line *&, int &, Tri *&);
  void extract_from_stlfile(char *, int, int,
                            std::map<std::tuple<double, double, double, int>, int> *, int &, int &,
                            Point *&, int &, Tri *&);

  void connectivity2d_global(int, int, Line *, Connect2d *&, int **&, int **&);
  int connectivity3d_global(int, int, Tri *, Connect3d *&, int **&, int **&, int **&, int **&,
                            int **&, int **&);
};

}    // namespace LAMMPS_NS

#endif
