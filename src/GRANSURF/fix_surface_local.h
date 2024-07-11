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

FixStyle(surface/local,FixSurfaceLocal)

#else

#ifndef LMP_FIX_SURFACE_LOCAL_H
#define LMP_FIX_SURFACE_LOCAL_H

#include <stdio.h>
#include "fix.h"
#include "my_pool_chunk.h"

namespace LAMMPS_NS {

class FixSurfaceLocal : public Fix {
 public:

  // 2d/3d connectivity

  struct Connect2d {      // line connectivity
    int np1,np2;          // # of other lines connected to pts 1,2
    tagint *neigh_p1;     // IDs of other lines connected to pt1
    tagint *neigh_p2;     // ditto for pt2
    int flags;            // future flags for end pt coupling
    int ilocal;           // local index of line particle
  };

  struct Connect3d {      // tri connectivity
    int ne1,ne2,ne3;      // # of other tris connected to edges 1,2,3
    int nc1,nc2,nc3;      // # of tris connected to corner pts 1,2,3
    tagint *neigh_e1;     // IDs of other tris connected to 1-2 edge
    tagint *neigh_e2;     // ditto for 2-3 edge
    tagint *neigh_e3;     // ditto for 3-1 edge
    tagint *neigh_c1;     // IDs of other tris connected to corner pt 1
    tagint *neigh_c2;     // ditto for corner pt 2
    tagint *neigh_c3;     // ditto for corner pt 3
    int flags;            // future flags for edge and corner pt coupling
    int indexe1,indexe2,indexe3;   // pool indices of neigh_e123 chunks
    int indexc1,indexc2,indexc3;   // pool indices of neigh_c123 chunks
    int ilocal;           // local index of triangle particle
  };

  int *index;                       // per-atom index into connect 2d/3d vecs
  Connect2d *connect2d;             // 2d connection info
  Connect3d *connect3d;             // 3d connection info
  MyPoolChunk<tagint> *tcp;         // allocator for 2d or 3d connectivity

  // size of local/ghost connection info vectors

  int nlocal_connect,nghost_connect,nmax_connect;

  FixSurfaceLocal(class LAMMPS *, int, char **);
  virtual ~FixSurfaceLocal();
  int setmask() override;
  void setup_pre_neighbor() override;
  void pre_neighbor() override;

  void grow_arrays(int) override;
  void grow_connect();
  void copy_arrays(int, int, int) override;
  void set_arrays(int) override;
  void clear_bonus() override;

  int pack_border(int, int *, double *) override;
  int unpack_border(int, int, double *) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;

  double memory_usage() override;

 private:
  int dimension,mode;
  char *idmol;
  class AtomVecLine *avec_line;
  class AtomVecTri *avec_tri;

  // data structs for extracting surfs from molecule files

  struct Point {
    double x[3];
  };

  struct Line {
    int mol,type;           // molID and type of the element
    int p1,p2;              // indices of points in line segment
                            // rhand rule: Z x (p2-p1) = outward normal
  };

  struct Tri {
    int mol,type;           // modID and type of the element
    int p1,p2,p3;           // indices of points in triangle
                            // rhand rule: (p2-p1) x (p3-p1) = outward normal
  };

  Point *points;              // global list of points
  Line *lines;                // global list of lines
  Tri *tris;                  // global list of tris
  int npoints,nlines,ntris;   // count of each

  Connect2d *connect2dall;    // global connectivity info
  Connect3d *connect3dall;
  int **clist;                // ragged 2d array for global corner pt lists

  // data structs for binning end pts of lines or tris from data file
  // these are surfs procs already own as line or tri style atoms

  double **endpts;         // current end pts of lines I own
                           // Nlocal x 4 array for local atoms
  double **corners;        // current corner pts of tris I own
                           // Nlocal x 9 array for local atoms

  struct OnePt {               // one end/corner point of iline/itri in a bin
    int iatom;                 // local index of the line/tri in atom list
    int iconnect;              // local index of the line/tri in connect list
    int ptwhich;               // 1/2 for two end pts of line
                               // 1/2/3 for three corner pts of tri
  };

  OnePt *pts;                      // list of pts in all bins
  int *bincount;                   // # of pts per bin
  int *binfirst;                   // index into pts of first pt in bin
  int nbins,nbinx,nbiny,nbinz;     // # of total bins and in each dim
  double invbinx,invbiny,invbinz;  // inverse of bin size
  double binlo[3],binhi[3];        // bounds of region that is binned

  double epssq;                // distance tolerance for end pts
                               // from different lines to be connected

  int nmatch;                  // # of line connections
  int nmatch1,nmatch2;         // # of tri connections
  int errormatch;              // # of errors with line connectivity
  int errormatch1,errormatch2; // # of errors with tri connectivity

  int vecflag;            // 0/1 whether tri matching should also
                          // store variable-length vecs of corner connections

  // static variable for ring communication callback to access class data
  // callback functions for ring communication

  static FixSurfaceLocal *fptr;
  static void linematch(int, char *);
  static void trimatch(int, char *);

  // private methods

  void copy_connect(int, int);

  void connectivity2d_local();
  void calculate_endpts(int);
  int pt2bin2d(double *);
  int overlap2bin2d(double *, double, int *);

  void connectivity3d_local();
  void calculate_corners(int);
  int pt2bin3d(double *);
  int overlap2bin3d(double *, double, int *);

  void extract_from_molecules(char *);
  void connectivity2d_global();
  void connectivity3d_global();
  void assign2d();
  void assign3d();
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
