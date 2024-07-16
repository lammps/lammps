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
    int np1,np2;          // # of lines connected to pts 1,2 (including self)
    tagint *neigh_p1;     // IDs of all lines connected to pt1 (if np1 > 1)
    tagint *neigh_p2;     // ditto for pt2
    int flags;            // future flags for end pt coupling
    int indexp1,indexp2;  // pool indices of neigh_p12 chunks
    int ilocal;           // local index of line particle
  };

  struct Connect3d {      // tri connectivity
    int ne1,ne2,ne3;      // # of tris connected to edges 1,2,3 (including self)
    int nc1,nc2,nc3;      // # of tris connected to corner pts 1,2,3 (including self)
    tagint *neigh_e1;     // IDs of all tris connected to 1-2 edge (if ne1 > 1)
    tagint *neigh_e2;     // ditto for 2-3 edge
    tagint *neigh_e3;     // ditto for 3-1 edge
    tagint *neigh_c1;     // IDs of all tris connected to corner pt 1 (if nc1 > 1)
    tagint *neigh_c2;     // ditto for corner pt 2
    tagint *neigh_c3;     // ditto for corner pt 3
    int flags;            // future flags for edge and corner pt coupling
    int indexe1,indexe2,indexe3;   // pool indices of neigh_e123 chunks
    int indexc1,indexc2,indexc3;   // pool indices of neigh_c123 chunks
    int ilocal;           // local index of triangle particle
  };

  // NOTE: can some of these be private ?

  int *cindex;            // per-atom index into connect 2d/3d vecs, -1 if none
  Connect2d *connect2d;         // 2d connection info
  Connect3d *connect3d;         // 3d connection info
  MyPoolChunk<tagint> *tcp;     // allocator for 2d/3d connectivity vecs

  // size of local/ghost connection info vectors

  int nlocal_connect,nghost_connect,nmax_connect;

  FixSurfaceLocal(class LAMMPS *, int, char **);
  virtual ~FixSurfaceLocal();
  int setmask() override;
  void post_constructor() override;
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
  char *sourceID;

  class AtomVecLine *avec_line;
  class AtomVecTri *avec_tri;

  // data structs for calculating global connectivity of line/tri particles
  // only used by Rvous comm during setup
  
  struct InRvous {
    int proc, ibin, ilocal, ipoint;
    tagint atomID;
    double x[3];
  };
  
  struct OutRvous {
    int ilocal, ipoint;
    tagint atomID;
  };

  double bboxlo[3],bboxhi[3];    // bounding box around all lines/tris

  // data structs for extracting surfs from molecule template or STL file
  // these are global data structs for all surfs, including connectivity
  // only used during setup, then deleted
  
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

  int **plist;                // ragged 2d array for global line end pt lists
  int **elist;                // ragged 2d array for global tri edge lists
  int **clist;                // ragged 2d array for global tri corner pt lists

  // data structs for binning end pts of lines or tris from data file
  // these are for local surfs which procs already own as line or tri style atoms
  // only used during setup, then deleted

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
  static int point_match(int, char *, int &, int *&, char *&, void *);

  // private methods

  void connectivity2d_local();
  void connectivity3d_local();
  void calculate_endpts();
  int overlap_bins_2d(double *, double, int *);
  void calculate_corners();
  int overlap_bins_3d(double *, double, int *);

  int check_exist();
  void extract_from_molecules(char *);
  void extract_from_stlfile(char *);
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
