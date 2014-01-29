/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(voronoi/atom,ComputeVoronoi)

#else

#ifndef LMP_COMPUTE_VORONOI_H
#define LMP_COMPUTE_VORONOI_H

#include "compute.h"
#include "voro++.hh"

namespace LAMMPS_NS {

class ComputeVoronoi : public Compute {
 public:
  ComputeVoronoi(class LAMMPS *, int, char **);
  ~ComputeVoronoi();
  void init();
  void compute_peratom();
  void compute_vector();
  double memory_usage();

  int pack_comm(int, int *, double *, int, int *);
  void unpack_comm(int, int, double *);

 private:
  void processCell(voro::voronoicell_neighbor&, int);

  int nmax, rmax, maxedge, sgroupbit;
  char *radstr;
  double fthresh, ethresh;
  double **voro;
  double *edge, *sendvector, *rfield;
  enum { VOROSURF_NONE, VOROSURF_ALL, VOROSURF_GROUP } surface;
  bool onlyGroup;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Missing atom style variable for radical voronoi tesselation radius.

UNDOCUMENTED

E: Missing group name after keyword 'surface'.

UNDOCUMENTED

E: Could not find compute/voronoi surface group ID

UNDOCUMENTED

E: Missing maximum edge count after keyword 'edge_histo'.

UNDOCUMENTED

E: Missing minimum face area after keyword 'face_threshold'.

UNDOCUMENTED

E: Missing minimum edge length after keyword 'edge_threshold'.

UNDOCUMENTED

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

W: More than one compute voronoi/atom command

It is not efficient to use compute voronoi/atom more than once.

E: Variable name for voronoi radius set does not exist

UNDOCUMENTED

E: Variable for voronoi radius is not atom style

UNDOCUMENTED

U: Compute voronoi/atom not allowed for triclinic boxes

This is a current restriction of this command.

*/
