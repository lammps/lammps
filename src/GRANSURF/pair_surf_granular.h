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

PairStyle(surf/granular,PairSurfGranular)

#else

#ifndef LMP_PAIR_SURF_GRANULAR_H
#define LMP_PAIR_SURF_GRANULAR_H

#include "pair_granular.h"
#include "fix_surface_local.h"

namespace LAMMPS_NS {

class PairSurfGranular : public PairGranular {
 public:
  PairSurfGranular(class LAMMPS *);
  ~PairSurfGranular() override;
  void compute(int, int) override;
  void init_style() override;
  double memory_usage() override;

 protected:
  int surfmoveflag;

  int style;
  int emax;                // allocated size of endpt list
  double **endpts;         // current end pts of each line
                           // Nall x 4 array for local + ghost atoms

  int cmax;                // allocated size of corners
  double **corners;        // current corner pts and norm of each tri
                           // Nall x 12 array for local + ghost atoms

  // ptr to AtomVec for bonus info

  class AtomVecLine *avecline;
  class AtomVecTri *avectri;

  // line connectivity info for owned and ghost lines

  class FixSurfaceLocal *fsl;              // ptr to surface/local fix
  FixSurfaceLocal::Connect2d *connect2d;   // ptr to connectivity info
  FixSurfaceLocal::Connect3d *connect3d;   // ptr to connectivity info
  MyPoolChunk<int> *tcp;                   // allocator for connectivity info

  // lines

  void calculate_endpts();
  int overlap_sphere_line(int, int, double *, double *, double &);
  int endpt_neigh_check(int, int, int);

  // tri

  void calculate_corners();
  void corners2norm(double *, double *);
  int overlap_sphere_tri(int, int, double *, double *, double &);
  int nearest_point_line(double *, double *, double *, double *);
  int edge_neigh_check(int, int, int);
  int corner_neigh_check(int, int, int);
};

}    // namespace LAMMPS_NS

#endif
#endif
