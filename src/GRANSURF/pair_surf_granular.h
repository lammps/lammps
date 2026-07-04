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
PairStyle(surf/granular,PairSurfGranular)
// clang-format on
#else

#ifndef LMP_PAIR_SURF_GRANULAR_H
#define LMP_PAIR_SURF_GRANULAR_H

#include "pair_granular.h"

#include "fix_surface.h"
#include "fix_surface_local.h"

#include <map>
#include <unordered_set>
#include <vector>

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
  int emax;           // allocated size of endpt list
  double **endpts;    // current end pts of each line
                      // Nall x 6 array for local + ghost atoms

  int cmax;            // allocated size of corners
  double **corners;    // current corner pts and norm of each tri
                       // Nall x 12 array for local + ghost atoms

  // ptr to AtomVec for bonus info

  class AtomVecLine *avecline;
  class AtomVecTri *avectri;

  // line connectivity info for owned and ghost lines

  class FixSurfaceLocal *fsl;               // ptr to surface/local fix
  FixSurfaceLocal::Connect2d *connect2d;    // ptr to connectivity info
  FixSurfaceLocal::Connect3d *connect3d;    // ptr to connectivity info
  MyPoolChunk<int> *tcp;                    // allocator for connectivity info

  std::vector<FixSurface::ContactSurf> contact_surfs;
  std::map<int, int> contacts_map;

  // arrays from fix surface/local

  int *atom2connect;    // per-atom index into connect 2d/3d vecs, -1 if none

  // lines and tris

  void calculate_endpts();
  void calculate_corners();
  void corners2norm(double *, double *);

  // contact processing and force calculation

  void prewalk_connections2d();
  void prewalk_connections3d();
  void walk_connections2d(std::vector<int> &, std::unordered_set<int> &);
  void walk_connections3d(std::vector<int> &, std::unordered_set<int> &);
  double calculate_2d_forces(std::vector<int> &);
  double calculate_3d_forces(std::vector<int> &);
  void calculate_3d_edge_force(int, double *, double *, double *, double *);
  double dist_away_from_edge(double *, double *, double *, double *, double *);
  double dist_away_from_2_tris(double *, double *, double *);
};

}    // namespace LAMMPS_NS

#endif
#endif
