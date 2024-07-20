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

PairStyle(tri/gran/hooke/history,PairTriGranHookeHistory)

#else

#ifndef LMP_PAIR_TRI_GRAN_HOOKE_HISTORY_H
#define LMP_PAIR_TRI_GRAN_HOOKE_HISTORY_H

#include "pair.h"
#include "fix_surface_local.h"
#include "my_pool_chunk.h"

namespace LAMMPS_NS {

class PairTriGranHookeHistory : public Pair {
 public:
  PairTriGranHookeHistory(class LAMMPS *);
  virtual ~PairTriGranHookeHistory();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  void reset_dt();
  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);
  double memory_usage();

 protected:
  double kn,kt,gamman,gammat,xmu;
  int dampflag,limit_damping;
  double dt;
  int freeze_group_bit;
  int history;

  int neighprev;
  double *onerad_dynamic,*onerad_frozen;
  double *maxrad_dynamic,*maxrad_frozen;

  int size_history;
  int surfmoveflag;

  class FixDummy *fix_dummy;
  class FixNeighHistory *fix_history;

  int cmax;                // allocated size of corners
  double **corners;        // current corner pts and norm of each tri
                           // Nall x 12 array for local + ghost atoms

  // ptr to AtomVec for bonus tri info

  class AtomVecTri *avec;

  // tri connectivity info for owned and ghost tris

  class FixSurfaceLocal *fsl;              // ptr to surface/local fix
  FixSurfaceLocal::Connect3d *connect3d;   // ptr to connectivity info
  MyPoolChunk<int> *tcp;                   // allocator for connectivity info

  // storage of rigid body masses for use in granular interactions

  class Fix *fix_rigid;    // ptr to rigid body fix, NULL if none
  double *mass_rigid;      // rigid mass for owned+ghost atoms
  int nmax;                // allocated size of mass_rigid

  // local methods

  void allocate();
  void calculate_corners();
  void corners2norm(double *, double *);
  int overlap_sphere_tri(int, int, double *, double *, double &);
  int nearest_point_line(double *, double *, double *, double *);
  int edge_neigh_check(int, int, int);
  int corner_neigh_check(int, int, int);
};

}

#endif
#endif
