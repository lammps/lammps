/* -*- c++ -*- ----------------------------------------------------------
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

ComputeStyle(snap,ComputeSnap)

#else

#ifndef LMP_COMPUTE_SNAP_H
#define LMP_COMPUTE_SNAP_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeSnap : public Compute {
 public:
  ComputeSnap(class LAMMPS *, int, char **);
  ~ComputeSnap();
  void init();
  void init_list(int, class NeighList *);
  void compute_array();
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  double memory_usage();

 private:
  int natoms, nmax, size_peratom;
  int ncoeff, nperdim, yoffset, zoffset;
  int virialoffset, ndims_peratom;
  double **cutsq;
  class NeighList *list;
  double **snap;
  double **snap_peratom;
  double rcutfac;
  double *radelem;
  double *wjelem;
  class SNA* snaptr;
  double cutmax;
  int quadraticflag;

  Compute *c_pe;
  Compute *c_virial;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute snap requires a pair style be defined

Self-explanatory.

E: Compute snap cutoff is longer than pairwise cutoff

UNDOCUMENTED

W: More than one compute snad/atom

Self-explanatory.

*/
