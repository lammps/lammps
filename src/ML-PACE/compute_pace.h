/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(pace,ComputePACE);
// clang-format on
#else

#ifndef LMP_COMPUTE_PACE_H
#define LMP_COMPUTE_PACE_H

#include "compute.h"
#include "ace-evaluator/ace_c_basis.h"
#include "ace-evaluator/ace_evaluator.h"
#include "ace-evaluator/ace_abstract_basis.h"

namespace LAMMPS_NS {

class ComputePACE : public Compute {
 public:
  ComputePACE(class LAMMPS *, int, char **);
  ~ComputePACE();
  void init();
  void init_list(int, class NeighList *);
  void compute_array();
  double memory_usage();

 private:
  int natoms, nmax, size_peratom, lastcol;
  int ncoeff, nvalues, nperdim, yoffset, zoffset;
  int ndims_peratom, ndims_force, ndims_virial;
  double **cutsq;
  class NeighList *list;
  double **pace, **paceall;
  double **pace_peratom;
  double rcutfac;
  int *map;    // map types to [0,nelements)
  int nelements, chemflag;
  int bikflag, bik_rows, dgradflag, dgrad_rows;
  double *cg;
  class ACECTildeEvaluator *ace;
  class ACECTildeBasisSet *basis_set;
  double cutmax;
  Compute *c_pe;
  Compute *c_virial;

  void dbdotr_compute();
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute pace requires a pair style be defined

Self-explanatory.

E: Compute pace cutoff is longer than pairwise cutoff

UNDOCUMENTED

W: More than one compute pace/atom

Self-explanatory.

*/
