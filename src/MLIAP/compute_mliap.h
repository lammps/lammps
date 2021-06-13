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
ComputeStyle(mliap,ComputeMLIAP);
// clang-format on
#else

#ifndef LMP_COMPUTE_MLIAP_H
#define LMP_COMPUTE_MLIAP_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeMLIAP : public Compute {
 public:
  ComputeMLIAP(class LAMMPS *, int, char **);
  ~ComputeMLIAP();
  void init();
  void init_list(int, class NeighList *);
  void compute_array();
  void generate_neigharrays();
  void grow_neigharrays();
  double memory_usage();

 private:
  double **mliaparray, **mliaparrayall;
  class NeighList *list;
  int *map;            // map types to [0,nelements)
  int ndescriptors;    // number of descriptors
  int nparams;         // number of model parameters per element
  int nelements;
  int gradgradflag;    // 1 for graddesc, 0 for gamma
  class MLIAPModel *model;
  class MLIAPDescriptor *descriptor;
  class MLIAPData *data;

  Compute *c_pe;
  Compute *c_virial;
  std::string id_virial;

  int lastcol;

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

E: Compute snap requires a pair style be defined

Self-explanatory.

E: Compute snap cutoff is longer than pairwise cutoff

UNDOCUMENTED

W: More than one compute snad/atom

Self-explanatory.

*/
