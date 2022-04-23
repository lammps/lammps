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
  ~ComputeMLIAP() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void compute_array() override;
  void generate_neigharrays();
  void grow_neigharrays();
  double memory_usage() override;

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
