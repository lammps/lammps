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
ComputeStyle(reduce,ComputeReduce);
// clang-format on
#else

#ifndef LMP_COMPUTE_REDUCE_H
#define LMP_COMPUTE_REDUCE_H

#include "compute.h"    // IWYU pragma: export

namespace LAMMPS_NS {

class ComputeReduce : public Compute {
 public:
  enum { SUM, SUMSQ, SUMABS, MINN, MAXX, AVE, AVESQ, AVEABS };
  enum { PERATOM, LOCAL };

  ComputeReduce(class LAMMPS *, int, char **);
  ~ComputeReduce() override;
  void init() override;
  double compute_scalar() override;
  void compute_vector() override;
  double memory_usage() override;

 protected:
  int me;
  int mode, nvalues;
  int *which, *argindex, *flavor, *value2index;
  char **ids;
  double *onevec;
  int *replace, *indices, *owner;
  int index;
  char *idregion;
  class Region *region;
  int maxatom;
  double *varatom;

  struct Pair {
    double value;
    int proc;
  };
  Pair pairme, pairall;

  virtual double compute_one(int, int);
  virtual bigint count(int);
  void combine(double &, double, int);
};

}    // namespace LAMMPS_NS

#endif
#endif
