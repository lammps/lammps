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
ComputeStyle(property/local,ComputePropertyLocal);
// clang-format on
#else

#ifndef LMP_COMPUTE_PROPERTY_LOCAL_H
#define LMP_COMPUTE_PROPERTY_LOCAL_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputePropertyLocal : public Compute {
 public:
  ComputePropertyLocal(class LAMMPS *, int, char **);
  ~ComputePropertyLocal() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void compute_local() override;
  double memory_usage() override;

 private:
  int nvalues, kindflag, cutstyle;

  int nmax;
  double *vlocal;
  double **alocal;
  double *buf;

  class NeighList *list;

  int ncount;
  int **indices;

  int count_pairs(int, int);
  int count_bonds(int);
  int count_angles(int);
  int count_dihedrals(int);
  int count_impropers(int);
  void reallocate(int);

  typedef void (ComputePropertyLocal::*FnPtrPack)(int);
  FnPtrPack *pack_choice;    // ptrs to pack functions

  void pack_patom1(int);
  void pack_patom2(int);
  void pack_ptype1(int);
  void pack_ptype2(int);

  void pack_batom1(int);
  void pack_batom2(int);
  void pack_btype(int);

  void pack_aatom1(int);
  void pack_aatom2(int);
  void pack_aatom3(int);
  void pack_atype(int);

  void pack_datom1(int);
  void pack_datom2(int);
  void pack_datom3(int);
  void pack_datom4(int);
  void pack_dtype(int);

  void pack_iatom1(int);
  void pack_iatom2(int);
  void pack_iatom3(int);
  void pack_iatom4(int);
  void pack_itype(int);
};

}    // namespace LAMMPS_NS

#endif
#endif
