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

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(snad/atom,ComputeSNADAtom);
// clang-format on
#else

#ifndef LMP_COMPUTE_SNAD_ATOM_H
#define LMP_COMPUTE_SNAD_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeSNADAtom : public Compute {
 public:
  ComputeSNADAtom(class LAMMPS *, int, char **);
  ~ComputeSNADAtom() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void compute_peratom() override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  double memory_usage() override;

 private:
  int nmax;
  int ncoeff, nvalues, yoffset, zoffset;
  double **cutsq;
  class NeighList *list;
  double **snad;
  double rcutfac;
  double *radelem;
  double *wjelem;
  int *map;    // map types to [0,nelements)
  int nelements, chemflag;
  int switchinnerflag;
  double *sinnerelem;
  double *dinnerelem;
  class SNA *snaptr;
  double cutmax;
  int quadraticflag;
};

}    // namespace LAMMPS_NS

#endif
#endif
