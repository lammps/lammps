/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(sna/grid,ComputeSNAGrid);
// clang-format on
#else

#ifndef LMP_COMPUTE_SNA_GRID_H
#define LMP_COMPUTE_SNA_GRID_H

#include "compute_grid.h"

namespace LAMMPS_NS {

class ComputeSNAGrid : public ComputeGrid {
 public:
  ComputeSNAGrid(class LAMMPS *, int, char **);
  ~ComputeSNAGrid() override;
  void init() override;
  void compute_array() override;
  double memory_usage() override;

 private:
  int ncoeff;
  double **cutsq;
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
