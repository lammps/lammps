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
ComputeStyle(esp/grid/local, ComputeESPGridLocal)
// clang-format on
#else

#ifndef LMP_COMPUTE_ESP_GRID_LOCAL_H
#define LMP_COMPUTE_ESP_GRID_LOCAL_H

#include "compute_grid_local.h"

namespace LAMMPS_NS {

class Grid3d;

class ComputeESPGridLocal : public ComputeGridLocal {
 public:
  ComputeESPGridLocal(class LAMMPS *, int, char **);
  ~ComputeESPGridLocal() override;

  void init() override;
  void compute_local() override;

  void *extract_reference();

 protected:
  int reaxflag;
  Grid3d *grid;
  double ***esp_ref;
  double *bcut_acks2;
  int nx, ny, nz;
  double xlo, ylo, zlo, spacing;
};

}

#endif
#endif
