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
ComputeStyle(esp/grid, ComputeESPGrid)
// clang-format on
#else

#ifndef LMP_COMPUTE_ESP_GRID_H
#define LMP_COMPUTE_ESP_GRID_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeESPGrid : public Compute {
 public:
  ComputeESPGrid(class LAMMPS *, int, char **);
  ~ComputeESPGrid() override;

  void init() override;
  double compute_scalar() override;
  void *extract_reference();

 private:
  double spacing;
  int nx, ny, nz;
  double xlo, ylo, zlo;

  double ***esp_ref;   // to be filled externally
  double *bcut_acks2;  // per-type ACKS2 outer cutoff
  int reaxflag;

  inline double weight(double r, double rcut) const;
};

}  // namespace LAMMPS_NS

#endif
#endif

