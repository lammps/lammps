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

#ifdef DUMP_CLASS
// clang-format off
DumpStyle(grid/vtk,DumpGridVTK);
// clang-format on
#else

#ifndef LMP_DUMP_GRID_VTK_H
#define LMP_DUMP_GRID_VTK_H

#include "dump_grid.h"

namespace LAMMPS_NS {

class DumpGridVTK : public DumpGrid {
 public:
  DumpGridVTK(class LAMMPS *, int, char **);
  ~DumpGridVTK() override;

 protected:
  int mode;
  double *xcoord, *ycoord, *zcoord;

  // methods

  void init_style() override;
  void write_header(bigint) override;
  void write_data(int, double *) override;
  void write_footer() override;

  void xyz_grid();
};

}    // namespace LAMMPS_NS

#endif
#endif
