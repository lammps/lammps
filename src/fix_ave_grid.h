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

#ifdef FIX_CLASS
// clang-format off
FixStyle(ave/grid,FixAveGrid);
// clang-format on
#else

#ifndef LMP_FIX_AVE_GRID_H
#define LMP_FIX_AVE_GRID_H

#include "fix.h"

namespace LAMMPS_NS {

class FixAveGrid : public Fix {
 public:
  FixAveGrid(class LAMMPS *, int, char **);
  ~FixAveGrid() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void end_of_step() override;

  void reset_grid() override;

  int get_grid_by_name(char *, int &) override;
  void *get_grid_by_index(int) override;
  int get_griddata_by_name(int, char *, int &) override;
  void *get_griddata_by_index(int) override;
  double memory_usage() override;

 private:
  int nxgrid,nygrid,nzgrid;
  int nvalues;
  int nrepeat, irepeat;
  bigint nvalid, nvalid_last;
  int modeatom,modegrid;
  int normflag,scaleflag,ave,nwindow;
  int dimension;

  int *which, *argindex;
  char **ids;
  int *value2index, *value2grid, *value2data;

  class Grid2d *grid2d;
  class Grid3d *grid3d;
  int ngrid_buf1, ngrid_buf2;
  double *grid_buf1, *grid_buf2;

  int nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in;
  int nxlo_out,nxhi_out,nylo_out,nyhi_out,nzlo_out,nzhi_out;
  int ngridout;
  double shift;

  double **vec2d,***vec3d;
  double ***array2d,****array3d;
  double **count2d,***count3d;

  int **bin;
  int maxatom;

  double *vresult;
  int maxvar;

  void atom2grid();
  void grid2grid();
  void zero_grid();
  bigint nextvalid();
};

}    // namespace LAMMPS_NS

#endif
#endif
