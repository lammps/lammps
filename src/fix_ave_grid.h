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

  void pack_reverse_grid(int, void *, int, int *) override;
  void unpack_reverse_grid(int, void *, int, int *) override;
  void pack_remap_grid(int, void *, int, int *) override;
  void unpack_remap_grid(int, void *, int, int *) override;

  void reset_grid() override;

  int get_grid_by_name(const std::string &, int &) override;
  void *get_grid_by_index(int) override;
  int get_griddata_by_name(int, const std::string &, int &) override;
  void *get_griddata_by_index(int) override;

  double memory_usage() override;

 private:
  int nxgrid,nygrid,nzgrid;
  int nvalues;
  int nrepeat, irepeat;
  bigint nvalid, nvalid_last;
  int modeatom, modegrid;
  int discardflag, normflag, aveflag, nwindow;

  double maxdist;
  int running_count;
  int window_count,window_oldest,window_newest;

  int biasflag;
  char *id_bias;
  class Compute *tbias;    // ptr to additional bias compute
  double adof,cdof;

  int dimension,triclinic;

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

  struct GridData {
    double **vec2d,***vec3d;
    double ***array2d,****array3d;
    double **count2d,***count3d;
  };

  GridData *grid_output;
  GridData *grid_sample,*grid_nfreq,*grid_running;
  GridData **grid_window;

  // old grid data for remap operation

  class Grid2d *grid2d_previous;
  class Grid3d *grid3d_previous;
  int nxlo_out_previous,nylo_out_previous,nzlo_out_previous;
  GridData *grid_sample_previous,*grid_nfreq_previous,*grid_running_previous;
  GridData **grid_window_previous;

  int **bin;
  int *skip;
  int maxatom;

  double *vresult;
  int maxvar;

  void atom2grid();
  void grid2grid();

  void normalize_atom(int, GridData *);
  void normalize_grid(int, GridData *);
  void normalize_count(int, GridData *);

  void allocate_grid();
  GridData *allocate_one_grid();
  GridData *clone_one_grid(GridData *);
  void deallocate_one_grid(GridData *, int, int, int);

  int pack_one_grid(GridData *, int, double *);
  int unpack_one_grid(double *, GridData *, int);

  double size_grid(GridData *);
  void zero_grid(GridData *);
  void copy_grid(GridData *, GridData *);
  void add_grid(GridData *, GridData *);
  void subtract_grid(GridData *, GridData *);
  void output_grid(GridData *);

  bigint nextvalid();
};

}    // namespace LAMMPS_NS

#endif
#endif
