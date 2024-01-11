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
FixStyle(cmap,FixCMAP);
// clang-format on
#else

#ifndef LMP_FIX_CMAP_H
#define LMP_FIX_CMAP_H

#include "fix.h"
namespace LAMMPS_NS {

class FixCMAP : public Fix {
 public:
  FixCMAP(class LAMMPS *, int, char **);
  ~FixCMAP() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void setup_pre_neighbor() override;
  void setup_pre_reverse(int, int) override;
  void min_setup(int) override;
  void pre_neighbor() override;
  void pre_reverse(int, int) override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  void min_post_force(int) override;
  double compute_scalar() override;

  void read_data_header(char *) override;
  void read_data_section(char *, int, char *, tagint) override;
  bigint read_data_skip_lines(char *) override;
  void write_data_header(FILE *, int) override;
  void write_data_section_size(int, int &, int &) override;
  void write_data_section_pack(int, double **) override;
  void write_data_section_keyword(int, FILE *) override;
  void write_data_section(int, FILE *, int, double **, int) override;

  void write_restart(FILE *) override;
  void restart(char *) override;
  int pack_restart(int, double *) override;
  void unpack_restart(int, int) override;
  int size_restart(int) override;
  int maxsize_restart() override;

  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  void set_arrays(int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;

  double memory_usage() override;

 private:
  int eflag_caller;
  int ctype, ilevel_respa;
  int ncrosstermtypes, crossterm_per_atom, maxcrossterm;
  int ncrosstermlist;
  bigint ncmap;

  int **crosstermlist;

  int nmax_previous;
  int *num_crossterm;
  int **crossterm_type;
  tagint **crossterm_atom1, **crossterm_atom2, **crossterm_atom3;
  tagint **crossterm_atom4, **crossterm_atom5;

  double E, dEdPhi, dEdPsi;
  double ecmap;
  double fcmap[4], cij[4][4];
  double *g_axis;

  // CMAP grid points obtained from external file

  double ***cmapgrid;

  // partial derivatives and cross-derivatives of the grid data

  double ***d1cmapgrid, ***d2cmapgrid, ***d12cmapgrid;

  // read map grid data

  void read_grid_map(char *);

  // read in CMAP cross terms from LAMMPS data file

  void read_cmap_data(int, char *);

  // pre-compute the partial and cross-derivatives of map grid points

  void set_map_derivatives(double **, double **, double **, double **);

  // cubic spline interpolation functions for derivatives of map grid points

  void spline(double *, double *, int);
  void spl_interpolate(double, double *, double *, double &, double &);

  // calculate dihedral angles

  double dihedral_angle_atan2(double, double, double, double, double, double, double, double,
                              double, double);

  // calculate bicubic interpolation coefficient matrix c_ij

  void bc_coeff(double *, double *, double *, double *);

  // perform bicubic interpolation at point of interest

  void bc_interpol(double, double, int, int, double *, double *, double *, double *);
};
}    // namespace LAMMPS_NS

#endif
#endif
