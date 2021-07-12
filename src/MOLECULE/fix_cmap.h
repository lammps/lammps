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
  ~FixCMAP();
  int setmask();
  void init();
  void setup(int);
  void setup_pre_neighbor();
  void setup_pre_reverse(int, int);
  void min_setup(int);
  void pre_neighbor();
  void pre_reverse(int, int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);
  double compute_scalar();

  void read_data_header(char *);
  void read_data_section(char *, int, char *, tagint);
  bigint read_data_skip_lines(char *);
  void write_data_header(FILE *, int);
  void write_data_section_size(int, int &, int &);
  void write_data_section_pack(int, double **);
  void write_data_section_keyword(int, FILE *);
  void write_data_section(int, FILE *, int, double **, int);

  void write_restart(FILE *);
  void restart(char *);
  int pack_restart(int, double *);
  void unpack_restart(int, int);
  int size_restart(int);
  int maxsize_restart();

  void grow_arrays(int);
  void copy_arrays(int, int, int);
  void set_arrays(int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);

  double memory_usage();

 private:
  int nprocs, me;
  int newton_bond, eflag_caller;
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

/* ERROR/WARNING messages:

E: Illegal ... command

UNDOCUMENTED

E: CMAP atoms %d %d %d %d %d missing on proc %d at step %ld

UNDOCUMENTED

E: Invalid CMAP crossterm_type

UNDOCUMENTED

E: Cannot open fix cmap file %s

UNDOCUMENTED

E: CMAP: atan2 function cannot take 2 zero arguments

UNDOCUMENTED

E: Invalid read data header line for fix cmap

UNDOCUMENTED

E: Incorrect %s format in data file

UNDOCUMENTED

E: Too many CMAP crossterms for one atom

UNDOCUMENTED

*/
