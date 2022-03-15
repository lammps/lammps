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
FixStyle(amoeba/pitorsion,FixAmoebaPiTorsion);
// clang-format on
#else

#ifndef LMP_FIX_AMOEBA_PITORSION_H
#define LMP_FIX_AMOEBA_PITORSION_H

#include "fix.h"
namespace LAMMPS_NS {

class FixAmoebaPiTorsion : public Fix {
 public:
  FixAmoebaPiTorsion(class LAMMPS *, int, char **);
  ~FixAmoebaPiTorsion();
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
  int eflag_caller;
  int ilevel_respa;
  int disable;
  bigint npitorsions;
  int npitorsion_types;
  double epitorsion;
  double onesixth;

  double *kpit;

  // per-atom data for pitorsions stored with each owned atom

  int *num_pitorsion;
  int **pitorsion_type;
  tagint **pitorsion_atom1, **pitorsion_atom2, **pitorsion_atom3;
  tagint **pitorsion_atom4, **pitorsion_atom5, **pitorsion_atom6;

  // previous max atoms on this proc before grow() is called

  int nmax_previous;

  // list of all pitorsions to compute on this proc

  int npitorsion_list;
  int max_pitorsion_list;
  int **pitorsion_list;
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
