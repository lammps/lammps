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
  ~FixAmoebaPiTorsion() override;
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
