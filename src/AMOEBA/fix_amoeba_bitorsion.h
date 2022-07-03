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
FixStyle(amoeba/bitorsion,FixAmoebaBiTorsion);
// clang-format on
#else

#ifndef LMP_FIX_AMOEBA_BITORSION_H
#define LMP_FIX_AMOEBA_BITORSION_H

#include "fix.h"
namespace LAMMPS_NS {

class FixAmoebaBiTorsion : public Fix {
 public:
  FixAmoebaBiTorsion(class LAMMPS *, int, char **);
  ~FixAmoebaBiTorsion();
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
  int pack_border(int, int *, double *) override;
  int unpack_border(int, int, double *) override;
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);

  double memory_usage();

 private:
  int nprocs, me;
  int eflag_caller;
  int ilevel_respa;
  int disable;
  bigint nbitorsions;      // total count of all bitorsions in system
  double ebitorsion;
  double onefifth;

  // per-atom data for bitorsions stored with each owned atom

  int *num_bitorsion;
  int **bitorsion_type;
  tagint **bitorsion_atom1, **bitorsion_atom2, **bitorsion_atom3;
  tagint **bitorsion_atom4, **bitorsion_atom5;

  // previous max atoms on this proc before grow() is called

  int nmax_previous;

  // list of all bitorsions to compute on this proc

  int nbitorsion_list;
  int max_bitorsion_list;
  int **bitorsion_list;

  // BiTorsion grid and spline data

  int nbitypes;
  int *nxgrid, *nygrid;
  double **ttx, **tty, **tbf;
  double **tbx, **tby, **tbxy;

  // data from PairAmoeba

  class Pair *pair;
  int *amtype, *atomic_num;

  // local methods

  void read_grid_data(char *);
  void create_splines();
  void nspline(int, double *, double *, double *, double *, double *, double *, double *, double *,
               double *);
  void cspline(int, double *, double *, double *, double *, double *, double *, double *, double *,
               double *, double *);
  void cytsy(int, double *, double *, double *, double *, double *, int &);
  void cytsyp(int, double *, double *, double *, int &);
  void cytsys(int, double *, double *, double *, double *, double *);

  void chkttor(int, int, int, double &, double &, double &);
  void bcuint1(double *, double *, double *, double *, double, double, double, double, double,
               double, double &, double &, double &);
  void bcucof(double *, double *, double *, double *, double, double, double[][4]);
};

}    // namespace LAMMPS_NS

#endif
#endif
