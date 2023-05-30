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

#ifdef ATOM_CLASS
// clang-format off
AtomStyle(line,AtomVecLine);
// clang-format on
#else

#ifndef LMP_ATOM_VEC_LINE_H
#define LMP_ATOM_VEC_LINE_H

#include "atom_vec.h"

namespace LAMMPS_NS {

class AtomVecLine : public AtomVec {
 public:
  struct Bonus {
    double length, theta;
    int ilocal;
  };
  struct Bonus *bonus;

  AtomVecLine(class LAMMPS *);
  ~AtomVecLine() override;
  void init() override;

  void grow_pointers() override;
  void copy_bonus(int, int, int) override;
  void clear_bonus() override;
  int pack_comm_bonus(int, int *, double *) override;
  void unpack_comm_bonus(int, int, double *) override;
  int pack_border_bonus(int, int *, double *) override;
  int unpack_border_bonus(int, int, double *) override;
  int pack_exchange_bonus(int, double *) override;
  int unpack_exchange_bonus(int, double *) override;
  int size_restart_bonus() override;
  int pack_restart_bonus(int, double *) override;
  int unpack_restart_bonus(int, double *) override;
  void data_atom_bonus(int, const std::vector<std::string> &) override;
  double memory_usage_bonus() override;

  void create_atom_post(int) override;
  void data_atom_post(int) override;
  void pack_data_pre(int) override;
  void pack_data_post(int) override;

  int pack_data_bonus(double *, int) override;
  void write_data_bonus(FILE *, int, double *, int) override;

  // unique to AtomVecLine

  void set_length(int, double);

  int nlocal_bonus;

 private:
  int *line;
  double *radius, *rmass;
  double **omega;

  int nghost_bonus, nmax_bonus;
  int line_flag;
  double rmass_one;

  void grow_bonus();
  void copy_bonus_all(int, int);
  // void consistency_check(int, char *);
};

}    // namespace LAMMPS_NS

#endif
#endif
