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

#ifdef ATOM_CLASS
// clang-format off
AtomStyle(tri,AtomVecTri);
// clang-format on
#else

#ifndef LMP_ATOM_VEC_TRI_H
#define LMP_ATOM_VEC_TRI_H

#include "atom_vec.h"

namespace LAMMPS_NS {

class AtomVecTri : public AtomVec {
 public:
  struct Bonus {
    double quat[4];
    double c1[3], c2[3], c3[3];
    double inertia[3];
    int ilocal;
  };
  struct Bonus *bonus;

  AtomVecTri(class LAMMPS *);
  ~AtomVecTri() override;
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

  // unique to AtomVecTri

  void set_equilateral(int, double);

  int nlocal_bonus;

 private:
  int *tri;
  double *radius, *rmass;
  double **omega, **angmom;

  int nghost_bonus, nmax_bonus;
  int tri_flag;
  double rmass_one;

  void grow_bonus();
  void copy_bonus_all(int, int);
};

}    // namespace LAMMPS_NS

#endif
#endif
