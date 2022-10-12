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
AtomStyle(body,AtomVecBody);
// clang-format on
#else

#ifndef LMP_ATOM_VEC_BODY_H
#define LMP_ATOM_VEC_BODY_H

#include "atom_vec.h"

namespace LAMMPS_NS {

class AtomVecBody : public AtomVec {
 public:
  class Body *bptr;

  struct Bonus {
    double quat[4];
    double inertia[3];
    int ninteger, ndouble;
    int iindex, dindex;
    int *ivalue;
    double *dvalue;
    int ilocal;
  };
  struct Bonus *bonus;

  AtomVecBody(class LAMMPS *);
  ~AtomVecBody() override;
  void process_args(int, char **) override;

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
  void data_body(int, int, int, int *, double *) override;
  double memory_usage_bonus() override;

  void create_atom_post(int) override;
  void data_atom_post(int) override;
  void pack_data_pre(int) override;
  void pack_data_post(int) override;

  int pack_data_bonus(double *, int) override;
  void write_data_bonus(FILE *, int, double *, int) override;

  // methods used by other classes to query/set body info

  double radius_body(int, int, int *, double *);
  void set_quat(int, double *);

  int nlocal_bonus;

 private:
  int *body;
  double *rmass, *radius;
  double **angmom;

  int nghost_bonus, nmax_bonus;
  int intdoubleratio;    // sizeof(double) / sizeof(int)
  int body_flag;

  MyPoolChunk<int> *icp;
  MyPoolChunk<double> *dcp;

  void grow_bonus();
  void copy_bonus_all(int, int);
  // check(int);
};

}    // namespace LAMMPS_NS

#endif
#endif
