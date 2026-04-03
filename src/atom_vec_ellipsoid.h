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
AtomStyle(ellipsoid,AtomVecEllipsoid);
// clang-format on
#else

#ifndef LMP_ATOM_VEC_ELLIPSOID_H
#define LMP_ATOM_VEC_ELLIPSOID_H

#include "atom_vec.h"

namespace LAMMPS_NS {

class AtomVecEllipsoid : virtual public AtomVec {
 public:
  enum BlockType {
    ELLIPSOID = 0, // n1 = n2 = 2
    N1_EQUAL_N2 = 1, // n1 = n2 > 2
    GENERAL = 2, // n1 != n2, n1 > 2, n2 > 2
  };
  struct Bonus {
    double shape[3];
    double quat[4];
    int ilocal;
  };
  struct Bonus *bonus;

  struct BonusSuper : public Bonus {
    double block[2];
    double inertia[3];
    BlockType type;
  };
  struct BonusSuper *bonus_super;

  AtomVecEllipsoid(class LAMMPS *);
  ~AtomVecEllipsoid() override;

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

  void read_data_general_to_restricted(int, int) override;
  void write_data_restricted_to_general() override;
  void write_data_restore_restricted() override;

  // unique to AtomVecEllipsoid

  void set_shape(int, double, double, double);
  void set_block(int, double, double);

  int nlocal_bonus;

 protected:
  int *ellipsoid;
  double *radius, *rmass;
  double **angmom;
  double **quat_hold;

  int nghost_bonus, nmax_bonus;
  int ellipsoid_flag;
  double rmass_one;

  virtual void grow_bonus();
  void copy_bonus_all(int, int);

  static BlockType determine_type(double *);
  static double radius_ellipsoid(double *, double *, BlockType);
  static void inertia_ellipsoid_principal(double *, double, double *,
                                   double *block, BlockType);


  template <bool is_super>
  int pack_comm_bonus_templated(int, int *, double *);

  template <bool is_super>
  void unpack_comm_bonus_templated(int, int, double *);

  template <bool is_super>
  int pack_border_bonus_templated(int, int *, double *);

  template <bool is_super>
  int unpack_border_bonus_templated(int, int, double *);

  void process_args(int, char **) override;


};

}    // namespace LAMMPS_NS

#endif
#endif
