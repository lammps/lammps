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
AtomStyle(hybrid,AtomVecHybrid);
// clang-format on
#else

#ifndef LMP_ATOM_VEC_HYBRID_H
#define LMP_ATOM_VEC_HYBRID_H

#include "atom_vec.h"

namespace LAMMPS_NS {

class AtomVecHybrid : virtual public AtomVec {
 public:
  int nstyles;
  class AtomVec **styles;
  char **keywords;

  AtomVecHybrid(class LAMMPS *);
  ~AtomVecHybrid() override;
  void process_args(int, char **) override;
  void init() override;

  void grow_pointers() override;
  void force_clear(int, size_t) override;
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
  double memory_usage_bonus() override;

  void pack_restart_pre(int) override;
  void pack_restart_post(int) override;
  void unpack_restart_init(int) override;
  void create_atom_post(int) override;
  void data_atom_post(int) override;

  void data_bonds_post(int, int, tagint, tagint, tagint) override;

  void pack_data_pre(int) override;
  void pack_data_post(int) override;

  int pack_data_bonus(double *, int) override;
  void write_data_bonus(FILE *, int, double *, int) override;

  int property_atom(const std::string &) override;
  void pack_property_atom(int, double *, int, int) override;

 private:
  int nstyles_bonus;
  class AtomVec **styles_bonus;

  void merge_fields(std::vector<std::string> &, const std::vector<std::string> &, int,
                    std::vector<std::string> &);
};

}    // namespace LAMMPS_NS

#endif
#endif
