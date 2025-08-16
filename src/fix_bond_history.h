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
FixStyle(BOND_HISTORY,FixBondHistory);
// clang-format on
#else

#ifndef LMP_FIX_BOND_HISTORY_H
#define LMP_FIX_BOND_HISTORY_H

#include "fix.h"

#include <map>
#include <utility>

namespace LAMMPS_NS {

class FixBondHistory : public Fix {
 public:
  FixBondHistory(class LAMMPS *, int, char **);
  ~FixBondHistory() override;
  int setmask() override;
  void post_constructor() override;
  void setup_post_neighbor() override;
  void setup_pre_exchange() override;
  void post_neighbor() override;
  void pre_exchange() override;
  double memory_usage() override;
  void write_restart(FILE *fp) override;
  void restart(char *buf) override;
  void set_arrays(int) override;

  void update_atom_value(int, int, int, double);
  double get_atom_value(int, int, int);
  int get_ndata() const { return ndata; }

  // methods to reorder/delete elements of atom->bond_atom
  void delete_history(int, int);
  void shift_history(int, int, int);
  void cache_history(int, int);
  void check_cache(int, int);
  void clear_cache();

  // methods for bond style hybrid
  void compress_history();
  void uncompress_history();

  // if data is temporarily stored while the bond_atom array
  // is being reordered, use map of vectors with pairs for keys
  // to enable quick look up
  std::map<std::pair<tagint, tagint>, std::vector<double>> cached_histories;

  int *setflag;    // Set by BondBPM, which bond types are used
  double **bondstore;
  int stored_flag;
  int ndata;

 protected:
  void allocate();

  int hybrid_flag;
  int nbondlist_orig;
  int *bondtype_orig;
  double **bondstore_comp;
  double **bondstore_orig;

  int update_flag;    // Flag whether history values can evolve
  int updated_bond_flag;
  int nbond, maxbond;
  int index;
  char *id_fix;
  char *id_array;
};

}    // namespace LAMMPS_NS

#endif
#endif
