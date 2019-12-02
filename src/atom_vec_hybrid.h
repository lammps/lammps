/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef ATOM_CLASS

AtomStyle(hybrid,AtomVecHybrid)

#else

#ifndef LMP_ATOM_VEC_HYBRID_H
#define LMP_ATOM_VEC_HYBRID_H

#include "atom_vec.h"

namespace LAMMPS_NS {

class AtomVecHybrid : public AtomVec {
 public:
  int nstyles;
  class AtomVec **styles;
  char **keywords;

  AtomVecHybrid(class LAMMPS *);
  ~AtomVecHybrid();
  void process_args(int, char **);
  void init();
  void force_clear(int, size_t);

  void copy_bonus(int, int, int);
  void clear_bonus() {}
  int pack_comm_bonus(int, int *, double *);
  void unpack_comm_bonus(int, int, double *);
  int pack_border_bonus(int, int *, double *);
  int unpack_border_bonus(int, int, double *);
  int pack_exchange_bonus(int, double *);
  int unpack_exchange_bonus(int, double *);
  int size_restart_bonus();
  int pack_restart_bonus(int, double *);
  int unpack_restart_bonus(int, double *);
  bigint memory_usage_bonus();

  void pack_restart_pre(int);
  void pack_restart_post(int);
  void unpack_restart_init(int);
  void create_atom_post(int);
  void data_atom_post(int);
  void pack_data_pre(int);
  void pack_data_post(int);

  //void create_atom_post(int);
  //void data_atom_post(int);
  //void pack_data_pre(int);
  //void pack_data_post(int);

  int property_atom(char *);
  void pack_property_atom(int, double *, int, int);

 private:
  int nallstyles;
  char **allstyles;

  struct FieldStrings {
    char **fstr;
  };
  FieldStrings *fieldstrings;

  int nstyles_bonus;
  class AtomVec **styles_bonus;

  char *merge_fields(int, char *, int, char *&);
  void build_styles();
  int known_style(char *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Atom style hybrid cannot have hybrid as an argument

Self-explanatory.

E: Atom style hybrid cannot use same atom style twice

Self-explanatory.

E: Cannot mix molecular and molecule template atom styles

Self-explanatory.

E: Per-processor system is too big

The number of owned atoms plus ghost atoms on a single
processor must fit in 32-bit integer.

E: Invalid atom type in Atoms section of data file

Atom types must range from 1 to specified # of types.

*/
