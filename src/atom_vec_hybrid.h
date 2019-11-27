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
  int property_atom(char *);
  void pack_property_atom(int, double *, int, int);

 private:
  int nallstyles;
  char **allstyles;

  void concatenate_fields();
  void concatenate(char *&, char *);
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
