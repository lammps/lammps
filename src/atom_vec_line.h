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
  ~AtomVecLine();
  void init();

  void grow_pointers();
  void copy_bonus(int, int, int);
  void clear_bonus();
  int pack_comm_bonus(int, int *, double *);
  void unpack_comm_bonus(int, int, double *);
  int pack_border_bonus(int, int *, double *);
  int unpack_border_bonus(int, int, double *);
  int pack_exchange_bonus(int, double *);
  int unpack_exchange_bonus(int, double *);
  int size_restart_bonus();
  int pack_restart_bonus(int, double *);
  int unpack_restart_bonus(int, double *);
  void data_atom_bonus(int, char **);
  double memory_usage_bonus();

  void create_atom_post(int);
  void data_atom_post(int);
  void pack_data_pre(int);
  void pack_data_post(int);

  int pack_data_bonus(double *, int);
  void write_data_bonus(FILE *, int, double *, int);

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

/* ERROR/WARNING messages:

E: Atom_style line can only be used in 2d simulations

Self-explanatory.

E: Per-processor system is too big

The number of owned atoms plus ghost atoms on a single
processor must fit in 32-bit integer.

E: Invalid atom type in Atoms section of data file

Atom types must range from 1 to specified # of types.

E: Invalid density in Atoms section of data file

Density value cannot be <= 0.0.

E: Assigning line parameters to non-line atom

Self-explanatory.

E: Inconsistent line segment in data file

The end points of the line segment are not equal distances from the
center point which is the atom coordinate.

*/
