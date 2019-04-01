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
/* ----------------------------------------------------------------------
   Electronic stopping power
   Contributing authors: K. Avchaciov and T. Metspalu
   Information: k.avchachov@gmail.com
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(elstop,FixElstop)

#else

#ifndef LMP_FIX_ELSTOP_H
#define LMP_FIX_ELSTOP_H

#include "fix.h"


namespace LAMMPS_NS {

class FixElstop : public Fix {
 public:
  FixElstop(class LAMMPS *, int, char **);
  ~FixElstop();
  int setmask();
  void init();
  void post_force(int);
  void init_list(int, class NeighList *);
  double compute_scalar();

 private:
  void read_table(const char *);
  void grow_table();

  double Ecut;               // cutoff energy
  double SeLoss, SeLoss_all; // electronic energy loss
  int SeLoss_sync_flag;      // sync done since last change?

  int maxlines;           // max number of lines in table
  int table_entries;      // number of table entries actually read
  double **elstop_ranges; // [ 0][i]: energies
                          // [>0][i]: stopping powers per type

  int iregion; // region index if used, else -1
  int minneigh; // minimum number of neighbors

  class NeighList *list;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Region ID for fix elstop does not exist

Self-explanatory.

E: Atom kinetic energy too high for fix elstop

The group given in the fix elstop command includes an atom that has
a kinetic energy higher than the largest energy in the elstop table.
Reconsider whether the table is physically applicable to your system.

E: Cannot open stopping range table ...

The file containing the elstop table could not be opened. Chck the
given path and the file's permissions.

E: fix elstop: Invalid table line

A line in the elstop table file contained too many or too few columns.

E: fix elstop: Energies must be in ascending order

The first column in the elstop table must be sorted from the smallest
energy to the largest.

E: Did not find any data in elstop table file

Parsing the elstop table file produced no lines that were identifiable
as energies/stopping powers. Most likely the file is empty or contains
only comments.

*/
