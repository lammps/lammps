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

#ifdef FIX_CLASS
// clang-format off
FixStyle(property/atom,FixPropertyAtom);
// clang-format on
#else

#ifndef LMP_FIX_PROPERTY_ATOM_H
#define LMP_FIX_PROPERTY_ATOM_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPropertyAtom : public Fix {
 public:
  FixPropertyAtom(class LAMMPS *, int, char **);
  virtual ~FixPropertyAtom();
  int setmask();
  void init();

  void read_data_section(char *, int, char *, tagint);
  bigint read_data_skip_lines(char *);
  void write_data_section_size(int, int &, int &);
  void write_data_section_pack(int, double **);
  void write_data_section_keyword(int, FILE *);
  void write_data_section(int, FILE *, int, double **, int);

  virtual void grow_arrays(int);
  void copy_arrays(int, int, int);
  int pack_border(int, int *, double *);
  int unpack_border(int, int, double *);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  int pack_restart(int, double *);
  void unpack_restart(int, int);
  int size_restart(int);
  int maxsize_restart();
  double memory_usage();

 protected:
  int nvalue, border;
  int molecule_flag, q_flag, rmass_flag;    // flags for specific fields
  int *styles;                              // style of each value, see enum
  int *index;                               // indices into atom custom data structs
  int *cols;                                // columns per value, for arrays
  char *astyle;                             // atom style at instantiation

  int values_peratom;    // # of values per atom, including multiple for arrays
  int nmax_old;          // length of peratom arrays the last time they grew
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix property/atom mol when atom_style already has molecule attribute

Self-explanatory.

E: Fix property/atom cannot specify mol twice

Self-explanatory.

E: Fix property/atom q when atom_style already has charge attribute

Self-explanatory.

E: Fix property/atom cannot specify q twice

Self-explanatory.

E: Fix property/atom rmass when atom_style already has rmass attribute

UNDOCUMENTED

E: Fix property/atom cannot specify rmass twice

UNDOCUMENTED

E: Fix property/atom vector name already exists

The name for an integer or floating-point vector must be unique.

W: Fix property/atom mol or charge or rmass w/out ghost communication

UNDOCUMENTED

E: Atom style was redefined after using fix property/atom

This is not allowed.

E: Incorrect %s format in data file

A section of the data file being read by fix property/atom does
not have the correct number of values per line.

E: Too few lines in %s section of data file

Self-explanatory.

E: Invalid atom ID in %s section of data file

An atom in a section of the data file being read by fix property/atom
has an invalid atom ID that is <= 0 or > the maximum existing atom ID.

U: Fix property/atom mol or charge w/out ghost communication

A model typically needs these properties defined for ghost atoms.

*/
