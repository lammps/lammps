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

#ifdef FIX_CLASS

FixStyle(property/atom,FixPropertyAtom)

#else

#ifndef LMP_FIX_PROPERTY_ATOM_H
#define LMP_FIX_PROPERTY_ATOM_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPropertyAtom : public Fix {
 public:
  FixPropertyAtom(class LAMMPS *, int, char **);
  ~FixPropertyAtom();
  int setmask();
  void init();

  void read_data_section(char *, int, char *);
  bigint read_data_skip_lines(char *);
  void write_data_section_size(int, int &, int &);
  void write_data_section_pack(int, double **);
  void write_data_section_keyword(int, FILE *);
  void write_data_section(int, FILE *, int, double **, int);

  void grow_arrays(int);
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

 private:
  int nvalue,border;
  int molecule_flag,q_flag;
  int *style,*index;
  char *astyle;

  int nmax_old;         // length of peratom arrays the last time they grew
};

}

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

E: Fix property/atom vector name already exists

The name for an integer or floating-point vector must be unique.

E: Atom style was redefined after using fix property/atom

This is not allowed.

E: Incorrect %s format in data file

A section of the data file being read by fix property/atom does
not have the correct number of values per line.

E: Invalid atom ID in %s section of data file

An atom in a section of the data file being read by fix property/atom
has an invalid atom ID that is <= 0 or > the maximum existing atom ID.

*/
