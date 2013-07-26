/* ----------------------------------------------------------------------
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
  int nvalue,border,molecule_flag;
  int *style,*index;
  int nmax_old;         // length of peratom arrays the last time they grew
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
