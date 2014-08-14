/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   Contributing author: Pierre de Buyl (KU Leuven)
------------------------------------------------------------------------- */

#ifdef DUMP_CLASS

DumpStyle(h5md,DumpH5MD)

#else

#ifndef LMP_DUMP_H5MD_H
#define LMP_DUMP_H5MD_H

#include "dump.h"
#include "ch5md.h"

namespace LAMMPS_NS {

class DumpH5MD : public Dump {
 public:
  DumpH5MD(class LAMMPS *, int, char**);
  virtual ~DumpH5MD();

 private:
  int natoms,ntotal;
  int nevery_save;
  int unwrap_flag;            // 1 if atom coords are unwrapped, 0 if no
  double *coords;
  h5md_file datafile;
  h5md_particles_group particles_data;

  void init_style();
  int modify_param(int, char **);
  void openfile();
  void write_header(bigint);
  void pack(tagint *);
  void write_data(int, double *);

  void write_frame();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Dump h5md requires sorting by atom ID

Use the dump_modify sort command to enable this.

E: Cannot use variable every setting for dump xtc

The format of this file requires snapshots at regular intervals.

E: Cannot change dump_modify every for dump xtc

The frequency of writing dump xtc snapshots cannot be changed.

*/
