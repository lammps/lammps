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

#ifdef DUMP_CLASS

DumpStyle(dcd,DumpDCD)

#else

#ifndef LMP_DUMP_DCD_H
#define LMP_DUMP_DCD_H

#include "dump.h"

namespace LAMMPS_NS {

class DumpDCD : public Dump {
 public:
  DumpDCD(LAMMPS *, int, char**);
  ~DumpDCD();

 private:
  int natoms,ntotal;
  int headerflag,nevery_save,nframes;

  float *coords,*xf,*yf,*zf;
  int unwrap_flag;            // 1 if atom coords are unwrapped, 0 if no

  void init_style();
  void openfile();
  void write_header(bigint);
  void pack(tagint *);
  void write_data(int, double *);
  int modify_param(int, char **);
  bigint memory_usage();

  void write_frame();
  void write_dcd_header(const char *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Invalid dump dcd filename

Filenames used with the dump dcd style cannot be binary or compressed
or cause multiple files to be written.

E: Too many atoms for dump dcd

The system size must fit in a 32-bit integer to use this dump
style.

E: Dump dcd requires sorting by atom ID

Use the dump_modify sort command to enable this.

E: Cannot use variable every setting for dump dcd

The format of DCD dump files requires snapshots be output
at a constant frequency.

E: Cannot change dump_modify every for dump dcd

The frequency of writing dump dcd snapshots cannot be changed.

E: Cannot open dump file

The output file for the dump command cannot be opened.  Check that the
path and name are correct.

E: Dump dcd of non-matching # of atoms

Every snapshot written by dump dcd must contain the same # of atoms.

E: Too big a timestep for dump dcd

The timestep must fit in a 32-bit integer to use this dump style.

*/
