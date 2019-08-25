/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#ifdef DUMP_CLASS

DumpStyle(molfile,DumpMolfile)

#else

#ifndef LMP_DUMP_MOLFILE_H
#define LMP_DUMP_MOLFILE_H

#include "dump.h"

namespace LAMMPS_NS {

class DumpMolfile : public Dump {
 public:
  DumpMolfile(LAMMPS *, int, char**);
  virtual ~DumpMolfile();
  virtual void write();

 protected:
  class MolfileInterface *mf; //< handles low-level I/O
  // per-atom data
  float *coords, *vels, *masses, *charges, *radiuses;
  int   *types, *molids;
  char **typenames;

  int natoms,me,ntotal,ntypes;
  int need_structure;
  int unwrap_flag;   // 1 if writing unwrapped atom coords, 0 if not
  int velocity_flag; // 1 if writing velocities, 0 if not
  int topology_flag; // 1 if writing topology data, 0 if not
  float cell[6];     // cell parameters: A, B, C, alpha, beta, gamma

  virtual void init_style();
  virtual int modify_param(int, char **);
  virtual void write_header(bigint) {};
  virtual void pack(tagint *);
  virtual void write_data(int, double *);
  virtual bigint memory_usage();
  virtual void openfile();
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

The format of Molfile dump files requires snapshots be output
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
