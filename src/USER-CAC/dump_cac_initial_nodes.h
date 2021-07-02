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

DumpStyle(cac/initial/nodes,DumpCACInitialNodes)

#else

#ifndef LMP_DUMP_CACINITIALNODES_H
#define LMP_DUMP_CACINITIALNODES_H

#include "dump.h"

namespace LAMMPS_NS {

class DumpCACInitialNodes : public Dump {
 public:
  DumpCACInitialNodes(class LAMMPS *, int, char**);
  virtual ~DumpCACInitialNodes();

 protected:
  int ntypes;
  char **typenames;
  bigint ptimestep;
  int total_node_count, total_element_count;
  void init_style();
  void write_header(bigint);
  void pack(tagint *);
  int convert_string(int, double *);
  int nodes_per_element,maxpoly;
  void write_data(int, double *);
  int count();
  int modify_param(int, char **);

  typedef void (DumpCACInitialNodes::*FnPtrWrite)(int, double *);
  FnPtrWrite write_choice;              // ptr to write data functions
  void write_string(int, double *);
  void write_lines(int, double *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Invalid dump xyz filename

Filenames used with the dump xyz style cannot be binary or cause files
to be written by each processor.

E: CAC dump styles require a CAC atom style

Self-explanatory

E: Dump modify element names do not match atom types

Number of element names must equal number of atom types.

*/
