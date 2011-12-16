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

#ifdef DUMP_CLASS

DumpStyle(atom,DumpAtom)

#else

#ifndef LMP_DUMP_ATOM_H
#define LMP_DUMP_ATOM_H

#include "dump.h"

namespace LAMMPS_NS {

class DumpAtom : public Dump {
 public:
  DumpAtom(LAMMPS *, int, char**);

 private:
  int scale_flag;            // 1 if atom coords are scaled, 0 if no
  int image_flag;            // 1 if append box count to atom coords, 0 if no

  char *columns;             // column labels

  void init_style();
  int modify_param(int, char **);
  void write_header(bigint);
  int count();
  void pack(int *);
  void write_data(int, double *);

  typedef void (DumpAtom::*FnPtrHeader)(bigint);
  FnPtrHeader header_choice;           // ptr to write header functions
  void header_binary(bigint);
  void header_binary_triclinic(bigint);
  void header_item(bigint);
  void header_item_triclinic(bigint);

  typedef void (DumpAtom::*FnPtrPack)(int *);
  FnPtrPack pack_choice;               // ptr to pack functions
  void pack_scale_image(int *);
  void pack_scale_noimage(int *);
  void pack_noscale_image(int *);
  void pack_noscale_noimage(int *);
  void pack_scale_image_triclinic(int *);
  void pack_scale_noimage_triclinic(int *);

  typedef void (DumpAtom::*FnPtrData)(int, double *);
  FnPtrData write_choice;              // ptr to write data functions
  void write_binary(int, double *);
  void write_image(int, double *);
  void write_noimage(int, double *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
