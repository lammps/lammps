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

#ifndef DUMP_ATOM_H
#define DUMP_ATOM_H

#include "dump.h"

namespace LAMMPS_NS {

class DumpAtom : public Dump {
 public:
  DumpAtom(LAMMPS *, int, char**);
  void init();

 private:
  int scale_flag;            // 1 if atom coords are scaled, 0 if no
  int image_flag;            // 1 if append box count to atom coords, 0 if no

  char *columns;             // column labels

  int modify_param(int, char **);
  void write_header(int);
  int count();
  int pack();
  void write_data(int, double *);

  typedef void (DumpAtom::*FnPtrHeader)(int);
  FnPtrHeader header_choice;           // ptr to write header functions
  void header_binary(int);
  void header_binary_triclinic(int);
  void header_item(int);
  void header_item_triclinic(int);

  typedef int (DumpAtom::*FnPtrPack)();
  FnPtrPack pack_choice;               // ptr to pack functions
  int pack_scale_image();
  int pack_scale_noimage();
  int pack_noscale_image();
  int pack_noscale_noimage();
  int pack_scale_image_triclinic();
  int pack_scale_noimage_triclinic();

  typedef void (DumpAtom::*FnPtrData)(int, double *);
  FnPtrData write_choice;              // ptr to write data functions
  void write_binary(int, double *);
  void write_image(int, double *);
  void write_noimage(int, double *);
};

}

#endif
