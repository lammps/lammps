/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef DUMP_CLASS
// clang-format off
DumpStyle(atom,DumpAtom);
// clang-format on
#else

#ifndef LMP_DUMP_ATOM_H
#define LMP_DUMP_ATOM_H

#include "dump.h"

namespace LAMMPS_NS {

class DumpAtom : public Dump {
 public:
  DumpAtom(LAMMPS *, int, char **);

  const char *MAGIC_STRING = "DUMPATOM";
  const int FORMAT_REVISION = 0x0002;
  const int ENDIAN = 0x0001;

 protected:
  int scale_flag;        // 1 if atom coords are scaled, 0 if no
  int image_flag;        // 1 if append box count to atom coords, 0 if no
  int triclinic_general; // 1 if output box & coords for general triclinic, 0 if no

  std::string columns;    // column labels

  void init_style() override;
  int modify_param(int, char **) override;
  void write_header(bigint) override;
  void pack(tagint *) override;
  int convert_string(int, double *) override;
  void write_data(int, double *) override;

  void header_format_binary();
  void header_unit_style_binary();
  void header_time_binary();
  void header_columns_binary();
  void format_magic_string_binary();
  void format_endian_binary();
  void format_revision_binary();

  typedef void (DumpAtom::*FnPtrHeader)(bigint);
  FnPtrHeader header_choice;    // ptr to write header functions
  void header_binary(bigint);
  void header_binary_triclinic(bigint);
  void header_binary_triclinic_general(bigint);
  void header_item(bigint);
  void header_item_triclinic(bigint);
  void header_item_triclinic_general(bigint);

  typedef void (DumpAtom::*FnPtrPack)(tagint *);
  FnPtrPack pack_choice;    // ptr to pack functions
  void pack_scale_image(tagint *);
  void pack_scale_noimage(tagint *);
  void pack_noscale_image(tagint *);
  void pack_noscale_noimage(tagint *);
  void pack_scale_image_triclinic(tagint *);
  void pack_scale_noimage_triclinic(tagint *);
  void pack_noscale_image_triclinic_general(tagint *);
  void pack_noscale_noimage_triclinic_general(tagint *);

  typedef int (DumpAtom::*FnPtrConvert)(int, double *);
  FnPtrConvert convert_choice;    // ptr to convert data functions
  int convert_image(int, double *);
  int convert_noimage(int, double *);

  typedef void (DumpAtom::*FnPtrWrite)(int, double *);
  FnPtrWrite write_choice;    // ptr to write data functions
  void write_binary(int, double *);
  void write_string(int, double *);
  void write_lines_image(int, double *);
  void write_lines_noimage(int, double *);
};

}    // namespace LAMMPS_NS

#endif
#endif
