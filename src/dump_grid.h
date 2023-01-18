/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef DUMP_CLASS
// clang-format off
DumpStyle(grid,DumpGrid);
// clang-format on
#else

#ifndef LMP_DUMP_GRID_H
#define LMP_DUMP_GRID_H

#include "dump.h"

namespace LAMMPS_NS {

class DumpGrid : public Dump {
 public:
  DumpGrid(class LAMMPS *, int, char **);
  ~DumpGrid() override;

  const char *MAGIC_STRING = "DUMPGRID";
  const int FORMAT_REVISION = 0x0002;
  const int ENDIAN = 0x0001;

 protected:
  int nevery;        // dump frequency for output
  char *idregion;    // region ID, nullptr if no region

  int expand;        // flag for whether field args were expanded
  char **earg;       // field names with wildcard expansion
  int nargnew;       // size of earg
                     //
  int *vtype;        // type of each vector (INT, DOUBLE)
  char **vformat;    // format string for each vector element
                     //
  char *columns;     // column labels
  char *columns_default;

  int dimension;

  int nxgrid, nygrid, nzgrid;    // global grid size

  int nfield;       // # of keywords listed by user
  int ioptional;    // index of start of optional args

  // per field info
  int *field2index;     // which compute/fix
  int *field2source;    // COMPUTE or FIX
  int *field2grid;      // index of grid within compute/fix
  int *field2data;      // index of data within compute/fix
  int *argindex;        // index into compute,fix,custom per-atom data
                        // 0 for per-atom vector, 1-N for cols of per-atom array

  int ncompute;                            // # of Computes accessed by dump
  char **id_compute;                       // their IDs
  std::vector<class Compute *> compute;    // list of ptrs to the Computes

  int nfix;                        // # of Fixes used by dump
  char **id_fix;                   // their IDs
  std::vector<class Fix *> fix;    // list of ptrs to the Fixes

  int nxlo_in, nxhi_in;    // bounds of this proc's portion of grids
  int nylo_in, nyhi_in;
  int nzlo_in, nzhi_in;

  // methods

  void init_style() override;
  void write_header(bigint) override;
  int count() override;
  void pack(tagint *) override;
  int convert_string(int, double *) override;
  void write_data(int, double *) override;
  double memory_usage() override;

  int parse_fields(int, char **);
  int add_compute(const std::string &, class Compute *);
  int add_fix(const std::string &, class Fix *);
  int modify_param(int, char **) override;

  void header_format_binary();
  void header_unit_style_binary();
  void header_time_binary();
  void header_columns_binary();

  typedef void (DumpGrid::*FnPtrHeader)(bigint);
  FnPtrHeader header_choice;    // ptr to write header functions
  void header_binary(bigint);
  void header_binary_triclinic(bigint);
  void header_item(bigint);
  void header_item_triclinic(bigint);

  void format_magic_string_binary();
  void format_endian_binary();
  void format_revision_binary();

  typedef void (DumpGrid::*FnPtrWrite)(int, double *);
  FnPtrWrite write_choice;    // ptr to write data functions
  void write_binary(int, double *);
  void write_string(int, double *);
  void write_lines(int, double *);

  // customize by adding a method prototype

  typedef void (DumpGrid::*FnPtrPack)(int);
  FnPtrPack *pack_choice;    // ptrs to pack functions

  void pack_grid2d(int);
  void pack_grid3d(int);
};

}    // namespace LAMMPS_NS

#endif
#endif
