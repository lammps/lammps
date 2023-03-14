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
DumpStyle(custom,DumpCustom);
// clang-format on
#else

#ifndef LMP_DUMP_CUSTOM_H
#define LMP_DUMP_CUSTOM_H

#include "dump.h"

namespace LAMMPS_NS {

class DumpCustom : public Dump {
 public:
  DumpCustom(class LAMMPS *, int, char **);
  ~DumpCustom() override;

  const char *MAGIC_STRING = "DUMPCUSTOM";
  const int FORMAT_REVISION = 0x0002;
  const int ENDIAN = 0x0001;

 protected:
  int nevery;        // dump frequency for output
  char *idregion;    // region ID, nullptr if no region

  int nthresh;                        // # of defined thresholds
  int nthreshlast;                    // # of defined thresholds with value = LAST
                                      //
  int *thresh_array;                  // array to threshold on for each nthresh
  int *thresh_op;                     // threshold operation for each nthresh
  double *thresh_value;               // threshold value for each nthresh
  int *thresh_last;                   // for threshold value = LAST,
                                      // index into thresh_fix
                                      // -1 if not LAST, value is numeric
                                      //
  class FixStoreAtom **thresh_fix;    // stores values for each threshold LAST
  char **thresh_fixID;                // IDs of thresh_fixes
  int *thresh_first;                  // 1 the first time a FixStore values accessed

  int expand;        // flag for whether field args were expanded
  char **earg;       // field names with wildcard expansion
  int nargnew;       // size of earg
                     //
  int *vtype;        // type of each vector (INT, DOUBLE)
  char **vformat;    // format string for each vector element
                     //
  char *columns;     // column labels
  char *columns_default;
  //
  int nchoose;        // # of selected atoms
  int maxlocal;       // size of atom selection and variable arrays
  int *choose;        // local indices of selected atoms
  double *dchoose;    // value for each atom to threshold against
  int *clist;         // compressed list of indices of selected atoms

  int nfield;                 // # of keywords listed by user
  int ioptional;              // index of start of optional args
                              //
  int *field2index;           // which compute,fix,variable,custom calcs this field
  int *argindex;              // index into compute,fix,custom per-atom data
                              // 0 for per-atom vector, 1-N for cols of per-atom array
                              //
  int ncompute;               // # of Computes accessed by dump
  char **id_compute;          // their IDs
  class Compute **compute;    // list of ptrs to the Computes
                              //
  int nfix;                   // # of Fixes used by dump
  char **id_fix;              // their IDs
  class Fix **fix;            // list of ptrs to the Fixes
                              //
  int nvariable;              // # of Variables used by dump
  char **id_variable;         // their names
  int *variable;              // list of Variable indices in Variable class
  double **vbuf;              // local storage for variable evaluation
                              //
  int ncustom;                // # of Custom atom properties used by dump
  char **id_custom;           // their names
  int *custom;                // list of Custom indices in Atom class
  int *custom_flag;           // list of IVEC,DVEC,IARRAY,DARRAY styles
                              //
  int ntypes;                 // # of atom types
  char **typenames;           // array of element names for each type

  // private methods

  void init_style() override;
  void write_header(bigint) override;
  int count() override;
  void pack(tagint *) override;
  int convert_string(int, double *) override;
  void write_data(int, double *) override;
  double memory_usage() override;

  int parse_fields(int, char **);
  int add_compute(const char *);
  int add_fix(const char *);
  int add_variable(const char *);
  int add_custom(const char *, int);
  int modify_param(int, char **) override;

  void header_format_binary();
  void header_unit_style_binary();
  void header_time_binary();
  void header_columns_binary();
  void format_magic_string_binary();
  void format_endian_binary();
  void format_revision_binary();

  typedef void (DumpCustom::*FnPtrHeader)(bigint);
  FnPtrHeader header_choice;    // ptr to write header functions
  void header_binary(bigint);
  void header_binary_triclinic(bigint);
  void header_item(bigint);
  void header_item_triclinic(bigint);

  typedef void (DumpCustom::*FnPtrWrite)(int, double *);
  FnPtrWrite write_choice;    // ptr to write data functions
  void write_binary(int, double *);
  void write_string(int, double *);
  void write_lines(int, double *);

  // customize by adding a method prototype

  typedef void (DumpCustom::*FnPtrPack)(int);
  FnPtrPack *pack_choice;    // ptrs to pack functions

  void pack_compute(int);
  void pack_fix(int);
  void pack_variable(int);
  void pack_custom(int);

  void pack_id(int);
  void pack_molecule(int);
  void pack_proc(int);
  void pack_procp1(int);
  void pack_type(int);
  void pack_mass(int);

  void pack_x(int);
  void pack_y(int);
  void pack_z(int);
  void pack_xs(int);
  void pack_ys(int);
  void pack_zs(int);
  void pack_xs_triclinic(int);
  void pack_ys_triclinic(int);
  void pack_zs_triclinic(int);
  void pack_xu(int);
  void pack_yu(int);
  void pack_zu(int);
  void pack_xu_triclinic(int);
  void pack_yu_triclinic(int);
  void pack_zu_triclinic(int);
  void pack_xsu(int);
  void pack_ysu(int);
  void pack_zsu(int);
  void pack_xsu_triclinic(int);
  void pack_ysu_triclinic(int);
  void pack_zsu_triclinic(int);
  void pack_ix(int);
  void pack_iy(int);
  void pack_iz(int);

  void pack_vx(int);
  void pack_vy(int);
  void pack_vz(int);
  void pack_fx(int);
  void pack_fy(int);
  void pack_fz(int);
  void pack_q(int);
  void pack_mux(int);
  void pack_muy(int);
  void pack_muz(int);
  void pack_mu(int);
  void pack_radius(int);
  void pack_diameter(int);

  void pack_omegax(int);
  void pack_omegay(int);
  void pack_omegaz(int);
  void pack_angmomx(int);
  void pack_angmomy(int);
  void pack_angmomz(int);
  void pack_tqx(int);
  void pack_tqy(int);
  void pack_tqz(int);
};

}    // namespace LAMMPS_NS

#endif
#endif
