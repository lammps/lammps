/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef DUMP_CUSTOM_H
#define DUMP_CUSTOM_H

#include "dump.h"

class DumpCustom : public Dump {
 public:
  DumpCustom(int, char **);
  ~DumpCustom();
  void init();
  int memory_usage();

 private:
  int iregion;               // -1 if no region, else which region
  int nthresh;               // # of defined threshholds
  int *thresh_array;         // array to threshhhold on for each nthresh
  int *thresh_op;            // threshhold operation for each nthresh
  double *thresh_value;      // threshhold value for each nthresh

  int nmine;                 // # of lines I am dumping
  int *vtype;                // type of each vector (INT, DOUBLE)
  char **vformat;            // format string for each vector element

  int maxlocal;              // size of choose and local-compute arrays
  int *choose;               // 1 if output this atom, 0 if no
  double *dchoose;           // value for each atom to threshhold against

  int nevery;

  int fixflag;               // 1 if this dump invokes any fixes, 0 if not
  int fixflag_epair;         // 1 if this fix is invoked, 0 if not
  int fixflag_centro;
  int fixflag_stress;
  int ifix_epair;            // which fix it is (0 to nfixes-1)
  int ifix_centro;
  int ifix_stress;

  int computeflag;           // 1 if this dump invokes any local comp, 0 if not
  int computeflag_etotal;    // 1 if this local comp is invoked, 0 if not
  int computeflag_ke;

  // local computation arrays

  double *etotal;
  double *ke;

  // private methods

  int modify_param(int, char **);
  void write_header(int);
  int count();
  int pack();
  void write_data(int, double *);

  void fix_create(char *);
  void fix_delete(char *);
  int fix_match(char *);

  typedef void (DumpCustom::*FnPtrHeader)(int);
  FnPtrHeader header_choice;           // ptr to write header functions
  void header_binary(int);
  void header_item(int);

  typedef void (DumpCustom::*FnPtrData)(int, double *);
  FnPtrData write_choice;              // ptr to write data functions
  void write_binary(int, double *);
  void write_text(int, double *);

  // customize by adding a method prototype

  typedef void (DumpCustom::*FnPtrPack)(int);
  FnPtrPack *pack_choice;              // ptrs to pack functions
  void pack_tag(int);
  void pack_molecule(int);
  void pack_type(int);
  void pack_x(int);
  void pack_y(int);
  void pack_z(int);
  void pack_xs(int);
  void pack_ys(int);
  void pack_zs(int);
  void pack_xu(int);
  void pack_yu(int);
  void pack_zu(int);
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
  void pack_tqx(int);
  void pack_tqy(int);
  void pack_tqz(int);
  void pack_etotal(int);
  void pack_ke(int);
  void pack_epair(int);
  void pack_centro(int);
  void pack_sxx(int);
  void pack_syy(int);
  void pack_szz(int);
  void pack_sxy(int);
  void pack_sxz(int);
  void pack_syz(int);

  // local computation methods

  void compute_etotal();
  void compute_ke();
};

#endif
