/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_BOND_BPM_H
#define LMP_BOND_BPM_H

#include "bond.h"

namespace LAMMPS_NS {

class BondBPM : public Bond {
 public:
  BondBPM(class LAMMPS *);
  virtual ~BondBPM();
  virtual void compute(int, int) = 0;
  virtual void coeff(int, char **) = 0;
  virtual void init_style();
  void settings(int, char **);  
  double equilibrium_distance(int);
  void write_restart(FILE *){};
  void read_restart(FILE *){};
  void write_data(FILE *) {};
  double single(int, double, int, int, double &) = 0;

 protected:
  double r0_max_estimate;
  double max_stretch;
  
  char *id_fix_store_local, *id_fix_prop_atom;
  class FixStoreLocal *fix_store_local;
  class FixUpdateSpecialBonds *fix_update_special_bonds;
  class FixBondHistory *fix_bond_history;

  void process_broken(int, int);  
  typedef void (BondBPM::*FnPtrPack)(int,int,int);
  FnPtrPack *pack_choice;              // ptrs to pack functions
  double *output_data;

  int prop_atom_flag, nvalues;
  int index_x_ref, index_y_ref, index_z_ref;

  void pack_id1(int,int,int);
  void pack_id2(int,int,int);
  void pack_time(int,int,int);
  void pack_x(int,int,int);
  void pack_y(int,int,int);
  void pack_z(int,int,int);
  void pack_x_ref(int,int,int);
  void pack_y_ref(int,int,int);
  void pack_z_ref(int,int,int);  
};

}    // namespace LAMMPS_NS

#endif

/* ERROR/WARNING messages:

E: Cannot find fix store/local

Fix id cannot be found.

E: Illegal bond_style command

Self-explanatory.

E: Bond style bpm/rotational must include at least one value to output

Must include at least one bond property to store in fix store/local

E: Bond style bpm/rotational cannot be used with 3,4-body interactions

No angle, dihedral, or improper styles can be defined when using
bond style bpm/rotational.

E: Bond style bpm/rotational cannot be used with atom style template

This bond style can change the bond topology which is not
allowed with this atom style.

*/
