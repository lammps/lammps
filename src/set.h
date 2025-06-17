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

#ifdef COMMAND_CLASS
// clang-format off
CommandStyle(set,Set);
// clang-format on
#else

#ifndef LMP_SET_H
#define LMP_SET_H

#include "command.h"

namespace LAMMPS_NS {

class Set : public Command {
 public:
  Set(class LAMMPS *lmp);
  ~Set() override;

  void command(int, char **) override;

  void process_args(int, int, char **);
  void selection(int);
  void invoke_actions();

 private:
  int caller;    // SETCOMMAND or FIXSET

  // params for atom selection

  int style;
  char *id;
  int nlo, nhi;
  bigint nlobig, nhibig;
  int groupbit;
  class Region *region;

  // one Action = one keyword/value pair

  struct Action {
    int keyword;
    int argindex;
    int count_select, count_action;
    int varflag;
    int varflag1, varflag2, varflag3, varflag4;
    int ivar1, ivar2, ivar3, ivar4;
    int ivalue1, ivalue2, ivalue3, ivalue4, ivalue5, ivalue6;
    tagint tvalue1;
    bigint bvalue1;
    double dvalue1, dvalue2, dvalue3, dvalue4;
  };

  int naction, maxaction;
  Action *actions;

  typedef void (Set::*FnPtrPack)(Action *);
  FnPtrPack *invoke_choice;    // list of ptrs to invoke functions

  double *vec1, *vec2, *vec3, *vec4;             // storage for evaluated peratom variables
  int varflag;                                   // 1 if any action uses a per-atom variable
  int varflag1, varflag2, varflag3, varflag4;    // 1 if any action uses this variable
  int maxvariable;                               // size of peratom variable vectors

  int *select;      // flag for selected atoms
  int maxselect;    // size of select vector

  int count_select;    // count of selected atoms on this proc
  int count_action;    // count of actions on this proc, only if different than selected

  // private functions

  void varparse(const char *, int, Action *);
  void setrandom(int, Action *);
  void topology(int, Action *);

  // customize by adding a process method

  void process_angle(int &, int, char **, Action *);
  void process_angmom(int &, int, char **, Action *);
  void process_bond(int &, int, char **, Action *);
  void process_cc(int &, int, char **, Action *);
  void process_charge(int &, int, char **, Action *);
  void process_density(int &, int, char **, Action *);
  void process_diameter(int &, int, char **, Action *);
  void process_dihedral(int &, int, char **, Action *);
  void process_dipole(int &, int, char **, Action *);
  void process_dipole_random(int &, int, char **, Action *);
  void process_dpd_theta(int &, int, char **, Action *);
  void process_edpd_cv(int &, int, char **, Action *);
  void process_edpd_temp(int &, int, char **, Action *);
  void process_epsilon(int &, int, char **, Action *);
  void process_image(int &, int, char **, Action *);
  void process_improper(int &, int, char **, Action *);
  void process_length(int &, int, char **, Action *);
  void process_mass(int &, int, char **, Action *);
  void process_mol(int &, int, char **, Action *);
  void process_omega(int &, int, char **, Action *);
  void process_quat(int &, int, char **, Action *);
  void process_quat_random(int &, int, char **, Action *);
  void process_radius_election(int &, int, char **, Action *);
  void process_shape(int &, int, char **, Action *);
  void process_smd_contact_radius(int &, int, char **, Action *);
  void process_smd_mass_density(int &, int, char **, Action *);
  void process_sph_cv(int &, int, char **, Action *);
  void process_sph_e(int &, int, char **, Action *);
  void process_sph_rho(int &, int, char **, Action *);
  void process_spin_atom(int &, int, char **, Action *);
  void process_spin_atom_random(int &, int, char **, Action *);
  void process_spin_electron(int &, int, char **, Action *);
  void process_temperature(int &, int, char **, Action *);
  void process_theta(int &, int, char **, Action *);
  void process_theta_random(int &, int, char **, Action *);
  void process_tri(int &, int, char **, Action *);
  void process_type(int &, int, char **, Action *);
  void process_type_fraction(int &, int, char **, Action *);
  void process_type_ratio(int &, int, char **, Action *);
  void process_type_subset(int &, int, char **, Action *);
  void process_volume(int &, int, char **, Action *);
  void process_vx(int &, int, char **, Action *);
  void process_vy(int &, int, char **, Action *);
  void process_vz(int &, int, char **, Action *);
  void process_x(int &, int, char **, Action *);
  void process_y(int &, int, char **, Action *);
  void process_z(int &, int, char **, Action *);

  void process_custom(int &, int, char **, Action *);

  // customize by adding an invoke method

  void invoke_angle(Action *);
  void invoke_angmom(Action *);
  void invoke_bond(Action *);
  void invoke_cc(Action *);
  void invoke_charge(Action *);
  void invoke_density(Action *);
  void invoke_diameter(Action *);
  void invoke_dihedral(Action *);
  void invoke_dipole(Action *);
  void invoke_dipole_random(Action *);
  void invoke_dpd_theta(Action *);
  void invoke_edpd_cv(Action *);
  void invoke_edpd_temp(Action *);
  void invoke_epsilon(Action *);
  void invoke_image(Action *);
  void invoke_improper(Action *);
  void invoke_length(Action *);
  void invoke_mass(Action *);
  void invoke_mol(Action *);
  void invoke_omega(Action *);
  void invoke_quat(Action *);
  void invoke_quat_random(Action *);
  void invoke_radius_election(Action *);
  void invoke_shape(Action *);
  void invoke_smd_contact_radius(Action *);
  void invoke_smd_mass_density(Action *);
  void invoke_sph_cv(Action *);
  void invoke_sph_e(Action *);
  void invoke_sph_rho(Action *);
  void invoke_spin_atom(Action *);
  void invoke_spin_atom_random(Action *);
  void invoke_spin_electron(Action *);
  void invoke_temperature(Action *);
  void invoke_theta(Action *);
  void invoke_theta_random(Action *);
  void invoke_tri(Action *);
  void invoke_type(Action *);
  void invoke_type_fraction(Action *);
  void invoke_type_ratio(Action *);
  void invoke_type_subset(Action *);
  void invoke_volume(Action *);
  void invoke_vx(Action *);
  void invoke_vy(Action *);
  void invoke_vz(Action *);
  void invoke_x(Action *);
  void invoke_y(Action *);
  void invoke_z(Action *);

  void invoke_custom(Action *);
};

}    // namespace LAMMPS_NS

#endif
#endif
