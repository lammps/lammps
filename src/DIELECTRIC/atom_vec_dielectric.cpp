// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "atom_vec_dielectric.h"

#include "atom.h"
#include "citeme.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

static const char cite_user_dielectric_package[] =
  "DIELECTRIC package:\n\n"
  "@Article{TrungCPC19,\n"
  " author = {Trung Dac Nguyen, Honghao Li, Debarshee Bagchi,"
  " Francisco J. Solis, Monica Olvera de la Cruz,\n"
  " title = {Incorporating surface polarization effects into large-scale"
  " coarse-grained Molecular Dynamics simulation},\n"
  " journal = {Comp.~Phys.~Comm.},\n"
  " year =    2019,\n"
  " volume =  241,\n"
  " pages =   {80--91}\n"
  "}\n\n"
  ;

/* ---------------------------------------------------------------------- */

AtomVecDielectric::AtomVecDielectric(LAMMPS *lmp) : AtomVec(lmp)
{
  if (lmp->citeme) lmp->citeme->add(cite_user_dielectric_package);

  molecular = Atom::MOLECULAR;
  bonds_allow = angles_allow = dihedrals_allow = impropers_allow = 1;
  mass_type = PER_TYPE;

  atom->molecule_flag = atom->q_flag = atom->mu_flag = 1;

  // strings with peratom variables to include in each AtomVec method
  // strings cannot contain fields in corresponding AtomVec default strings
  // order of fields in a string does not matter
  // except: fields_data_atom & fields_data_vel must match data file

  fields_grow = (char *)
    "q molecule num_bond bond_type bond_atom "
    "num_angle angle_type angle_atom1 angle_atom2 angle_atom3 "
    "num_dihedral dihedral_type dihedral_atom1 dihedral_atom2 "
    "dihedral_atom3 dihedral_atom4 "
    "num_improper improper_type improper_atom1 improper_atom2 "
    "improper_atom3 improper_atom4 "
    "nspecial special "
    "mu area ed em epsilon curvature q_unscaled";
  fields_copy = (char *)
    "q molecule num_bond bond_type bond_atom "
    "num_angle angle_type angle_atom1 angle_atom2 angle_atom3 "
    "num_dihedral dihedral_type dihedral_atom1 dihedral_atom2 "
    "dihedral_atom3 dihedral_atom4 "
    "num_improper improper_type improper_atom1 improper_atom2 "
    "improper_atom3 improper_atom4 "
    "nspecial special "
    "mu area ed em epsilon curvature q_unscaled";
  fields_comm = (char *) "q mu area ed em epsilon curvature q_unscaled";
  fields_comm_vel = (char *) "";
  fields_reverse = (char *) "";
  fields_border = (char *) "q molecule mu area ed em epsilon curvature q_unscaled";
  fields_border_vel = (char *) "q molecule mu area ed em epsilon curvature q_unscaled";
  fields_exchange = (char *)
    "q molecule num_bond bond_type bond_atom "
    "num_angle angle_type angle_atom1 angle_atom2 angle_atom3 "
    "num_dihedral dihedral_type dihedral_atom1 dihedral_atom2 "
    "dihedral_atom3 dihedral_atom4 "
    "num_improper improper_type improper_atom1 improper_atom2 "
    "improper_atom3 improper_atom4 "
    "nspecial special "
    "mu area ed em epsilon curvature q_unscaled";
  fields_restart = (char *)
    "q molecule num_bond bond_type bond_atom "
    "num_angle angle_type angle_atom1 angle_atom2 angle_atom3 "
    "num_dihedral dihedral_type dihedral_atom1 dihedral_atom2 "
    "dihedral_atom3 dihedral_atom4 "
    "num_improper improper_type improper_atom1 improper_atom2 "
    "improper_atom3 improper_atom4 "
    "mu area ed em epsilon curvature q_unscaled";
  fields_create = (char *)
    "q molecule num_bond num_angle num_dihedral num_improper nspecial "
    "mu area ed em epsilon curvature q_unscaled";
  fields_data_atom = (char *) "id molecule type q x "
    "mu3 area ed em epsilon curvature";
  fields_data_vel = (char *) "id v";

  setup_fields();

  bond_per_atom = angle_per_atom = dihedral_per_atom = improper_per_atom = 0;
}

/* ----------------------------------------------------------------------
   set local copies of all grow ptrs used by this class, except defaults
   needed in replicate when 2 atom classes exist and it calls pack_restart()
------------------------------------------------------------------------- */

void AtomVecDielectric::grow_pointers()
{
  num_bond = atom->num_bond;
  bond_type = atom->bond_type;
  num_angle = atom->num_angle;
  angle_type = atom->angle_type;
  num_dihedral = atom->num_dihedral;
  dihedral_type = atom->dihedral_type;
  num_improper = atom->num_improper;
  improper_type = atom->improper_type;
  nspecial = atom->nspecial;

  mu = atom->mu;
  area = atom->area;
  ed = atom->ed;
  em = atom->em;
  epsilon = atom->epsilon;
  curvature = atom->curvature;
  q_unscaled = atom->q_unscaled;
}

/* ----------------------------------------------------------------------
   initialize non-zero atom quantities
------------------------------------------------------------------------- */

void AtomVecDielectric::create_atom_post(int ilocal)
{
  area[ilocal] = 1.0;
  em[ilocal] = 1.0;
  epsilon[ilocal] = 1.0;
}

/* ----------------------------------------------------------------------
   modify what AtomVec::data_atom() just unpacked
   or initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecDielectric::data_atom_post(int ilocal)
{
  num_bond[ilocal] = 0;
  num_angle[ilocal] = 0;
  num_dihedral[ilocal] = 0;
  num_improper[ilocal] = 0;
  nspecial[ilocal][0] = 0;
  nspecial[ilocal][1] = 0;
  nspecial[ilocal][2] = 0;

  double* q = atom->q;
  q_unscaled[ilocal] = q[ilocal];
  q[ilocal] /= epsilon[ilocal];

  double *mu_one = mu[ilocal];
  mu_one[3] =
    sqrt(mu_one[0]*mu_one[0] + mu_one[1]*mu_one[1] + mu_one[2]*mu_one[2]);
}

/* ----------------------------------------------------------------------
   initialize other atom quantities after AtomVec::unpack_restart()
------------------------------------------------------------------------- */

void AtomVecDielectric::unpack_restart_init(int ilocal)
{
  nspecial[ilocal][0] = 0;
  nspecial[ilocal][1] = 0;
  nspecial[ilocal][2] = 0;
}

/* ----------------------------------------------------------------------
   assign an index to named atom property and return index
   return -1 if name is unknown to this atom style
------------------------------------------------------------------------- */

int AtomVecDielectric::property_atom(char *name)
{
  if (strcmp(name,"area") == 0) return 0;
  if (strcmp(name,"ed") == 0) return 1;
  if (strcmp(name,"em") == 0) return 2;
  if (strcmp(name,"epsilon") == 0) return 3;
  if (strcmp(name,"curvature") == 0) return 4;
  if (strcmp(name,"q_unscaled") == 0) return 5;
  return -1;
}

/* ----------------------------------------------------------------------
   pack per-atom data into buf for ComputePropertyAtom
   index maps to data specific to this atom style
------------------------------------------------------------------------- */

void AtomVecDielectric::pack_property_atom(int index, double *buf,
                                           int nvalues, int groupbit)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int n = 0;

  if (index == 0) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = area[i];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 1) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = ed[i];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 2) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = em[i];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 3) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = epsilon[i];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 4) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = curvature[i];
      else buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 5) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) buf[n] = q_unscaled[i];
      else buf[n] = 0.0;
      n += nvalues;
    }
  }
}
