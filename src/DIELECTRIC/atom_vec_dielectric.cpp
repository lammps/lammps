/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "atom_vec_dielectric.h"

#include "atom.h"
#include "citeme.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "pair.h"
#include "pair_hybrid.h"

#include <cmath>

using namespace LAMMPS_NS;

static const char cite_user_dielectric_package[] =
    "DIELECTRIC package: doi:10.1016/j.cpc.2019.03.006\n\n"
    "@Article{TrungCPC19,\n"
    " author = {Trung Dac Nguyen and Honghao Li and Debarshee Bagchi and"
    "   Francisco J. Solis and Olvera de la Cruz, Monica}\n"
    " title = {Incorporating Surface Polarization Effects Into Large-Scale\n"
    "   Coarse-Grained Molecular Dynamics Simulation},\n"
    " journal = {Comput.\\ Phys.\\ Commun.},\n"
    " year =    2019,\n"
    " volume =  241,\n"
    " pages =   {80--91}\n"
    "}\n\n";

/* ---------------------------------------------------------------------- */

AtomVecDielectric::AtomVecDielectric(LAMMPS *_lmp) : AtomVec(_lmp)
{
  if (lmp->citeme) lmp->citeme->add(cite_user_dielectric_package);

  molecular = Atom::MOLECULAR;
  bonds_allow = angles_allow = dihedrals_allow = impropers_allow = 1;
  mass_type = PER_TYPE;

  atom->molecule_flag = atom->q_flag = atom->mu_flag = 1;
  atom->dielectric_flag = 1;

  mu_hold = nullptr;

  // strings with peratom variables to include in each AtomVec method
  // strings cannot contain fields in corresponding AtomVec default strings
  // order of fields in a string does not matter
  // except: fields_data_atom & fields_data_vel must match data file

  // clang-format off
  fields_grow = {"q", "molecule", "num_bond", "bond_type", "bond_atom", "num_angle", "angle_type",
    "angle_atom1", "angle_atom2", "angle_atom3", "num_dihedral", "dihedral_type", "dihedral_atom1",
    "dihedral_atom2", "dihedral_atom3", "dihedral_atom4", "num_improper", "improper_type",
    "improper_atom1", "improper_atom2", "improper_atom3", "improper_atom4", "nspecial", "special",
    "mu", "area", "ed", "em", "epsilon", "curvature", "q_scaled"};
  fields_copy = {"q", "molecule", "num_bond", "bond_type", "bond_atom", "num_angle", "angle_type",
    "angle_atom1", "angle_atom2", "angle_atom3", "num_dihedral", "dihedral_type", "dihedral_atom1",
    "dihedral_atom2", "dihedral_atom3", "dihedral_atom4", "num_improper", "improper_type",
    "improper_atom1", "improper_atom2", "improper_atom3", "improper_atom4", "nspecial", "special",
    "mu", "area", "ed", "em", "epsilon", "curvature", "q_scaled"};
  fields_comm = {"q", "mu", "area", "ed", "em", "epsilon", "curvature", "q_scaled"};
  fields_border = {"q", "molecule", "mu", "area", "ed", "em", "epsilon", "curvature", "q_scaled"};
  fields_border_vel = {"q", "molecule", "mu", "area", "ed", "em", "epsilon", "curvature",
    "q_scaled"};
  fields_exchange = {"q", "molecule", "num_bond", "bond_type", "bond_atom", "num_angle",
    "angle_type", "angle_atom1", "angle_atom2", "angle_atom3", "num_dihedral", "dihedral_type",
    "dihedral_atom1", "dihedral_atom2", "dihedral_atom3", "dihedral_atom4", "num_improper",
    "improper_type", "improper_atom1", "improper_atom2", "improper_atom3", "improper_atom4",
    "nspecial", "special", "mu", "area", "ed", "em", "epsilon", "curvature", "q_scaled"};
  fields_restart = {"q", "molecule", "num_bond", "bond_type", "bond_atom", "num_angle",
    "angle_type", "angle_atom1", "angle_atom2", "angle_atom3", "num_dihedral", "dihedral_type",
    "dihedral_atom1", "dihedral_atom2", "dihedral_atom3", "dihedral_atom4", "num_improper",
    "improper_type", "improper_atom1", "improper_atom2", "improper_atom3", "improper_atom4",
    "mu", "area", "ed", "em", "epsilon", "curvature", "q_scaled"};
  fields_create = {"q", "molecule", "num_bond", "num_angle", "num_dihedral", "num_improper",
    "nspecial", "mu", "area", "ed", "em", "epsilon", "curvature", "q_scaled"};
  fields_data_atom = {"id", "molecule", "type", "q", "x", "mu3", "area", "ed", "em", "epsilon",
    "curvature"};
  fields_data_vel = {"id", "v"};
  // clang-format on

  setup_fields();
  bond_per_atom = angle_per_atom = dihedral_per_atom = improper_per_atom = 0;
}

/* ---------------------------------------------------------------------- */

void AtomVecDielectric::init()
{
  AtomVec::init();

  // since atom style dielectric modifies the charge q, it will produce incorrect results
  // with pair styles using coulomb without dielectric support.

  std::string pair_style(force->pair_style);
  if ((pair_style != "none") && (pair_style != "zero") &&
      !utils::strmatch(force->pair_style, "/dielectric")) {
    bool mismatch = false;
    if (utils::strmatch(force->pair_style, "^reaxff")) mismatch = true;
    if (utils::strmatch(force->pair_style, "^comb")) mismatch = true;
    if (utils::strmatch(force->pair_style, "coul")) mismatch = true;
    if (utils::strmatch(force->pair_style, "tip4p")) mismatch = true;
    if (utils::strmatch(force->pair_style, "dipole")) mismatch = true;

    if (utils::strmatch(force->pair_style, "^hybrid")) {
      auto hybrid = dynamic_cast<PairHybrid *>(force->pair);
      if (hybrid) {
        for (int i = 0; i < hybrid->nstyles; i++) {
          if (utils::strmatch(hybrid->keywords[i], "^reaxff")) mismatch = true;
          if (utils::strmatch(hybrid->keywords[i], "^comb")) mismatch = true;
          if (utils::strmatch(hybrid->keywords[i], "coul")) mismatch = true;
          if (utils::strmatch(hybrid->keywords[i], "tip4p")) mismatch = true;
          if (utils::strmatch(hybrid->keywords[i], "dipole")) mismatch = true;
        }
      }
    }
    if (mismatch)
      error->all(FLERR, "Pair style {} is not compatible with atom style {}", pair_style,
                 atom->get_style());
  }
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
  q_scaled = atom->q_scaled;
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

  double *q = atom->q;
  q_scaled[ilocal] = q[ilocal]/epsilon[ilocal];

  double *mu_one = mu[ilocal];
  mu_one[3] = sqrt(mu_one[0] * mu_one[0] + mu_one[1] * mu_one[1] + mu_one[2] * mu_one[2]);
}

/* ----------------------------------------------------------------------
   convert read_data file info from general to restricted triclinic
   parent class operates on data from Velocities section of data file
   child class operates on dipole moment mu
------------------------------------------------------------------------- */

void AtomVecDielectric::read_data_general_to_restricted(int nlocal_previous, int nlocal)
{
  AtomVec::read_data_general_to_restricted(nlocal_previous, nlocal);

  for (int i = nlocal_previous; i < nlocal; i++)
    domain->general_to_restricted_vector(mu[i]);
}

/* ----------------------------------------------------------------------
   convert info output by write_data from restricted to general triclinic
   parent class operates on x and data from Velocities section of data file
   child class operates on dipole momemt mu which has 4 values per atom
------------------------------------------------------------------------- */

void AtomVecDielectric::write_data_restricted_to_general()
{
  AtomVec::write_data_restricted_to_general();

  int nlocal = atom->nlocal;
  memory->create(mu_hold,nlocal,3,"atomvec:mu_hold");
    for (int i = 0; i < nlocal; i++) {
    memcpy(&mu_hold[i],&mu[i],3*sizeof(double));
    domain->restricted_to_general_vector(mu[i]);
  }
}

/* ----------------------------------------------------------------------
   restore info output by write_data to restricted triclinic
   original data is in "hold" arrays
   parent class operates on x and data from Velocities section of data file
   child class operates on dipole moment mu which has 4 values per atom
------------------------------------------------------------------------- */

void AtomVecDielectric::write_data_restore_restricted()
{
  AtomVec::write_data_restore_restricted();

  if (!mu_hold) return;

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++)
    memcpy(&mu[i],&mu_hold[i],3*sizeof(double));
  memory->destroy(mu_hold);
  mu_hold = nullptr;
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

int AtomVecDielectric::property_atom(const std::string &name)
{
  if (name == "area") return 0;
  if (name == "ed") return 1;
  if (name == "em") return 2;
  if (name == "epsilon") return 3;
  if (name == "curvature") return 4;
  if (name == "q_scaled") return 5;
  return -1;
}

/* ----------------------------------------------------------------------
   pack per-atom data into buf for ComputePropertyAtom
   index maps to data specific to this atom style
------------------------------------------------------------------------- */

void AtomVecDielectric::pack_property_atom(int index, double *buf, int nvalues, int groupbit)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int n = 0;

  if (index == 0) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit)
        buf[n] = area[i];
      else
        buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 1) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit)
        buf[n] = ed[i];
      else
        buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 2) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit)
        buf[n] = em[i];
      else
        buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 3) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit)
        buf[n] = epsilon[i];
      else
        buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 4) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit)
        buf[n] = curvature[i];
      else
        buf[n] = 0.0;
      n += nvalues;
    }
  } else if (index == 5) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit)
        buf[n] = q_scaled[i];
      else
        buf[n] = 0.0;
      n += nvalues;
    }
  }
}
