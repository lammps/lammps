/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "atom_vec_bpm_sphere.h"

#include "atom.h"
#include "error.h"
#include "fix.h"
#include "math_const.h"
#include "modify.h"

using namespace LAMMPS_NS;
using MathConst::MY_PI;

/* ---------------------------------------------------------------------- */

AtomVecBPMSphere::AtomVecBPMSphere(LAMMPS *_lmp) : AtomVec(_lmp)
{
  mass_type = PER_ATOM;
  molecular = Atom::MOLECULAR;
  bonds_allow = 1;
  radvary = 0;

  atom->molecule_flag = 1;
  atom->radius_flag = atom->rmass_flag = atom->omega_flag = atom->torque_flag = atom->quat_flag = 1;

  // strings with peratom variables to include in each AtomVec method
  // strings cannot contain fields in corresponding AtomVec default strings
  // order of fields in a string does not matter
  // except: fields_data_atom & fields_data_vel must match data file

  fields_grow = {"molecule", "num_bond", "bond_type", "bond_atom", "nspecial", "special",
                 "radius",   "rmass",    "omega",     "torque",    "quat"};
  fields_copy = {"molecule", "num_bond", "bond_type", "bond_atom", "nspecial",
                 "special",  "radius",   "rmass",     "omega",     "quat"};
  fields_comm_vel = {"omega", "quat"};
  fields_reverse = {"torque"};
  fields_border = {"molecule", "radius", "rmass"};
  fields_border_vel = {"molecule", "radius", "rmass", "omega", "quat"};
  fields_exchange = {"molecule", "num_bond", "bond_type", "bond_atom", "nspecial",
                     "special",  "radius",   "rmass",     "omega",     "quat"};
  fields_restart = {"molecule", "num_bond", "bond_type", "bond_atom",
                    "radius",   "rmass",    "omega",     "quat"};
  fields_create = {"molecule", "num_bond", "nspecial", "radius", "rmass", "omega", "quat"};
  fields_data_atom = {"id", "molecule", "type", "radius", "rmass", "x"};
  fields_data_vel = {"id", "v", "omega"};

  bond_per_atom = 0;
  bond_negative = nullptr;
}

/* ----------------------------------------------------------------------
   process sub-style args
   optional arg = 0/1 for static/dynamic particle radii
------------------------------------------------------------------------- */

void AtomVecBPMSphere::process_args(int narg, char **arg)
{
  if (narg != 0 && narg != 1) error->all(FLERR, "Illegal atom_style bpm/sphere command");

  radvary = 0;
  if (narg == 1) {
    radvary = utils::numeric(FLERR, arg[0], true, lmp);
    if (radvary < 0 || radvary > 1) error->all(FLERR, "Illegal atom_style bpm/sphere command");
  }

  // dynamic particle radius and mass must be communicated every step

  if (radvary) {
    fields_comm = {"radius", "rmass"};
    fields_comm_vel = {"radius", "rmass", "omega"};
  }

  // delay setting up of fields until now

  setup_fields();
}

/* ---------------------------------------------------------------------- */

void AtomVecBPMSphere::init()
{
  AtomVec::init();

  // check if optional radvary setting should have been set to 1

  if (radvary == 0)
    for (const auto &ifix : modify->get_fix_by_style("^adapt")) {
      if (ifix->diam_flag)
        error->all(FLERR, "Fix {} changes atom radii but atom_style bpm/sphere is not dynamic",
                   ifix->style);
    }
}

/* ----------------------------------------------------------------------
   set local copies of all grow ptrs used by this class, except defaults
   needed in replicate when 2 atom classes exist and it calls pack_restart()
------------------------------------------------------------------------- */

void AtomVecBPMSphere::grow_pointers()
{
  radius = atom->radius;
  rmass = atom->rmass;
  omega = atom->omega;
  quat = atom->quat;

  num_bond = atom->num_bond;
  bond_type = atom->bond_type;
  nspecial = atom->nspecial;
}

/* ----------------------------------------------------------------------
   initialize non-zero atom quantities
------------------------------------------------------------------------- */

void AtomVecBPMSphere::create_atom_post(int ilocal)
{
  radius[ilocal] = 0.5;
  rmass[ilocal] = 4.0 * MY_PI / 3.0 * 0.5 * 0.5 * 0.5;

  quat[ilocal][0] = 1.0;
  quat[ilocal][1] = 0.0;
  quat[ilocal][2] = 0.0;
  quat[ilocal][3] = 0.0;
}

/* ----------------------------------------------------------------------
   modify values for AtomVec::pack_restart() to pack
------------------------------------------------------------------------- */

void AtomVecBPMSphere::pack_restart_pre(int ilocal)
{
  // ensure bond_negative vector is needed length

  if (bond_per_atom < atom->bond_per_atom) {
    delete[] bond_negative;
    bond_per_atom = atom->bond_per_atom;
    bond_negative = new int[bond_per_atom];
  }

  // flip any negative types to positive and flag which ones

  any_bond_negative = 0;
  for (int m = 0; m < num_bond[ilocal]; m++) {
    if (bond_type[ilocal][m] < 0) {
      bond_negative[m] = 1;
      bond_type[ilocal][m] = -bond_type[ilocal][m];
      any_bond_negative = 1;
    } else
      bond_negative[m] = 0;
  }
}

/* ----------------------------------------------------------------------
   unmodify values packed by AtomVec::pack_restart()
------------------------------------------------------------------------- */

void AtomVecBPMSphere::pack_restart_post(int ilocal)
{
  // restore the flagged types to their negative values

  if (any_bond_negative) {
    for (int m = 0; m < num_bond[ilocal]; m++)
      if (bond_negative[m]) bond_type[ilocal][m] = -bond_type[ilocal][m];
  }
}

/* ----------------------------------------------------------------------
   initialize other atom quantities after AtomVec::unpack_restart()
------------------------------------------------------------------------- */

void AtomVecBPMSphere::unpack_restart_init(int ilocal)
{
  nspecial[ilocal][0] = 0;
  nspecial[ilocal][1] = 0;
  nspecial[ilocal][2] = 0;
}

/* ----------------------------------------------------------------------
   modify what AtomVec::data_atom() just unpacked
   or initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecBPMSphere::data_atom_post(int ilocal)
{
  radius_one = 0.5 * atom->radius[ilocal];
  radius[ilocal] = radius_one;
  if (radius_one > 0.0) rmass[ilocal] *= 4.0 * MY_PI / 3.0 * radius_one * radius_one * radius_one;

  if (rmass[ilocal] <= 0.0) error->one(FLERR, "Invalid density in Atoms section of data file");

  omega[ilocal][0] = 0.0;
  omega[ilocal][1] = 0.0;
  omega[ilocal][2] = 0.0;

  quat[ilocal][0] = 1.0;
  quat[ilocal][1] = 0.0;
  quat[ilocal][2] = 0.0;
  quat[ilocal][3] = 0.0;

  num_bond[ilocal] = 0;
  nspecial[ilocal][0] = 0;
  nspecial[ilocal][1] = 0;
  nspecial[ilocal][2] = 0;
}

/* ----------------------------------------------------------------------
   modify values for AtomVec::pack_data() to pack
------------------------------------------------------------------------- */

void AtomVecBPMSphere::pack_data_pre(int ilocal)
{
  radius_one = radius[ilocal];
  rmass_one = rmass[ilocal];

  radius[ilocal] *= 2.0;
  if (radius_one != 0.0)
    rmass[ilocal] = rmass_one / (4.0 * MY_PI / 3.0 * radius_one * radius_one * radius_one);
}

/* ----------------------------------------------------------------------
   unmodify values packed by AtomVec::pack_data()
------------------------------------------------------------------------- */

void AtomVecBPMSphere::pack_data_post(int ilocal)
{
  radius[ilocal] = radius_one;
  rmass[ilocal] = rmass_one;
}
