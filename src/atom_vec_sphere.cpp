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

#include "atom_vec_sphere.h"

#include "atom.h"
#include "error.h"
#include "fix.h"
#include "fix_adapt.h"
#include "math_const.h"
#include "modify.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

AtomVecSphere::AtomVecSphere(LAMMPS *lmp) : AtomVec(lmp)
{
  mass_type = PER_ATOM;
  molecular = Atom::ATOMIC;

  atom->sphere_flag = 1;
  atom->radius_flag = atom->rmass_flag = atom->omega_flag = atom->torque_flag = 1;

  // strings with peratom variables to include in each AtomVec method
  // strings cannot contain fields in corresponding AtomVec default strings
  // order of fields in a string does not matter
  // except: fields_data_atom & fields_data_vel must match data file

  fields_grow = {"radius", "rmass", "omega", "torque"};
  fields_copy = {"radius", "rmass", "omega"};
  fields_comm_vel = {"omega"};
  fields_reverse = {"torque"};
  fields_border = {"radius", "rmass"};
  fields_border_vel = {"radius", "rmass", "omega"};
  fields_exchange = {"radius", "rmass", "omega"};
  fields_restart = {"radius", "rmass", "omega"};
  fields_create = {"radius", "rmass", "omega"};
  fields_data_atom = {"id", "type", "radius", "rmass", "x"};
  fields_data_vel = {"id", "v", "omega"};
}

/* ----------------------------------------------------------------------
   process sub-style args
   optional arg = 0/1 for static/dynamic particle radii
------------------------------------------------------------------------- */

void AtomVecSphere::process_args(int narg, char **arg)
{
  if (narg != 0 && narg != 1) error->all(FLERR, "Illegal atom_style sphere command");

  radvary = 0;
  if (narg == 1) {
    radvary = utils::numeric(FLERR, arg[0], true, lmp);
    if (radvary < 0 || radvary > 1) error->all(FLERR, "Illegal atom_style sphere command");
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

void AtomVecSphere::init()
{
  AtomVec::init();

  // check if optional radvary setting should have been set to 1

  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style, "adapt") == 0) {
      auto fix = dynamic_cast<FixAdapt *>(modify->fix[i]);
      if (fix->diamflag && radvary == 0)
        error->all(FLERR, "Fix adapt changes particle radii but atom_style sphere is not dynamic");
    }
}

/* ----------------------------------------------------------------------
   set local copies of all grow ptrs used by this class, except defaults
   needed in replicate when 2 atom classes exist and it calls pack_restart()
------------------------------------------------------------------------- */

void AtomVecSphere::grow_pointers()
{
  radius = atom->radius;
  rmass = atom->rmass;
  omega = atom->omega;
}

/* ----------------------------------------------------------------------
   initialize non-zero atom quantities
------------------------------------------------------------------------- */

void AtomVecSphere::create_atom_post(int ilocal)
{
  radius[ilocal] = 0.5;
  rmass[ilocal] = 4.0 * MY_PI / 3.0 * 0.5 * 0.5 * 0.5;
}

/* ----------------------------------------------------------------------
   modify what AtomVec::data_atom() just unpacked
   or initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecSphere::data_atom_post(int ilocal)
{
  radius_one = 0.5 * atom->radius[ilocal];
  radius[ilocal] = radius_one;
  if (radius_one > 0.0) rmass[ilocal] *= 4.0 * MY_PI / 3.0 * radius_one * radius_one * radius_one;

  if (rmass[ilocal] <= 0.0) error->one(FLERR, "Invalid density in Atoms section of data file");

  omega[ilocal][0] = 0.0;
  omega[ilocal][1] = 0.0;
  omega[ilocal][2] = 0.0;
}

/* ----------------------------------------------------------------------
   modify values for AtomVec::pack_data() to pack
------------------------------------------------------------------------- */

void AtomVecSphere::pack_data_pre(int ilocal)
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

void AtomVecSphere::pack_data_post(int ilocal)
{
  radius[ilocal] = radius_one;
  rmass[ilocal] = rmass_one;
}
