/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributead under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "atom_vec_sphere.h"
#include <cstring>
#include "atom.h"
#include "modify.h"
#include "fix.h"
#include "fix_adapt.h"
#include "math_const.h"
#include "error.h"
#include "utils.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

AtomVecSphere::AtomVecSphere(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = 0;

  atom->sphere_flag = 1;
  atom->radius_flag = atom->rmass_flag = atom->omega_flag =
    atom->torque_flag = 1;

  // strings with peratom variables to include in each AtomVec method
  // strings cannot contain fields in corresponding AtomVec default strings
  // order of fields in the string does not matter
  //   except fields_data_atom and fields_data_vel which must match data file

  fields_grow = (char *) "radius rmass omega torque";
  fields_copy = (char *) "radius rmass omega";
  fields_comm = NULL;
  fields_comm_vel = (char *) "omega";
  fields_reverse = (char *) "torque";
  fields_border = (char *) "radius rmass";
  fields_border_vel = (char *) "radius rmass omega";
  fields_exchange = (char *) "radius rmass omega";
  fields_restart = (char *) "radius rmass omega";
  fields_create = (char *) "radius rmass omega";
  fields_data_atom = (char *) "id type radius rmass x";
  fields_data_vel = (char *) "omega";

  setup_fields();
}

/* ----------------------------------------------------------------------
   process sub-style args
   optional arg = 0/1 for static/dynamic particle radii
------------------------------------------------------------------------- */

void AtomVecSphere::process_args(int narg, char **arg)
{
  if (narg == 0) return;
  if (narg != 1) error->all(FLERR,"Illegal atom_style sphere command");

  radvary = utils::numeric(FLERR,arg[0],true,lmp);
  if (radvary < 0 || radvary > 1)
    error->all(FLERR,"Illegal atom_style sphere command");
  if (radvary == 0) return;

  // dynamic particle radius and mass must be communicated every step

  fields_comm = (char *) "radius rmass";
  fields_comm_vel = (char *) "radius rmass omega";
}

/* ---------------------------------------------------------------------- */

void AtomVecSphere::init()
{
  AtomVec::init();

  // check if optional radvary setting should have been set to 1

  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"adapt") == 0) {
      FixAdapt *fix = (FixAdapt *) modify->fix[i];
      if (fix->diamflag && radvary == 0)
        error->all(FLERR,"Fix adapt changes particle radii "
                   "but atom_style sphere is not dynamic");
    }
}

/* ----------------------------------------------------------------------
   initialize non-zero atom quantities
------------------------------------------------------------------------- */

void AtomVecSphere::create_atom_post(int ilocal)
{
  atom->radius[ilocal] = 0.5;
  atom->rmass[ilocal] = 4.0*MY_PI/3.0 * 0.5*0.5*0.5;
}

/* ----------------------------------------------------------------------
   modify what AtomVec::data_atom() just unpacked
   or initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecSphere::data_atom_post(int ilocal)
{
  double radius = 0.5 * atom->radius[ilocal];
  atom->radius[ilocal] = radius;
  if (radius > 0.0) 
    atom->rmass[ilocal] = 
      4.0*MY_PI/3.0 * radius*radius*radius * atom->rmass[ilocal];

  if (atom->rmass[ilocal] <= 0.0) 
    error->one(FLERR,"Invalid density in Atoms section of data file");
}

/* ----------------------------------------------------------------------
   modify values for AtomVec::pack_data() to pack
------------------------------------------------------------------------- */

void AtomVecSphere::pack_data_pre(int ilocal)
{ 
  radius = atom->radius[ilocal];
  rmass = atom->rmass[ilocal];

  atom->radius[ilocal] *= 2.0;
  if (radius == 0.0) 
    atom->rmass[ilocal] = rmass / (4.0*MY_PI/3.0 * radius*radius*radius);
}

/* ----------------------------------------------------------------------
   unmodify values packed by AtomVec::pack_data()
------------------------------------------------------------------------- */

void AtomVecSphere::pack_data_post(int ilocal)
{ 
  atom->radius[ilocal] = radius;
  atom->rmass[ilocal] = rmass;
}
