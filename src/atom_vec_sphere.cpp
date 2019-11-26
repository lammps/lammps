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

/* ---------------------------------------------------------------------- */

void AtomVecSphere::init()
{
  AtomVec::init();

  // set radvary if particle diameters are time-varying due to fix adapt
  // NOTE: change this to a atom_style sphere optional arg

  radvary = 0;
  //comm_x_only = 1;
  //size_forward = 3;

  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"adapt") == 0) {
      FixAdapt *fix = (FixAdapt *) modify->fix[i];
      if (fix->diamflag) {
        radvary = 1;
        comm_x_only = 0;
        size_forward = 5;
      }
    }

  //fields_comm = (char *) "radius rmass";
  //fields_comm_vel = (char *) "radius rmass";
}

/* ----------------------------------------------------------------------
   create one atom of itype at coord
   modify what default AtomVec::create_atom() just created
------------------------------------------------------------------------- */

void AtomVecSphere::create_atom(int itype, double *coord)
{
  AtomVec::create_atom(itype,coord);
  int ilocal = atom->nlocal-1;

  atom->radius[ilocal] = 0.5;
  atom->rmass[ilocal] = 4.0*MY_PI/3.0 * 0.5*0.5*0.5;
}

/* ----------------------------------------------------------------------
   unpack one line from Atoms section of data file
   modify what default AtomVec::data_atom() just unpacked
   or initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecSphere::data_atom(double *coord, imageint imagetmp, char **values)
{
  AtomVec::data_atom(coord,imagetmp,values);
  int ilocal = atom->nlocal-1;

  double radius = 0.5 * atom->radius[ilocal];
  atom->radius[ilocal] = radius;
  if (radius > 0.0) 
    atom->rmass[ilocal] = 
      4.0*MY_PI/3.0 * radius*radius*radius * atom->rmass[ilocal];

  if (atom->rmass[ilocal] <= 0.0) 
    error->one(FLERR,"Invalid mass in Atoms section of data file");
}
