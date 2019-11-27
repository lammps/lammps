/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: James Larentzos (U.S. Army Research Laboratory)
------------------------------------------------------------------------- */

#include "atom_vec_dpd.h"
#include "atom.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecDPD::AtomVecDPD(LAMMPS *lmp) : AtomVec(lmp)
{
  molecular = 0;
  mass_type = 1;

  atom->rho_flag = 1;
  atom->dpd_flag = 1;

  // strings with peratom variables to include in each AtomVec method
  // strings cannot contain fields in corresponding AtomVec default strings
  // order of fields in the string does not matter
  //   except fields_data_atom and fields_data_vel which must match data file

  fields_grow = (char *) "rho dpdTheta uCond uMech uChem uCG uCGnew duChem";
  fields_copy = (char *) "dpdTheta uCond uMech uChem uCG uCGnew";
  fields_comm = (char *) "dpdTheta uCond uMech uChem";
  fields_comm_vel = (char *) "dpdTheta uCond uMech uChem";
  fields_reverse = NULL;
  fields_border = (char *) "dpdTheta uCond uMech uChem uCG uCGnew";
  fields_border_vel = (char *) "dpdTheta uCond uMech uChem uCG uCGnew";
  fields_exchange = (char *) "dpdTheta uCond uMech uChem uCG uCGnew";
  fields_restart = (char *) "dpdTheta uCond uMech uChem";
  fields_create = (char *) "rho dpdTheta uCond uMech uChem uCG uCGnew duChem";
  fields_data_atom = (char *) "id type dpdTheta x";
  fields_data_vel = (char *) "omega";

  setup_fields();
}

/* ----------------------------------------------------------------------
   initialize other atom quantities after AtomVec::unpack_restart()
------------------------------------------------------------------------- */

void AtomVecDPD::unpack_restart_init(int ilocal)
{
  atom->uCG[ilocal] = 0.0;
  atom->uCGnew[ilocal] = 0.0;
}

/* ----------------------------------------------------------------------
   modify what AtomVec::data_atom() just unpacked
   or initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecDPD::data_atom_post(int ilocal)
{
  atom->rho[ilocal] = 0.0;
  atom->uCond[ilocal] = 0.0;
  atom->uMech[ilocal] = 0.0;
  atom->uChem[ilocal] = 0.0;
  atom->uCG[ilocal] = 0.0;
  atom->uCGnew[ilocal] = 0.0;

  if (atom->dpdTheta[ilocal] <= 0)
    error->one(FLERR,"Internal temperature in Atoms section of date file "
               "must be > zero");
}
