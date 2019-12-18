/* ----------------------------------------------------------------------
 *
 *                    *** Smooth Mach Dynamics ***
 *
 * This file is part of the USER-SMD package for LAMMPS.
 * Copyright (2014) Georg C. Ganzenmueller, georg.ganzenmueller@emi.fhg.de
 * Fraunhofer Ernst-Mach Institute for High-Speed Dynamics, EMI,
 * Eckerstrasse 4, D-79104 Freiburg i.Br, Germany.
 *
 * ----------------------------------------------------------------------- */

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

#include "atom_vec_smd.h"
#include <cstring>
#include "atom.h"
#include "error.h"

using namespace LAMMPS_NS;

#define NMAT_FULL 9
#define NMAT_SYMM 6

/* ---------------------------------------------------------------------- */

AtomVecSMD::AtomVecSMD(LAMMPS *lmp) : AtomVec(lmp) 
{
  molecular = 0;
  mass_type = 1;
  forceclearflag = 1;

  atom->smd_flag = 1;

  atom->radius_flag = 1;
  atom->rmass_flag = 1;
  atom->vfrac_flag = 1;
  atom->contact_radius_flag = 1;
  atom->molecule_flag = 1;
  atom->smd_data_9_flag = 1;
  atom->e_flag = 1;
  atom->vest_flag = 1;
  atom->smd_stress_flag = 1;
  atom->eff_plastic_strain_flag = 1;
  atom->x0_flag = 1;
  atom->damage_flag = 1;
  atom->eff_plastic_strain_rate_flag = 1;

  // strings with peratom variables to include in each AtomVec method
  // strings cannot contain fields in corresponding AtomVec default strings
  // order of fields in a string does not matter
  // except: fields_data_atom & fields_data_vel must match data file

  fields_grow = (char *) 
    "e de vfrac rmass x0 radius contact_radius molecule "
    "smd_data_9 vest smd_stress "
    "eff_plastic_strain eff_plastic_strain_rate damage";
  fields_copy = (char *) 
    "e vfrac rmass x0 radius contact_radius molecule "
    "eff_plastic_strain eff_plastic_strain_rate vest "
    "smd_data_9 smd_stress damage";
  fields_comm = (char *) "radius vfrac vest e";
  fields_comm_vel = (char *) "radius vfrac vest e";
  fields_reverse = (char *) "de";
  fields_border = (char *) 
    "x0 molecule radius rmass vfrac contact_radius e "
    "eff_plastic_strain smd_data_9 smd_stress";
  fields_border_vel = (char *) 
    "x0 molecule radius rmass vfrac contact_radius e "
    "eff_plastic_strain smd_data_9 smd_stress vest";
  fields_exchange = (char *) 
    "x0 molecule radius rmass vfrac contact_radius e "
    "eff_plastic_strain eff_plastic_strain_rate smd_data_9 smd_stress "
    "vest damage";
  fields_restart = (char *) 
    "x0 molecule radius rmass vfrac contact_radius e "
    "eff_plastic_strain eff_plastic_strain_rate smd_data_9 smd_stress "
    "vest damage";
  fields_create = (char *) 
    "x0 vest vfrac rmass radius contact_radius molecule e "
    "eff_plastic_strain eff_plastic_strain_rate smd_data_9 smd_stress damage";
  fields_data_atom = (char *) 
    "id type molecule vfrac rmass radius contact_radius x0 x";
  fields_data_vel = (char *) "id v";

  // set these array sizes based on defines

  atom->add_peratom_change_columns("smd_data_9",NMAT_FULL);
  atom->add_peratom_change_columns("smd_stress",NMAT_SYMM);

  setup_fields();
}

/* ----------------------------------------------------------------------
   set local copies of all grow ptrs used by this class, except defaults
   needed in replicate when 2 atom classes exist and it calls pack_restart()
------------------------------------------------------------------------- */

void AtomVecSMD::grow_pointers()
{
  e = atom->e;
  de = atom->de;
  vfrac = atom->vfrac;
  rmass = atom->rmass;
  x0 = atom->x0;
  radius = atom->radius;
  contact_radius = atom->contact_radius;
  molecule = atom->molecule;
  smd_data_9 = atom->smd_data_9;
  vest = atom->vest;
  smd_stress = atom->smd_stress;
  eff_plastic_strain = atom->eff_plastic_strain;
  eff_plastic_strain_rate = atom->eff_plastic_strain_rate;
  damage = atom->damage;
}

/* ----------------------------------------------------------------------
   clear extra forces starting at atom N
   nbytes = # of bytes to clear for a per-atom vector
   NOTE: does f need to be re-cleared?
------------------------------------------------------------------------- */

void AtomVecSMD::force_clear(int n, size_t nbytes) 
{
  memset(&de[n],0,nbytes);
  memset(&f[n][0],0,3*nbytes);
}

/* ----------------------------------------------------------------------
   initialize non-zero atom quantities
------------------------------------------------------------------------- */

void AtomVecSMD::create_atom_post(int ilocal)
{
  x0[ilocal][0] = x[ilocal][0];
  x0[ilocal][1] = x[ilocal][1];
  x0[ilocal][2] = x[ilocal][2];

  vfrac[ilocal] = 1.0;
  rmass[ilocal] = 1.0;
  radius[ilocal] = 0.5;
  contact_radius[ilocal] = 0.5;
  molecule[ilocal] = 1;
  
  smd_data_9[ilocal][0] = 1.0; // xx
  smd_data_9[ilocal][4] = 1.0; // yy
  smd_data_9[ilocal][8] = 1.0; // zz
}

/* ----------------------------------------------------------------------
   modify what AtomVec::data_atom() just unpacked
   or initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecSMD::data_atom_post(int ilocal)
{
  e[ilocal] = 0.0;
  x0[ilocal][0] = x[ilocal][0];
  x0[ilocal][1] = x[ilocal][1];
  x0[ilocal][2] = x[ilocal][2];

  vest[ilocal][0] = 0.0;
  vest[ilocal][1] = 0.0;
  vest[ilocal][2] = 0.0;

  damage[ilocal] = 0.0;

  eff_plastic_strain[ilocal] = 0.0;
  eff_plastic_strain_rate[ilocal] = 0.0;

  for (int k = 0; k < NMAT_FULL; k++)
    smd_data_9[ilocal][k] = 0.0;

  for (int k = 0; k < NMAT_SYMM; k++)
    smd_stress[ilocal][k] = 0.0;

  smd_data_9[ilocal][0] = 1.0; // xx
  smd_data_9[ilocal][4] = 1.0; // yy
  smd_data_9[ilocal][8] = 1.0; // zz
}
