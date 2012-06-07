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
   Common parameters for the CMM coarse grained MD potentials.
   Contributing author: Axel Kohlmeyer <akohlmey@gmail.com>
------------------------------------------------------------------------- */

#include "cg_cmm_parms.h"

#include "string.h"

using namespace LAMMPS_NS;

/* static constant class members */
const char * const CGCMMParms::cg_type_list[] = {"none", "lj9_6", "lj12_4", "lj12_6"};
const double CGCMMParms::cg_prefact[] = {0.0, 6.75, 2.59807621135332, 4.0};
const double CGCMMParms::cg_pow1[] = {0.0, 9.0, 12.0, 12.0};
const double CGCMMParms::cg_pow2[] = {0.0, 6.0,  4.0,  6.0};

/* ---------------------------------------------------------------------- */

int CGCMMParms::find_cg_type(const char *label)
{
  for (int i=0; i < NUM_CG_TYPES; ++i) {
    if (strcmp(label,cg_type_list[i]) == 0) {
      return i;
    }
  }
  return CG_NOT_SET;
}

/* ---------------------------------------------------------------------- */
