/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Contributing Author: David Nicholson (MIT)
------------------------------------------------------------------------- */


#include <cmath>
#include <cstdlib>
#include <cstring>
#include "dump_cfg.h"
#include "atom.h"
#include "domain.h"
#include "modify.h"
#include "compute.h"
#include "fix.h"
#include "error.h"
#include "uef_utils.h"
#include "dump_cfg_uef.h"
#include "fix_nh_uef.h"

using namespace LAMMPS_NS;

#define UNWRAPEXPAND 10.0
#define ONEFIELD 32
#define DELTA 1048576

/* ----------------------------------------------------------------------
 * base method is mostly fine, just need to find the FixNHUef
 * ----------------------------------------------------------------------*/
void DumpCFGUef::init_style()
{
  DumpCFG::init_style();

  // check to make sure the other uef fix is on
  // borrowed from Pieter's nvt/sllod code
  int i=0;
  for (i=0; i<modify->nfix; i++)
  {
    if (strcmp(modify->fix[i]->style,"nvt/uef")==0)
      break;
    if (strcmp(modify->fix[i]->style,"npt/uef")==0)
      break;
  }
  if (i==modify->nfix)
    error->all(FLERR,"Can't use dump cfg/uef without defining a fix nvt/npt/uef");
  ifix_uef=i;
}

/* ----------------------------------------------------------------------
 * this is really the only difference between the base class and this one.
 * since the output is in scaled coordinates, changing the simulation box
 * edges to the flow frame will put coordinates in the flow frame too.
 * ----------------------------------------------------------------------*/

void DumpCFGUef::write_header(bigint n)
{
  // set scale factor used by AtomEye for CFG viz
  // default = 1.0
  // for peridynamics, set to pre-computed PD scale factor
  //   so PD particles mimic C atoms
  // for unwrapped coords, set to UNWRAPEXPAND (10.0)
  //   so molecules are not split across periodic box boundaries

  double box[3][3],rot[3][3];
  ((FixNHUef*) modify->fix[ifix_uef])->get_box(box);
  ((FixNHUef*) modify->fix[ifix_uef])->get_rot(rot);
  // rot goes from "lab frame" to "upper triangular frame"
  // it's transpose takes the simulation box to the flow frame
  for (int i=0;i<3;i++)
    for(int j=i+1;j<3;j++)
    {
      double t=rot[i][j];
      rot[i][j]=rot[j][i];
      rot[j][i]=t;
    }
  UEF_utils::mul_m2(rot,box);


  double scale = 1.0;
  if (atom->peri_flag) scale = atom->pdscale;
  else if (unwrapflag == 1) scale = UNWRAPEXPAND;

  fprintf(fp,"Number of particles = " BIGINT_FORMAT "\n",n);
  fprintf(fp,"A = %g Angstrom (basic length-scale)\n",scale);
  // in box[][] columns are cell edges
  // in H0, rows are cell edges
  fprintf(fp,"H0(1,1) = %g A\n",box[0][0]);
  fprintf(fp,"H0(1,2) = %g A\n",box[1][0]);
  fprintf(fp,"H0(1,3) = %g A\n",box[2][0]);
  fprintf(fp,"H0(2,1) = %g A\n",box[0][1]);
  fprintf(fp,"H0(2,2) = %g A\n",box[1][1]);
  fprintf(fp,"H0(2,3) = %g A\n",box[2][1]);
  fprintf(fp,"H0(3,1) = %g A\n",box[0][2]);
  fprintf(fp,"H0(3,2) = %g A\n",box[1][2]);
  fprintf(fp,"H0(3,3) = %g A\n",box[2][2]);
  fprintf(fp,".NO_VELOCITY.\n");
  fprintf(fp,"entry_count = %d\n",nfield-2);
  for (int i = 0; i < nfield-5; i++)
    fprintf(fp,"auxiliary[%d] = %s\n",i,auxname[i]);
}
