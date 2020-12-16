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
   Contributing author: Sam Cameron (University of Bristol)
------------------------------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "fix_selfpropel.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "comm.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixSelfPropel::FixSelfPropel(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  thermo_virial = 1;
  virial_flag = 1;
  
  if (narg != 4)
    error->all(FLERR,"Illegal fix selfpropel command");

  selfpropulsionforce = utils::numeric(FLERR,arg[3],false,lmp);


}

/* ---------------------------------------------------------------------- */

int FixSelfPropel::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

FixSelfPropel::~FixSelfPropel()
{

}


/* ---------------------------------------------------------------------- */

void FixSelfPropel::init()
{

  if (!atom->mu_flag)
    error->all(FLERR,"Fix selfpropel requires atom attributes mu");
}

void FixSelfPropel::setup(int vflag)
{
  post_force(vflag);
}


/* ----------------------------------------------------------------------
   apply self-propulsion force, stolen from MISC/fix_efield.cpp 
------------------------------------------------------------------------- */

void FixSelfPropel::post_force(int vflag)
{
  double **f = atom->f;
  int *mask = atom->mask;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  // energy and virial setup

  if (vflag) v_setup(vflag);
  else evflag = 0;


  // fsum[0] = "potential energy" for added force
  // fsum[123] = extra force added to atoms


  double **x = atom->x;
  double **mu = atom->mu;
  double fx,fy,fz;
  double v[6];

  // constant activity parameter

  double unwrap[3];

  // charge interactions

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {

      fx = selfpropulsionforce*mu[i][0];
      fy = selfpropulsionforce*mu[i][1];
      fz = selfpropulsionforce*mu[i][2];
      f[i][0] += fx;
      f[i][1] += fy;
      f[i][2] += fz;
      
      domain->unmap(x[i],image[i],unwrap);
      
      if (evflag) {
	v[0] = fx*unwrap[0];
	v[1] = fy*unwrap[1];
	v[2] = fz*unwrap[2];
	v[3] = fx*unwrap[1];
	v[4] = fx*unwrap[2];
	v[5] = fy*unwrap[2];
	v_tally(i, v);
      }
    }
}
