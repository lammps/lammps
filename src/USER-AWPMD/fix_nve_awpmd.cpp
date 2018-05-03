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
   Contributing author: Ilya Valuev (JIHT, Moscow, Russia)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstring>
#include "fix_nve_awpmd.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"

#include "TCP/wpmd_split.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixNVEAwpmd::FixNVEAwpmd(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (!atom->wavepacket_flag)
    error->all(FLERR,"Fix nve/awpmd requires atom style wavepacket");
  //if (!atom->mass_type != 1)
   // error->all(FLERR,"Fix nve/awpmd requires per type mass");

  time_integrate = 1;
}

/* ---------------------------------------------------------------------- */

int FixNVEAwpmd::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= FINAL_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixNVEAwpmd::init()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;

  if (strstr(update->integrate_style,"respa"))
    step_respa = ((Respa *) update->integrate)->step;

  awpmd_pair=(PairAWPMDCut *)force->pair;
  awpmd_pair->wpmd->norm_needed=1;
}

/* ----------------------------------------------------------------------
   allow for only per-type  mass
------------------------------------------------------------------------- */

void FixNVEAwpmd::initial_integrate(int vflag)
{


  // update v,vr and x,radius of atoms in group

  double **x = atom->x;
  double *eradius = atom->eradius;
  double **v = atom->v;
  double *ervel = atom->ervel;
  double **f = atom->f;
  double *erforce = atom->erforce;
  double *vforce=atom->vforce;
  double *ervelforce=atom->ervelforce;

  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // x + dt * [v + 0.5 * dt * (f / m)];

  // simple Euler update
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      double dtfm = dtf / mass[type[i]];
      double dtfmr=dtfm;
      for(int j=0;j<3;j++){
        x[i][j] += dtv*vforce[3*i+j];
        v[i][j] += dtfm*f[i][j];
      }
      eradius[i]+= dtv*ervelforce[i];
      ervel[i] += dtfmr*erforce[i];
    }
  }

}

/* ---------------------------------------------------------------------- */

void FixNVEAwpmd::final_integrate(){}

/* ---------------------------------------------------------------------- */

void FixNVEAwpmd::initial_integrate_respa(int vflag, int ilevel, int iloop)
{
  dtv = step_respa[ilevel];
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;

  // innermost level - NVE update of v and x
  // all other levels - NVE update of v

  if (ilevel == 0) initial_integrate(vflag);
  else final_integrate();
}

/* ---------------------------------------------------------------------- */

void FixNVEAwpmd::final_integrate_respa(int ilevel, int iloop)
{
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;
  final_integrate();
}

/* ---------------------------------------------------------------------- */

void FixNVEAwpmd::reset_dt()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}
