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

#include "fix_enforce2d.h"
#include <cstring>
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "modify.h"
#include "respa.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixEnforce2D::FixEnforce2D(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  flist(NULL)
{
  if (narg != 3) error->all(FLERR,"Illegal fix enforce2d command");

  nfixlist = 0;
}

/* ---------------------------------------------------------------------- */

FixEnforce2D::~FixEnforce2D()
{
  if (copymode) return;

  delete [] flist;
}

/* ---------------------------------------------------------------------- */

int FixEnforce2D::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixEnforce2D::init()
{
  if (domain->dimension == 3)
    error->all(FLERR,"Cannot use fix enforce2d with 3d simulation");

  // list of fixes with enforce2d methods

  nfixlist = 0;
  for (int i = 0; i < modify->nfix; i++)
    if (modify->fix[i]->enforce2d_flag) nfixlist++;

  if (nfixlist) {
    int myindex = -1;
    delete [] flist;
    flist = new Fix*[nfixlist];
    nfixlist = 0;
    for (int i = 0; i < modify->nfix; i++) {
      if (modify->fix[i]->enforce2d_flag) {
        if (myindex < 0)
          flist[nfixlist++] = modify->fix[i];
        else {
          char msg[256];
          snprintf(msg,256,"Fix enforce2d must be defined after fix %s",modify->fix[i]->style);
          error->all(FLERR,msg);
        }
      }
      if (modify->fix[i] == this) myindex = i;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixEnforce2D::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    int nlevels_respa = ((Respa *) update->integrate)->nlevels;
    for (int ilevel = 0; ilevel < nlevels_respa; ilevel++) {
      ((Respa *) update->integrate)->copy_flevel_f(ilevel);
      post_force_respa(vflag,ilevel,0);
      ((Respa *) update->integrate)->copy_f_flevel(ilevel);
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixEnforce2D::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixEnforce2D::post_force(int /*vflag*/)
{
  double **v = atom->v;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      v[i][2] = 0.0;
      f[i][2] = 0.0;
    }

  // for systems with omega/angmom/torque, zero x and y components

  if (atom->omega_flag) {
    double **omega = atom->omega;
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        omega[i][0] = 0.0;
        omega[i][1] = 0.0;
      }
  }

  if (atom->angmom_flag) {
    double **angmom = atom->angmom;
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        angmom[i][0] = 0.0;
        angmom[i][1] = 0.0;
      }
  }

  if (atom->torque_flag) {
    double **torque = atom->torque;
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        torque[i][0] = 0.0;
        torque[i][1] = 0.0;
      }
  }

  // invoke other fixes that enforce 2d
  // fix rigid variants

  for (int m = 0; m < nfixlist; m++)
    flist[m]->enforce2d();
}

/* ---------------------------------------------------------------------- */

void FixEnforce2D::post_force_respa(int vflag, int /*ilevel*/, int /*iloop*/)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixEnforce2D::min_post_force(int vflag)
{
  post_force(vflag);
}
