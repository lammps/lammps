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

/* ----------------------------------------------------------------------
   Contributing authors: Frances Mackay, Santtu Ollila, Colin Denniston (UWO)

   Based on fix_momentum,
   Contributing author: Naveen Michaud-Agrawal (Johns Hopkins U)
------------------------------------------------------------------------- */

#include "fix_lb_momentum.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix_lb_fluid.h"
#include "group.h"
#include "modify.h"

#include "latboltz_const.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixLbMomentum::FixLbMomentum(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), fix_lb_fluid(nullptr)
{
  if (narg < 4) error->all(FLERR, "Illegal fix lb/momentum command");
  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  if (nevery <= 0) error->all(FLERR, "Illegal fix lb/momentum command");

  linear = 1;
  xflag = 1;
  yflag = 1;
  zflag = 1;

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "linear") == 0) {
      if (iarg + 4 > narg) error->all(FLERR, "Illegal fix lb/momentum command");
      linear = 1;
      xflag = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      yflag = utils::inumeric(FLERR, arg[iarg + 2], false, lmp);
      zflag = utils::inumeric(FLERR, arg[iarg + 3], false, lmp);
      iarg += 4;
    } else
      error->all(FLERR, "Illegal fix lb/momentum command");
  }

  if (linear == 0) error->all(FLERR, "Illegal fix lb/momentum command");

  if (linear)
    if (xflag < 0 || xflag > 1 || yflag < 0 || yflag > 1 || zflag < 0 || zflag > 1)
      error->all(FLERR, "Illegal fix lb/momentum command");
}

/* ---------------------------------------------------------------------- */

int FixLbMomentum::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixLbMomentum::init()
{
  // warn if 0 atoms in group
  if ((count = group->count(igroup)) == 0)
    error->warning(FLERR, "Fix lb/momentum group has no atoms: Only fluid momentum affected");

  auto ifix = modify->get_fix_by_style("lb/fluid");
  if (ifix.size() > 0) fix_lb_fluid = (FixLbFluid *) ifix[0];

  count ? masstotal = group->mass(igroup) : 0;
}

/* ---------------------------------------------------------------------- */

void FixLbMomentum::end_of_step()
{
  // particle com velcity
  double vcmp[3] = {0, 0, 0};
  if (count) group->vcm(igroup, masstotal, vcmp);

  // total fluid mass and momentum.
  double masslb, momentumlb[3];
  fix_lb_fluid->calc_mass_momentum(masslb, momentumlb);

  //Calculate the center of mass velocity of the particles + fluid.
  double vcmall[3];
  vcmall[0] = (masstotal * vcmp[0] + momentumlb[0]) / (masslb + masstotal);
  vcmall[1] = (masstotal * vcmp[1] + momentumlb[1]) / (masslb + masstotal);
  vcmall[2] = (masstotal * vcmp[2] + momentumlb[2]) / (masslb + masstotal);

  // adjust velocities to zero net linear momentum
  // only adjust a component if flag is set

  //Subtract vcmall from the particles.
  if (count) {
    double **v = atom->v;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (xflag) v[i][0] -= vcmall[0];
        if (yflag) v[i][1] -= vcmall[1];
        if (zflag) v[i][2] -= vcmall[2];
      }
  }

  //Subtract vcmall from the fluid.

  // convert vcmall to lattice units used in lb/fluid
  double dx_lb = fix_lb_fluid->dx_lb;
  double dt_lb = fix_lb_fluid->dt_lb;
  vcmall[0] *= dt_lb / dx_lb;
  vcmall[1] *= dt_lb / dx_lb;
  vcmall[2] *= dt_lb / dx_lb;

  double ucmx, ucmy, ucmz;
  ucmx = xflag ? vcmall[0] : 0.0;
  ucmy = yflag ? vcmall[1] : 0.0;
  ucmz = zflag ? vcmall[2] : 0.0;

  int numvel = fix_lb_fluid->numvel;
  double etacov[14];    // This is for changes and only moments up to 13 are affected
  double rho;
  int subNbx = fix_lb_fluid->subNbx;
  int subNby = fix_lb_fluid->subNby;
  int subNbz = fix_lb_fluid->subNbz;
  double ***density_lb = fix_lb_fluid->density_lb;
  double ****f_lb = fix_lb_fluid->f_lb;
  double ****u_lb = fix_lb_fluid->u_lb;

  etacov[0] = etacov[10] = etacov[11] = etacov[12] = etacov[13] = 0.0;

  for (int i = 0; i < subNbx; i++)
    for (int j = 0; j < subNby; j++)
      for (int k = 0; k < subNbz; k++) {
        rho = density_lb[i][j][k];

        etacov[1] = rho * ucmx;
        etacov[2] = rho * ucmy;
        etacov[3] = rho * ucmz;
        etacov[4] = rho * (2.0 * u_lb[i][j][k][0] * ucmx - ucmx * ucmx);
        etacov[5] = rho * (2.0 * u_lb[i][j][k][1] * ucmy - ucmy * ucmy);
        etacov[6] = rho * (2.0 * u_lb[i][j][k][2] * ucmz - ucmz * ucmz);
        etacov[7] = rho * (u_lb[i][j][k][0] * ucmy + u_lb[i][j][k][1] * ucmx - ucmx * ucmy);
        etacov[8] = rho * (u_lb[i][j][k][1] * ucmz + u_lb[i][j][k][2] * ucmy - ucmy * ucmz);
        etacov[9] = rho * (u_lb[i][j][k][0] * ucmz + u_lb[i][j][k][2] * ucmx - ucmx * ucmz);

        if (numvel == 15) {
          etacov[13] = rho *
              (u_lb[i][j][k][0] * u_lb[i][j][k][1] * ucmz +
               u_lb[i][j][k][0] * ucmy * u_lb[i][j][k][2] - u_lb[i][j][k][0] * ucmy * ucmz +
               ucmx * u_lb[i][j][k][1] * u_lb[i][j][k][2] - ucmx * u_lb[i][j][k][1] * ucmz -
               ucmx * ucmy * u_lb[i][j][k][2] + ucmx * ucmy * ucmz);

          for (int l = 0; l < 15; l++)
            for (int ii = 1; ii < 14; ii++)
              f_lb[i][j][k][l] -= w_lb15[l] * mg_lb15[ii][l] * etacov[ii] * Ng_lb15[ii];
        } else    // 19-velocity model
          for (int l = 0; l < 19; l++)
            for (int ii = 1; ii < 14; ii++)
              f_lb[i][j][k][l] -= w_lb19[l] * mg_lb19[ii][l] * etacov[ii] * Ng_lb19[ii];

        if (xflag) u_lb[i][j][k][0] -= vcmall[0];
        if (yflag) u_lb[i][j][k][1] -= vcmall[1];
        if (zflag) u_lb[i][j][k][2] -= vcmall[2];
      }
}
