/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Philipp Kloza (University of Cambridge)
                        pak37@cam.ac.uk
------------------------------------------------------------------------- */

#include "pair_mesocnt_viscous.h"

#include "atom.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

void PairMesoCNTViscous::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  mesolj();
  viscous();

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairMesoCNTViscous::coeff(int narg, char **arg)
{
  if (narg < 6) error->all(FLERR, "Incorrect args for pair coefficients");
  read_file(arg[2]);

  visc = utils::numeric(FLERR, arg[3], false, lmp);
  visc_cutoff = utils::numeric(FLERR, arg[4], false, lmp);
  visc_cutoffsq = visc_cutoff * visc_cutoff;

  nend_types = narg - 5;

  if (!allocated) allocate();

  // end atom types
  for (int i = 5; i < narg; i++) end_types[i - 5] = utils::inumeric(FLERR, arg[i], false, lmp);

  // units, eV to energy unit conversion
  ang = force->angstrom;
  ang_inv = 1.0 / ang;
  if (strcmp(update->unit_style, "lj") == 0)
    error->all(FLERR, "Pair style mesocnt does not support lj units");
  else if (strcmp(update->unit_style, "real") == 0)
    eunit = 23.06054966;
  else if (strcmp(update->unit_style, "metal") == 0)
    eunit = 1.0;
  else if (strcmp(update->unit_style, "si") == 0)
    eunit = 1.6021765e-19;
  else if (strcmp(update->unit_style, "cgs") == 0)
    eunit = 1.6021765e-12;
  else if (strcmp(update->unit_style, "electron") == 0)
    eunit = 3.674932248e-2;
  else if (strcmp(update->unit_style, "micro") == 0)
    eunit = 1.6021765e-4;
  else if (strcmp(update->unit_style, "nano") == 0)
    eunit = 1.6021765e2;
  else
    error->all(FLERR, "Pair style mesocnt does not recognize this units style");
  funit = eunit * ang_inv;

  // potential variables
  sig = sig_ang * ang;
  r = r_ang * ang;
  rsq = r * r;
  d = 2.0 * r;
  d_ang = 2.0 * r_ang;
  rc = 3.0 * sig;
  cutoff = rc + d;
  cutoffsq = cutoff * cutoff;
  cutoff_ang = cutoff * ang_inv;
  cutoffsq_ang = cutoff_ang * cutoff_ang;
  comega = 0.275 * (1.0 - 1.0 / (1.0 + 0.59 * r_ang));
  ctheta = 0.35 + 0.0226 * (r_ang - 6.785);

  // compute spline coefficients
  spline_coeff(uinf_data, uinf_coeff, delh_uinf, uinf_points);
  spline_coeff(gamma_data, gamma_coeff, delh_gamma, gamma_points);
  spline_coeff(phi_data, phi_coeff, delh_phi, delpsi_phi, phi_points);
  spline_coeff(usemi_data, usemi_coeff, delh_usemi, delxi_usemi, usemi_points);

  memory->destroy(uinf_data);
  memory->destroy(gamma_data);
  memory->destroy(phi_data);
  memory->destroy(usemi_data);

  int ntypes = atom->ntypes;
  for (int i = 1; i <= ntypes; i++)
    for (int j = i; j <= ntypes; j++) setflag[i][j] = 1;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairMesoCNTViscous::init_one(int /* i */, int /* j */)
{
  return (cutoff > visc_cutoff) ? cutoff : visc_cutoff;
}


/* ----------------------------------------------------------------------
   calculates viscous damping for all interacting segments within
   visc_cutoff
------------------------------------------------------------------------- */

void PairMesoCNTViscous::viscous()
{
  int i, j, ii, jj, inum, jnum;
  double xtmp, ytmp, ztmp, delx, dely, delz;
  double vxtmp, vytmp, vztmp;
  double fx, fy, fz;
  double rsq;
  int *ilist, *jlist, *numneigh, **firstneigh;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *type = atom->type;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  
  for (ii = 0; ii < inum; i++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    vxtmp = v[i][0];
    vytmp = v[i][1];
    vztmp = v[i][2];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      
      if (rsq < visc_cutoffsq) {
        fx = visc * (v[j][0] - vxtmp);
        fy = visc * (v[j][1] - vytmp);
        fz = visc * (v[j][2] - vztmp);
        
        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;
        
        f[j][0] -= fx;
        f[j][1] -= fy;
        f[j][2] -= fz;
      }
    }
  }
}
