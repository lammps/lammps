/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Wengen Ouyang (Wuhan University)
   e-mail: w.g.ouyang at gmail dot com

   This is a full version of the potential described in
   [Ouyang et al., J. Chem. Theory Comput. 17, 7215-7223 (2021)]
------------------------------------------------------------------------- */

#include "pair_saip_metal.h"

#include "atom.h"
#include "citeme.h"
#include "error.h"
#include "force.h"
#include "interlayer_taper.h"
#include "neigh_list.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace InterLayer;

static const char cite_saip[] =
    "saip/metal potential: doi:10.1021/acs.jctc.1c00622\n\n"
    "@Article{Ouyang2021\n"
    " author = {W. Ouyang and O. Hod and R. Guerra},\n"
    " title = {Registry-Dependent Potential for Interfaces of Gold with Graphitic Systems},\n"
    " journal = {J.~Chem.\\ Theory Comput.},\n"
    " volume =  17,\n"
    " pages =   {7215--7223}\n"
    " year =    2021,\n"
    "}\n\n";

/* ---------------------------------------------------------------------- */

PairSAIPMETAL::PairSAIPMETAL(LAMMPS *lmp) : PairILPGrapheneHBN(lmp)
{
  variant = SAIP_METAL;
  single_enable = 0;
  if (lmp->citeme) lmp->citeme->add(cite_saip);
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairSAIPMETAL::settings(int narg, char **arg)
{
  if (narg < 1 || narg > 2) error->all(FLERR, "Illegal pair_style command");
  if (!utils::strmatch(force->pair_style, "^hybrid/overlay"))
    error->all(FLERR, "Pair style saip/metal must be used as sub-style with hybrid/overlay");

  cut_global = utils::numeric(FLERR, arg[0], false, lmp);
  if (narg == 2) tap_flag = utils::numeric(FLERR, arg[1], false, lmp);
}

/* ----------------------------------------------------------------------
   Repulsive forces and energy
------------------------------------------------------------------------- */

void PairSAIPMETAL::calc_FRep(int eflag, int /* vflag */)
{
  int i, j, ii, jj, inum, jnum, itype, jtype, k, kk;
  double prodnorm1, fkcx, fkcy, fkcz, filp;
  double xtmp, ytmp, ztmp, delx, dely, delz, evdwl, fpair, fpair1;
  double rsq, r, Rcut, rhosq1, exp0, exp1, Tap, dTap, Vilp;
  double frho1, Erep, fsum, rdsq1;
  int *ilist, *jlist, *numneigh, **firstneigh;
  int *ILP_neighs_i;

  evdwl = 0.0;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  double dprodnorm1[3] = {0.0, 0.0, 0.0};
  double fp1[3] = {0.0, 0.0, 0.0};
  double fprod1[3] = {0.0, 0.0, 0.0};
  double delki[3] = {0.0, 0.0, 0.0};
  double fk[3] = {0.0, 0.0, 0.0};

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  //calculate exp(-lambda*(r-z0))*[epsilon/2 + f(rho_ij)]
  // loop over neighbors of owned atoms
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = type[j];

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;

      // only include the interaction between different layers
      if (rsq < cutsq[itype][jtype] && atom->molecule[i] != atom->molecule[j]) {

        int iparam_ij = elem2param[map[itype]][map[jtype]];
        Param &p = params[iparam_ij];

        r = sqrt(rsq);
        // turn on/off taper function
        if (tap_flag) {
          Rcut = sqrt(cutsq[itype][jtype]);
          Tap = calc_Tap(r, Rcut);
          dTap = calc_dTap(r, Rcut);
        } else {
          Tap = 1.0;
          dTap = 0.0;
        }
        // for atoms in bulk materials
        if (strcmp(elements[map[itype]], "C") != 0 && strcmp(elements[map[itype]], "H") != 0 &&
            strcmp(elements[map[itype]], "B") != 0 && strcmp(elements[map[itype]], "N") != 0) {
          // Set ni along rij
          exp0 = exp(-p.lambda * (r - p.z0));
          frho1 = 1.0 * p.C;
          Erep = 0.5 * p.epsilon + frho1;
          Vilp = exp0 * Erep;
          // derivatives
          fpair = p.lambda * exp0 / r * Erep;
          filp = fpair * Tap - Vilp * dTap / r;
          fkcx = delx * filp;
          fkcy = dely * filp;
          fkcz = delz * filp;

          f[i][0] += fkcx;
          f[i][1] += fkcy;
          f[i][2] += fkcz;
          f[j][0] -= fkcx;
          f[j][1] -= fkcy;
          f[j][2] -= fkcz;

          if (eflag) pvector[1] += evdwl = Tap * Vilp;
          if (evflag) ev_tally(i, j, nlocal, newton_pair, evdwl, 0.0, filp, delx, dely, delz);
        } else {    // for atoms in 2D materials
          // Calculate the transverse distance
          prodnorm1 = normal[i][0] * delx + normal[i][1] * dely + normal[i][2] * delz;
          rhosq1 = rsq - prodnorm1 * prodnorm1;    // rho_ij
          rdsq1 = rhosq1 * p.delta2inv;            // (rho_ij/delta)^2

          // store exponents
          exp0 = exp(-p.lambda * (r - p.z0));
          exp1 = exp(-rdsq1);

          frho1 = exp1 * p.C;
          Erep = 0.5 * p.epsilon + frho1;
          Vilp = exp0 * Erep;

          // derivatives
          fpair = p.lambda * exp0 / r * Erep;
          fpair1 = 2.0 * exp0 * frho1 * p.delta2inv;
          fsum = fpair + fpair1;
          // derivatives of the product of rij and ni, the result is a vector
          dprodnorm1[0] =
              dnormdri[0][0][i] * delx + dnormdri[1][0][i] * dely + dnormdri[2][0][i] * delz;
          dprodnorm1[1] =
              dnormdri[0][1][i] * delx + dnormdri[1][1][i] * dely + dnormdri[2][1][i] * delz;
          dprodnorm1[2] =
              dnormdri[0][2][i] * delx + dnormdri[1][2][i] * dely + dnormdri[2][2][i] * delz;
          fp1[0] = prodnorm1 * normal[i][0] * fpair1;
          fp1[1] = prodnorm1 * normal[i][1] * fpair1;
          fp1[2] = prodnorm1 * normal[i][2] * fpair1;
          fprod1[0] = prodnorm1 * dprodnorm1[0] * fpair1;
          fprod1[1] = prodnorm1 * dprodnorm1[1] * fpair1;
          fprod1[2] = prodnorm1 * dprodnorm1[2] * fpair1;

          fkcx = (delx * fsum - fp1[0]) * Tap - Vilp * dTap * delx / r;
          fkcy = (dely * fsum - fp1[1]) * Tap - Vilp * dTap * dely / r;
          fkcz = (delz * fsum - fp1[2]) * Tap - Vilp * dTap * delz / r;

          f[i][0] += fkcx - fprod1[0] * Tap;
          f[i][1] += fkcy - fprod1[1] * Tap;
          f[i][2] += fkcz - fprod1[2] * Tap;
          f[j][0] -= fkcx;
          f[j][1] -= fkcy;
          f[j][2] -= fkcz;

          // calculate the forces acted on the neighbors of atom i from atom j
          ILP_neighs_i = ILP_firstneigh[i];
          for (kk = 0; kk < ILP_numneigh[i]; kk++) {
            k = ILP_neighs_i[kk];
            if (k == i) continue;
            // derivatives of the product of rij and ni respect to rk, k=0,1,2, where atom k is the neighbors of atom i
            dprodnorm1[0] = dnormal[0][0][kk][i] * delx + dnormal[1][0][kk][i] * dely +
                dnormal[2][0][kk][i] * delz;
            dprodnorm1[1] = dnormal[0][1][kk][i] * delx + dnormal[1][1][kk][i] * dely +
                dnormal[2][1][kk][i] * delz;
            dprodnorm1[2] = dnormal[0][2][kk][i] * delx + dnormal[1][2][kk][i] * dely +
                dnormal[2][2][kk][i] * delz;
            fk[0] = (-prodnorm1 * dprodnorm1[0] * fpair1) * Tap;
            fk[1] = (-prodnorm1 * dprodnorm1[1] * fpair1) * Tap;
            fk[2] = (-prodnorm1 * dprodnorm1[2] * fpair1) * Tap;
            f[k][0] += fk[0];
            f[k][1] += fk[1];
            f[k][2] += fk[2];
            delki[0] = x[k][0] - x[i][0];
            delki[1] = x[k][1] - x[i][1];
            delki[2] = x[k][2] - x[i][2];
            if (evflag)
              ev_tally_xyz(k, j, nlocal, newton_pair, 0.0, 0.0, fk[0], fk[1], fk[2], delki[0],
                           delki[1], delki[2]);
          }

          if (eflag) pvector[1] += evdwl = Tap * Vilp;
          if (evflag)
            ev_tally_xyz(i, j, nlocal, newton_pair, evdwl, 0.0, fkcx, fkcy, fkcz, delx, dely, delz);
        }    // end of speration atoms for bulk materials
      }      // end of rsq < cutoff
    }        // loop over jj
  }          // loop over ii
}
