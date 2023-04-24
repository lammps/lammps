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
   Contributing authors: Xiaowng Zhou (Sandia)
------------------------------------------------------------------------- */

#include "pair_eam_he.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairEAMHE::PairEAMHE(LAMMPS *lmp) : PairEAM(lmp), PairEAMFS(lmp)
{
  he_flag = 1;
}

void PairEAMHE::compute(int eflag, int vflag)
{
  int i, j, ii, jj, m, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, evdwl, fpair;
  double rsq, r, p, rhoip, rhojp, z2, z2p, recip, phip, psip, phi;
  double *coeff;
  int *ilist, *jlist, *numneigh, **firstneigh;

  evdwl = 0.0;
  ev_init(eflag, vflag);

  // grow energy and fp arrays if necessary
  // need to be atom->nmax in length

  if (atom->nmax > nmax) {
    memory->destroy(rho);
    memory->destroy(fp);
    memory->destroy(numforce);
    nmax = atom->nmax;
    memory->create(rho, nmax, "pair:rho");
    memory->create(fp, nmax, "pair:fp");
    memory->create(numforce, nmax, "pair:numforce");
  }

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // zero out density

  if (newton_pair) {
    for (i = 0; i < nall; i++) rho[i] = 0.0;
  } else
    for (i = 0; i < nlocal; i++) rho[i] = 0.0;

  // rho = density at each atom
  // loop over neighbors of my atoms

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

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;

      if (rsq < cutforcesq) {
        jtype = type[j];
        p = sqrt(rsq) * rdr + 1.0;
        m = static_cast<int>(p);
        m = MIN(m, nr - 1);
        p -= m;
        p = MIN(p, 1.0);
        coeff = rhor_spline[type2rhor[jtype][itype]][m];
        rho[i] += ((coeff[3] * p + coeff[4]) * p + coeff[5]) * p + coeff[6];
        if (newton_pair || j < nlocal) {
          coeff = rhor_spline[type2rhor[itype][jtype]][m];
          rho[j] += ((coeff[3] * p + coeff[4]) * p + coeff[5]) * p + coeff[6];
        }
      }
    }
  }

  // communicate and sum densities

  if (newton_pair) comm->reverse_comm(this);

  // fp = derivative of embedding energy at each atom
  // phi = embedding energy at each atom
  // if rho > rhomax or rho < rhomin (i.e., table is exceeded),
  //   add linear term to conserve energy

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    p = (rho[i] - rhomin) * rdrho + 1.0;
    m = static_cast<int>(p);
    if (m < 2) {
      m = 2;
    } else if (m > nrho - 1) {
      m = nrho - 1;
    }
    p -= m;
    if (p < -1.0) {
      p = -1.0;
    } else if (p > 1.0) {
      p = 1.0;
    }
    coeff = frho_spline[type2frho[type[i]]][m];
    fp[i] = (coeff[0] * p + coeff[1]) * p + coeff[2];
    if (eflag) {
      phi = ((coeff[3] * p + coeff[4]) * p + coeff[5]) * p + coeff[6];
      if (rho[i] < rhomin) {
        phi += fp[i] * (rho[i] - rhomin);
      } else if (rho[i] > rhomax) {
        phi += fp[i] * (rho[i] - rhomax);
      }
      phi *= scale[type[i]][type[i]];
      if (eflag_global) eng_vdwl += phi;
      if (eflag_atom) eatom[i] += phi;
    }
  }

  // communicate derivative of embedding function

  comm->forward_comm(this);
  embedstep = update->ntimestep;

  // compute forces on each atom
  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];

    jlist = firstneigh[i];
    jnum = numneigh[i];
    numforce[i] = 0;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;

      if (rsq < cutforcesq) {
        ++numforce[i];
        jtype = type[j];
        r = sqrt(rsq);
        p = r * rdr + 1.0;
        m = static_cast<int>(p);
        m = MIN(m, nr - 1);
        p -= m;
        p = MIN(p, 1.0);

        // rhoip = derivative of (density at atom j due to atom i)
        // rhojp = derivative of (density at atom i due to atom j)
        // phi = pair potential energy
        // phip = phi'
        // z2 = phi * r
        // z2p = (phi * r)' = (phi' r) + phi
        // psip needs both fp[i] and fp[j] terms since r_ij appears in two
        //   terms of embed eng: Fi(sum rho_ij) and Fj(sum rho_ji)
        //   hence embed' = Fi(sum rho_ij) rhojp + Fj(sum rho_ji) rhoip
        // scale factor can be applied by thermodynamic integration

        coeff = rhor_spline[type2rhor[itype][jtype]][m];
        rhoip = (coeff[0] * p + coeff[1]) * p + coeff[2];
        coeff = rhor_spline[type2rhor[jtype][itype]][m];
        rhojp = (coeff[0] * p + coeff[1]) * p + coeff[2];
        coeff = z2r_spline[type2z2r[itype][jtype]][m];
        z2p = (coeff[0] * p + coeff[1]) * p + coeff[2];
        z2 = ((coeff[3] * p + coeff[4]) * p + coeff[5]) * p + coeff[6];

        recip = 1.0 / r;
        phi = z2 * recip;
        phip = z2p * recip - phi * recip;
        psip = fp[i] * rhojp + fp[j] * rhoip + phip;
        fpair = -scale[itype][jtype] * psip * recip;

        f[i][0] += delx * fpair;
        f[i][1] += dely * fpair;
        f[i][2] += delz * fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx * fpair;
          f[j][1] -= dely * fpair;
          f[j][2] -= delz * fpair;
        }

        if (eflag) evdwl = scale[itype][jtype] * phi;
        if (evflag) ev_tally(i, j, nlocal, newton_pair, evdwl, 0.0, fpair, delx, dely, delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}
