// clang-format off
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

#include "pair_sph_idealgas.h"
#include <cmath>
#include "atom.h"
#include "force.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "domain.h"


using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairSPHIdealGas::PairSPHIdealGas(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;
}

/* ---------------------------------------------------------------------- */

PairSPHIdealGas::~PairSPHIdealGas() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(viscosity);
  }
}

/* ---------------------------------------------------------------------- */

void PairSPHIdealGas::compute(int eflag, int vflag) {
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, fpair;

  int *ilist, *jlist, *numneigh, **firstneigh;
  double vxtmp, vytmp, vztmp, imass, jmass, fi, fj, fvisc, h, ih, ihsq;
  double rsq, wfd, delVdotDelR, mu, deltaE, ci, cj;

  ev_init(eflag, vflag);

  double **v = atom->vest;
  double **x = atom->x;
  double **f = atom->f;
  double *rho = atom->rho;
  double *mass = atom->mass;
  double *desph = atom->desph;
  double *esph = atom->esph;
  double *drho = atom->drho;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    vxtmp = v[i][0];
    vytmp = v[i][1];
    vztmp = v[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    imass = mass[itype];

    fi = 0.4 * esph[i] / imass / rho[i]; // ideal gas EOS; this expression is fi = pressure / rho^2
    ci = sqrt(0.4*esph[i]/imass); // speed of sound with heat capacity ratio gamma=1.4

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];
      jmass = mass[jtype];

      if (rsq < cutsq[itype][jtype]) {
        h = cut[itype][jtype];
        ih = 1. / h;
        ihsq = ih * ih;

        wfd = h - sqrt(rsq);
        if (domain->dimension == 3) {
          // Lucy Kernel, 3d
          // Note that wfd, the derivative of the weight function with respect to r,
          // is lacking a factor of r.
          // The missing factor of r is recovered by
          // (1) using delV . delX instead of delV . (delX/r) and
          // (2) using f[i][0] += delx * fpair instead of f[i][0] += (delx/r) * fpair
          wfd = -25.066903536973515383e0 * wfd * wfd * ihsq * ihsq * ihsq * ih;
        } else {
          // Lucy Kernel, 2d
          wfd = -19.098593171027440292e0 * wfd * wfd * ihsq * ihsq * ihsq;
        }

        fj = 0.4 * esph[j] / jmass / rho[j];

        // dot product of velocity delta and distance vector
        delVdotDelR = delx * (vxtmp - v[j][0]) + dely * (vytmp - v[j][1])
            + delz * (vztmp - v[j][2]);

        // artificial viscosity (Monaghan 1992)
        if (delVdotDelR < 0.) {
          cj = sqrt(0.4*esph[j]/jmass);
          mu = h * delVdotDelR / (rsq + 0.01 * h * h);
          fvisc = -viscosity[itype][jtype] * (ci + cj) * mu / (rho[i] + rho[j]);
        } else {
          fvisc = 0.;
        }

        // total pair force & thermal energy increment
        fpair = -imass * jmass * (fi + fj + fvisc) * wfd;
        deltaE = -0.5 * fpair * delVdotDelR;

        f[i][0] += delx * fpair;
        f[i][1] += dely * fpair;
        f[i][2] += delz * fpair;

        // and change in density
        drho[i] += jmass * delVdotDelR * wfd;

        // change in thermal energy
        desph[i] += deltaE;

        if (newton_pair || j < nlocal) {
          f[j][0] -= delx * fpair;
          f[j][1] -= dely * fpair;
          f[j][2] -= delz * fpair;
          desph[j] += deltaE;
          drho[j] += imass * delVdotDelR * wfd;
        }

        if (evflag)
          ev_tally(i, j, nlocal, newton_pair, 0.0, 0.0, fpair, delx, dely, delz);

        if (evflag)
          ev_tally(i, j, nlocal, newton_pair, 0.0, 0.0, fpair, delx, dely,
              delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
 allocate all arrays
 ------------------------------------------------------------------------- */

void PairSPHIdealGas::allocate() {
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");

  memory->create(cut, n + 1, n + 1, "pair:cut");
  memory->create(viscosity, n + 1, n + 1, "pair:viscosity");
}

/* ----------------------------------------------------------------------
 global settings
 ------------------------------------------------------------------------- */

void PairSPHIdealGas::settings(int narg, char **/*arg*/) {
  if (narg != 0)
    error->all(FLERR,
        "Illegal number of arguments for pair_style sph/idealgas");
}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairSPHIdealGas::coeff(int narg, char **arg) {
  if (narg != 4)
    error->all(FLERR,"Incorrect number of args for pair_style sph/idealgas coefficients");
  if (!allocated)
    allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR,arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR,arg[1], 1, atom->ntypes, jlo, jhi, error);

  double viscosity_one = utils::numeric(FLERR,arg[2],false,lmp);
  double cut_one = utils::numeric(FLERR,arg[3],false,lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      viscosity[i][j] = viscosity_one;
      //printf("setting cut[%d][%d] = %f\n", i, j, cut_one);
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0)
    error->all(FLERR,"Incorrect args for pair sph/idealgas coefficients");
}

/* ----------------------------------------------------------------------
 init for one type pair i,j and corresponding j,i
 ------------------------------------------------------------------------- */

double PairSPHIdealGas::init_one(int i, int j) {

  if (setflag[i][j] == 0) {
      error->all(FLERR,"All pair sph/idealgas coeffs are not set");
  }

  cut[j][i] = cut[i][j];

  return cut[i][j];
}

/* ---------------------------------------------------------------------- */

double PairSPHIdealGas::single(int /*i*/, int /*j*/, int /*itype*/, int /*jtype*/,
    double /*rsq*/, double /*factor_coul*/, double /*factor_lj*/, double &fforce) {
  fforce = 0.0;

  return 0.0;
}
