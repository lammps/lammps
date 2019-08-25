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
   Contributing author:
      Morteza Jalalvand (IASBS)  jalalvand.m AT gmail.com

    references: Espanol and Revenga, Phys Rev E 67, 026705 (2003)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdlib>
#include "pair_sdpd_taitwater_isothermal.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "domain.h"
#include "update.h"
#ifndef USE_ZEST
#include "random_mars.h"
#endif

using namespace LAMMPS_NS;

static const double sqrt_2_inv = std::sqrt(0.5);

/* ---------------------------------------------------------------------- */

PairSDPDTaitwaterIsothermal::PairSDPDTaitwaterIsothermal (LAMMPS *lmp)
: Pair (lmp) {
  restartinfo = 0;
  single_enable =0;
}

/* ---------------------------------------------------------------------- */

PairSDPDTaitwaterIsothermal::~PairSDPDTaitwaterIsothermal () {
  if (allocated) {
    memory->destroy (setflag);
    memory->destroy (cutsq);

    memory->destroy (cut);
    memory->destroy (rho0);
    memory->destroy (soundspeed);
    memory->destroy (B);
  }
}

/* ---------------------------------------------------------------------- */

void PairSDPDTaitwaterIsothermal::compute (int eflag, int vflag) {
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, fpair;

  int *ilist, *jlist, *numneigh, **firstneigh;
  double vxtmp, vytmp, vztmp, imass, jmass, fi, fj, fvisc;
  double h, ih, ihsq, velx, vely, velz;
  double rsq, tmp, wfd, delVdotDelR;
  double prefactor, wiener[3][3], f_random[3];

  ev_init(eflag, vflag);

  double **v = atom->vest;
  double **x = atom->x;
  double **f = atom->f;
  double *rho = atom->rho;
  double *mass = atom->mass;
  double *drho = atom->drho;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  int dimension = domain->dimension;
  double dtinv = 1.0 / update->dt;
  double kBoltzmann = force->boltz;

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

    // compute pressure of atom i with Tait EOS
    tmp = rho[i] / rho0[itype];
    fi = tmp * tmp * tmp;
    fi = B[itype] * (fi * fi * tmp - 1.0) / (rho[i] * rho[i]);

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
        ih = 1.0 / h;
        ihsq = ih * ih;

        double r = sqrt (rsq);
        wfd = h - r;
        if (dimension == 3) {
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

        // compute pressure  of atom j with Tait EOS
        tmp = rho[j] / rho0[jtype];
        fj = tmp * tmp * tmp;
        fj = B[jtype] * (fj * fj * tmp - 1.0) / (rho[j] * rho[j]);

        velx=vxtmp - v[j][0];
        vely=vytmp - v[j][1];
        velz=vztmp - v[j][2];

        // dot product of velocity delta and distance vector
        delVdotDelR = delx * velx + dely * vely + delz * velz;

        // Espanol Viscosity (Espanol, 2003)

        fvisc = (5. / 3.) * viscosity * imass * jmass * wfd / (rho[i]*rho[j]);

        // total pair force
        fpair = -imass * jmass * (fi + fj) * wfd;

        // random force calculation
        // independent increments of a Wiener process matrix
#ifdef USE_ZEST
        wiener[0][0] = gaussian (generator);
        wiener[1][1] = gaussian (generator);
        wiener[2][2] = gaussian (generator);

        wiener[0][1] = wiener[1][0] = sqrt_2_inv * gaussian (generator);
        wiener[0][2] = wiener[2][0] = sqrt_2_inv * gaussian (generator);
        wiener[1][2] = wiener[2][1] = sqrt_2_inv * gaussian (generator);
#else
        wiener[0][0] = random->gaussian ();
        wiener[1][1] = random->gaussian ();
        wiener[2][2] = random->gaussian ();

        wiener[0][1] = wiener[1][0] = sqrt_2_inv * random->gaussian ();
        wiener[0][2] = wiener[2][0] = sqrt_2_inv * random->gaussian ();
        wiener[1][2] = wiener[2][1] = sqrt_2_inv * random->gaussian ();
#endif

        prefactor = sqrt (-4. * kBoltzmann*temperature * fvisc * dtinv) / r;

        f_random[0] = prefactor * (wiener[0][0]*delx + wiener[0][1]*dely + wiener[0][2]*delz);
        f_random[1] = prefactor * (wiener[1][0]*delx + wiener[1][1]*dely + wiener[1][2]*delz);
        f_random[2] = prefactor * (wiener[2][0]*delx + wiener[2][1]*dely + wiener[2][2]*delz);

        f[i][0] += delx * fpair + (velx + delx * delVdotDelR / rsq) * fvisc + f_random[0];
        f[i][1] += dely * fpair + (vely + dely * delVdotDelR / rsq) * fvisc + f_random[1];
        f[i][2] += delz * fpair + (velz + delz * delVdotDelR / rsq) * fvisc + f_random[2];

        // and change in density
        drho[i] += jmass * delVdotDelR * wfd;

        if (newton_pair || j < nlocal) {
          f[j][0] -= delx * fpair + (velx + delx * delVdotDelR / rsq) * fvisc + f_random[0];
          f[j][1] -= dely * fpair + (vely + dely * delVdotDelR / rsq) * fvisc + f_random[1];
          f[j][2] -= delz * fpair + (velz + delz * delVdotDelR / rsq) * fvisc + f_random[2];
          drho[j] += imass * delVdotDelR * wfd;
        }

        if (evflag)
          ev_tally (i, j, nlocal, newton_pair, 0.0, 0.0, fpair, delx, dely, delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute ();
}

/* ----------------------------------------------------------------------
 allocate all arrays
 ------------------------------------------------------------------------- */

void PairSDPDTaitwaterIsothermal::allocate () {
  allocated = 1;
  int n = atom->ntypes;

  memory->create (setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create (cutsq, n + 1, n + 1, "pair:cutsq");

  memory->create (rho0, n + 1, "pair:rho0");
  memory->create (soundspeed, n + 1, "pair:soundspeed");
  memory->create (B, n + 1, "pair:B");
  memory->create (cut, n + 1, n + 1, "pair:cut");
}

/* ----------------------------------------------------------------------
 global settings
 ------------------------------------------------------------------------- */

void PairSDPDTaitwaterIsothermal::settings (int narg, char **arg) {
  if (narg != 2 && narg != 3)
    error->all (FLERR, "Illegal number of arguments for "
                "pair_style sdpd/taitwater/morris/isothermal");

  temperature = force->numeric (FLERR, arg[0]);
  viscosity = force->numeric (FLERR, arg[1]);

  if (temperature <= 0) error->all (FLERR, "Temperature must be positive");
  if (viscosity <= 0) error->all (FLERR, "Viscosity must be positive");

  // seed is immune to underflow/overflow because it is unsigned
  seed = comm->nprocs + comm->me + atom->nlocal;
  if (narg == 3) seed += force->inumeric (FLERR, arg[2]);
#ifdef USE_ZEST
  generator.seed (seed);
#else
  random = new RanMars (lmp, seed);
#endif
}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairSDPDTaitwaterIsothermal::coeff (int narg, char **arg) {
  if (narg != 5)
    error->all (FLERR, "Incorrect args for pair_style "
                "sph/taitwater/morris coefficients");

  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  force->bounds (FLERR, arg[0], atom->ntypes, ilo, ihi);
  force->bounds (FLERR, arg[1], atom->ntypes, jlo, jhi);

  double rho0_one = force->numeric (FLERR,arg[2]);
  double soundspeed_one = force->numeric (FLERR,arg[3]);
  double cut_one = force->numeric (FLERR,arg[4]);
  double B_one = soundspeed_one * soundspeed_one * rho0_one / 7.0;

  if (rho0_one <= 0) error->all (FLERR, "Density must be positive");
  if (soundspeed_one <= 0) error->all (FLERR, "Sound speed must be positive");
  if (cut_one <= 0) error->all (FLERR, "Cutoff must be positive");

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    rho0[i] = rho0_one;
    soundspeed[i] = soundspeed_one;
    B[i] = B_one;
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      cut[i][j] = cut_one;

      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0)
    error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
 init for one type pair i,j and corresponding j,i
 ------------------------------------------------------------------------- */

double PairSDPDTaitwaterIsothermal::init_one (int i, int j) {
  if (setflag[i][j] == 0)
    error->all(FLERR,"Not all pair sph/taitwater/morris coeffs are set");

  cut[j][i] = cut[i][j];

  return cut[i][j];
}

