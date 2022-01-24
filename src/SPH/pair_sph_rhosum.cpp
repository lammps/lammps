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

#include "pair_sph_rhosum.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairSPHRhoSum::PairSPHRhoSum(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;

  // set comm size needed by this Pair

  comm_forward = 1;
  first = 1;
}

/* ---------------------------------------------------------------------- */

PairSPHRhoSum::~PairSPHRhoSum() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
  }
}

/* ----------------------------------------------------------------------
 init specific to this pair style
 ------------------------------------------------------------------------- */

void PairSPHRhoSum::init_style() {
  // need a full neighbor list
  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ---------------------------------------------------------------------- */

void PairSPHRhoSum::compute(int eflag, int vflag) {
  int i, j, ii, jj, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz;
  double rsq, imass, h, ih, ihsq;
  int *jlist;
  double wf;
  // neighbor list variables
  int inum, *ilist, *numneigh, **firstneigh;

  ev_init(eflag, vflag);

  double **x = atom->x;
  double *rho = atom->rho;
  int *type = atom->type;
  double *mass = atom->mass;

  // check consistency of pair coefficients

  if (first) {
    for (i = 1; i <= atom->ntypes; i++) {
      for (j = 1; i <= atom->ntypes; i++) {
        if (cutsq[i][j] > 0.0) {
          if (!setflag[i][i] || !setflag[j][j]) {
            if (comm->me == 0) {
              printf(
                  "SPH particle types %d and %d interact, but not all of their single particle properties are set.\n",
                  i, j);
            }
          }
        }
      }
    }
    first = 0;
  }

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // recompute density
  // we use a full neighborlist here

  if (nstep != 0) {
    if ((update->ntimestep % nstep) == 0) {

      // initialize density with self-contribution,
      for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];
        itype = type[i];
        imass = mass[itype];

        h = cut[itype][itype];
        if (domain->dimension == 3) {
          /*
          // Lucy kernel, 3d
          wf = 2.0889086280811262819e0 / (h * h * h);
          */

          // quadric kernel, 3d
          wf = 2.1541870227086614782 / (h * h * h);
        } else {
          /*
          // Lucy kernel, 2d
          wf = 1.5915494309189533576e0 / (h * h);
          */

          // quadric kernel, 2d
          wf = 1.5915494309189533576e0 / (h * h);
        }

        rho[i] = imass * wf;
      }

      // add density at each atom via kernel function overlap
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

          if (rsq < cutsq[itype][jtype]) {
            h = cut[itype][jtype];
            ih = 1.0 / h;
            ihsq = ih * ih;

            if (domain->dimension == 3) {
              /*
              // Lucy kernel, 3d
              r = sqrt(rsq);
              wf = (h - r) * ihsq;
              wf =  2.0889086280811262819e0 * (h + 3. * r) * wf * wf * wf * ih;
              */

              // quadric kernel, 3d
              wf = 1.0 - rsq * ihsq;
              wf = wf * wf;
              wf = wf * wf;
              wf = 2.1541870227086614782e0 * wf * ihsq * ih;
            } else {
              // Lucy kernel, 2d
              //r = sqrt(rsq);
              //wf = (h - r) * ihsq;
              //wf = 1.5915494309189533576e0 * (h + 3. * r) * wf * wf * wf;

              // quadric kernel, 2d
              wf = 1.0 - rsq * ihsq;
              wf = wf * wf;
              wf = wf * wf;
              wf = 1.5915494309189533576e0 * wf * ihsq;
            }

            rho[i] += mass[jtype] * wf;
          }

        }
      }
    }
  }

  // communicate densities
  comm->forward_comm_pair(this);
}

/* ----------------------------------------------------------------------
 allocate all arrays
 ------------------------------------------------------------------------- */

void PairSPHRhoSum::allocate() {
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");

  memory->create(cut, n + 1, n + 1, "pair:cut");
}

/* ----------------------------------------------------------------------
 global settings
 ------------------------------------------------------------------------- */

void PairSPHRhoSum::settings(int narg, char **arg) {
  if (narg != 1)
    error->all(FLERR,
        "Illegal number of arguments for pair_style sph/rhosum");
  nstep = utils::inumeric(FLERR,arg[0],false,lmp);
}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairSPHRhoSum::coeff(int narg, char **arg) {
  if (narg != 3)
    error->all(FLERR,"Incorrect number of args for sph/rhosum coefficients");
  if (!allocated)
    allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR,arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR,arg[1], 1, atom->ntypes, jlo, jhi, error);

  double cut_one = utils::numeric(FLERR,arg[2],false,lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      //printf("setting cut[%d][%d] = %f\n", i, j, cut_one);
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

double PairSPHRhoSum::init_one(int i, int j) {
  if (setflag[i][j] == 0) {
    error->all(FLERR,"All pair sph/rhosum coeffs are not set");
  }

  cut[j][i] = cut[i][j];

  return cut[i][j];
}

/* ---------------------------------------------------------------------- */

double PairSPHRhoSum::single(int /*i*/, int /*j*/, int /*itype*/, int /*jtype*/, double /*rsq*/,
    double /*factor_coul*/, double /*factor_lj*/, double &fforce) {
  fforce = 0.0;

  return 0.0;
}

/* ---------------------------------------------------------------------- */

int PairSPHRhoSum::pack_forward_comm(int n, int *list, double *buf,
                                     int /*pbc_flag*/, int * /*pbc*/) {
  int i, j, m;
  double *rho = atom->rho;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = rho[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairSPHRhoSum::unpack_forward_comm(int n, int first, double *buf) {
  int i, m, last;
  double *rho = atom->rho;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++)
    rho[i] = buf[m++];
}
