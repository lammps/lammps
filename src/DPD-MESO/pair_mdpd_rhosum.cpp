// clang-format off
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

/*-----------------------------------------------------------------------
   This is a Child Class PairMDPD for taking care of density summation
   before the force calculation.
   The code uses 3D Lucy kernel, it can be modified for other kernels.

   Contributing author: Zhen Li (Clemson University)
------------------------------------------------------------------------- */

#include "pair_mdpd_rhosum.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"

#include <cmath>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairMDPDRhoSum::PairMDPDRhoSum(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;

  // set comm size needed by this Pair

  comm_forward = 1;
  first = 1;
}

/* ---------------------------------------------------------------------- */

PairMDPDRhoSum::~PairMDPDRhoSum() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
  }
}

/* ----------------------------------------------------------------------
 init specific to this pair style
 ------------------------------------------------------------------------- */

void PairMDPDRhoSum::init_style()
{
  if (!atom->rho_flag)
    error->all(FLERR,"Pair style mdpd/rhosum requires atom attribute rho");

  // need a full neighbor list
  neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ---------------------------------------------------------------------- */

void PairMDPDRhoSum::compute(int eflag, int vflag) {
  int i, j, ii, jj, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz;
  double r, rsq, h, ih, ihsq;
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
    for (i = 1; i <= atom->ntypes; i++)
    for (j = 1; i <= atom->ntypes; i++)
    if (cutsq[i][j] > 0.0)
    if (!setflag[i][i] || !setflag[j][j])
    if (comm->me == 0)
      printf("mDPD particle types %d and %d interact, but not all of their single particle properties are set.\n", i, j);

    first = 0;
  }


  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // recompute density
  // we use a full neighborlist here

  // initialize density with zero self-contribution,
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];

    h = cut[itype][itype];
    // Lucy kernel, 3d
    wf = 2.0889086280811262819e0 / (h * h * h);
    rho[i] = 0;
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

        // Lucy kernel, 3d
        r = sqrt(rsq);
        wf = (h - r) * ihsq;
        wf =  2.0889086280811262819e0 * (h + 3. * r) * wf * wf * wf * ih;
        rho[i] += mass[jtype]*wf;
      }
    }
  }

  // communicate densities
  comm->forward_comm(this);
}

/* ----------------------------------------------------------------------
 allocate all arrays
 ------------------------------------------------------------------------- */

void PairMDPDRhoSum::allocate() {
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

void PairMDPDRhoSum::settings(int narg, char **/*arg*/) {
  if (narg != 0)
    error->all(FLERR,"Illegal number of setting arguments for pair_style mdpd/rhosum");
}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairMDPDRhoSum::coeff(int narg, char **arg) {
  if (narg != 3)
    error->all(FLERR,"Incorrect number of args for mdpd/rhosum coefficients");
  if (!allocated)
    allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR,arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR,arg[1], 1, atom->ntypes, jlo, jhi, error);

  double cut_one = utils::numeric(FLERR,arg[2],false,lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
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

double PairMDPDRhoSum::init_one(int i, int j) {
  if (setflag[i][j] == 0) {
    error->all(FLERR,"All pair mdpd/rhosum coeffs are not set");
  }

  cut[j][i] = cut[i][j];

  return cut[i][j];
}

/* ---------------------------------------------------------------------- */

double PairMDPDRhoSum::single(int /*i*/, int /*j*/, int /*itype*/, int /*jtype*/, double /*rsq*/,
    double /*factor_coul*/, double /*factor_lj*/, double &fforce) {
  fforce = 0.0;

  return 0.0;
}

/* ---------------------------------------------------------------------- */

int PairMDPDRhoSum::pack_forward_comm(int n, int *list, double *buf,
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

void PairMDPDRhoSum::unpack_forward_comm(int n, int first, double *buf) {
  int i, m, last;
  double *rho = atom->rho;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++)
    rho[i] = buf[m++];
}
