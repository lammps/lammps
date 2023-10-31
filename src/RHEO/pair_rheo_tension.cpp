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
   Contributing authors:
   Joel Clemmer (SNL)
----------------------------------------------------------------------- */

#include "pair_rheo_tension.h"

#include "atom.h"
#include "comm.h"
#include "compute_rheo_kernel.h"
#include "domain.h"
#include "error.h"
#include "fix_rheo.h"
#include "force.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "update.h"
#include "utils.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace RHEO_NS;
using namespace MathExtra;

static constexpr double EPSILON = 1e-2;

/* ---------------------------------------------------------------------- */

PairRHEOTension::PairRHEOTension(LAMMPS *lmp) :
  Pair(lmp), compute_kernel(nullptr), fix_rheo(nullptr)
{
  restartinfo = 0;
  single_enable = 0;

  comm_forward = 3;
  comm_reverse = 3;
}

/* ---------------------------------------------------------------------- */

PairRHEOTension::~PairRHEOTension()
{
  // Remove custom property if it exists
  int tmp1, tmp2, index;

  index = atom->find_custom("rheo_c_tension", tmp1, tmp2);
  if (index != -1) atom->remove_custom(index, 1, 0);

  index = atom->find_custom("rheo_n_tension", tmp1, tmp2);
  if (index != -1) atom->remove_custom(index, 1, 3);

  if (allocated) {
    memory->destroy(alpha);
    memory->destroy(setflag);
    memory->destroy(cutsq);
  }
}

/* ---------------------------------------------------------------------- */

void PairRHEOTension::compute(int eflag, int vflag)
{
  int i, j, a, b, ii, jj, inum, jnum, itype, jtype;
  int fluidi, fluidj;
  double xtmp, ytmp, ztmp, w, wp;
  double rhoi, rhoj, voli, volj;
  double *dWij, *dWji;
  double dx[3], ft[3];

  int *ilist, *jlist, *numneigh, **firstneigh;
  double imass, jmass, rsq, r, rinv;

  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  int dim = domain->dimension;

  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  double *rho = atom->rho;
  double *mass = atom->mass;
  double *special_lj = force->special_lj;
  int *type = atom->type;
  int *status = atom->status;
  tagint *tag = atom->tag;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
/*
  int nmax = atom->nmax;
  if (nmax_store <= nmax) {
    memory->grow(ct, nmax, "atom:rheo_c_tension");
    memory->grow(nnt_tension, nmax, 3, "atom:rheo_n_tension");
    nmax_store = atom->nmax;
  }

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    imass = mass[itype];
    rhoi = rho[i];
    voli = imass / rhoi;
    fluidi = !(status[i] & PHASECHECK);

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      dx[0] = xtmp - x[j][0];
      dx[1] = ytmp - x[j][1];
      dx[2] = ztmp - x[j][2];
      rsq = lensq3(dx);
      jtype = type[j];

      if (rsq > hsq) continue;

      r = sqrt(rsq);
      rinv = 1 / r;

      jmass = mass[jtype];
      rhoj = rho[j];
      volj = jmass / rhoj;
      fluidj = !(status[j] & PHASECHECK);

      wp = compute_kernel->calc_dw(i, j, dx[0], dx[1], dx[2],r);
      dWij = compute_kernel->dWij;
      dWji = compute_kernel->dWji;

      f[i][0] += ft[0];
      f[i][1] += ft[1];
      f[i][2] += ft[2];

      if (evflag) // Does not account for unbalanced forces
        ev_tally_xyz(i, j, nlocal, newton_pair, 0.0, 0.0, ft[0], ft[1], ft[2], dx[0], dx[1], dx[2]);

      if (newton_pair || j < nlocal) {

        f[j][0] -= ft[0];
        f[j][1] -= ft[1];
        f[j][2] -= ft[2];
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();

  comm->reverse_comm(this);
  comm->forward_comm(this);
*/
}

/* ----------------------------------------------------------------------
 allocate all arrays
 ------------------------------------------------------------------------- */

void PairRHEOTension::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(alpha, n + 1, n + 1, "pair:alpha");
  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
}

/* ----------------------------------------------------------------------
 global settings
 ------------------------------------------------------------------------- */

void PairRHEOTension::settings(int narg, char **arg)
{
}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairRHEOTension::coeff(int narg, char **arg)
{
  if (narg != 3)
    error->all(FLERR,"Incorrect number of args for pair_style rheo coefficients");
  if (!allocated)
    allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi,error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi,error);

  double alpha_one = utils::numeric(FLERR, arg[2], false, lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = 0; j <= atom->ntypes; j++) {
      alpha[i][j] = alpha_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0)
    error->all(FLERR,"Incorrect args for pair rheo/tension coefficients");
}

/* ----------------------------------------------------------------------
 setup specific to this pair style
 ------------------------------------------------------------------------- */

void PairRHEOTension::setup()
{
  auto fixes = modify->get_fix_by_style("rheo");
  if (fixes.size() == 0) error->all(FLERR, "Need to define fix rheo to use pair rheo");
  fix_rheo = dynamic_cast<FixRHEO *>(fixes[0]);
  /*
  compute_kernel = fix_rheo->compute_kernel;
  compute_grad = fix_rheo->compute_grad;
  compute_interface = fix_rheo->compute_interface;
  h = fix_rheo->h;
  csq = fix_rheo->csq;
  rho0 = fix_rheo->rho0;

  hsq = h * h;
  hinv = 1.0 / h;
  hinv3 = hinv * 3.0;
  */
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairRHEOTension::init_style()
{
  neighbor->add_request(this);


  // Create c_tension arrays n_tension arrays if they don't already exist
  // Create a custom atom property so it works with compute property/atom
  // Do not create grow callback as there's no reason to copy/exchange data
  // Manually grow if nmax_store exceeded
  // For B and gradC, create a local array since they are unlikely to be printed

  int tmp1, tmp2;
  int index = atom->find_custom("rheo_c_tension", tmp1, tmp2);
  if (index == -1)  index = atom->add_custom("rheo_c_tension", 1, 0);
  ct = atom->dvector[index];

  index = atom->find_custom("rheo_n_tension", tmp1, tmp2);
  if (index == -1)  index = atom->add_custom("rheo_n_tension", 1, 3);
  nt = atom->darray[index];

  nmax_store = atom->nmax;
}

/* ----------------------------------------------------------------------
 init for one type pair i,j and corresponding j,i
 ------------------------------------------------------------------------- */

double PairRHEOTension::init_one(int i, int j)
{
  if (setflag[i][j] == 0)
    error->all(FLERR,"All pair rheo/tension coeffs are not set");

  alpha[j][i] = alpha[i][j];

  return h;
}


/* ---------------------------------------------------------------------- */

int PairRHEOTension::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)
{
  /*
  int i,j,k,m;
  m = 0;
  double *rho = atom->rho;

  for (i = 0; i < n; i++) {
    j = list[i];
    if (comm_stage == 0) {
      buf[m++] = fp_store[j][0];
      buf[m++] = fp_store[j][1];
      buf[m++] = fp_store[j][2];
    } else {
      buf[m++] = chi[j];
      buf[m++] = rho[j];
    }
  }
  return m;
  */
}

/* ---------------------------------------------------------------------- */

void PairRHEOTension::unpack_forward_comm(int n, int first, double *buf)
{
  /*
  int i, k, m, last;
  double *rho = atom->rho;
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (comm_stage == 0) {
      fp_store[i][0] = buf[m++];
      fp_store[i][1] = buf[m++];
      fp_store[i][2] = buf[m++];
    } else {
      chi[i] = buf[m++];
      rho[i] = buf[m++];
    }
  }
  */
}


/* ---------------------------------------------------------------------- */

int PairRHEOTension::pack_reverse_comm(int n, int first, double *buf)
{
  /*
  int i, k, m, last;
  double **fp_store = compute_interface->fp_store;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = fp_store[i][0];
    buf[m++] = fp_store[i][1];
    buf[m++] = fp_store[i][2];
  }

  return m;
  */
}

/* ---------------------------------------------------------------------- */

void PairRHEOTension::unpack_reverse_comm(int n, int *list, double *buf)
{
  /*
  int i, j, k, m;
  double **fp_store = compute_interface->fp_store;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    fp_store[j][0] += buf[m++];
    fp_store[j][1] += buf[m++];
    fp_store[j][2] += buf[m++];
  }
  */
}
