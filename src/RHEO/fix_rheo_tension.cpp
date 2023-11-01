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
   Joel Clemmer (SNL), Thomas O'Connor (CMU), Eric Palermo (CMU)
----------------------------------------------------------------------- */

#include "fix_rheo_tension.h"

#include "atom.h"
#include "comm.h"
#include "compute_rheo_kernel.h"
#include "compute_rheo_interface.h"
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
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixRHEOTension::FixRHEOTension(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), compute_kernel(nullptr), compute_interface(nullptr), fix_rheo(nullptr)
{
  if (narg != 4) error->all(FLERR,"Illegal fix command");
  alpha = utils::numeric(FLERR,arg[3],false,lmp);

  comm_forward = 3;
  comm_reverse = 3;

  // Create cgrad, n, and divr arrays as custom atom properties,
  //   can print with compute property/atom
  //   no grow callback as there's no reason to copy/exchange data, manually grow
  // For norm, create a local array since they are unlikely to be printed

  int tmp1, tmp2;
  index_cgradt = atom->find_custom("cgrad_rheo_tension", tmp1, tmp2);
  if (index_cgradt == -1)  index_cgradt = atom->add_custom("cgrad_rheo_tension", 1, 3);
  cgradt = atom->darray[index_cgradt];

  index_nt = atom->find_custom("n_rheo_tension", tmp1, tmp2);
  if (index_nt == -1)  index_nt = atom->add_custom("n_rheo_tension", 1, 3);
  nt = atom->darray[index_nt];

  index_divnt = atom->find_custom("divn_rheo_tension", tmp1, tmp2);
  if (index_divnt == -1) index_divnt = atom->add_custom("divn_rheo_tension", 1, 0);
  divnt = atom->dvector[index_divnt];

  norm = nullptr;
  nmax_store = 0;
}

/* ---------------------------------------------------------------------- */

FixRHEOTension::~FixRHEOTension()
{
  // Remove custom property if it exists
  int tmp1, tmp2, index;

  index = atom->find_custom("cgrad_rheo_tension", tmp1, tmp2);
  if (index != -1) atom->remove_custom(index, 1, 3);

  index = atom->find_custom("n_rheo_tension", tmp1, tmp2);
  if (index != -1) atom->remove_custom(index, 1, 3);

  index = atom->find_custom("divn_rheo_tension", tmp1, tmp2);
  if (index != -1) atom->remove_custom(index, 1, 0);

  memory->destroy(norm);
}

/* ---------------------------------------------------------------------- */

int FixRHEOTension::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixRHEOTension::init()
{
  auto fixes = modify->get_fix_by_style("^rheo$");
  if (fixes.size() == 0) error->all(FLERR, "Need to define fix rheo to use fix rheo/tension");
  fix_rheo = dynamic_cast<FixRHEO *>(fixes[0]);

  compute_kernel = fix_rheo->compute_kernel;
  compute_interface = fix_rheo->compute_interface;
  interface_flag = fix_rheo->interface_flag;
  h = fix_rheo->h;
  rho0 = fix_rheo->rho0;

  hsq = h * h;

  neighbor->add_request(this, NeighConst::REQ_DEFAULT);
}

/* ---------------------------------------------------------------------- */

void FixRHEOTension::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}


/* ---------------------------------------------------------------------- */

void FixRHEOTension::setup(int vflag)
{
  // Grow and populate arrays
  post_force(vflag);
}

/* ----------------------------------------------------------------------
  Calculate and apply tension forces
------------------------------------------------------------------------- */

void FixRHEOTension::post_force(int vflag)
{
  int i, j, a, ii, jj, inum, jnum, itype, jtype;
  int fluidi, fluidj;
  double xtmp, ytmp, ztmp, w, wp, c;
  double rhoi, rhoj, Voli, Volj;
  double *dWij, *dWji;
  double dx[3], ft[3];

  int *ilist, *jlist, *numneigh, **firstneigh;
  double imass, jmass, rsq, r, rinv;

  int nlocal = atom->nlocal;
  int newton = force->newton;
  int dim = domain->dimension;

  v_init(vflag);

  double **x = atom->x;
  double **f = atom->f;
  double *rho = atom->rho;
  double *mass = atom->mass;
  imageint *image = atom->image;
  int *type = atom->type;
  int *status = atom->status;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  if (nmax_store <= atom->nmax)
    grow_arrays(atom->nmax);

  for (i = 0; i < nlocal+atom->nghost; i++) {
    cgradt[i][0] = 0.0;
    cgradt[i][1] = 0.0;
    cgradt[i][2] = 0.0;
    norm[i] = 0.0;
    divnt[i] = 0.0;
  }

  // Calculate color gradient
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    fluidi = !(status[i] & PHASECHECK);
    jlist = firstneigh[i];
    jnum = numneigh[i];
    imass = mass[itype];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      dx[0] = xtmp - x[j][0];
      dx[1] = ytmp - x[j][1];
      dx[2] = ztmp - x[j][2];
      rsq = lensq3(dx);

      if (rsq > hsq) continue;

      fluidj = !(status[j] & PHASECHECK);
      jtype = type[j];
      r = sqrt(rsq);
      rinv = 1 / r;

      rhoi = rho[i];
      rhoj = rho[j];

      // Add corrections for walls
      if (interface_flag) {
        if (fluidi && (!fluidj)) {
          rhoj = compute_interface->correct_rho(j, i);
        } else if ((!fluidi) && fluidj) {
          rhoi = compute_interface->correct_rho(i, j);
        } else if ((!fluidi) && (!fluidj)) {
          rhoi = rho0;
          rhoj = rho0;
        }
      }

      Voli = mass[itype] / rhoi;
      Volj = mass[jtype] / rhoj;

      wp = compute_kernel->calc_dw(i, j, dx[0], dx[1], dx[2],r);
      dWij = compute_kernel->dWij;
      dWji = compute_kernel->dWji;

      c = 0;
      if (itype == jtype) c += rhoi;
      c /= (rhoi + rhoj);

      for (a = 0; a < 3; a++) {
        cgradt[i][a] -= c * Volj * dWij[a];
        if (newton || j < nlocal)
          cgradt[j][a] -= c * Voli * dWji[a];
      }
    }
  }

  comm_stage = 0;
  comm_reverse = 3;
  if (newton) comm->reverse_comm(this);

  // Calculate normal direction
  double minv;
  for (i = 0; i < nlocal; i++) {
    minv = 1.0 / sqrt(cgradt[i][0] * cgradt[i][0] + cgradt[i][1] * cgradt[i][1] + cgradt[i][2] * cgradt[i][2]);

    for (a = 0; a < 3; a++)
      nt[i][a] = cgradt[i][a] * minv;
  }

  comm->forward_comm(this);

  // Calculate divergence
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    fluidi = !(status[i] & PHASECHECK);
    jlist = firstneigh[i];
    jnum = numneigh[i];
    imass = mass[itype];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      dx[0] = xtmp - x[j][0];
      dx[1] = ytmp - x[j][1];
      dx[2] = ztmp - x[j][2];
      rsq = lensq3(dx);

      if (rsq > hsq) continue;

      fluidj = !(status[j] & PHASECHECK);
      jtype = type[j];
      r = sqrt(rsq);
      rinv = 1 / r;

      rhoi = rho[i];
      rhoj = rho[j];

      // Add corrections for walls
      if (interface_flag) {
        if (fluidi && (!fluidj)) {
          rhoj = compute_interface->correct_rho(j, i);
        } else if ((!fluidi) && fluidj) {
          rhoi = compute_interface->correct_rho(i, j);
        } else if ((!fluidi) && (!fluidj)) {
          rhoi = rho0;
          rhoj = rho0;
        }
      }

      Voli = mass[itype] / rhoi;
      Volj = mass[jtype] / rhoj;

      wp = compute_kernel->calc_dw(i, j, dx[0], dx[1], dx[2],r);
      dWij = compute_kernel->dWij;
      dWji = compute_kernel->dWji;

      for (a = 0; a < 3; a++) {
        divnt[i] -= nt[i][a] * Volj * dWij[a];
        norm[i] -= dx[a] * Volj * dWij[a];
        if (newton || j < nlocal) {
          divnt[j] -= nt[j][a] * Voli * dWji[a];
          norm[j] += dx[a] * Voli * dWji[a];
        }
      }
    }
  }

  comm_stage = 1;
  comm_reverse = 2;
  if (newton) comm->reverse_comm(this);

  // Skip forces if it's setup
  if (update->setupflag) return;

  // apply force
  int prefactor;
  double unwrap[3];
  double v[6];
  for (i = 0; i < nlocal; i++) {
    itype = type[i];
    divnt[i] /= norm[i];

    prefactor *= -alpha * divnt[i] / mass[itype];

    for (a = 0; a < 3; a++)
      f[i][a] += prefactor * cgradt[i][a];

    if (evflag) {
      domain->unmap(x[i], image[i], unwrap);
      v[0] = prefactor * cgradt[i][0] * unwrap[0];
      v[1] = prefactor * cgradt[i][1] * unwrap[1];
      v[2] = prefactor * cgradt[i][2] * unwrap[2];
      v[3] = prefactor * cgradt[i][0] * unwrap[1];
      v[4] = prefactor * cgradt[i][0] * unwrap[2];
      v[5] = prefactor * cgradt[i][1] * unwrap[2];
      v_tally(i, v);
    }
  }



  if (evflag) {

        }

}


/* ---------------------------------------------------------------------- */

int FixRHEOTension::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)
{
  int i, j, a, m;
  m = 0;

  for (i = 0; i < n; i++) {
    j = list[i];
    for (a = 0; a < 3; a++)
      buf[m++] = nt[j][a];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixRHEOTension::unpack_forward_comm(int n, int first, double *buf)
{
  int i, a, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++)
    for (a = 0; a < 3; a++)
      nt[i][a] = buf[m++];
}


/* ---------------------------------------------------------------------- */

int FixRHEOTension::pack_reverse_comm(int n, int first, double *buf)
{
  int i, a, m, last;

  m = 0;
  last = first + n;
  if (comm_stage == 0)
    for (i = first; i < last; i++)
      for (a = 0; a < 3; a++)
        buf[m++] = cgradt[i][a];
  else
    for (i = first; i < last; i++) {
      buf[m++] = norm[i];
      buf[m++] = divnt[i];
    }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixRHEOTension::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i, j, a, m;

  m = 0;
  if (comm_stage == 0)
    for (i = 0; i < n; i++) {
      j = list[i];
      for (a = 0; a < 3; a++)
        cgradt[j][a] += buf[m++];
    }
  else
    for (i = 0; i < n; i++) {
      j = list[i];
      norm[j] += buf[m++];
      divnt[j] += buf[m++];
    }
}

/* ---------------------------------------------------------------------- */

void FixRHEOTension::grow_arrays(int nmax)
{
  // Grow atom variables and reassign pointers
  memory->grow(atom->darray[index_cgradt], nmax, 3, "atom:rheo_cgradt");
  memory->grow(atom->darray[index_nt], nmax, 3, "atom:rheo_nt");
  memory->grow(atom->dvector[index_divnt], nmax, "atom:rheo_divnt");

  cgradt = atom->darray[index_cgradt];
  nt = atom->darray[index_nt];
  divnt = atom->dvector[index_divnt];

  // Grow local variables
  memory->grow(norm, nmax, "rheo/tension:norm");

  nmax_store = atom->nmax;
}