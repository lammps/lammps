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

// Todo:
// add citations
// remove (or fix) pairwise forces on undercoordinated atoms
// add option for vacuum tension (Frustenau 2020?)

#include "fix_rheo_tension.h"

#include "atom.h"
#include "comm.h"
#include "compute_rheo_kernel.h"
#include "compute_rheo_interface.h"
#include "compute_rheo_vshift.h"
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
  Fix(lmp, narg, arg), compute_kernel(nullptr), compute_interface(nullptr), compute_vshift(nullptr), fix_rheo(nullptr), rho0(nullptr)
{
  if (narg != 8) error->all(FLERR,"Illegal fix command");
  alpha = utils::numeric(FLERR, arg[3], false, lmp);
  beta = utils::numeric(FLERR, arg[4], false, lmp);
  wmin = utils::numeric(FLERR, arg[5], false, lmp);
  cmin = utils::numeric(FLERR, arg[6], false, lmp);
  vshift_strength = utils::numeric(FLERR, arg[7], false, lmp);

  comm_forward = 3;
  comm_reverse = 3;

  // Create cgrad, n, and divr arrays as custom atom properties,
  //   can print with compute property/atom
  //   no grow callback as there's no reason to copy/exchange data, manually grow
  // For norm, create a local array since they are unlikely to be printed

  int tmp1, tmp2;
  index_ct = atom->find_custom("c_rheo_tension", tmp1, tmp2);
  if (index_ct == -1)  index_ct = atom->add_custom("c_rheo_tension", 1, 0);
  ct = atom->dvector[index_ct];

  index_cgradt = atom->find_custom("cgrad_rheo_tension", tmp1, tmp2);
  if (index_cgradt == -1)  index_cgradt = atom->add_custom("cgrad_rheo_tension", 1, 3);
  cgradt = atom->darray[index_cgradt];

  index_nt = atom->find_custom("n_rheo_tension", tmp1, tmp2);
  if (index_nt == -1)  index_nt = atom->add_custom("n_rheo_tension", 1, 3);
  nt = atom->darray[index_nt];

  index_divnt = atom->find_custom("divn_rheo_tension", tmp1, tmp2);
  if (index_divnt == -1) index_divnt = atom->add_custom("divn_rheo_tension", 1, 0);
  divnt = atom->dvector[index_divnt];

  index_wsame = atom->find_custom("wsame_rheo_tension", tmp1, tmp2);
  if (index_wsame == -1) index_wsame = atom->add_custom("wsame_rheo_tension", 1, 0);
  wsame = atom->dvector[index_wsame];

  index_ft = atom->find_custom("f_rheo_tension", tmp1, tmp2);
  if (index_ft == -1)  index_ft = atom->add_custom("f_rheo_tension", 1, 3);
  ft = atom->darray[index_ft];

  norm = nullptr;
  nmax_store = 0;
}

/* ---------------------------------------------------------------------- */

FixRHEOTension::~FixRHEOTension()
{
  // Remove custom property if it exists
  int tmp1, tmp2, index;

  index = atom->find_custom("c_rheo_tension", tmp1, tmp2);
  if (index != -1) atom->remove_custom(index, 1, 0);

  index = atom->find_custom("cgrad_rheo_tension", tmp1, tmp2);
  if (index != -1) atom->remove_custom(index, 1, 3);

  index = atom->find_custom("n_rheo_tension", tmp1, tmp2);
  if (index != -1) atom->remove_custom(index, 1, 3);

  index = atom->find_custom("divn_rheo_tension", tmp1, tmp2);
  if (index != -1) atom->remove_custom(index, 1, 0);

  index = atom->find_custom("wsame_rheo_tension", tmp1, tmp2);
  if (index != -1) atom->remove_custom(index, 1, 0);

  index = atom->find_custom("f_rheo_tension", tmp1, tmp2);
  if (index != -1) atom->remove_custom(index, 1, 3);

  memory->destroy(norm);
}

/* ---------------------------------------------------------------------- */

int FixRHEOTension::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
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
  compute_vshift = fix_rheo->compute_vshift;
  interface_flag = fix_rheo->interface_flag;
  shift_flag = fix_rheo->shift_flag;
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
  // Grow and populate arrays for dump files
  if (nmax_store <= atom->nmax)
    grow_arrays(atom->nmax);

  size_t nbytes = nmax_store * sizeof(double);
  memset(&ct[0], 0, nbytes);
  memset(&norm[0], 0, nbytes);
  memset(&wsame[0], 0, nbytes);
  memset(&divnt[0], 0, nbytes);
  memset(&cgradt[0][0], 0, 3 * nbytes);
  memset(&ft[0][0], 0, 3 * nbytes);
  memset(&nt[0][0], 0, 3 * nbytes);
}

/* ----------------------------------------------------------------------
  Calculate and apply tension forces
------------------------------------------------------------------------- */

void FixRHEOTension::pre_force(int vflag)
{
  int i, j, a, ii, jj, inum, jnum, itype, jtype;
  int fluidi, fluidj;
  double xtmp, ytmp, ztmp, w, wp, ctmp;
  double rhoi, rhoj, Voli, Volj;
  double *dWij, *dWji;
  double dx[3];

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

  size_t nbytes = nmax_store * sizeof(double);
  memset(&ct[0], 0, nbytes);
  memset(&norm[0], 0, nbytes);
  memset(&wsame[0], 0, nbytes);
  memset(&divnt[0], 0, nbytes);
  memset(&cgradt[0][0], 0, 3 * nbytes);
  memset(&ft[0][0], 0, 3 * nbytes);

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

      rhoi = rho[i];
      rhoj = rho[j];

      // Add corrections for walls
      if (interface_flag) {
        if (fluidi && (!fluidj)) {
          rhoj = compute_interface->correct_rho(j, i);
        } else if ((!fluidi) && fluidj) {
          rhoi = compute_interface->correct_rho(i, j);
        } else if ((!fluidi) && (!fluidj)) {
          rhoi = rho0[itype];
          rhoj = rho0[jtype];
        }
      }

      Voli = mass[itype] / rhoi;
      Volj = mass[jtype] / rhoj;

      w = compute_kernel->calc_w(i, j, dx[0], dx[1], dx[2], r);

      if (itype != jtype) ctmp = 1;
      else ctmp = 0;

      ct[i] += ctmp * Volj * w;
      if (newton || j < nlocal)
        ct[j] += ctmp * Voli * w;
    }
  }

  comm_stage = 0;
  comm_reverse = 1;
  if (newton) comm->reverse_comm(this);

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

      rhoi = rho[i];
      rhoj = rho[j];

      // Add corrections for walls
      if (interface_flag) {
        if (fluidi && (!fluidj)) {
          rhoj = compute_interface->correct_rho(j, i);
        } else if ((!fluidi) && fluidj) {
          rhoi = compute_interface->correct_rho(i, j);
        } else if ((!fluidi) && (!fluidj)) {
          rhoi = rho0[itype];
          rhoj = rho0[jtype];
        }
      }

      Voli = mass[itype] / rhoi;
      Volj = mass[jtype] / rhoj;

      wp = compute_kernel->calc_dw(i, j, dx[0], dx[1], dx[2], r);
      dWij = compute_kernel->dWij;
      dWji = compute_kernel->dWji;

      //c = 0;
      //if (itype != jtype) c += rhoi;
      //c /= (rhoi + rhoj);

      if (itype != jtype) ctmp = 1;
      else ctmp = 0;

      for (a = 0; a < dim; a++) {
        cgradt[i][a] -= ctmp * Volj * dWij[a];
        if (newton || j < nlocal)
          cgradt[j][a] -= ctmp * Voli * dWji[a];
      }
    }
  }

  comm_stage = 1;
  comm_reverse = 3;
  if (newton) comm->reverse_comm(this);

  // Calculate normal direction
  double minv;
  for (i = 0; i < nlocal; i++) {
    minv = cgradt[i][0] * cgradt[i][0] + cgradt[i][1] * cgradt[i][1];
    if (dim == 3) minv += cgradt[i][2] * cgradt[i][2];
    minv = sqrt(minv);
    if (minv != 0) minv = 1 / minv;

    for (a = 0; a < dim; a++)
      nt[i][a] = cgradt[i][a] * minv;
  }

  comm_forward = 3;
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
          rhoi = rho0[itype];
          rhoj = rho0[jtype];
        }
      }

      Voli = mass[itype] / rhoi;
      Volj = mass[jtype] / rhoj;

      w = compute_kernel->calc_w(i, j, dx[0], dx[1], dx[2], r);
      wp = compute_kernel->calc_dw(i, j, dx[0], dx[1], dx[2], r);
      dWij = compute_kernel->dWij;
      dWji = compute_kernel->dWji;

      for (a = 0; a < dim; a++) {
        if (itype != jtype) {
          divnt[i] -= (nt[i][a] + nt[j][a]) * Volj * dWij[a];
        } else {
          divnt[i] -= (nt[i][a] - nt[j][a]) * Volj * dWij[a];
          wsame[i] += w * r;
        }
        norm[i] -= dx[a] * Volj * dWij[a];
        if (newton || j < nlocal) {
          if (itype != jtype) {
            divnt[j] -= (nt[j][a] + nt[i][a]) * Voli * dWji[a];
          } else {
            divnt[j] -= (nt[j][a] - nt[i][a]) * Voli * dWji[a];
            wsame[j] += w * r;
          }
          norm[j] += dx[a] * Voli * dWji[a];
        }
      }
    }
  }

  comm_stage = 2;
  comm_reverse = 3;
  if (newton) comm->reverse_comm(this);

  comm_forward = 1;
  comm->forward_comm(this);

  // Skip forces if it's setup
  if (update->setupflag) return;

  // apply force, remove normal vshift

  double **vshift;
  if (shift_flag)
    vshift = compute_vshift->vshift;
  double nx, ny, nz, vx, vy, vz, dot;
  double wmin_inv, weight, prefactor, unwrap[3], v[6];

  if (wmin > 0) wmin_inv = 1.0 / wmin;
  else wmin_inv = 0.0;

  for (i = 0; i < nlocal; i++) {

    if (wsame[i] < wmin) continue;

    weight = MAX(1.0, wsame[i] * wmin_inv);
    itype = type[i];

    if (norm[i] != 0)
      divnt[i] *= dim * norm[i];
    else
      divnt[i] = 0.0;

    // Tension force from Adami, Hu, Adams 2010
    prefactor = -alpha * divnt[i] * weight;
    for (a = 0; a < dim; a++) {
      f[i][a] += prefactor * cgradt[i][a];
      ft[i][a] += prefactor * cgradt[i][a];
    }

    // remove normal shifting component for interfacial particles
    // Based on Yang, Rakhsha, Hu, & Negrut 2022
    if (shift_flag && (vshift_strength != 1.0)) {
      if (ct[i] > cmin) {
        nx = nt[i][0];
        ny = nt[i][1];
        vx = vshift[i][0];
        vy = vshift[i][1];

        dot = nx * vx + ny * vy;
        if (dim == 3) {
          nz = nt[i][2];
          vz = vshift[i][2];
          dot += nz * vz;
        }

        // Allowing shifting into the bulk
        //if (dot > 0.0) continue;

        vshift[i][0] -= (1.0 - vshift_strength) * nx * dot;
        vshift[i][1] -= (1.0 - vshift_strength) * ny * dot;
        if (dim == 3) {
          vshift[i][2] -= (1.0 - vshift_strength) * nz * dot;
        }
      }
    }

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

  // If there is no lower limit, apply optional pairwise forces
  // This is totally ad hoc, needs some work
  // Attempts to deal with stray single particles
  if (wmin <= 0 || beta == 0.0) return;

  int newton_pair = force->newton_pair;
  double fpair, wi, wj;
  double cut_two_thirds = 2.0 * h / 3.0;
  double cut_five_sixths = 5.0 * h / 6.0;
  double cut_sixth_sq = (h / 6.0) * (h / 6.0);
  double cut_third_sq = (h / 3.0) * (h / 3.0);
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    wi = MAX(MIN(1.0, (wmin - wsame[i]) * wmin_inv), 0.0);

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      if (wsame[i] >= wmin && wsame[j] >= wmin) continue;

      dx[0] = xtmp - x[j][0];
      dx[1] = ytmp - x[j][1];
      dx[2] = ztmp - x[j][2];
      rsq = lensq3(dx);

      if (rsq > hsq) continue;

      r = sqrt(rsq);
      jtype = type[j];

      if (itype == jtype) {
        fpair = (r - cut_two_thirds);
        fpair *= fpair;
        fpair -= cut_third_sq;
      } else {
        //fpair = 0.0;

        if (r > (0.5*cut_two_thirds)) continue;
        fpair = (r - cut_two_thirds);
        fpair *= fpair;
        fpair -= cut_third_sq;

        //if (r > cut_two_thirds) continue;
        //fpair = (r - cut_five_sixths);
        //fpair *= fpair;
        //fpair -= cut_sixth_sq;

        //fpair = (h - r) * 0.66666666666666;
      }

      wj = MAX(MIN(1.0, (wmin - wsame[j]) * wmin_inv), 0.0);
      rinv = 1.0 / r;
      fpair *= MAX(wi, wj) * beta * rinv;

      f[i][0] += dx[0] * fpair;
      f[i][1] += dx[1] * fpair;
      f[i][2] += dx[2] * fpair;

      if (newton_pair || j < nlocal) {
        f[j][0] -= dx[0] * fpair;
        f[j][1] -= dx[1] * fpair;
        f[j][2] -= dx[2] * fpair;
      }

      if (evflag) {
        // In progress
      }
    }
  }
}


/* ---------------------------------------------------------------------- */

int FixRHEOTension::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)
{
  int i, j, a, m;
  m = 0;

  if (comm_stage == 1)
    for (i = 0; i < n; i++) {
      j = list[i];
      for (a = 0; a < 3; a++)
        buf[m++] = nt[j][a];
    }
  else if (comm_stage == 2)
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = wsame[j];
    }

  return m;
}

/* ---------------------------------------------------------------------- */

void FixRHEOTension::unpack_forward_comm(int n, int first, double *buf)
{
  int i, a, m, last;

  m = 0;
  last = first + n;
  if (comm_stage == 1)
    for (i = first; i < last; i++)
      for (a = 0; a < 3; a++)
        nt[i][a] = buf[m++];
  else if (comm_stage == 2)
    for (i = first; i < last; i++)
      wsame[i] = buf[m++];
}


/* ---------------------------------------------------------------------- */

int FixRHEOTension::pack_reverse_comm(int n, int first, double *buf)
{
  int i, a, m, last;

  m = 0;
  last = first + n;
  if (comm_stage == 0)
    for (i = first; i < last; i++)
      buf[m++] = ct[i];
  else if (comm_stage == 1)
    for (i = first; i < last; i++)
      for (a = 0; a < 3; a++)
        buf[m++] = cgradt[i][a];
  else if (comm_stage == 2)
    for (i = first; i < last; i++) {
      buf[m++] = norm[i];
      buf[m++] = divnt[i];
      buf[m++] = wsame[i];
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
      ct[j] += buf[m++];
    }
  else if (comm_stage == 1)
    for (i = 0; i < n; i++) {
      j = list[i];
      for (a = 0; a < 3; a++)
        cgradt[j][a] += buf[m++];
    }
  else if (comm_stage == 2)
    for (i = 0; i < n; i++) {
      j = list[i];
      norm[j] += buf[m++];
      divnt[j] += buf[m++];
      wsame[j] += buf[m++];
    }
}

/* ---------------------------------------------------------------------- */

void FixRHEOTension::grow_arrays(int nmax)
{
  // Grow atom variables and reassign pointers
  memory->grow(atom->dvector[index_ct], nmax, "atom:rheo_ct");
  memory->grow(atom->darray[index_cgradt], nmax, 3, "atom:rheo_cgradt");
  memory->grow(atom->darray[index_nt], nmax, 3, "atom:rheo_nt");
  memory->grow(atom->dvector[index_divnt], nmax, "atom:rheo_divnt");
  memory->grow(atom->dvector[index_wsame], nmax, "atom:rheo_wsame");
  memory->grow(atom->darray[index_ft], nmax, 3, "atom:rheo_ft");

  ct = atom->dvector[index_ct];
  cgradt = atom->darray[index_cgradt];
  nt = atom->darray[index_nt];
  divnt = atom->dvector[index_divnt];
  wsame = atom->dvector[index_wsame];
  ft = atom->darray[index_ft];

  // Grow local variables
  memory->grow(norm, nmax, "rheo/tension:norm");

  nmax_store = atom->nmax;
}