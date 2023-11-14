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

#include "compute_rheo_interface.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "compute_rheo_kernel.h"
#include "error.h"
#include "force.h"
#include "fix_rheo.h"
#include "fix_rheo_pressure.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace RHEO_NS;

static constexpr double EPSILON = 1e-1;

/* ---------------------------------------------------------------------- */

ComputeRHEOInterface::ComputeRHEOInterface(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg), fix_rheo(nullptr), compute_kernel(nullptr), fp_store(nullptr),
  rho0(nullptr), norm(nullptr), normwf(nullptr), chi(nullptr), id_fix_pa(nullptr)
{
  if (narg != 3) error->all(FLERR,"Illegal compute rheo/interface command");

  comm_forward = 3;
  comm_reverse = 4;

  nmax_store = atom->nmax;
  memory->create(chi, nmax_store, "rheo:chi");
  memory->create(norm, nmax_store, "rheo/interface:norm");
  memory->create(normwf, nmax_store, "rheo/interface:normwf");

  // For fp_store, create an instance of fix property atom
  // Need restarts + exchanging with neighbors since it needs to persist
  // between timesteps (fix property atom will handle callbacks)

  int tmp1, tmp2;
  int index = atom->find_custom("fp_store", tmp1, tmp2);
  if (index == -1) {
    id_fix_pa = utils::strdup(id + std::string("_fix_property_atom"));
    modify->add_fix(fmt::format("{} all property/atom d2_fp_store 3", id_fix_pa));
    index = atom->find_custom("fp_store", tmp1, tmp2);
  }
  fp_store = atom->darray[index];
}

/* ---------------------------------------------------------------------- */

ComputeRHEOInterface::~ComputeRHEOInterface()
{
  if (id_fix_pa && modify->nfix) modify->delete_fix(id_fix_pa);
  delete[] id_fix_pa;
  memory->destroy(chi);
  memory->destroy(norm);
  memory->destroy(normwf);
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOInterface::init()
{
  compute_kernel = fix_rheo->compute_kernel;
  rho0 = fix_rheo->rho0;
  cut = fix_rheo->cut;
  cutsq = cut * cut;
  wall_max = sqrt(3.0) / 12.0 * cut;

  auto fixes = modify->get_fix_by_style("rheo/pressure");
  fix_pressure = dynamic_cast<FixRHEOPressure *>(fixes[0]);

  neighbor->add_request(this, NeighConst::REQ_DEFAULT);
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOInterface::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOInterface::compute_peratom()
{
  int i, j, ii, jj, jnum, itype, jtype, fluidi, fluidj, status_match;
  double xtmp, ytmp, ztmp, delx, dely, delz, rsq, w, dot;

  int inum, *ilist, *jlist, *numneigh, **firstneigh;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  double **x = atom->x;
  int *type = atom->type;
  int newton = force->newton;
  int *status = atom->status;
  double *rho = atom->rho;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  if (atom->nmax > nmax_store) {
    nmax_store = atom->nmax;
    memory->grow(norm, nmax_store, "rheo/interface:norm");
    memory->grow(normwf, nmax_store, "rheo/interface:normwf");
    memory->grow(chi, nmax_store, "rheo:chi");
  }

  for (i = 0; i < nall; i++) {
    if (status[i] & PHASECHECK) rho[i] = 0.0;
    normwf[i] = 0.0;
    norm[i] = 0.0;
    chi[i] = 0.0;
  }

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    fluidi = !(status[i] & PHASECHECK);
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;

      if (rsq < cutsq) {
        jtype = type[j];
        fluidj = !(status[j] & PHASECHECK);
        w = compute_kernel->calc_w_quintic(i, j, delx, dely, delz, sqrt(rsq));

        status_match = 0;
        norm[i] += w;
        if ((fluidi && fluidj) || ((!fluidi) && (!fluidj)))
          status_match = 1;

        if (status_match) {
          chi[i] += w;
        } else {
          if (!fluidi) {
            dot = (-fp_store[j][0] + fp_store[i][0]) * delx;
            dot += (-fp_store[j][1] + fp_store[i][1]) * dely;
            dot += (-fp_store[j][2] + fp_store[i][2]) * delz;

            rho[i] += w * (fix_pressure->calc_pressure(rho[j], jtype) - rho[j] * dot);
            normwf[i] += w;
          }
        }

        if (newton || j < nlocal) {
          norm[j] += w;
          if (status_match) {
            chi[j] += w;
          } else {
            if (!fluidj) {
              dot = (-fp_store[i][0] + fp_store[j][0]) * delx;
              dot += (-fp_store[i][1] + fp_store[j][1]) * dely;
              dot += (-fp_store[i][2] + fp_store[j][2]) * delz;

              rho[j] += w * (fix_pressure->calc_pressure(rho[i], itype) + rho[i] * dot);
              normwf[j] += w;
            }
          }
        }
      }
    }
  }

  if (newton) comm->reverse_comm(this);

  for (i = 0; i < nlocal; i++) {
    if (norm[i] != 0.0) chi[i] /= norm[i];

    // Recalculate rho for non-fluid particles
    if (status[i] & PHASECHECK) {
      if (normwf[i] != 0.0) {
        // Stores rho for solid particles 1+Pw in Adami Adams 2012
        rho[i] = MAX(EPSILON, fix_pressure->calc_rho(rho[i] / normwf[i], type[i]));
      } else {
        rho[i] = rho0[itype];
      }
    }
  }

  comm_stage = 1;
  comm_forward = 2;
  comm->forward_comm(this);
}

/* ---------------------------------------------------------------------- */

int ComputeRHEOInterface::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)
{
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
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOInterface::unpack_forward_comm(int n, int first, double *buf)
{
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
}

/* ---------------------------------------------------------------------- */

int ComputeRHEOInterface::pack_reverse_comm(int n, int first, double *buf)
{
  int i,k,m,last;
  double *rho = atom->rho;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = norm[i];
    buf[m++] = chi[i];
    buf[m++] = normwf[i];
    buf[m++] = rho[i];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOInterface::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,k,j,m;
  double *rho = atom->rho;
  int *status = atom->status;
  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    norm[j] += buf[m++];
    chi[j] += buf[m++];
    if (status[j] & PHASECHECK){
      normwf[j] += buf[m++];
      rho[j] += buf[m++];
    } else {
      m++;
      m++;
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOInterface::correct_v(double *vi, double *vj, int i, int j)
{
  double wall_prefactor, wall_denom, wall_numer;

  wall_numer = 2.0 * cut * (chi[i] - 0.5);
  if (wall_numer < 0) wall_numer = 0;
  wall_denom = 2.0 * cut * (chi[j] - 0.5);
  if (wall_denom < wall_max) wall_denom = wall_max;

  wall_prefactor = wall_numer / wall_denom;

  vi[0] = (vi[0] - vj[0]) * wall_prefactor + vi[0];
  vi[1] = (vi[1] - vj[1]) * wall_prefactor + vi[1];
  vi[2] = (vi[2] - vj[2]) * wall_prefactor + vi[2];
}

/* ---------------------------------------------------------------------- */

double ComputeRHEOInterface::correct_rho(int i, int j)
{
  // i is wall, j is fluid
  //In future may depend on atom type j's pressure equation
  return atom->rho[i];
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOInterface::store_forces()
{
  double minv;
  int *type = atom->type;
  int *mask = atom->mask;
  double *mass = atom->mass;
  double **f = atom->f;

  // When this is called, fp_store stores the pressure force
  // After this method, fp_store instead stores non-pressure forces
  // and is also normalized by the particles mass
  // If forces are overwritten by a fix, there are no pressure forces
  // so just normalize
  auto fixlist = modify->get_fix_by_style("setforce");
  if (fixlist.size() != 0) {
    for (const auto &fix : fixlist) {
      for (int i = 0; i < atom->nlocal; i++) {
        minv = 1.0 / mass[type[i]];
        if (mask[i] & fix->groupbit) {
          fp_store[i][0] = f[i][0] * minv;
          fp_store[i][1] = f[i][1] * minv;
          fp_store[i][2] = f[i][2] * minv;
        } else {
          fp_store[i][0] = (f[i][0] - fp_store[i][0]) * minv;
          fp_store[i][1] = (f[i][1] - fp_store[i][1]) * minv;
          fp_store[i][2] = (f[i][2] - fp_store[i][2]) * minv;
        }
      }
    }
  } else {
    for (int i = 0; i < atom->nlocal; i++) {
      minv = 1.0 / mass[type[i]];
      fp_store[i][0] = (f[i][0] - fp_store[i][0]) * minv;
      fp_store[i][1] = (f[i][1] - fp_store[i][1]) * minv;
      fp_store[i][2] = (f[i][2] - fp_store[i][2]) * minv;
    }
  }

  // Forward comm forces
  comm_forward = 3;
  comm_stage = 0;
  comm->forward_comm(this);
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeRHEOInterface::memory_usage()
{
  double bytes = 3 * nmax_store * sizeof(double);
  return bytes;
}

