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
#include "compute_rheo_kernel.h"
#include "domain.h"
#include "error.h"
#include "fix_rheo.h"
#include "fix_rheo_pressure.h"
#include "force.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace RHEO_NS;
using namespace MathExtra;

static constexpr double EPSILON = 1e-1;

/* ---------------------------------------------------------------------- */

ComputeRHEOInterface::ComputeRHEOInterface(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), chi(nullptr), fp_store(nullptr), fix_rheo(nullptr), rho0(nullptr),
    norm(nullptr), normwf(nullptr), id_fix_pa(nullptr), list(nullptr), compute_kernel(nullptr),
    fix_pressure(nullptr)

{
  if (narg != 3) error->all(FLERR, "Illegal compute rheo/interface command");

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
  int a, i, j, ii, jj, jnum, itype, jtype, fluidi, fluidj, status_match;
  double xtmp, ytmp, ztmp, rsq, w, dot, dx[3];

  int inum, *ilist, *jlist, *numneigh, **firstneigh;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  double **x = atom->x;
  int *type = atom->type;
  int newton = force->newton;
  int *status = atom->rheo_status;
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

      dx[0] = xtmp - x[j][0];
      dx[1] = ytmp - x[j][1];
      dx[2] = ztmp - x[j][2];
      rsq = lensq3(dx);

      if (rsq < cutsq) {
        jtype = type[j];
        fluidj = !(status[j] & PHASECHECK);
        w = compute_kernel->calc_w_quintic(sqrt(rsq));

        norm[i] += w;

        status_match = 0;
        if ((fluidi && fluidj) || ((!fluidi) && (!fluidj))) status_match = 1;
        if (status_match) {
          chi[i] += w;
        } else {
          if (!fluidi) {
            dot = 0;
            for (a = 0; a < 3; a++) dot += (-fp_store[j][a] + fp_store[i][a]) * dx[a];

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
              dot = 0;
              for (a = 0; a < 3; a++) dot += (-fp_store[i][a] + fp_store[j][a]) * dx[a];

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

int ComputeRHEOInterface::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/,
                                            int * /*pbc*/)
{
  int m = 0;
  double *rho = atom->rho;

  for (int i = 0; i < n; i++) {
    int j = list[i];
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
  double *rho = atom->rho;
  int m = 0;
  int last = first + n;
  for (int i = first; i < last; i++) {
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
  double *rho = atom->rho;
  int m = 0;
  int last = first + n;
  for (int i = first; i < last; i++) {
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
  double *rho = atom->rho;
  int *status = atom->rheo_status;
  int m = 0;
  for (int i = 0; i < n; i++) {
    int j = list[i];
    norm[j] += buf[m++];
    chi[j] += buf[m++];
    if (status[j] & PHASECHECK) {
      normwf[j] += buf[m++];
      rho[j] += buf[m++];
    } else {
      m++;
      m++;
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOInterface::correct_v(double *v_solid, double *v_fluid, int i_solid, int i_fluid)
{
  double wall_prefactor, wall_denom, wall_numer;

  wall_numer = MAX(2.0 * cut * (chi[i_solid] - 0.5), 0.0);
  wall_denom = MAX(2.0 * cut * (chi[i_fluid] - 0.5), wall_max);

  wall_prefactor = wall_numer / wall_denom;

  v_solid[0] = (v_solid[0] - v_fluid[0]) * wall_prefactor;
  v_solid[1] = (v_solid[1] - v_fluid[1]) * wall_prefactor;
  v_solid[2] = (v_solid[2] - v_fluid[2]) * wall_prefactor;
}

/* ---------------------------------------------------------------------- */

double ComputeRHEOInterface::correct_rho(int i_solid)
{
  int itype = atom->type[i_solid];
  return MAX(rho0[itype], atom->rho[i_solid]);
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOInterface::store_forces()
{
  double minv;
  int *type = atom->type;
  int *mask = atom->mask;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
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
        if (rmass)
          minv = 1.0 / rmass[i];
        else
          minv = 1.0 / mass[type[i]];
        if (mask[i] & fix->groupbit)
          for (int a = 0; a < 3; a++) fp_store[i][a] = f[i][a] * minv;
        else
          for (int a = 0; a < 3; a++) fp_store[i][a] = (f[i][a] - fp_store[i][a]) * minv;
      }
    }
  } else {
    if (rmass) {
      for (int i = 0; i < atom->nlocal; i++) {
        minv = 1.0 / rmass[i];
        for (int a = 0; a < 3; a++) fp_store[i][a] = (f[i][a] - fp_store[i][a]) * minv;
      }
    } else {
      for (int i = 0; i < atom->nlocal; i++) {
        minv = 1.0 / mass[type[i]];
        for (int a = 0; a < 3; a++) fp_store[i][a] = (f[i][a] - fp_store[i][a]) * minv;
      }
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
