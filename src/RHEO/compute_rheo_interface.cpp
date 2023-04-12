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
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"

#include <cmath>

using namespace LAMMPS_NS;

#define EPSILON 1e-1

/* ---------------------------------------------------------------------- */

ComputeRHEOInterface::ComputeRHEOInterface(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg), fix_rheo(nullptr), compute_kernel(nullptr), fx_m_norm(nullptr),
  norm(nullptr), normwf(nullptr), chi(nullptr), f_pressure(nullptr), id_fix_pa(nullptr)
{
  if (narg != 3) error->all(FLERR,"Illegal compute rheo/interface command");

  nmax = 0;

  comm_forward = 3;
  comm_reverse = 4;
}

/* ---------------------------------------------------------------------- */

ComputeRHEOInterface::~ComputeRHEOInterface()
{
  // Remove custom property if it exists
  int tmp1, tmp2, index;
  index = atom->find_custom("rheo_chi", tmp1, tmp2);
  if (index != -1) atom->remove_custom(index, 1, 0);

  if (id_fix_pa && modify->nfix) modify->delete_fix(id_fix_pa);

  memory->destroy(fx_m_norm);
  memory->destroy(norm);
  memory->destroy(normwf);
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOInterface::init()
{
  compute_kernel = fix_rheo->compute_kernel;
  cut = fix_rheo->cut;
  cutsq = cut * cut;
  wall_max = sqrt(3.0) / 12.0 * cut;

  // Create chi array if it doesn't already exist
  // Create a custom atom property so it works with compute property/atom
  // Do not create grow callback as there's no reason to copy/exchange data
  // Manually grow if nmax_old exceeded

  int create_flag = 0;
  int tmp1, tmp2;
  int nmax = atom->nmax;
  int index = atom->find_custom("rheo_chi", tmp1, tmp2);
  if (index == -1) {
    index = atom->add_custom("rheo_chi", 1, 0);
    nmax_old = nmax;
  }
  chi = atom->dvector[index];

  // For fpressure, go ahead and create an instance of fix property atom
  // Need restarts + exchanging with neighbors since it needs to persist
  // between timesteps (fix property atom will handle callbacks)

  index = atom->find_custom("rheo_pressure", tmp1, tmp2);
  if (index == -1) {
    id_fix_pa = utils::strdup(id + std::string("_fix_property_atom"));
    modify->add_fix(fmt::format("{} all property/atom d2_f_pressure 3", id_fix_pa)));
    index = atom->find_custom("rheo_pressure", tmp1, tmp2);
  }
  f_pressure = atom->darray[index];

  // need an occasional half neighbor list
  neighbor->add_request(this, NeighConst::REQ_HALF);
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOInterface::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

// Left off here
void ComputeRHEOInterface::compute_peratom()
{
  int i, j, ii, jj, jnum, itype, jtype, phase_match;
  double xtmp, ytmp, ztmp, delx, dely, delz, rsq;
  int *jlist;
  double w;

  // neighbor list ariables
  int inum, *ilist, *numneigh, **firstneigh;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  double **x = atom->x;
  int *type = atom->type;
  int newton = force->newton;
  int *phase = atom->phase;
  double *rho = atom->rho;
  double *fx = atom->dvector[index_fx];
  double *fy = atom->dvector[index_fy];
  double *fz = atom->dvector[index_fz];
  double mi, mj;

  //Declare mass pointer to calculate acceleration from force
  double *mass = atom->mass;

  cs2 = 1.0;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    fix_chi->grow_arrays(nmax);
    chi = fix_chi->vstore;
    memory->destroy(norm);
    memory->destroy(normwf);
    memory->create(norm, nmax, "RHEO/chi:norm");
    memory->create(normwf, nmax, "RHEO/chi:normwf");
  }

  for (i = 0; i < nall; i++) {
    if (phase[i] > FixRHEO::FLUID_MAX) rho[i] = 0.0;
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
    jlist = firstneigh[i];
    jnum = numneigh[i];
    mi = mass[type[i]];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;

      if (rsq < cutsq) {
        jtype = type[j];
        mj = mass[type[j]];
        w = compute_kernel->calc_w_quintic(i, j, delx, dely, delz, sqrt(rsq));

        phase_match = 0;
        norm[i] += w;
        if ((phase[i] <= FixRHEO::FLUID_MAX && phase[j] <= FixRHEO::FLUID_MAX)
            || (phase[i] > FixRHEO::FLUID_MAX && phase[j] > FixRHEO::FLUID_MAX)) {
          phase_match = 1;
        }

        if (phase_match) {
          chi[i] += w;
        } else {
          if (phase[i] > FixRHEO::FLUID_MAX) {
            //speed of sound and rho0 assumed to = 1 (units in density not pressure)
            //In general, rho is calculated using the force vector on the wall particle
            //fx stores f-fp
            rho[i] += w*(cs2*(rho[j] - 1.0) - rho[j]*((-fx[j]/mj+fx[i]/mi)*delx + (-fy[j]/mj+fy[i]/mi)*dely + (-fz[j]/mj+fz[i]/mi)*delz));
            //For the specific taste case whre force on wall particles = 0
            //rho[i] += w*(1.0*(rho[j] - 1.0) + rho[j]*(1e-3*delx));
            normwf[i] += w;
          }
        }

        if (newton || j < nlocal) {
          norm[j] += w;
          if (phase_match) {
            chi[j] += w;
          } else {
            if (phase[j] > FixRHEO::FLUID_MAX) {
              rho[j] += w*(cs2*(rho[i] - 1.0) + rho[i]*((-fx[i]/mi+fx[j]/mj)*delx + (-fy[i]/mi+fy[j]/mj)*dely + (-fz[i]/mi+fz[j]/mj)*delz));
              normwf[j] += w;
            }
          }
        }
      }
    }
  }

  if (newton) comm->reverse_comm_compute(this);

  for (i = 0; i < nlocal; i++) {
    if (norm[i] != 0.0) chi[i] /= norm[i];
    if (normwf[i] != 0.0) { // Only if it's a wall particle
      rho[i] = 1.0 + (rho[i] / normwf[i])/cs2; // Stores rho for solid particles 1+Pw in Adami Adams 2012
      if (rho[i] < EPSILON) rho[i] = EPSILON;
    }

    if (normwf[i] == 0.0 && phase[i] > FixRHEO::FLUID_MAX) rho[i] = 1.0;
  }

  comm_stage = 1;
  comm_forward = 2;
  comm->forward_comm_compute(this);
}

/* ---------------------------------------------------------------------- */

int ComputeRHEOInterface::pack_forward_comm(int n, int *list, double *buf,
                                        int /*pbc_flag*/, int * /*pbc*/)
{
  int i,j,k,m;
  m = 0;
  double *rho = atom->rho;
  double *fx = atom->dvector[index_fx];
  double *fy = atom->dvector[index_fy];
  double *fz = atom->dvector[index_fz];

  for (i = 0; i < n; i++) {
    j = list[i];
    if (comm_stage == 0) {
      buf[m++] = fx[j];
      buf[m++] = fy[j];
      buf[m++] = fz[j];
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
  double *fx = atom->dvector[index_fx];
  double *fy = atom->dvector[index_fy];
  double *fz = atom->dvector[index_fz];

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (comm_stage == 0) {
      fx[i] = buf[m++];
      fy[i] = buf[m++];
      fz[i] = buf[m++];
    } else {
      chi[i] = buf[m++];
      rho[i] = buf[m++]; // Won't do anything for fluids
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
  int *phase = atom->phase;
  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    norm[j] += buf[m++];
    chi[j] += buf[m++];
    if (phase[j] > FixRHEO::FLUID_MAX){
      normwf[j] += buf[m++];
      rho[j] += buf[m++];
    } else {
      m++;
      m++;
    }
  }
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOInterface::correct_v(double *vi, double *vj, double *vi_out, int i, int j)
{
  double wall_prefactor, wall_denom, wall_numer;

  wall_numer = 2.0*cut*(chi[i]-0.5);
  if (wall_numer < 0) wall_numer = 0;
  wall_denom = 2.0*cut*(chi[j]-0.5);
  if (wall_denom < wall_max) wall_denom = wall_max;

  wall_prefactor = wall_numer / wall_denom;

  vi_out[0] = (vi[0]-vj[0])*wall_prefactor + vi[0];
  vi_out[1] = (vi[1]-vj[1])*wall_prefactor + vi[1];
  vi_out[2] = (vi[2]-vj[2])*wall_prefactor + vi[2];
}

/* ---------------------------------------------------------------------- */

double ComputeRHEOInterface::correct_rho(int i, int j) // i is wall, j is fluid
{
  //In future may depend on atom type j's pressure equation
  return atom->rho[i];
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOInterface::store_forces()
{
  double *fx = atom->dvector[index_fx];
  double *fy = atom->dvector[index_fy];
  double *fz = atom->dvector[index_fz];
  double **f = atom->f;
  double **fp = atom->fp;
  int *mask = atom->mask;

  int flag;
  int ifix = modify->find_fix_by_style("setforce");
  if (ifix != -1) {
    for (int i = 0; i < atom->nlocal; i++) {
      if (mask[i] & modify->fix[ifix]->groupbit) {
        fx[i] = f[i][0];
        fy[i] = f[i][1];
        fz[i] = f[i][2];
      } else {
        fx[i] = f[i][0] - fp[i][0];
        fy[i] = f[i][1] - fp[i][1];
        fz[i] = f[i][2] - fp[i][2];
      }
    }
  } else {
    for (int i = 0; i < atom->nlocal; i++) {
      fx[i] = f[i][0] - fp[i][0];
      fy[i] = f[i][1] - fp[i][1];
      fz[i] = f[i][2] - fp[i][2];
    }
  }

  //Forward comm forces -note only needed here b/c property atom will forward otherwise
  comm_forward = 3;
  comm_stage = 0;
  comm->forward_comm_compute(this);
}


/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeRHEOInterface::memory_usage()
{
  double bytes = nmax * sizeof(double);
  return bytes;
}

