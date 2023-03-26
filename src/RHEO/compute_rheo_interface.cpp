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

#include "fix_rheo.h"
#include "compute_rheo_kernel.h"
#include "fix_store.h"
#include "fix.h"
#include <cmath>
#include <cstring>
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "pair.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define EPSILON 1e-1

/* ---------------------------------------------------------------------- */

ComputeRHEOInterface::ComputeRHEOInterface(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  id_fix_chi(nullptr)
{
  if (narg != 4) error->all(FLERR,"Illegal compute RHEO/chi command");

  cut = utils::numeric(FLERR,arg[3],false,lmp);
  cutsq = cut*cut;

  wall_max = sqrt(3)/12.0*cut;

  nmax = 0;

  comm_forward = 3;
  comm_reverse = 4;

  fix_chi = nullptr;
  norm = nullptr;
  normwf = nullptr;

  // new id = fix-ID + FIX_STORE_ATTRIBUTE
  // new fix group = group for this fix

  id_fix_chi = nullptr;
  std::string fixcmd = id + std::string("_chi");
  id_fix_chi = new char[fixcmd.size()+1];
  strcpy(id_fix_chi,fixcmd.c_str());
  fixcmd += fmt::format(" all STORE peratom 0 {}", 1);
  modify->add_fix(fixcmd);
  fix_chi = (FixStore *) modify->fix[modify->nfix-1];
  chi = fix_chi->vstore;
}

/* ---------------------------------------------------------------------- */

ComputeRHEOInterface::~ComputeRHEOInterface()
{
  if (id_fix_chi && modify->nfix) modify->delete_fix(id_fix_chi);
  if (modify->nfix) modify->delete_fix("PROPERTY_ATOM_COMP_RHEO_SOLIDS");

  memory->destroy(norm);
  memory->destroy(normwf);
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOInterface::init()
{
  // need an occasional full neighbor list
  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->half = 1;
  neighbor->requests[irequest]->full = 0;
  //neighbor->requests[irequest]->occasional = 1; //Anticipate needing regulalry

  int icompute = modify->find_compute("rheo_kernel");
  if (icompute == -1) error->all(FLERR, "Using compute/RHEO/chi without compute/RHEO/kernel");

  compute_kernel = ((ComputeRHEOKernel *) modify->compute[icompute]);

  //Store persistent per atom quantities - need to be exchanged
  char **fixarg = new char*[8];
  fixarg[0] = (char *) "PROPERTY_ATOM_COMP_RHEO_SOLIDS";
  fixarg[1] = (char *) "all";
  fixarg[2] = (char *) "property/atom";
  fixarg[3] = (char *) "d_fx";
  fixarg[4] = (char *) "d_fy";
  fixarg[5] = (char *) "d_fz";
  fixarg[6] = (char *) "ghost";
  fixarg[7] = (char *) "yes";
  modify->add_fix(8,fixarg,1);
  delete [] fixarg;

  int temp_flag;
  index_fx = atom->find_custom("fx", temp_flag);
  if ((index_fx < 0) || (temp_flag != 1))
      error->all(FLERR, "Compute rheo/solids can't find fix property/atom fx");
  index_fy = atom->find_custom("fy", temp_flag);
  if ((index_fy < 0) || (temp_flag != 1))
      error->all(FLERR, "Compute rheo/solids can't find fix property/atom fy");
  index_fz = atom->find_custom("fz", temp_flag);
  if ((index_fz < 0) || (temp_flag != 1))
      error->all(FLERR, "Compute rheo/solids can't find fix property/atom fz");
}

/* ---------------------------------------------------------------------- */

void ComputeRHEOInterface::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

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

