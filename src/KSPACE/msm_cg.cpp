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
   Contributing authors: Paul Crozier, Stan Moore, Stephen Bond, (all SNL)
------------------------------------------------------------------------- */

#include "lmptype.h"
#include "mpi.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "atom.h"
#include "gridcomm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "neighbor.h"
#include "memory.h"
#include "msm_cg.h"

#include "math_const.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define OFFSET 16384
#define SMALLQ 0.00001

enum{REVERSE_RHO,REVERSE_AD,REVERSE_AD_PERATOM};
enum{FORWARD_RHO,FORWARD_AD,FORWARD_AD_PERATOM};

/* ---------------------------------------------------------------------- */

MSMCG::MSMCG(LAMMPS *lmp, int narg, char **arg) : MSM(lmp, narg, arg)
{
  if ((narg < 1) || (narg > 2))
    error->all(FLERR,"Illegal kspace_style msm/cg command");

  triclinic_support = 0;

  if (narg == 2) smallq = fabs(force->numeric(FLERR,arg[1]));
  else smallq = SMALLQ;

  num_charged = -1;
  is_charged = NULL;
}

/* ----------------------------------------------------------------------
   free all memory
------------------------------------------------------------------------- */

MSMCG::~MSMCG()
{
  memory->destroy(is_charged);
}

/* ----------------------------------------------------------------------
   compute the MSM long-range force, energy, virial
------------------------------------------------------------------------- */

void MSMCG::compute(int eflag, int vflag)
{
  if (scalar_pressure_flag)
    error->all(FLERR,"Must use 'kspace_modify pressure/scalar no' with "
      "kspace_style msm/cg");

  const double * const q = atom->q;
  const int nlocal = atom->nlocal;
  int i,j,n;

  // set energy/virial flags

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = evflag_atom = eflag_global = vflag_global =
    eflag_atom = vflag_atom = eflag_either = vflag_either = 0;

  // invoke allocate_peratom() if needed for first time

  if (vflag_atom && !peratom_allocate_flag) {
    allocate_peratom();
    cg_peratom_all->ghost_notify();
    cg_peratom_all->setup();
    for (int n=0; n<levels; n++) {
      if (!active_flag[n]) continue;
      cg_peratom[n]->ghost_notify();
      cg_peratom[n]->setup();
    }
    peratom_allocate_flag = 1;
  }

  // extend size of per-atom arrays if necessary

  if (nlocal > nmax) {
    memory->destroy(part2grid);
    memory->destroy(is_charged);
    nmax = atom->nmax;
    memory->create(part2grid,nmax,3,"msm:part2grid");
    memory->create(is_charged,nmax,"msm/cg:is_charged");
  }

  // one time setup message

  if (num_charged < 0) {
    bigint charged_all, charged_num;
    double charged_frac, charged_fmax, charged_fmin;

    num_charged=0;
    for (i=0; i < nlocal; ++i)
      if (fabs(q[i]) > smallq)
        ++num_charged;

    // get fraction of charged particles per domain

    if (nlocal > 0)
      charged_frac = static_cast<double>(num_charged) * 100.0
                   / static_cast<double>(nlocal);
    else
      charged_frac = 0.0;

    MPI_Reduce(&charged_frac,&charged_fmax,1,MPI_DOUBLE,MPI_MAX,0,world);
    MPI_Reduce(&charged_frac,&charged_fmin,1,MPI_DOUBLE,MPI_MIN,0,world);

    // get fraction of charged particles overall

    charged_num = num_charged;
    MPI_Reduce(&charged_num,&charged_all,1,MPI_LMP_BIGINT,MPI_SUM,0,world);
    charged_frac = static_cast<double>(charged_all) * 100.0
                   / static_cast<double>(atom->natoms);

    if (me == 0) {
      if (screen)
        fprintf(screen,
                "  MSM/cg optimization cutoff: %g\n"
                "  Total charged atoms: %.1f%%\n"
                "  Min/max charged atoms/proc: %.1f%% %.1f%%\n",
                smallq,charged_frac,charged_fmin,charged_fmax);
      if (logfile)
        fprintf(logfile,
                "  MSM/cg optimization cutoff: %g\n"
                "  Total charged atoms: %.1f%%\n"
                "  Min/max charged atoms/proc: %.1f%% %.1f%%\n",
                smallq,charged_frac,charged_fmin,charged_fmax);
    }
  }

  // only need to rebuild this list after a neighbor list update
  if (neighbor->ago == 0) {
    num_charged = 0;
    for (i = 0; i < nlocal; ++i) {
      if (fabs(q[i]) > smallq) {
        is_charged[num_charged] = i;
        ++num_charged;
      }
    }
  }

  // find grid points for all my particles
  // map my particle charge onto my local 3d density grid (aninterpolation)

  particle_map();
  make_rho();

  // all procs reverse communicate charge density values from their ghost grid points
  //   to fully sum contribution in their 3d grid

  current_level = 0;
  cg_all->reverse_comm(this,REVERSE_RHO);

  // forward communicate charge density values to fill ghost grid points
  // compute direct sum interaction and then restrict to coarser grid

  for (int n=0; n<=levels-2; n++) {
    if (!active_flag[n]) continue;
    current_level = n;
    cg[n]->forward_comm(this,FORWARD_RHO);

    direct(n);
    restriction(n);
  }


  // compute direct interation for top grid level for nonperiodic
  //   and for second from top grid level for periodic

  if (active_flag[levels-1]) {
    if (domain->nonperiodic) {
      current_level = levels-1;
      cg[levels-1]->forward_comm(this,FORWARD_RHO);
      direct_top(levels-1);
      cg[levels-1]->reverse_comm(this,REVERSE_AD);
      if (vflag_atom)
        cg_peratom[levels-1]->reverse_comm(this,REVERSE_AD_PERATOM);
    } else {
      // Here using MPI_Allreduce is cheaper than using commgrid
      grid_swap_forward(levels-1,qgrid[levels-1]);
      direct(levels-1);
      grid_swap_reverse(levels-1,egrid[levels-1]);
      current_level = levels-1;
      if (vflag_atom)
        cg_peratom[levels-1]->reverse_comm(this,REVERSE_AD_PERATOM);
    }
  }

  // prolongate energy/virial from coarser grid to finer grid
  // reverse communicate from ghost grid points to get full sum

  for (int n=levels-2; n>=0; n--) {
    if (!active_flag[n]) continue;
    prolongation(n);

    current_level = n;
    cg[n]->reverse_comm(this,REVERSE_AD);

    // extra per-atom virial communication

    if (vflag_atom)
      cg_peratom[n]->reverse_comm(this,REVERSE_AD_PERATOM);
  }

  // all procs communicate E-field values
  // to fill ghost cells surrounding their 3d bricks

  current_level = 0;
  cg_all->forward_comm(this,FORWARD_AD);

  // extra per-atom energy/virial communication

  if (vflag_atom)
    cg_peratom_all->forward_comm(this,FORWARD_AD_PERATOM);

  // calculate the force on my particles (interpolation)

  fieldforce();

  // calculate the per-atom energy/virial for my particles

  if (evflag_atom) fieldforce_peratom();

  // update qsum and qsqsum, if needed

  if (eflag_global || eflag_atom) {
    if (qsum_update_flag || (atom->natoms != natoms_original)) {
      qsum_qsq(0);
      natoms_original = atom->natoms;
    }
  }

  // sum global energy across procs and add in self-energy term

  const double qscale = force->qqrd2e * scale;

  if (eflag_global) {
    double energy_all;
    MPI_Allreduce(&energy,&energy_all,1,MPI_DOUBLE,MPI_SUM,world);
    energy = energy_all;

    double e_self = qsqsum*gamma(0.0)/cutoff;
    energy -= e_self;
    energy *= 0.5*qscale;
  }

  // total long-range virial

  if (vflag_global) {
    double virial_all[6];
    MPI_Allreduce(virial,virial_all,6,MPI_DOUBLE,MPI_SUM,world);
    for (i = 0; i < 6; i++) virial[i] = 0.5*qscale*virial_all[i];
  }

  // per-atom energy/virial
  // energy includes self-energy correction

  if (evflag_atom) {
    const double qs = 0.5*qscale;

    if (eflag_atom) {
      const double sf = gamma(0.0)/cutoff;
      for (j = 0; j < num_charged; j++) {
        i = is_charged[j];
        eatom[i] -= q[i]*q[i]*sf;
        eatom[i] *= qs;
      }
    }

    if (vflag_atom) {
      for (n = 0; n < num_charged; n++) {
        i = is_charged[n];
        for (j = 0; j < 6; j++)
          vatom[i][j] *= qs;
      }
    }
  }

}

/* ----------------------------------------------------------------------
   find center grid pt for each of my particles
   check that full stencil for the particle will fit in my 3d brick
   store central grid pt indices in part2grid array
------------------------------------------------------------------------- */

void MSMCG::particle_map()
{
  const double * const * const x = atom->x;

  int flag = 0;
  int i;

  if (!isfinite(boxlo[0]) || !isfinite(boxlo[1]) || !isfinite(boxlo[2]))
    error->one(FLERR,"Non-numeric box dimensions - simulation unstable");

  for (int j = 0; j < num_charged; j++) {
    i = is_charged[j];

    // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
    // current particle coord can be outside global and local box
    // add/subtract OFFSET to avoid int(-0.75) = 0 when want it to be -1

    const int nx=static_cast<int>((x[i][0]-boxlo[0])*delxinv[0]+OFFSET)-OFFSET;
    const int ny=static_cast<int>((x[i][1]-boxlo[1])*delyinv[0]+OFFSET)-OFFSET;
    const int nz=static_cast<int>((x[i][2]-boxlo[2])*delzinv[0]+OFFSET)-OFFSET;

    part2grid[i][0] = nx;
    part2grid[i][1] = ny;
    part2grid[i][2] = nz;

    // check that entire stencil around nx,ny,nz will fit in my 3d brick

    if (nx+nlower < nxlo_out[0] || nx+nupper > nxhi_out[0] ||
        ny+nlower < nylo_out[0] || ny+nupper > nyhi_out[0] ||
        nz+nlower < nzlo_out[0] || nz+nupper > nzhi_out[0])
      flag = 1;
  }

  if (flag) error->one(FLERR,"Out of range atoms - cannot compute MSM");
}

/* ----------------------------------------------------------------------
   create discretized "density" on section of global grid due to my particles
   density(x,y,z) = charge "density" at grid points of my 3d brick
   (nxlo:nxhi,nylo:nyhi,nzlo:nzhi) is extent of my brick (including ghosts)
   in global grid
------------------------------------------------------------------------- */

void MSMCG::make_rho()
{
  const double * const q = atom->q;
  const double * const * const x = atom->x;

  // clear 3d density array

  double * const * const * const qgridn = qgrid[0];

  memset(&(qgridn[nzlo_out[0]][nylo_out[0]][nxlo_out[0]]),0,ngrid[0]*sizeof(double));

  double dx,dy,dz,x0,y0,z0;
  int i,j,l,m,n,nx,ny,nz,mx,my,mz;

  // loop over my charges, add their contribution to nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt

  for (j = 0; j < num_charged; j++) {
    i = is_charged[j];

    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx - (x[i][0]-boxlo[0])*delxinv[0];
    dy = ny - (x[i][1]-boxlo[1])*delyinv[0];
    dz = nz - (x[i][2]-boxlo[2])*delzinv[0];

    compute_phis(dx,dy,dz);

    z0 = q[i];
    for (n = nlower; n <= nupper; n++) {
      mz = n+nz;
      y0 = z0*phi1d[2][n];
      for (m = nlower; m <= nupper; m++) {
        my = m+ny;
        x0 = y0*phi1d[1][m];
        for (l = nlower; l <= nupper; l++) {
          mx = l+nx;
          qgridn[mz][my][mx] += x0*phi1d[0][l];
        }
      }
    }
  }

}

/* ----------------------------------------------------------------------
   interpolate from grid to get force on my particles
------------------------------------------------------------------------- */

void MSMCG::fieldforce()
{

  const double * const * const * const egridn = egrid[0];
  const double * const * const x = atom->x;
  double * const * const f = atom->f;
  const double * const q = atom->q;

  int i,j,l,m,n,nx,ny,nz,mx,my,mz;
  double dx,dy,dz;
  double phi_x,phi_y,phi_z;
  double dphi_x,dphi_y,dphi_z;
  double ekx,eky,ekz;


  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of E-field on particle

  for (j = 0; j < num_charged; j++) {
    i = is_charged[j];
    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx - (x[i][0]-boxlo[0])*delxinv[0];
    dy = ny - (x[i][1]-boxlo[1])*delyinv[0];
    dz = nz - (x[i][2]-boxlo[2])*delzinv[0];

    compute_phis_and_dphis(dx,dy,dz);

    ekx = eky = ekz = 0.0;
    for (n = nlower; n <= nupper; n++) {
      mz = n+nz;
      phi_z = phi1d[2][n];
      dphi_z = dphi1d[2][n];
      for (m = nlower; m <= nupper; m++) {
        my = m+ny;
        phi_y = phi1d[1][m];
        dphi_y = dphi1d[1][m];
        for (l = nlower; l <= nupper; l++) {
          mx = l+nx;
          phi_x = phi1d[0][l];
          dphi_x = dphi1d[0][l];
          ekx += dphi_x*phi_y*phi_z*egridn[mz][my][mx];
          eky += phi_x*dphi_y*phi_z*egridn[mz][my][mx];
          ekz += phi_x*phi_y*dphi_z*egridn[mz][my][mx];
        }
      }
    }

    ekx *= delxinv[0];
    eky *= delyinv[0];
    ekz *= delzinv[0];

    // convert E-field to force

    const double qfactor = force->qqrd2e*scale*q[i];
    f[i][0] += qfactor*ekx;
    f[i][1] += qfactor*eky;
    f[i][2] += qfactor*ekz;
  }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get per-atom energy/virial
------------------------------------------------------------------------- */

void MSMCG::fieldforce_peratom()
{
  const double * const q = atom->q;
  const double * const * const x = atom->x;

  double ***egridn = egrid[0];

  double ***v0gridn = v0grid[0];
  double ***v1gridn = v1grid[0];
  double ***v2gridn = v2grid[0];
  double ***v3gridn = v3grid[0];
  double ***v4gridn = v4grid[0];
  double ***v5gridn = v5grid[0];

  int i,j,l,m,n,nx,ny,nz,mx,my,mz;
  double dx,dy,dz,x0,y0,z0;
  double u,v0,v1,v2,v3,v4,v5;

  // loop over my charges, interpolate from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt

  for (j = 0; j < num_charged; j++) {
    i = is_charged[j];
    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx - (x[i][0]-boxlo[0])*delxinv[0];
    dy = ny - (x[i][1]-boxlo[1])*delyinv[0];
    dz = nz - (x[i][2]-boxlo[2])*delzinv[0];

    compute_phis_and_dphis(dx,dy,dz);

    u = v0 = v1 = v2 = v3 = v4 = v5 = 0.0;
    for (n = nlower; n <= nupper; n++) {
      mz = n+nz;
      z0 = phi1d[2][n];
      for (m = nlower; m <= nupper; m++) {
        my = m+ny;
        y0 = z0*phi1d[1][m];
        for (l = nlower; l <= nupper; l++) {
          mx = l+nx;
          x0 = y0*phi1d[0][l];
          if (eflag_atom) u += x0*egridn[mz][my][mx];
          if (vflag_atom) {
            v0 += x0*v0gridn[mz][my][mx];
            v1 += x0*v1gridn[mz][my][mx];
            v2 += x0*v2gridn[mz][my][mx];
            v3 += x0*v3gridn[mz][my][mx];
            v4 += x0*v4gridn[mz][my][mx];
            v5 += x0*v5gridn[mz][my][mx];
          }
        }
      }
    }

    if (eflag_atom) eatom[i] += q[i]*u;
    if (vflag_atom) {
      vatom[i][0] += q[i]*v0;
      vatom[i][1] += q[i]*v1;
      vatom[i][2] += q[i]*v2;
      vatom[i][3] += q[i]*v3;
      vatom[i][4] += q[i]*v4;
      vatom[i][5] += q[i]*v5;
    }
  }
}


double MSMCG::memory_usage()
{
  double bytes = MSM::memory_usage();
  bytes += nmax * sizeof(int);
  return bytes;
}
