// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Trung Nguyen (Northwestern)
------------------------------------------------------------------------- */

#include "msm_dielectric.h"

#include "atom.h"
#include "atom_vec_dielectric.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "gridcomm.h"
#include "math_const.h"
#include "memory.h"

using namespace LAMMPS_NS;
using namespace MathConst;

enum{REVERSE_RHO,REVERSE_AD,REVERSE_AD_PERATOM};
enum{FORWARD_RHO,FORWARD_AD,FORWARD_AD_PERATOM};
/* ---------------------------------------------------------------------- */

MSMDielectric::MSMDielectric(LAMMPS *_lmp) : MSM(_lmp)
{
  efield = nullptr;
  phi = nullptr;
}

/* ----------------------------------------------------------------------
   free all memory
------------------------------------------------------------------------- */

MSMDielectric::~MSMDielectric()
{
  memory->destroy(efield);
  memory->destroy(phi);
}

/* ----------------------------------------------------------------------
   called once before run
------------------------------------------------------------------------- */

void MSMDielectric::init()
{
  MSM::init();

  avec = dynamic_cast<AtomVecDielectric *>(atom->style_match("dielectric"));
  if (!avec) error->all(FLERR,"msm/dielectric requires atom style dielectric");
}

/* ----------------------------------------------------------------------
   compute the MSMDielectric long-range force, energy, virial
------------------------------------------------------------------------- */

void MSMDielectric::compute(int eflag, int vflag)
{
  int i,j;

  // set energy/virial flags

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = evflag_atom = eflag_global = vflag_global =
    eflag_atom = vflag_atom = eflag_either = vflag_either = 0;

  if (scalar_pressure_flag && vflag_either) {
    if (vflag_atom)
      error->all(FLERR,"Must use 'kspace_modify pressure/scalar no' to obtain "
        "per-atom virial with kspace_style msm/dielectric");

    // must switch on global energy computation if not already on

    if (eflag == 0 || eflag == 2) {
      eflag++;
      ev_setup(eflag,vflag);
    }
  }

  // if atom count has changed, update qsum and qsqsum

  if (atom->natoms != natoms_original) {
    qsum_qsq();
    natoms_original = atom->natoms;
  }

  // return if there are no charges

  if (qsqsum == 0.0) return;

  // invoke allocate_peratom() if needed for first time

  if (vflag_atom && !peratom_allocate_flag) allocate_peratom();

  // convert atoms from box to lamda coords

  if (triclinic)
    domain->x2lamda(atom->nlocal);

  // extend size of per-atom arrays if necessary

  if (atom->nmax > nmax) {
    memory->destroy(part2grid);
    memory->destroy(efield);
    nmax = atom->nmax;
    memory->create(part2grid,nmax,3,"msm:part2grid");
    memory->create(efield,nmax,3,"msm:efield");
  }

  // find grid points for all my particles
  // map my particle charge onto my local 3d density grid (aninterpolation)

  particle_map();
  make_rho();

  // all procs reverse communicate charge density values from
  // their ghost grid points
  // to fully sum contribution in their 3d grid

  current_level = 0;
  gcall->reverse_comm(GridComm::KSPACE,this,1,sizeof(double),REVERSE_RHO,
                      gcall_buf1,gcall_buf2,MPI_DOUBLE);

  // forward communicate charge density values to fill ghost grid points
  // compute direct sum interaction and then restrict to coarser grid

  for (int n=0; n<=levels-2; n++) {
    if (!active_flag[n]) continue;
    current_level = n;
    gc[n]->forward_comm(GridComm::KSPACE,this,1,sizeof(double),FORWARD_RHO,
                        gc_buf1[n],gc_buf2[n],MPI_DOUBLE);
    direct(n);
    restriction(n);
  }

  // compute direct interation for top grid level for non-periodic
  //   and for second from top grid level for periodic

  if (active_flag[levels-1]) {
    if (domain->nonperiodic) {
      current_level = levels-1;
      gc[levels-1]->
        forward_comm(GridComm::KSPACE,this,1,sizeof(double),FORWARD_RHO,
                     gc_buf1[levels-1],gc_buf2[levels-1],MPI_DOUBLE);
      direct_top(levels-1);
      gc[levels-1]->
        reverse_comm(GridComm::KSPACE,this,1,sizeof(double),REVERSE_AD,
                     gc_buf1[levels-1],gc_buf2[levels-1],MPI_DOUBLE);
      if (vflag_atom)
        gc[levels-1]->
          reverse_comm(GridComm::KSPACE,this,6,sizeof(double),REVERSE_AD_PERATOM,
                       gc_buf1[levels-1],gc_buf2[levels-1],MPI_DOUBLE);

    } else {
      // Here using MPI_Allreduce is cheaper than using commgrid
      grid_swap_forward(levels-1,qgrid[levels-1]);
      direct(levels-1);
      grid_swap_reverse(levels-1,egrid[levels-1]);
      current_level = levels-1;
      if (vflag_atom)
        gc[levels-1]->
          reverse_comm(GridComm::KSPACE,this,6,sizeof(double),REVERSE_AD_PERATOM,
                       gc_buf1[levels-1],gc_buf2[levels-1],MPI_DOUBLE);
    }
  }

  // prolongate energy/virial from coarser grid to finer grid
  // reverse communicate from ghost grid points to get full sum

  for (int n=levels-2; n>=0; n--) {
    if (!active_flag[n]) continue;
    prolongation(n);

    current_level = n;
    gc[n]->reverse_comm(GridComm::KSPACE,this,1,sizeof(double),REVERSE_AD,
                        gc_buf1[n],gc_buf2[n],MPI_DOUBLE);

    // extra per-atom virial communication

    if (vflag_atom)
      gc[n]->reverse_comm(GridComm::KSPACE,this,6,sizeof(double),
                          REVERSE_AD_PERATOM,gc_buf1[n],gc_buf2[n],MPI_DOUBLE);
  }

  // all procs communicate E-field values
  // to fill ghost cells surrounding their 3d bricks

  current_level = 0;
  gcall->forward_comm(GridComm::KSPACE,this,1,sizeof(double),FORWARD_AD,
                      gcall_buf1,gcall_buf2,MPI_DOUBLE);

  // extra per-atom energy/virial communication

  if (vflag_atom)
    gcall->forward_comm(GridComm::KSPACE,this,6,sizeof(double),FORWARD_AD_PERATOM,
                        gcall_buf1,gcall_buf2,MPI_DOUBLE);

  // calculate the force on my particles (interpolation)

  fieldforce();

  // calculate the per-atom energy/virial for my particles

  if (evflag_atom) fieldforce_peratom();

  // sum global energy across procs and add in self-energy term

  const double qscale = qqrd2e * scale;

  if (eflag_global) {
    double energy_all;
    MPI_Allreduce(&energy,&energy_all,1,MPI_DOUBLE,MPI_SUM,world);
    energy = energy_all;

    double e_self = qsqsum*gamma(0.0)/cutoff;
    energy -= e_self;
    energy *= 0.5*qscale;
  }

  // total long-range virial

  if (vflag_global && !scalar_pressure_flag) {
    double virial_all[6];
    MPI_Allreduce(virial,virial_all,6,MPI_DOUBLE,MPI_SUM,world);
    for (i = 0; i < 6; i++) virial[i] = 0.5*qscale*virial_all[i];
  }

  // fast compute of scalar pressure (if requested)

  if (scalar_pressure_flag && vflag_global)
    for (i = 0; i < 3; i++) virial[i] = energy/3.0;

  // per-atom energy/virial
  // energy includes self-energy correction

  if (evflag_atom) {
    double *q = atom->q;
    int nlocal = atom->nlocal;

    if (eflag_atom) {
      for (i = 0; i < nlocal; i++) {
        eatom[i] -= q[i]*q[i]*gamma(0.0)/cutoff;
        eatom[i] *= 0.5*qscale;
      }
    }

    if (vflag_atom) {
      for (i = 0; i < nlocal; i++)
        for (j = 0; j < 6; j++) vatom[i][j] *= 0.5*qscale;
    }
  }

  // convert atoms back from lamda to box coords

  if (triclinic) domain->lamda2x(atom->nlocal);
}

/* ----------------------------------------------------------------------
   interpolate from grid to get force on my particles
------------------------------------------------------------------------- */

void MSMDielectric::fieldforce()
{
  double ***egridn = egrid[0];

  int i,l,m,n,nx,ny,nz,mx,my,mz;
  double dx,dy,dz;
  double phi_x,phi_y,phi_z,u;
  double dphi_x,dphi_y,dphi_z;
  double ekx,eky,ekz,etmp;


  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of E-field on particle

  double *q = atom->q;
  double **x = atom->x;
  double **f = atom->f;
  double *eps = atom->epsilon;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx - (x[i][0]-boxlo[0])*delxinv[0];
    dy = ny - (x[i][1]-boxlo[1])*delyinv[0];
    dz = nz - (x[i][2]-boxlo[2])*delzinv[0];

    compute_phis_and_dphis(dx,dy,dz);

    u = ekx = eky = ekz = 0.0;
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
          etmp = egridn[mz][my][mx];
          u += phi_z*phi_y*phi_x*etmp;
          ekx += dphi_x*phi_y*phi_z*etmp;
          eky += phi_x*dphi_y*phi_z*etmp;
          ekz += phi_x*phi_y*dphi_z*etmp;
        }
      }
    }

    ekx *= delxinv[0];
    eky *= delyinv[0];
    ekz *= delzinv[0];

    // electrical potential

    phi[i] = u;

    // effectively divide by length for a triclinic system

    if (triclinic) {
      double tmp[3];
      tmp[0] = ekx;
      tmp[1] = eky;
      tmp[2] = ekz;
      x2lamdaT(&tmp[0],&tmp[0]);
      ekx = tmp[0];
      eky = tmp[1];
      ekz = tmp[2];
    }

    // convert E-field to force
    const double efactor = scale * eps[i];
    efield[i][0] = efactor*ekx;
    efield[i][1] = efactor*eky;
    efield[i][2] = efactor*ekz;

    const double qfactor = qqrd2e*scale*q[i];
    f[i][0] += qfactor*ekx;
    f[i][1] += qfactor*eky;
    f[i][2] += qfactor*ekz;
  }
}
