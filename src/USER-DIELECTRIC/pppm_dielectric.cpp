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
   Contributing authors: Trung Nguyen (Northwestern)
     point-dipoles by Stan Moore (SNL)
------------------------------------------------------------------------- */

#include <mpi.h>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "pppm_dielectric.h"
#include "atom.h"
#include "atom_vec_dielectric.h"
#include "comm.h"
#include "gridcomm.h"
#include "neighbor.h"
#include "force.h"
#include "pair.h"
#include "domain.h"
#include "fft3d_wrap.h"
#include "remap_wrap.h"
#include "memory.h"
#include "error.h"

#include "math_const.h"
#include "math_special.h"

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathSpecial;

#define SMALL 0.00001

#ifdef FFT_SINGLE
#define ZEROF 0.0f
#define ONEF  1.0f
#else
#define ZEROF 0.0
#define ONEF  1.0
#endif

enum{REVERSE_RHO,REVERSE_MU};
enum{FORWARD_IK,FORWARD_MU,FORWARD_AD,FORWARD_IK_PERATOM,FORWARD_MU_PERATOM,FORWARD_AD_PERATOM};

/* ---------------------------------------------------------------------- */

PPPMDielectric::PPPMDielectric(LAMMPS *lmp) : PPPM(lmp), 
  cg_mu(NULL), densityx_brick_dipole(NULL), densityy_brick_dipole(NULL), densityz_brick_dipole(NULL),
  u_brick_dipole(NULL), ux_brick_dipole(NULL), uy_brick_dipole(NULL), uz_brick_dipole(NULL),
  vdxx_brick_dipole(NULL), vdxy_brick_dipole(NULL),  vdyy_brick_dipole(NULL), vdxz_brick_dipole(NULL),
  vdyz_brick_dipole(NULL), vdzz_brick_dipole(NULL), work3(NULL), work4(NULL),
  densityx_fft_dipole(NULL), densityy_fft_dipole(NULL), densityz_fft_dipole(NULL)
{
  dipoleflag = 0; // turned off for now, until dipole works
  group_group_enable = 0;

  mu_flag = 0;

  efield = NULL;
  phi = NULL;
  potflag = 0;

  avec = (AtomVecDielectric *) atom->style_match("dielectric");
  if (!avec) error->all(FLERR,"pppm/dielectric requires atom style dielectric");
}

/* ---------------------------------------------------------------------- */

PPPMDielectric::~PPPMDielectric()
{
  memory->destroy(efield);
  memory->destroy(phi);
}

/* ----------------------------------------------------------------------
   called once before run
------------------------------------------------------------------------- */

void PPPMDielectric::init()
{
  PPPM::init();

  if (mu_flag) musum_musq();

  if (mu_flag) {
    cg_mu->ghost_notify();
    cg_mu->setup();
  }
}

/* ----------------------------------------------------------------------
   reset local grid arrays and communication stencils
   called by fix balance b/c it changed sizes of processor sub-domains
------------------------------------------------------------------------- */

void PPPMDielectric::setup_grid()
{
  PPPM::setup_grid();

  if (mu_flag) {
    cg_mu->ghost_notify();
    if (overlap_allowed == 0 && cg_mu->ghost_overlap())
      error->all(FLERR,"PPPMDielectric grid stencil extends "
                 "beyond nearest neighbor processor");
    cg_mu->setup();
  }

  // pre-compute volume-dependent coeffs

  PPPM::setup();
}

/* ----------------------------------------------------------------------
   compute the PPPMDielectric long-range force, energy, virial
------------------------------------------------------------------------- */

void PPPMDielectric::compute(int eflag, int vflag)
{
  int i,j;

  // set energy/virial flags
  // invoke allocate_peratom() if needed for first time

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = evflag_atom = eflag_global = vflag_global =
         eflag_atom = vflag_atom = 0;

  if (potflag) evflag_atom = 1;

  if (evflag_atom && !peratom_allocate_flag) {
    allocate_peratom();
    cg_peratom->ghost_notify();
    cg_peratom->setup();
  }

  // if atom count has changed, update qsum and qsqsum

  if (atom->natoms != natoms_original) {
    qsum_qsq();
    musum_musq();
    natoms_original = atom->natoms;
  }

  // return if there are no charges or dipoles

  //if (qsqsum == 0.0 && musqsum == 0.0) return;

  // convert atoms from box to lamda coords

  if (triclinic == 0) boxlo = domain->boxlo;
  else {
    boxlo = domain->boxlo_lamda;
    domain->x2lamda(atom->nlocal);
  }

  // extend size of per-atom arrays if necessary

  if (atom->nmax > nmax) {
    memory->destroy(part2grid);
    memory->destroy(efield);
    memory->destroy(phi);
    nmax = atom->nmax;
    memory->create(part2grid,nmax,3,"pppm/dielectric:part2grid");
    memory->create(efield,nmax,3,"pppm/dielectric:efield");
    memory->create(phi,nmax,"pppm/dielectric:phi");
  }

  // find grid points for all my particles
  // map my particle charge onto my local 3d density grid

  particle_map();
  make_rho();

  if (mu_flag)
    make_rho_dipole();

  // all procs communicate density values from their ghost cells
  //   to fully sum contribution in their 3d bricks
  // remap from 3d decomposition to FFT decomposition

  cg->reverse_comm(this,REVERSE_RHO);
  brick2fft();

  if (mu_flag) {
    cg_mu->reverse_comm(this,REVERSE_MU);
  	brick2fft_dipole();
  }

  // compute potential gradient on my FFT grid and
  //   portion of e_long on this proc's FFT grid
  // return gradients (electric fields) in 3d brick decomposition
  // also performs per-atom calculations via poisson_peratom()

  poisson();

  if (mu_flag)
    poisson_ik_dipole();

  // all procs communicate E-field values
  // to fill ghost cells surrounding their 3d bricks

  if (differentiation_flag == 1) cg->forward_comm(this,FORWARD_AD);
  else cg->forward_comm(this,FORWARD_IK);

  if (mu_flag)
    cg_mu->forward_comm(this,FORWARD_MU);

  // extra per-atom energy/virial communication

  if (evflag_atom) {
    if (differentiation_flag == 1 && vflag_atom)
      cg_peratom->forward_comm(this,FORWARD_AD_PERATOM);
    else if (differentiation_flag == 0)
      cg_peratom->forward_comm(this,FORWARD_IK_PERATOM);
  }

  // calculate the force on my particles

  fieldforce();

  if (mu_flag)
    fieldforce_ik_dipole();

  // extra per-atom energy/virial communication

  if (evflag_atom) fieldforce_peratom();

  // sum global energy across procs and add in volume-dependent term

  const double qscale = qqrd2e * scale;
  const double g3 = g_ewald*g_ewald*g_ewald;

  if (eflag_global) {
    double energy_all;
    MPI_Allreduce(&energy,&energy_all,1,MPI_DOUBLE,MPI_SUM,world);
    energy = energy_all;

    energy *= 0.5*volume;
    energy -= g_ewald*qsqsum/MY_PIS +
      MY_PI2*qsum*qsum / (g_ewald*g_ewald*volume);
    if (mu_flag)
      energy -= musqsum*qqrd2e*2.0*g3/3.0/MY_PIS;
    energy *= qscale;
  }

  // sum global virial across procs

  if (vflag_global) {
    double virial_all[6];
    MPI_Allreduce(virial,virial_all,6,MPI_DOUBLE,MPI_SUM,world);
    for (i = 0; i < 6; i++) virial[i] = 0.5*qscale*volume*virial_all[i];
  }

  // per-atom energy/virial
  // energy includes self-energy correction
  // ntotal accounts for TIP4P tallying eatom/vatom for ghost atoms

  if (evflag_atom) {
    double *q = atom->q;
    double **mu = atom->mu;
    int nlocal = atom->nlocal;
    int ntotal = nlocal;
    if (tip4pflag) ntotal += atom->nghost;

    if (eflag_atom) {
      for (i = 0; i < nlocal; i++) {
        eatom[i] *= 0.5;
        eatom[i] -= g_ewald*q[i]*q[i]/MY_PIS + MY_PI2*q[i]*qsum /
          (g_ewald*g_ewald*volume);

      if (mu_flag)
        eatom[i] -= (mu[i][0]*mu[i][0] + mu[i][1]*mu[i][1] + mu[i][2]*mu[i][2])*qqrd2e*2.0*g3/3.0/MY_PIS;

        eatom[i] *= qscale;
      }
      for (i = nlocal; i < ntotal; i++) eatom[i] *= 0.5*qscale;
    }

    if (vflag_atom) {
      for (i = 0; i < ntotal; i++)
        for (j = 0; j < 6; j++) vatom[i][j] *= 0.5*qscale;
    }
  }

  // 2d slab correction

  if (slabflag == 1) slabcorr();

  // convert atoms back from lamda to box coords

  if (triclinic) domain->lamda2x(atom->nlocal);
}

/* ----------------------------------------------------------------------
   allocate memory that depends on # of K-vectors and order
------------------------------------------------------------------------- */

void PPPMDielectric::allocate()
{
  PPPM::allocate();

  if (mu_flag) {
    memory->create3d_offset(densityx_brick_dipole,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                            nxlo_out,nxhi_out,"pppm:densityx_brick_dipole");
    memory->create3d_offset(densityy_brick_dipole,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                            nxlo_out,nxhi_out,"pppm:densityy_brick_dipole");
    memory->create3d_offset(densityz_brick_dipole,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                            nxlo_out,nxhi_out,"pppm:densityz_brick_dipole");

    memory->create(densityx_fft_dipole,nfft_both,"pppm:densityy_fft_dipole");
    memory->create(densityy_fft_dipole,nfft_both,"pppm:densityy_fft_dipole");
    memory->create(densityz_fft_dipole,nfft_both,"pppm:densityz_fft_dipole");

    memory->create(work3,2*nfft_both,"pppm:work3");
    memory->create(work4,2*nfft_both,"pppm:work4");

    memory->create3d_offset(u_brick_dipole,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                            nxlo_out,nxhi_out,"pppm:u_brick_dipole");

    memory->create3d_offset(ux_brick_dipole,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                            nxlo_out,nxhi_out,"pppm:ux_brick_dipole");
    memory->create3d_offset(uy_brick_dipole,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                            nxlo_out,nxhi_out,"pppm:uy_brick_dipole");
    memory->create3d_offset(uz_brick_dipole,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                            nxlo_out,nxhi_out,"pppm:uz_brick_dipole");

    memory->create3d_offset(vdxx_brick_dipole,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                            nxlo_out,nxhi_out,"pppm:vdxx_brick_dipole");
    memory->create3d_offset(vdxy_brick_dipole,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                            nxlo_out,nxhi_out,"pppm:vdxy_brick_dipole");
    memory->create3d_offset(vdyy_brick_dipole,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                            nxlo_out,nxhi_out,"pppm:vdyy_brick_dipole");
    memory->create3d_offset(vdxz_brick_dipole,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                            nxlo_out,nxhi_out,"pppm:vdxz_brick_dipole");
    memory->create3d_offset(vdyz_brick_dipole,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                            nxlo_out,nxhi_out,"pppm:vdyz_brick_dipole");
    memory->create3d_offset(vdzz_brick_dipole,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                            nxlo_out,nxhi_out,"pppm:vdzz_brick_dipole");

    int (*procneigh)[2] = comm->procneigh;
    cg_mu = new GridComm(lmp,world,9,3,
                         nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in,
                         nxlo_out,nxhi_out,nylo_out,nyhi_out,nzlo_out,nzhi_out,
                         procneigh[0][0],procneigh[0][1],procneigh[1][0],
                         procneigh[1][1],procneigh[2][0],procneigh[2][1]);
  }
}

/* ----------------------------------------------------------------------
   deallocate memory that depends on # of K-vectors and order
------------------------------------------------------------------------- */

void PPPMDielectric::deallocate()
{
  PPPM::deallocate();

  if (mu_flag) {
    memory->destroy3d_offset(densityx_brick_dipole,nzlo_out,nylo_out,nxlo_out);
    memory->destroy3d_offset(densityy_brick_dipole,nzlo_out,nylo_out,nxlo_out);
    memory->destroy3d_offset(densityz_brick_dipole,nzlo_out,nylo_out,nxlo_out);

    memory->destroy3d_offset(u_brick_dipole,nzlo_out,nylo_out,nxlo_out);

    memory->destroy3d_offset(ux_brick_dipole,nzlo_out,nylo_out,nxlo_out);
    memory->destroy3d_offset(uy_brick_dipole,nzlo_out,nylo_out,nxlo_out);
    memory->destroy3d_offset(uz_brick_dipole,nzlo_out,nylo_out,nxlo_out);

    memory->destroy3d_offset(vdxx_brick_dipole,nzlo_out,nylo_out,nxlo_out);
    memory->destroy3d_offset(vdxy_brick_dipole,nzlo_out,nylo_out,nxlo_out);
    memory->destroy3d_offset(vdyy_brick_dipole,nzlo_out,nylo_out,nxlo_out);
    memory->destroy3d_offset(vdxz_brick_dipole,nzlo_out,nylo_out,nxlo_out);
    memory->destroy3d_offset(vdyz_brick_dipole,nzlo_out,nylo_out,nxlo_out);
    memory->destroy3d_offset(vdzz_brick_dipole,nzlo_out,nylo_out,nxlo_out);

    memory->destroy(densityx_fft_dipole);
    memory->destroy(densityy_fft_dipole);
    memory->destroy(densityz_fft_dipole);

    memory->destroy(work3);
    memory->destroy(work4);

    delete cg_mu;
  }
}

/* ----------------------------------------------------------------------
   create discretized "density" on section of global grid due to my particles
   density(x,y,z) = charge "density" at grid points of my 3d brick
   (nxlo:nxhi,nylo:nyhi,nzlo:nzhi) is extent of my brick (including ghosts)
   in global grid
------------------------------------------------------------------------- */

void PPPMDielectric::make_rho_dipole()
{
  int l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz;
  FFT_SCALAR x0,y0,z0;
  FFT_SCALAR x1,y1,z1;
  FFT_SCALAR x2,y2,z2;

  // clear 3d density array

  memset(&(densityx_brick_dipole[nzlo_out][nylo_out][nxlo_out]),0,
         ngrid*sizeof(FFT_SCALAR));
  memset(&(densityy_brick_dipole[nzlo_out][nylo_out][nxlo_out]),0,
         ngrid*sizeof(FFT_SCALAR));
  memset(&(densityz_brick_dipole[nzlo_out][nylo_out][nxlo_out]),0,
         ngrid*sizeof(FFT_SCALAR));

  // loop over my charges, add their contribution to nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt

  double **mu = atom->mu;
  double **x = atom->x;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {

    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv;
    dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv;
    dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv;

    compute_rho1d(dx,dy,dz);

    z0 = delvolinv * mu[i][0];
    z1 = delvolinv * mu[i][1];
    z2 = delvolinv * mu[i][2];
    for (n = nlower; n <= nupper; n++) {
      mz = n+nz;
      y0 = z0*rho1d[2][n];
      y1 = z1*rho1d[2][n];
      y2 = z2*rho1d[2][n];
      for (m = nlower; m <= nupper; m++) {
        my = m+ny;
        x0 = y0*rho1d[1][m];
        x1 = y1*rho1d[1][m];
        x2 = y2*rho1d[1][m];
        for (l = nlower; l <= nupper; l++) {
          mx = l+nx;
          densityx_brick_dipole[mz][my][mx] += x0*rho1d[0][l];
          densityy_brick_dipole[mz][my][mx] += x1*rho1d[0][l];
          densityz_brick_dipole[mz][my][mx] += x2*rho1d[0][l];
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   remap density from 3d brick decomposition to FFT decomposition
------------------------------------------------------------------------- */

void PPPMDielectric::brick2fft_dipole()
{
  int n,ix,iy,iz;

  // copy grabs inner portion of density from 3d brick
  // remap could be done as pre-stage of FFT,
  //   but this works optimally on only double values, not complex values

  n = 0;
  for (iz = nzlo_in; iz <= nzhi_in; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
        densityx_fft_dipole[n] = densityx_brick_dipole[iz][iy][ix];
        densityy_fft_dipole[n] = densityy_brick_dipole[iz][iy][ix];
        densityz_fft_dipole[n] = densityz_brick_dipole[iz][iy][ix];
        n++;
      }

  remap->perform(densityx_fft_dipole,densityx_fft_dipole,work1);
  remap->perform(densityy_fft_dipole,densityy_fft_dipole,work1);
  remap->perform(densityz_fft_dipole,densityz_fft_dipole,work1);
}

/* ----------------------------------------------------------------------
   FFT-based Poisson solver for ik
------------------------------------------------------------------------- */

void PPPMDielectric::poisson_ik_dipole()
{
  int i,j,k,n,ii;
  double eng;
  double wreal,wimg;

  // transform dipole density (r -> k)

  n = 0;
  for (i = 0; i < nfft; i++) {
    work1[n] = densityx_fft_dipole[i];
    work1[n+1] = ZEROF;
    work2[n] = densityy_fft_dipole[i];
    work2[n+1] = ZEROF;
    work3[n] = densityz_fft_dipole[i];
    work3[n+1] = ZEROF;
    n += 2;
  }

  fft1->compute(work1,work1,1);
  fft1->compute(work2,work2,1);
  fft1->compute(work3,work3,1);

  // global energy and virial contribution

  double scaleinv = 1.0/(nx_pppm*ny_pppm*nz_pppm);
  double s2 = scaleinv*scaleinv;

  if (eflag_global || vflag_global) {
    if (vflag_global) {
      n = 0;
      ii = 0;
      for (k = nzlo_fft; k <= nzhi_fft; k++)
        for (j = nylo_fft; j <= nyhi_fft; j++)
          for (i = nxlo_fft; i <= nxhi_fft; i++) {
            wreal = (work1[n]*fkx[i] + work2[n]*fky[j] + work3[n]*fkz[k]);
            wimg = (work1[n+1]*fkx[i] + work2[n+1]*fky[j] + work3[n+1]*fkz[k]);
            eng = s2 * greensfn[ii] * (wreal*wreal + wimg*wimg);
            //eng = s2 * greensfn[ii] * (wreal+wimg)*(wreal+wimg);
            for (int jj = 0; jj < 6; jj++) virial[jj] += eng*vg[ii][jj];
            if (eflag_global) energy += eng;
            ii++;
            n += 2;
          }
    } else {
      n = 0;
      ii = 0;
      for (k = nzlo_fft; k <= nzhi_fft; k++)
        for (j = nylo_fft; j <= nyhi_fft; j++)
          for (i = nxlo_fft; i <= nxhi_fft; i++) {
            wreal = (work1[n]*fkx[i] + work2[n]*fky[j] + work3[n]*fkz[k]);
            wimg = (work1[n+1]*fkx[i] + work2[n+1]*fky[j] + work3[n+1]*fkz[k]);
            energy +=
            //s2 * greensfn[ii] * (wreal+wimg)*(wreal+wimg);
            s2 * greensfn[ii] * (wreal*wreal + wimg*wimg);
	          ii++;
            n += 2;
          }
    }
  }

  // scale by 1/total-grid-pts to get rho(k)
  // multiply by Green's function to get V(k)

  n = 0;
  for (i = 0; i < nfft; i++) {
    work1[n]   *= scaleinv * greensfn[i];
    work1[n+1] *= scaleinv * greensfn[i];
    work2[n]   *= scaleinv * greensfn[i];
    work2[n+1] *= scaleinv * greensfn[i];
    work3[n]   *= scaleinv * greensfn[i];
    work3[n+1] *= scaleinv * greensfn[i];
    n += 2;
  }

  // triclinic system

  /*if (triclinic) {
    poisson_ik_triclinic();
    return;
  }*/

  // compute electric potential
  // FFT leaves data in 3d brick decomposition

  n = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        work4[n] = work1[n]*fkx[i] + work2[n]*fky[j] + work3[n]*fkz[k];
        work4[n+1] = work1[n+1]*fkx[i] + work2[n+1]*fky[j] + work3[n+1]*fkz[k];
        n += 2;
      }

  fft2->compute(work4,work4,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        u_brick_dipole[k][j][i] = work4[n];
        n += 2;
      }

  // Ex

  n = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        work4[n] = fkx[i]*(work1[n]*fkx[i] + work2[n]*fky[j] + work3[n]*fkz[k]);
        work4[n+1] = fkx[i]*(work1[n+1]*fkx[i] + work2[n+1]*fky[j] + work3[n+1]*fkz[k]);
        n += 2;
      }

  fft2->compute(work4,work4,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        ux_brick_dipole[k][j][i] = work4[n];
        n += 2;
      }

  // Ey

  n = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        work4[n] = fky[j]*(work1[n]*fkx[i] + work2[n]*fky[j] + work3[n]*fkz[k]);
        work4[n+1] = fky[j]*(work1[n+1]*fkx[i] + work2[n+1]*fky[j] + work3[n+1]*fkz[k]);
        n += 2;
      }

  fft2->compute(work4,work4,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        uy_brick_dipole[k][j][i] = work4[n];
        n += 2;
      }

  // Ez

  n = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        work4[n] = fkz[k]*(work1[n]*fkx[i] + work2[n]*fky[j] + work3[n]*fkz[k]);
        work4[n+1] = fkz[k]*(work1[n+1]*fkx[i] + work2[n+1]*fky[j] + work3[n+1]*fkz[k]);
        n += 2;
      }

  fft2->compute(work4,work4,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        uz_brick_dipole[k][j][i] = work4[n];
        n += 2;
      }

  // Vxx

  n = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        work4[n] = fkx[i]*fkx[i]*(work1[n+1]*fkx[i] + work2[n+1]*fky[j] + work3[n+1]*fkz[k]);
        work4[n+1] = -fkx[i]*fkx[i]*(work1[n]*fkx[i] + work2[n]*fky[j] + work3[n]*fkz[k]);
        n += 2;
      }

  fft2->compute(work4,work4,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        vdxx_brick_dipole[k][j][i] = work4[n];
        n += 2;
      }

  // Vyy

  n = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        work4[n] = fky[j]*fky[j]*(work1[n+1]*fkx[i] + work2[n+1]*fky[j] + work3[n+1]*fkz[k]);
        work4[n+1] = -fky[j]*fky[j]*(work1[n]*fkx[i] + work2[n]*fky[j] + work3[n]*fkz[k]);
        n += 2;
      }

  fft2->compute(work4,work4,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        vdyy_brick_dipole[k][j][i] = work4[n];
        n += 2;
      }

  // Vzz

  n = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        work4[n] = fkz[k]*fkz[k]*(work1[n+1]*fkx[i] + work2[n+1]*fky[j] + work3[n+1]*fkz[k]);
        work4[n+1] = -fkz[k]*fkz[k]*(work1[n]*fkx[i] + work2[n]*fky[j] + work3[n]*fkz[k]);
        n += 2;
      }

  fft2->compute(work4,work4,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        vdzz_brick_dipole[k][j][i] = work4[n];
        n += 2;
      }

  // Vxy

  n = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        work4[n] = fkx[i]*fky[j]*(work1[n+1]*fkx[i] + work2[n+1]*fky[j] + work3[n+1]*fkz[k]);
        work4[n+1] = -fkx[i]*fky[j]*(work1[n]*fkx[i] + work2[n]*fky[j] + work3[n]*fkz[k]);
        n += 2;
      }

  fft2->compute(work4,work4,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        vdxy_brick_dipole[k][j][i] = work4[n];
        n += 2;
      }

  // Vxz

  n = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        work4[n] = fkx[i]*fkz[k]*(work1[n+1]*fkx[i] + work2[n+1]*fky[j] + work3[n+1]*fkz[k]);
        work4[n+1] = -fkx[i]*fkz[k]*(work1[n]*fkx[i] + work2[n]*fky[j] + work3[n]*fkz[k]);
        n += 2;
      }

  fft2->compute(work4,work4,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        vdxz_brick_dipole[k][j][i] = work4[n];
        n += 2;
      }

  // Vyz

  n = 0;
  for (k = nzlo_fft; k <= nzhi_fft; k++)
    for (j = nylo_fft; j <= nyhi_fft; j++)
      for (i = nxlo_fft; i <= nxhi_fft; i++) {
        work4[n] = fky[j]*fkz[k]*(work1[n+1]*fkx[i] + work2[n+1]*fky[j] + work3[n+1]*fkz[k]);
        work4[n+1] = -fky[j]*fkz[k]*(work1[n]*fkx[i] + work2[n]*fky[j] + work3[n]*fkz[k]);
        n += 2;
      }

  fft2->compute(work4,work4,-1);

  n = 0;
  for (k = nzlo_in; k <= nzhi_in; k++)
    for (j = nylo_in; j <= nyhi_in; j++)
      for (i = nxlo_in; i <= nxhi_in; i++) {
        vdyz_brick_dipole[k][j][i] = work4[n];
        n += 2;
      }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get electric field & force on my particles for ik
------------------------------------------------------------------------- */

void PPPMDielectric::fieldforce_ik()
{
  int i,l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz,x0,y0,z0;
  FFT_SCALAR ekx,eky,ekz,u;

  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of E-field on particle

  double *q = atom->q;
  double **x = atom->x;
  double **f = atom->f;
  double *eps = avec->epsilon;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv;
    dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv;
    dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv;

    compute_rho1d(dx,dy,dz);

    u = ekx = eky = ekz = ZEROF;
    for (n = nlower; n <= nupper; n++) {
      mz = n+nz;
      z0 = rho1d[2][n];
      for (m = nlower; m <= nupper; m++) {
        my = m+ny;
        y0 = z0*rho1d[1][m];
        for (l = nlower; l <= nupper; l++) {
          mx = l+nx;
          x0 = y0*rho1d[0][l];
          if (potflag) u += x0*u_brick[mz][my][mx];
          ekx -= x0*vdx_brick[mz][my][mx];
          eky -= x0*vdy_brick[mz][my][mx];
          ekz -= x0*vdz_brick[mz][my][mx];
        }
      }
    }

    // electrostatic potential

    if (potflag) phi[i] = u;

    // convert E-field to force
    const double efactor = scale * eps[i];
    efield[i][0] = efactor*ekx;
    efield[i][1] = efactor*eky;
    efield[i][2] = efactor*ekz;

    const double qfactor = qqrd2e * scale * q[i] * eps[i];
    f[i][0] += qfactor*ekx;
    f[i][1] += qfactor*eky;
    if (slabflag != 2) f[i][2] += qfactor*ekz;
  }
}


/* ----------------------------------------------------------------------
   interpolate from grid to get electric field & force on my particles for ik
------------------------------------------------------------------------- */

void PPPMDielectric::fieldforce_ik_dipole()
{
  int i,l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz,u;
  FFT_SCALAR x0,y0,z0;
  FFT_SCALAR ex,ey,ez;
  FFT_SCALAR vxx,vyy,vzz,vxy,vxz,vyz;

  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt


  double **mu = atom->mu;
  double **x = atom->x;
  double **f = atom->f;
  double **t = atom->torque;

  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv;
    dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv;
    dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv;

    compute_rho1d(dx,dy,dz);

    u = ex = ey = ez = ZEROF;
    vxx = vyy = vzz = vxy = vxz = vyz = ZEROF;
    for (n = nlower; n <= nupper; n++) {
      mz = n+nz;
      z0 = rho1d[2][n];
      for (m = nlower; m <= nupper; m++) {
        my = m+ny;
        y0 = z0*rho1d[1][m];
        for (l = nlower; l <= nupper; l++) {
          mx = l+nx;
          x0 = y0*rho1d[0][l];
          u += x0*u_brick_dipole[mz][my][mx];
          ex -= x0*ux_brick_dipole[mz][my][mx];
          ey -= x0*uy_brick_dipole[mz][my][mx];
          ez -= x0*uz_brick_dipole[mz][my][mx];
          vxx -= x0*vdxx_brick_dipole[mz][my][mx];
          vyy -= x0*vdyy_brick_dipole[mz][my][mx];
          vzz -= x0*vdzz_brick_dipole[mz][my][mx];
          vxy -= x0*vdxy_brick_dipole[mz][my][mx];
          vxz -= x0*vdxz_brick_dipole[mz][my][mx];
          vyz -= x0*vdyz_brick_dipole[mz][my][mx];
        }
      }
    }

    // electrical potential due to dipoles

    if (potflag) phi[i] = u;

    // convert E-field to torque

    const double mufactor = qqrd2e * scale;
    f[i][0] += mufactor*(vxx*mu[i][0] + vxy*mu[i][1] + vxz*mu[i][2]);
    f[i][1] += mufactor*(vxy*mu[i][0] + vyy*mu[i][1] + vyz*mu[i][2]);
    f[i][2] += mufactor*(vxz*mu[i][0] + vyz*mu[i][1] + vzz*mu[i][2]);

    t[i][0] += mufactor*(mu[i][1]*ez - mu[i][2]*ey);
    t[i][1] += mufactor*(mu[i][2]*ex - mu[i][0]*ez);
    t[i][2] += mufactor*(mu[i][0]*ey - mu[i][1]*ex);
  }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get electric field & force on my particles for ad
------------------------------------------------------------------------- */

void PPPMDielectric::fieldforce_ad()
{
  int i,l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz;
  FFT_SCALAR ekx,eky,ekz,u;

  double s1,s2,s3;
  double sf = 0.0;
  double *prd;

  prd = domain->prd;
  double xprd = prd[0];
  double yprd = prd[1];
  double zprd = prd[2];

  double hx_inv = nx_pppm/xprd;
  double hy_inv = ny_pppm/yprd;
  double hz_inv = nz_pppm/zprd;

  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of E-field on particle

  double *q = atom->q;
  double **x = atom->x;
  double **f = atom->f;
  double *eps = avec->epsilon;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv;
    dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv;
    dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv;

    compute_rho1d(dx,dy,dz);
    compute_drho1d(dx,dy,dz);

    u = ekx = eky = ekz = ZEROF;
    for (n = nlower; n <= nupper; n++) {
      mz = n+nz;
      for (m = nlower; m <= nupper; m++) {
        my = m+ny;
        for (l = nlower; l <= nupper; l++) {
          mx = l+nx;
          u += rho1d[0][l]*rho1d[1][m]*rho1d[2][n]*u_brick[mz][my][mx];
          ekx += drho1d[0][l]*rho1d[1][m]*rho1d[2][n]*u_brick[mz][my][mx];
          eky += rho1d[0][l]*drho1d[1][m]*rho1d[2][n]*u_brick[mz][my][mx];
          ekz += rho1d[0][l]*rho1d[1][m]*drho1d[2][n]*u_brick[mz][my][mx];
        }
      }
    }
    ekx *= hx_inv;
    eky *= hy_inv;
    ekz *= hz_inv;

    // electrical potential

    if (potflag) phi[i] = u;

    // convert E-field to force and substract self forces

    const double qfactor = qqrd2e * scale;
    double qtmp = eps[i]*q[i];

    s1 = x[i][0]*hx_inv;
    s2 = x[i][1]*hy_inv;
    s3 = x[i][2]*hz_inv;
    sf = sf_coeff[0]*sin(2*MY_PI*s1);
    sf += sf_coeff[1]*sin(4*MY_PI*s1);
    sf *= 2*qtmp*qtmp;
    f[i][0] += qfactor*(ekx*qtmp - sf);
    if (qtmp != 0) efield[i][0] = qfactor*(ekx - sf/qtmp);
    else efield[i][0] = qfactor*ekx;

    sf = sf_coeff[2]*sin(2*MY_PI*s2);
    sf += sf_coeff[3]*sin(4*MY_PI*s2);
    sf *= 2*qtmp*qtmp;
    f[i][1] += qfactor*(eky*qtmp - sf);
    if (qtmp != 0) efield[i][1] = qfactor*(eky - sf/qtmp);
    else efield[i][1] = qfactor*eky;

    sf = sf_coeff[4]*sin(2*MY_PI*s3);
    sf += sf_coeff[5]*sin(4*MY_PI*s3);
    sf *= 2*qtmp*qtmp;
    if (slabflag != 2) {
      f[i][2] += qfactor*(ekz*qtmp - sf);
      if (qtmp != 0) efield[i][2] = qfactor*(ekz - sf/qtmp);
      else efield[i][2] = qfactor*ekz;
    }

  }
}

/* ----------------------------------------------------------------------
   pack own values to buf to send to another proc
------------------------------------------------------------------------- */

void PPPMDielectric::pack_forward(int flag, FFT_SCALAR *buf, int nlist, int *list)
{
  int n = 0;

  if (flag == FORWARD_IK) {
    FFT_SCALAR *xsrc = &vdx_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *ysrc = &vdy_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *zsrc = &vdz_brick[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++) {
      buf[n++] = xsrc[list[i]];
      buf[n++] = ysrc[list[i]];
      buf[n++] = zsrc[list[i]];
    }
  } else if (flag == FORWARD_MU) {
    FFT_SCALAR *src_ux = &ux_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *src_uy = &uy_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *src_uz = &uz_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *src_vxx = &vdxx_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *src_vyy = &vdyy_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *src_vzz = &vdzz_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *src_vxy = &vdxy_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *src_vxz = &vdxz_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *src_vyz = &vdyz_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++) {
      buf[n++] = src_ux[list[i]];
      buf[n++] = src_uy[list[i]];
      buf[n++] = src_uz[list[i]];
      buf[n++] = src_vxx[list[i]];
      buf[n++] = src_vyy[list[i]];
      buf[n++] = src_vzz[list[i]];
      buf[n++] = src_vxy[list[i]];
      buf[n++] = src_vxz[list[i]];
      buf[n++] = src_vyz[list[i]];
    }
  } else if (flag == FORWARD_AD) {
    FFT_SCALAR *src = &u_brick[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++)
      buf[i] = src[list[i]];
  } else if (flag == FORWARD_IK_PERATOM) {
    FFT_SCALAR *esrc = &u_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v0src = &v0_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v1src = &v1_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v2src = &v2_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v3src = &v3_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v4src = &v4_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v5src = &v5_brick[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++) {
      if (eflag_atom) buf[n++] = esrc[list[i]];
      if (vflag_atom) {
        buf[n++] = v0src[list[i]];
        buf[n++] = v1src[list[i]];
        buf[n++] = v2src[list[i]];
        buf[n++] = v3src[list[i]];
        buf[n++] = v4src[list[i]];
        buf[n++] = v5src[list[i]];
      }
    }
  } else if (flag == FORWARD_AD_PERATOM) {
    FFT_SCALAR *v0src = &v0_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v1src = &v1_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v2src = &v2_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v3src = &v3_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v4src = &v4_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v5src = &v5_brick[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++) {
      buf[n++] = v0src[list[i]];
      buf[n++] = v1src[list[i]];
      buf[n++] = v2src[list[i]];
      buf[n++] = v3src[list[i]];
      buf[n++] = v4src[list[i]];
      buf[n++] = v5src[list[i]];
    }
  }
}

/* ----------------------------------------------------------------------
   unpack another proc's own values from buf and set own ghost values
------------------------------------------------------------------------- */

void PPPMDielectric::unpack_forward(int flag, FFT_SCALAR *buf, int nlist, int *list)
{
  int n = 0;

  if (flag == FORWARD_IK) {
    FFT_SCALAR *xdest = &vdx_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *ydest = &vdy_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *zdest = &vdz_brick[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++) {
      xdest[list[i]] = buf[n++];
      ydest[list[i]] = buf[n++];
      zdest[list[i]] = buf[n++];
    }
  } else if (flag == FORWARD_MU) {
    FFT_SCALAR *dest_ux = &ux_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *dest_uy = &uy_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *dest_uz = &uz_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *dest_vxx = &vdxx_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *dest_vyy = &vdyy_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *dest_vzz = &vdzz_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *dest_vxy = &vdxy_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *dest_vxz = &vdxz_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *dest_vyz = &vdyz_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++) {
      dest_ux[list[i]] = buf[n++];
      dest_uy[list[i]] = buf[n++];
      dest_uz[list[i]] = buf[n++];
      dest_vxx[list[i]] = buf[n++];
      dest_vyy[list[i]] = buf[n++];
      dest_vzz[list[i]] = buf[n++];
      dest_vxy[list[i]] = buf[n++];
      dest_vxz[list[i]] = buf[n++];
      dest_vyz[list[i]] = buf[n++];
    }
  } else if (flag == FORWARD_AD) {
    FFT_SCALAR *dest = &u_brick[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++)
      dest[list[i]] = buf[i];
  } else if (flag == FORWARD_IK_PERATOM) {
    FFT_SCALAR *esrc = &u_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v0src = &v0_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v1src = &v1_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v2src = &v2_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v3src = &v3_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v4src = &v4_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v5src = &v5_brick[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++) {
      if (eflag_atom) esrc[list[i]] = buf[n++];
      if (vflag_atom) {
        v0src[list[i]] = buf[n++];
        v1src[list[i]] = buf[n++];
        v2src[list[i]] = buf[n++];
        v3src[list[i]] = buf[n++];
        v4src[list[i]] = buf[n++];
        v5src[list[i]] = buf[n++];
      }
    }
  } else if (flag == FORWARD_AD_PERATOM) {
    FFT_SCALAR *v0src = &v0_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v1src = &v1_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v2src = &v2_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v3src = &v3_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v4src = &v4_brick[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *v5src = &v5_brick[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++) {
      v0src[list[i]] = buf[n++];
      v1src[list[i]] = buf[n++];
      v2src[list[i]] = buf[n++];
      v3src[list[i]] = buf[n++];
      v4src[list[i]] = buf[n++];
      v5src[list[i]] = buf[n++];
    }
  }
}

/* ----------------------------------------------------------------------
   pack ghost values into buf to send to another proc
------------------------------------------------------------------------- */

void PPPMDielectric::pack_reverse(int flag, FFT_SCALAR *buf, int nlist, int *list)
{
  int n = 0;
  if (flag == REVERSE_RHO) {
    FFT_SCALAR *src = &density_brick[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++)
      buf[n++] = src[list[i]];
  } else if (flag == REVERSE_MU) {
    FFT_SCALAR *src_mu0 = &densityx_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *src_mu1 = &densityy_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *src_mu2 = &densityz_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++) {
      buf[n++] = src_mu0[list[i]];
      buf[n++] = src_mu1[list[i]];
      buf[n++] = src_mu2[list[i]];
    }
  }
}

/* ----------------------------------------------------------------------
   unpack another proc's ghost values from buf and add to own values
------------------------------------------------------------------------- */

void PPPMDielectric::unpack_reverse(int flag, FFT_SCALAR *buf, int nlist, int *list)
{
  int n = 0;
  if (flag == REVERSE_RHO) {
    FFT_SCALAR *dest = &density_brick[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++)
      dest[list[i]] += buf[n++];
  } else if (flag == REVERSE_MU) {
    FFT_SCALAR *dest_mu0 = &densityx_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *dest_mu1 = &densityy_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    FFT_SCALAR *dest_mu2 = &densityz_brick_dipole[nzlo_out][nylo_out][nxlo_out];
    for (int i = 0; i < nlist; i++) {
      dest_mu0[list[i]] += buf[n++];
      dest_mu1[list[i]] += buf[n++];
      dest_mu2[list[i]] += buf[n++];
    }
  }
}

/* ----------------------------------------------------------------------
   Slab-geometry correction term to dampen inter-slab interactions between
   periodically repeating slabs.  Yields good approximation to 2D Ewald if
   adequate empty space is left between repeating slabs (J. Chem. Phys.
   111, 3155).  Slabs defined here to be parallel to the xy plane. Also
   extended to non-neutral systems (J. Chem. Phys. 131, 094107).
------------------------------------------------------------------------- */

void PPPMDielectric::slabcorr()
{
  // compute local contribution to global dipole moment

  double *q = atom->q;
  double **x = atom->x;
  double *eps = avec->epsilon;
  double zprd = domain->zprd;
  int nlocal = atom->nlocal;

  double dipole = 0.0;
  for (int i = 0; i < nlocal; i++) dipole += q[i]*x[i][2];

  if (mu_flag) {
    double **mu = atom->mu;
    for (int i = 0; i < nlocal; i++) dipole += mu[i][2];
  }

  // sum local contributions to get global dipole moment

  double dipole_all;
  MPI_Allreduce(&dipole,&dipole_all,1,MPI_DOUBLE,MPI_SUM,world);

  // need to make non-neutral systems and/or
  //  per-atom energy translationally invariant

  double dipole_r2 = 0.0;
  if (eflag_atom || fabs(qsum) > SMALL) {

    if (mu_flag)
      error->all(FLERR,"Cannot (yet) use kspace slab correction with "
        "long-range dipoles and non-neutral systems or per-atom energy");

    for (int i = 0; i < nlocal; i++)
      dipole_r2 += q[i]*x[i][2]*x[i][2];

    // sum local contributions

    double tmp;
    MPI_Allreduce(&dipole_r2,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
    dipole_r2 = tmp;
  }

  // compute corrections

  const double e_slabcorr = MY_2PI*(dipole_all*dipole_all -
    qsum*dipole_r2 - qsum*qsum*zprd*zprd/12.0)/volume;
  const double qscale = qqrd2e * scale;

  if (eflag_global) energy += qscale * e_slabcorr;

  // per-atom energy

  if (eflag_atom) {
    double efact = qscale * MY_2PI/volume;
    for (int i = 0; i < nlocal; i++)
      eatom[i] += efact * eps[i]*q[i]*(x[i][2]*dipole_all - 0.5*(dipole_r2 +
        qsum*x[i][2]*x[i][2]) - qsum*zprd*zprd/12.0);
  }

  // add on force corrections

  double ffact = qscale * (-4.0*MY_PI/volume);
  double **f = atom->f;

  for (int i = 0; i < nlocal; i++) {
    f[i][2] += ffact * eps[i]*q[i]*(dipole_all - qsum*x[i][2]);
    efield[i][2] += ffact * eps[i]*(dipole_all - qsum*x[i][2]);
  }

  // add on torque corrections

  if (mu_flag && atom->torque) {
    double **mu = atom->mu;
    double **torque = atom->torque;
    for (int i = 0; i < nlocal; i++) {
      torque[i][0] += ffact * dipole_all * mu[i][1];
      torque[i][1] += -ffact * dipole_all * mu[i][0];
    }
  }
}

/* ----------------------------------------------------------------------
   compute qsum,qsqsum,q2 and give error/warning if not charge neutral
   called initially, when particle count changes, when charges are changed
------------------------------------------------------------------------- */

void PPPMDielectric::qsum_qsq()
{
  const double * const q = atom->q;
  const double * const eps = avec->epsilon;
  const int nlocal = atom->nlocal;
  double qsum_local(0.0), qsqsum_local(0.0);

#if defined(_OPENMP)
#pragma omp parallel for default(none) reduction(+:qsum_local,qsqsum_local)
#endif
  for (int i = 0; i < nlocal; i++) {
    double qtmp = eps[i]*q[i];
    qsum_local += qtmp;
    qsqsum_local += qtmp*qtmp;
  }

  MPI_Allreduce(&qsum_local,&qsum,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&qsqsum_local,&qsqsum,1,MPI_DOUBLE,MPI_SUM,world);

  q2 = qsqsum * force->qqrd2e;
}

/* ----------------------------------------------------------------------
   memory usage of local arrays
------------------------------------------------------------------------- */

double PPPMDielectric::memory_usage()
{
  double bytes = nmax*3 * sizeof(double);
  int nbrick = (nxhi_out-nxlo_out+1) * (nyhi_out-nylo_out+1) *
    (nzhi_out-nzlo_out+1);
  if (differentiation_flag == 1) {
    bytes += 2 * nbrick * sizeof(FFT_SCALAR);
  } else {
    bytes += 4 * nbrick * sizeof(FFT_SCALAR);
  }
  if (triclinic) bytes += 3 * nfft_both * sizeof(double);
  bytes += 6 * nfft_both * sizeof(double);
  bytes += nfft_both * sizeof(double);
  bytes += nfft_both*5 * sizeof(FFT_SCALAR);

  if (mu_flag) {
    bytes += 3 * nbrick * sizeof(FFT_SCALAR);
    //work?
  }

  if (peratom_allocate_flag)
    bytes += 6 * nbrick * sizeof(FFT_SCALAR);

  if (group_allocate_flag) {
    bytes += 2 * nbrick * sizeof(FFT_SCALAR);
    bytes += 2 * nfft_both * sizeof(FFT_SCALAR);;
  }

  if (cg) bytes += cg->memory_usage();

  if (mu_flag)
    bytes += cg_mu->memory_usage();

  return bytes;
}


/* ----------------------------------------------------------------------
   compute qsum,qsqsum,q2 and give error/warning if not charge neutral
   called initially, when particle count changes, when charges are changed
------------------------------------------------------------------------- */

void PPPMDielectric::musum_musq()
{
  double** mu = atom->mu;
  const int nlocal = atom->nlocal;
  double musum_local(0.0), musqsum_local(0.0);

  for (int i = 0; i < nlocal; i++) {
    musum_local += mu[i][0] + mu[i][1] + mu[i][2];
    musqsum_local += mu[i][0]*mu[i][0] + mu[i][1]*mu[i][1] + mu[i][2]*mu[i][2];
  }

  MPI_Allreduce(&musum_local,&musum,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&musqsum_local,&musqsum,1,MPI_DOUBLE,MPI_SUM,world);

  /*
  if ((qsqsum == 0.0) && (comm->me == 0) && warn_nocharge) {
    error->warning(FLERR,"Using kspace solver on system with no charge");
    warn_nocharge = 0;
  }*/

  mu2 = musqsum * force->qqrd2e;

  // not yet sure of the correction needed for non-neutral systems
  // so issue warning or error
  /*
  if (fabs(qsum) > SMALL) {
    char str[128];
    sprintf(str,"System is not charge neutral, net charge = %g",qsum);
    if (!warn_nonneutral) error->all(FLERR,str);
    if (warn_nonneutral == 1 && comm->me == 0) error->warning(FLERR,str);
    warn_nonneutral = 2;
  }*/
}

