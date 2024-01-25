// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Trung Nguyen (Northwestern)
------------------------------------------------------------------------- */

#include "pppm_disp_dielectric.h"

#include "atom.h"
#include "atom_vec_dielectric.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "grid3d.h"
#include "math_const.h"
#include "memory.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace MathConst;

static constexpr int MAXORDER =   7;
static constexpr int OFFSET = 16384;
static constexpr double SMALL = 0.00001;
static constexpr double LARGE = 10000.0;
static constexpr double EPS_HOC = 1.0e-7;

enum{REVERSE_RHO,REVERSE_RHO_GEOM,REVERSE_RHO_ARITH,REVERSE_RHO_NONE};
enum{FORWARD_IK,FORWARD_AD,FORWARD_IK_PERATOM,FORWARD_AD_PERATOM,
     FORWARD_IK_GEOM,FORWARD_AD_GEOM,
     FORWARD_IK_PERATOM_GEOM,FORWARD_AD_PERATOM_GEOM,
     FORWARD_IK_ARITH,FORWARD_AD_ARITH,
     FORWARD_IK_PERATOM_ARITH,FORWARD_AD_PERATOM_ARITH,
     FORWARD_IK_NONE,FORWARD_AD_NONE,FORWARD_IK_PERATOM_NONE,
     FORWARD_AD_PERATOM_NONE};

static constexpr FFT_SCALAR ZEROF = 0.0;
static constexpr FFT_SCALAR ONEF =  1.0;

/* ---------------------------------------------------------------------- */

PPPMDispDielectric::PPPMDispDielectric(LAMMPS *_lmp) : PPPMDisp(_lmp)
{
  dipoleflag = 0; // turned off for now, until dipole works
  group_group_enable = 0;

  mu_flag = 0;
  use_qscaled = true;

  // no warnings about non-neutral systems from qsum_qsq()
  warn_nonneutral = 2;

  efield = nullptr;
  phi = nullptr;
  potflag = 0;

  avec = dynamic_cast<AtomVecDielectric *>(atom->style_match("dielectric"));
  if (!avec) error->all(FLERR,"pppm/dielectric requires atom style dielectric");
}

/* ---------------------------------------------------------------------- */

PPPMDispDielectric::~PPPMDispDielectric()
{
  memory->destroy(efield);
  memory->destroy(phi);
}

/* ----------------------------------------------------------------------
   compute the PPPM long-range force, energy, virial
------------------------------------------------------------------------- */

void PPPMDispDielectric::compute(int eflag, int vflag)
{

  int i;

  // set energy/virial flags
  // invoke allocate_peratom() if needed for first time

  ev_init(eflag,vflag);

  if (evflag_atom && !peratom_allocate_flag) allocate_peratom();

  // convert atoms from box to lamda coords

  if (triclinic == 0) boxlo = domain->boxlo;
  else {
    boxlo = domain->boxlo_lamda;
    domain->x2lamda(atom->nlocal);
  }
  // extend size of per-atom arrays if necessary

  if (atom->nmax > nmax) {

    if (function[0]) {
      memory->destroy(part2grid);
      memory->destroy(efield);
      memory->destroy(phi);
    }
    if (function[1] + function[2] + function[3]) memory->destroy(part2grid_6);
    nmax = atom->nmax;
    if (function[0]) {
      memory->create(part2grid,nmax,3,"pppm/disp:part2grid");
      memory->create(efield,nmax,3,"pppm/disp:efield");
      memory->create(phi,nmax,"pppm/disp:phi");
    }
    if (function[1] + function[2] + function[3])
      memory->create(part2grid_6,nmax,3,"pppm/disp:part2grid_6");
  }

  energy = 0.0;
  energy_1 = 0.0;
  energy_6 = 0.0;
  if (vflag) for (i = 0; i < 6; i++) virial_6[i] = virial_1[i] = 0.0;

  // find grid points for all my particles
  // distribute partcles' charges/dispersion coefficients on the grid
  // communication between processors and remapping two fft
  // Solution of poissons equation in k-space and backtransformation
  // communication between processors
  // calculation of forces

  if (function[0]) {

    // perform calculations for coulomb interactions only

    particle_map_c(delxinv,delyinv,delzinv,shift,part2grid,nupper,nlower,
                   nxlo_out,nylo_out,nzlo_out,nxhi_out,nyhi_out,nzhi_out);

    make_rho_c();

    gc->reverse_comm(Grid3d::KSPACE,this,REVERSE_RHO,1,sizeof(FFT_SCALAR),
                     gc_buf1,gc_buf2,MPI_FFT_SCALAR);

    brick2fft(nxlo_in,nylo_in,nzlo_in,nxhi_in,nyhi_in,nzhi_in,
              density_brick,density_fft,work1,remap);

    if (differentiation_flag == 1) {
      poisson_ad(work1,work2,density_fft,fft1,fft2,
                 nx_pppm,ny_pppm,nz_pppm,nfft,
                 nxlo_fft,nylo_fft,nzlo_fft,nxhi_fft,nyhi_fft,nzhi_fft,
                 nxlo_in,nylo_in,nzlo_in,nxhi_in,nyhi_in,nzhi_in,
                 energy_1,greensfn,
                 virial_1,vg,vg2,
                 u_brick,v0_brick,v1_brick,v2_brick,v3_brick,v4_brick,v5_brick);

      gc->forward_comm(Grid3d::KSPACE,this,FORWARD_AD,1,sizeof(FFT_SCALAR),
                       gc_buf1,gc_buf2,MPI_FFT_SCALAR);

      fieldforce_c_ad();

      if (vflag_atom)
        gc->forward_comm(Grid3d::KSPACE,this,FORWARD_AD_PERATOM,6,sizeof(FFT_SCALAR),
                         gc_buf1,gc_buf2,MPI_FFT_SCALAR);

    } else {
      poisson_ik(work1,work2,density_fft,fft1,fft2,
                 nx_pppm,ny_pppm,nz_pppm,nfft,
                 nxlo_fft,nylo_fft,nzlo_fft,nxhi_fft,nyhi_fft,nzhi_fft,
                 nxlo_in,nylo_in,nzlo_in,nxhi_in,nyhi_in,nzhi_in,
                 energy_1,greensfn,
                 fkx,fky,fkz,fkx2,fky2,fkz2,
                 vdx_brick,vdy_brick,vdz_brick,virial_1,vg,vg2,
                 u_brick,v0_brick,v1_brick,v2_brick,v3_brick,v4_brick,v5_brick);

      gc->forward_comm(Grid3d::KSPACE,this,FORWARD_IK,3,sizeof(FFT_SCALAR),
                       gc_buf1,gc_buf2,MPI_FFT_SCALAR);

      fieldforce_c_ik();

      if (evflag_atom)
        gc->forward_comm(Grid3d::KSPACE,this,FORWARD_IK_PERATOM,7,sizeof(FFT_SCALAR),
                         gc_buf1,gc_buf2,MPI_FFT_SCALAR);
    }

    if (evflag_atom) fieldforce_c_peratom();
  }

  if (function[1]) {

    // perform calculations for geometric mixing

    particle_map(delxinv_6,delyinv_6,delzinv_6,shift_6,part2grid_6,
                 nupper_6,nlower_6,
                 nxlo_out_6,nylo_out_6,nzlo_out_6,
                 nxhi_out_6,nyhi_out_6,nzhi_out_6);

    make_rho_g();

    gc6->reverse_comm(Grid3d::KSPACE,this,REVERSE_RHO_GEOM,1,sizeof(FFT_SCALAR),
                      gc6_buf1,gc6_buf2,MPI_FFT_SCALAR);

    brick2fft(nxlo_in_6,nylo_in_6,nzlo_in_6,nxhi_in_6,nyhi_in_6,nzhi_in_6,
              density_brick_g,density_fft_g,work1_6,remap_6);

    if (differentiation_flag == 1) {
      poisson_ad(work1_6,work2_6,density_fft_g,fft1_6,fft2_6,
                 nx_pppm_6,ny_pppm_6,nz_pppm_6,nfft_6,
                 nxlo_fft_6,nylo_fft_6,nzlo_fft_6,nxhi_fft_6,nyhi_fft_6,nzhi_fft_6,
                 nxlo_in_6,nylo_in_6,nzlo_in_6,nxhi_in_6,nyhi_in_6,nzhi_in_6,
                 energy_6,greensfn_6,
                 virial_6,vg_6,vg2_6,
                 u_brick_g,v0_brick_g,v1_brick_g,v2_brick_g,
                 v3_brick_g,v4_brick_g,v5_brick_g);

      gc6->forward_comm(Grid3d::KSPACE,this,FORWARD_AD_GEOM,1,sizeof(FFT_SCALAR),
                        gc6_buf1,gc6_buf2,MPI_FFT_SCALAR);

      fieldforce_g_ad();

      if (vflag_atom)
        gc6->forward_comm(Grid3d::KSPACE,this,FORWARD_AD_PERATOM_GEOM,6,sizeof(FFT_SCALAR),
                          gc6_buf1,gc6_buf2,MPI_FFT_SCALAR);

    } else {
      poisson_ik(work1_6,work2_6,density_fft_g,fft1_6,fft2_6,
                 nx_pppm_6,ny_pppm_6,nz_pppm_6,nfft_6,
                 nxlo_fft_6,nylo_fft_6,nzlo_fft_6,nxhi_fft_6,nyhi_fft_6,nzhi_fft_6,
                 nxlo_in_6,nylo_in_6,nzlo_in_6,nxhi_in_6,nyhi_in_6,nzhi_in_6,
                 energy_6,greensfn_6,
                 fkx_6,fky_6,fkz_6,fkx2_6,fky2_6,fkz2_6,
                 vdx_brick_g,vdy_brick_g,vdz_brick_g,virial_6,vg_6,vg2_6,
                 u_brick_g,v0_brick_g,v1_brick_g,v2_brick_g,
                 v3_brick_g,v4_brick_g,v5_brick_g);

      gc6->forward_comm(Grid3d::KSPACE,this,FORWARD_IK_GEOM,3,sizeof(FFT_SCALAR),
                        gc6_buf1,gc6_buf2,MPI_FFT_SCALAR);

      fieldforce_g_ik();

      if (evflag_atom)
        gc6->forward_comm(Grid3d::KSPACE,this,FORWARD_IK_PERATOM_GEOM,7,sizeof(FFT_SCALAR),
                          gc6_buf1,gc6_buf2,MPI_FFT_SCALAR);
    }

    if (evflag_atom) fieldforce_g_peratom();
  }

  if (function[2]) {

    // perform calculations for arithmetic mixing

    particle_map(delxinv_6,delyinv_6,delzinv_6,shift_6,part2grid_6,
                 nupper_6,nlower_6,
                 nxlo_out_6,nylo_out_6,nzlo_out_6,
                 nxhi_out_6,nyhi_out_6,nzhi_out_6);

    make_rho_a();

    gc6->reverse_comm(Grid3d::KSPACE,this,REVERSE_RHO_ARITH,7,sizeof(FFT_SCALAR),
                      gc6_buf1,gc6_buf2,MPI_FFT_SCALAR);

    brick2fft_a();

    if (differentiation_flag == 1) {
      poisson_ad(work1_6,work2_6,density_fft_a3,fft1_6,fft2_6,
                 nx_pppm_6,ny_pppm_6,nz_pppm_6,nfft_6,
                 nxlo_fft_6,nylo_fft_6,nzlo_fft_6,nxhi_fft_6,nyhi_fft_6,nzhi_fft_6,
                 nxlo_in_6,nylo_in_6,nzlo_in_6,nxhi_in_6,nyhi_in_6,nzhi_in_6,
                 energy_6,greensfn_6,
                 virial_6,vg_6,vg2_6,
                 u_brick_a3,v0_brick_a3,v1_brick_a3,v2_brick_a3,
                 v3_brick_a3,v4_brick_a3,v5_brick_a3);
      poisson_2s_ad(density_fft_a0,density_fft_a6,
                    u_brick_a0,v0_brick_a0,v1_brick_a0,v2_brick_a0,
                    v3_brick_a0,v4_brick_a0,v5_brick_a0,
                    u_brick_a6,v0_brick_a6,v1_brick_a6,v2_brick_a6,
                    v3_brick_a6,v4_brick_a6,v5_brick_a6);
      poisson_2s_ad(density_fft_a1,density_fft_a5,
                    u_brick_a1,v0_brick_a1,v1_brick_a1,v2_brick_a1,
                    v3_brick_a1,v4_brick_a1,v5_brick_a1,
                    u_brick_a5,v0_brick_a5,v1_brick_a5,v2_brick_a5,
                    v3_brick_a5,v4_brick_a5,v5_brick_a5);
      poisson_2s_ad(density_fft_a2,density_fft_a4,
                    u_brick_a2,v0_brick_a2,v1_brick_a2,v2_brick_a2,
                    v3_brick_a2,v4_brick_a2,v5_brick_a2,
                    u_brick_a4,v0_brick_a4,v1_brick_a4,v2_brick_a4,
                    v3_brick_a4,v4_brick_a4,v5_brick_a4);

      gc6->forward_comm(Grid3d::KSPACE,this,FORWARD_AD_ARITH,7,sizeof(FFT_SCALAR),
                        gc6_buf1,gc6_buf2,MPI_FFT_SCALAR);

      fieldforce_a_ad();

      if (evflag_atom)
        gc6->forward_comm(Grid3d::KSPACE,this,FORWARD_AD_PERATOM_ARITH,42,sizeof(FFT_SCALAR),
                          gc6_buf1,gc6_buf2,MPI_FFT_SCALAR);

    }  else {
      poisson_ik(work1_6,work2_6,density_fft_a3,fft1_6,fft2_6,
                 nx_pppm_6,ny_pppm_6,nz_pppm_6,nfft_6,
                 nxlo_fft_6,nylo_fft_6,nzlo_fft_6,nxhi_fft_6,nyhi_fft_6,nzhi_fft_6,
                 nxlo_in_6,nylo_in_6,nzlo_in_6,nxhi_in_6,nyhi_in_6,nzhi_in_6,
                 energy_6,greensfn_6,
                 fkx_6,fky_6,fkz_6,fkx2_6,fky2_6,fkz2_6,
                 vdx_brick_a3,vdy_brick_a3,vdz_brick_a3,virial_6,vg_6,vg2_6,
                 u_brick_a3,v0_brick_a3,v1_brick_a3,v2_brick_a3,
                 v3_brick_a3,v4_brick_a3,v5_brick_a3);
      poisson_2s_ik(density_fft_a0,density_fft_a6,
                    vdx_brick_a0,vdy_brick_a0,vdz_brick_a0,
                    vdx_brick_a6,vdy_brick_a6,vdz_brick_a6,
                    u_brick_a0,v0_brick_a0,v1_brick_a0,v2_brick_a0,
                    v3_brick_a0,v4_brick_a0,v5_brick_a0,
                    u_brick_a6,v0_brick_a6,v1_brick_a6,v2_brick_a6,
                    v3_brick_a6,v4_brick_a6,v5_brick_a6);
      poisson_2s_ik(density_fft_a1,density_fft_a5,
                    vdx_brick_a1,vdy_brick_a1,vdz_brick_a1,
                    vdx_brick_a5,vdy_brick_a5,vdz_brick_a5,
                    u_brick_a1,v0_brick_a1,v1_brick_a1,v2_brick_a1,
                    v3_brick_a1,v4_brick_a1,v5_brick_a1,
                    u_brick_a5,v0_brick_a5,v1_brick_a5,v2_brick_a5,
                    v3_brick_a5,v4_brick_a5,v5_brick_a5);
      poisson_2s_ik(density_fft_a2,density_fft_a4,
                    vdx_brick_a2,vdy_brick_a2,vdz_brick_a2,
                    vdx_brick_a4,vdy_brick_a4,vdz_brick_a4,
                    u_brick_a2,v0_brick_a2,v1_brick_a2,v2_brick_a2,
                    v3_brick_a2,v4_brick_a2,v5_brick_a2,
                    u_brick_a4,v0_brick_a4,v1_brick_a4,v2_brick_a4,
                    v3_brick_a4,v4_brick_a4,v5_brick_a4);

      gc6->forward_comm(Grid3d::KSPACE,this,FORWARD_IK_ARITH,21,sizeof(FFT_SCALAR),
                        gc6_buf1,gc6_buf2,MPI_FFT_SCALAR);

      fieldforce_a_ik();

      if (evflag_atom)
        gc6->forward_comm(Grid3d::KSPACE,this,FORWARD_IK_PERATOM_ARITH,49,sizeof(FFT_SCALAR),
                          gc6_buf1,gc6_buf2,MPI_FFT_SCALAR);
    }

    if (evflag_atom) fieldforce_a_peratom();
  }

  if (function[3]) {

    // perform calculations if no mixing rule applies

    particle_map(delxinv_6,delyinv_6,delzinv_6,shift_6,part2grid_6,
                 nupper_6,nlower_6,
                 nxlo_out_6,nylo_out_6,nzlo_out_6,
                 nxhi_out_6,nyhi_out_6,nzhi_out_6);

    make_rho_none();

    gc6->reverse_comm(Grid3d::KSPACE,this,REVERSE_RHO_NONE,nsplit_alloc,sizeof(FFT_SCALAR),
                      gc6_buf1,gc6_buf2,MPI_FFT_SCALAR);

    brick2fft_none();

    if (differentiation_flag == 1) {
      int n = 0;
      for (int k = 0; k < nsplit_alloc/2; k++) {
        poisson_none_ad(n,n+1,density_fft_none[n],density_fft_none[n+1],
                        u_brick_none[n],u_brick_none[n+1],
                        v0_brick_none,v1_brick_none,v2_brick_none,
                        v3_brick_none,v4_brick_none,v5_brick_none);
        n += 2;
      }

      gc6->forward_comm(Grid3d::KSPACE,this,FORWARD_AD_NONE,1*nsplit_alloc,sizeof(FFT_SCALAR),
                        gc6_buf1,gc6_buf2,MPI_FFT_SCALAR);

      fieldforce_none_ad();

      if (vflag_atom)
        gc6->forward_comm(Grid3d::KSPACE,this,FORWARD_AD_PERATOM_NONE,6*nsplit_alloc,sizeof(FFT_SCALAR),
                          gc6_buf1,gc6_buf2,MPI_FFT_SCALAR);

    } else {
      int n = 0;
      for (int k = 0; k < nsplit_alloc/2; k++) {
        poisson_none_ik(n,n+1,density_fft_none[n],density_fft_none[n+1],
                        vdx_brick_none[n],vdy_brick_none[n],vdz_brick_none[n],
                        vdx_brick_none[n+1],vdy_brick_none[n+1],vdz_brick_none[n+1],
                        u_brick_none,v0_brick_none,v1_brick_none,v2_brick_none,
                        v3_brick_none,v4_brick_none,v5_brick_none);
        n += 2;
      }

      gc6->forward_comm(Grid3d::KSPACE,this,FORWARD_IK_NONE,3*nsplit_alloc,sizeof(FFT_SCALAR),
                        gc6_buf1,gc6_buf2,MPI_FFT_SCALAR);

      fieldforce_none_ik();

      if (evflag_atom)
        gc6->forward_comm(Grid3d::KSPACE,this,FORWARD_IK_PERATOM_NONE,7*nsplit_alloc,sizeof(FFT_SCALAR),
                          gc6_buf1,gc6_buf2,MPI_FFT_SCALAR);
    }

    if (evflag_atom) fieldforce_none_peratom();
  }

  // update qsum and qsqsum, if atom count has changed and energy needed

  if ((eflag_global || eflag_atom) && atom->natoms != natoms_original) {
    qsum_qsq();
    natoms_original = atom->natoms;
  }

  // recompute the average epsilon of all the atoms

  compute_ave_epsilon();

  // sum energy across procs and add in volume-dependent term

  const double qscale = force->qqrd2e * scale;

  if (eflag_global) {

    if (function[0]) {

      // switch to unscaled charges to find charge density

      use_qscaled = false;

      make_rho_c();

      gc->reverse_comm(Grid3d::KSPACE,this,REVERSE_RHO,1,sizeof(FFT_SCALAR),
                     gc_buf1,gc_buf2,MPI_FFT_SCALAR);

      brick2fft(nxlo_in,nylo_in,nzlo_in,nxhi_in,nyhi_in,nzhi_in,
                density_brick,density_fft,work1,remap);

      // compute electrostatic energy with the unscaled charges and average epsilon

      if (differentiation_flag == 1) {
        poisson_ad(work1,work2,density_fft,fft1,fft2,
                  nx_pppm,ny_pppm,nz_pppm,nfft,
                  nxlo_fft,nylo_fft,nzlo_fft,nxhi_fft,nyhi_fft,nzhi_fft,
                  nxlo_in,nylo_in,nzlo_in,nxhi_in,nyhi_in,nzhi_in,
                  energy_1,greensfn,
                  virial_1,vg,vg2,
                  u_brick,v0_brick,v1_brick,v2_brick,v3_brick,v4_brick,v5_brick);

      } else {
        poisson_ik(work1,work2,density_fft,fft1,fft2,
                  nx_pppm,ny_pppm,nz_pppm,nfft,
                  nxlo_fft,nylo_fft,nzlo_fft,nxhi_fft,nyhi_fft,nzhi_fft,
                  nxlo_in,nylo_in,nzlo_in,nxhi_in,nyhi_in,nzhi_in,
                  energy_1,greensfn,
                  fkx,fky,fkz,fkx2,fky2,fkz2,
                  vdx_brick,vdy_brick,vdz_brick,virial_1,vg,vg2,
                  u_brick,v0_brick,v1_brick,v2_brick,v3_brick,v4_brick,v5_brick);

        gc->forward_comm(Grid3d::KSPACE,this,FORWARD_IK,3,sizeof(FFT_SCALAR),
                        gc_buf1,gc_buf2,MPI_FFT_SCALAR);
      }
    }

    double energy_all;
    MPI_Allreduce(&energy_1,&energy_all,1,MPI_DOUBLE,MPI_SUM,world);
    energy_1 = energy_all;
    MPI_Allreduce(&energy_6,&energy_all,1,MPI_DOUBLE,MPI_SUM,world);
    energy_6 = energy_all;

    energy_1 *= 0.5*volume;
    energy_6 *= 0.5*volume;

    energy_1 -= g_ewald*qsqsum/MY_PIS +
      MY_PI2*qsum*qsum / (g_ewald*g_ewald*volume);
    energy_6 += - MY_PI*MY_PIS/(6*volume)*pow(g_ewald_6,3)*csumij +
      1.0/12.0*pow(g_ewald_6,6)*csum;
    energy_1 *= qscale;

    // revert to qscaled charges (for force in the next time step)

    use_qscaled = true;
  }

  // sum virial across procs

  if (vflag_global) {
    double virial_all[6];
    MPI_Allreduce(virial_1,virial_all,6,MPI_DOUBLE,MPI_SUM,world);
    for (i = 0; i < 6; i++) virial[i] = 0.5*qscale*volume*virial_all[i];
    MPI_Allreduce(virial_6,virial_all,6,MPI_DOUBLE,MPI_SUM,world);
    for (i = 0; i < 6; i++) virial[i] += 0.5*volume*virial_all[i];
    if (function[1]+function[2]+function[3]) {
      double a =  MY_PI*MY_PIS/(6*volume)*pow(g_ewald_6,3)*csumij;
      virial[0] -= a;
      virial[1] -= a;
      virial[2] -= a;
    }
  }

  if (eflag_atom) {
    if (function[0]) {
      double *q = atom->q;
      // coulomb self energy correction
      for (i = 0; i < atom->nlocal; i++) {
        eatom[i] -= qscale*g_ewald*q[i]*q[i]/MY_PIS +
          qscale*MY_PI2*q[i]*qsum / (g_ewald*g_ewald*volume);
      }
    }
    if (function[1] + function[2] + function[3]) {
      int tmp;
      for (i = 0; i < atom->nlocal; i++) {
        tmp = atom->type[i];
        eatom[i] += - MY_PI*MY_PIS/(6*volume)*pow(g_ewald_6,3)*csumi[tmp] +
                      1.0/12.0*pow(g_ewald_6,6)*cii[tmp];
      }
    }
  }

  if (vflag_atom) {
    if (function[1] + function[2] + function[3]) {
      int tmp;
      // dispersion self virial correction
      for (i = 0; i < atom->nlocal; i++) {
        tmp = atom->type[i];
        for (int n = 0; n < 3; n++)
          vatom[i][n] -= MY_PI*MY_PIS/(6*volume)*pow(g_ewald_6,3)*csumi[tmp];
      }
    }
  }

  // 2d slab correction

  if (slabflag) slabcorr(eflag);
  if (function[0]) energy += energy_1;
  if (function[1] + function[2] + function[3]) energy += energy_6;

  // convert atoms back from lamda to box coords

  if (triclinic) domain->lamda2x(atom->nlocal);
}

/* ----------------------------------------------------------------------
   compute the average dielectric constant of all the atoms
   NOTE: for dielectric use cases
------------------------------------------------------------------------- */

void PPPMDispDielectric::compute_ave_epsilon()
{
  const double * const epsilon = atom->epsilon;
  const int nlocal = atom->nlocal;
  double epsilon_local(0.0);

#if defined(_OPENMP)
#pragma omp parallel for default(shared) reduction(+:epsilon_local)
#endif
  for (int i = 0; i < nlocal; i++) {
    epsilon_local += epsilon[i];
  }

  MPI_Allreduce(&epsilon_local,&epsilon_ave,1,MPI_DOUBLE,MPI_SUM,world);
  epsilon_ave /= (double)atom->natoms;
}

/* ----------------------------------------------------------------------
   compute qsum,qsqsum,q2 and give error/warning if not charge neutral
   called initially, when particle count changes, when charges are changed
------------------------------------------------------------------------- */

void PPPMDispDielectric::qsum_qsq(int warning_flag)
{
  const double * const q = atom->q;
  const double * const epsilon = atom->epsilon;
  const int nlocal = atom->nlocal;
  double qsum_local(0.0), qsqsum_local(0.0), qsqsume_local(0.0);
  double qsqsume;

#if defined(_OPENMP)
#pragma omp parallel for default(shared) reduction(+:qsum_local,qsqsum_local)
#endif
  for (int i = 0; i < nlocal; i++) {
    qsum_local += q[i];
    qsqsum_local += q[i]*q[i];
    qsqsume_local += q[i]*q[i]/epsilon[i];
  }

  MPI_Allreduce(&qsum_local,&qsum,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&qsqsum_local,&qsqsum,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&qsqsume_local,&qsqsume,1,MPI_DOUBLE,MPI_SUM,world);

  if ((qsqsum == 0.0) && (comm->me == 0) && warn_nocharge && warning_flag) {
    error->warning(FLERR,"Using kspace solver on system with no charge");
    warn_nocharge = 0;
  }

  // q2 is used to compute the mesh spacing, here using qsqsume to match with regular pppm
  q2 = qsqsume * force->qqrd2e; //q2 = qsqsum * force->qqrd2e;

  // not yet sure of the correction needed for non-neutral systems
  // so issue warning or error

  if (fabs(qsum) > SMALL) {
    std::string message = fmt::format("System is not charge neutral, net "
                                      "charge = {:.8}",qsum);
    if (!warn_nonneutral) error->all(FLERR,message);
    if (warn_nonneutral == 1 && comm->me == 0) error->warning(FLERR,message);
    warn_nonneutral = 2;
  }
}

/* ----------------------------------------------------------------------
   create discretized "density" on section of global grid due to my particles
   density(x,y,z) = charge "density" at grid points of my 3d brick
   (nxlo:nxhi,nylo:nyhi,nzlo:nzhi) is extent of my brick (including ghosts)
   in global grid
------------------------------------------------------------------------- */

void PPPMDispDielectric::make_rho_c()
{
  int l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz,x0,y0,z0;

  // clear 3d density array

  memset(&(density_brick[nzlo_out][nylo_out][nxlo_out]),0,
         ngrid*sizeof(FFT_SCALAR));

  // loop over my charges, add their contribution to nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt

  double *q = atom->q_scaled;
  if (!use_qscaled) q = atom->q;
  double **x = atom->x;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv;
    dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv;
    dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv;

    compute_rho1d(dx,dy,dz, order, rho_coeff, rho1d);

    z0 = delvolinv * q[i];
    for (n = nlower; n <= nupper; n++) {
      mz = n+nz;
      y0 = z0*rho1d[2][n];
      for (m = nlower; m <= nupper; m++) {
        my = m+ny;
        x0 = y0*rho1d[1][m];
        for (l = nlower; l <= nupper; l++) {
          mx = l+nx;
          density_brick[mz][my][mx] += x0*rho1d[0][l];
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get electric field & force on my particles
   for ik scheme
------------------------------------------------------------------------- */

void PPPMDispDielectric::fieldforce_c_ik()
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
  double *eps = atom->epsilon;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv;
    dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv;
    dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv;

    compute_rho1d(dx,dy,dz, order, rho_coeff, rho1d);

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

    // convert E-field to force

    const double qfactor = force->qqrd2e * scale * q[i];
    f[i][0] += qfactor*ekx;
    f[i][1] += qfactor*eky;
    if (slabflag != 2) f[i][2] += qfactor*ekz;
  }
}
/* ----------------------------------------------------------------------
   interpolate from grid to get electric field & force on my particles
   for ad scheme
------------------------------------------------------------------------- */

void PPPMDispDielectric::fieldforce_c_ad()
{
  int i,l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz;
  FFT_SCALAR ekx,eky,ekz,u;
  double s1,s2,s3;
  double sf = 0.0;

  double *prd;

  if (triclinic == 0) prd = domain->prd;
  else prd = domain->prd_lamda;

  double xprd = prd[0];
  double yprd = prd[1];
  double zprd = prd[2];
  double zprd_slab = zprd*slab_volfactor;

  double hx_inv = nx_pppm/xprd;
  double hy_inv = ny_pppm/yprd;
  double hz_inv = nz_pppm/zprd_slab;

  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of E-field on particle

  double *q = atom->q;
  double **x = atom->x;
  double **f = atom->f;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv;
    dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv;
    dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv;

    compute_rho1d(dx,dy,dz, order, rho_coeff, rho1d);
    compute_drho1d(dx,dy,dz, order, drho_coeff, drho1d);

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

    s1 = x[i][0]*hx_inv;
    s2 = x[i][1]*hy_inv;
    s3 = x[i][2]*hz_inv;
    sf = sf_coeff[0]*sin(2*MY_PI*s1);
    sf += sf_coeff[1]*sin(4*MY_PI*s1);
    sf *= 2*q[i]*q[i];
    f[i][0] += qfactor*(ekx*q[i] - sf);

    sf = sf_coeff[2]*sin(2*MY_PI*s2);
    sf += sf_coeff[3]*sin(4*MY_PI*s2);
    sf *= 2*q[i]*q[i];
    f[i][1] += qfactor*(eky*q[i] - sf);

    sf = sf_coeff[4]*sin(2*MY_PI*s3);
    sf += sf_coeff[5]*sin(4*MY_PI*s3);
    sf *= 2*q[i]*q[i];
    if (slabflag != 2) f[i][2] += qfactor*(ekz*q[i] - sf);
  }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get electric field & force on my particles
------------------------------------------------------------------------- */

void PPPMDispDielectric::fieldforce_c_peratom()
{
  int i,l,m,n,nx,ny,nz,mx,my,mz;
  FFT_SCALAR dx,dy,dz,x0,y0,z0;
  FFT_SCALAR u_pa,v0,v1,v2,v3,v4,v5;

  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of E-field on particle

  double *q = atom->q;
  double **x = atom->x;

  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++) {
    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];
    dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv;
    dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv;
    dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv;

    compute_rho1d(dx,dy,dz, order, rho_coeff, rho1d);

    u_pa = v0 = v1 = v2 = v3 = v4 = v5 = ZEROF;
    for (n = nlower; n <= nupper; n++) {
      mz = n+nz;
      z0 = rho1d[2][n];
      for (m = nlower; m <= nupper; m++) {
        my = m+ny;
        y0 = z0*rho1d[1][m];
        for (l = nlower; l <= nupper; l++) {
          mx = l+nx;
          x0 = y0*rho1d[0][l];
          if (eflag_atom) u_pa += x0*u_brick[mz][my][mx];
          if (vflag_atom) {
            v0 += x0*v0_brick[mz][my][mx];
            v1 += x0*v1_brick[mz][my][mx];
            v2 += x0*v2_brick[mz][my][mx];
            v3 += x0*v3_brick[mz][my][mx];
            v4 += x0*v4_brick[mz][my][mx];
            v5 += x0*v5_brick[mz][my][mx];
          }
        }
      }
    }

    // electrostatic potential

    phi[i] = u_pa;

    // convert E-field to force

    const double qfactor = 0.5*force->qqrd2e * scale * q[i];

    if (eflag_atom) eatom[i] += u_pa*qfactor;
    if (vflag_atom) {
      vatom[i][0] += v0*qfactor;
      vatom[i][1] += v1*qfactor;
      vatom[i][2] += v2*qfactor;
      vatom[i][3] += v3*qfactor;
      vatom[i][4] += v4*qfactor;
      vatom[i][5] += v5*qfactor;
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

void PPPMDispDielectric::slabcorr(int /*eflag*/)
{
  // compute local contribution to global dipole moment

  double *q = atom->q;
  double **x = atom->x;
  double *eps = atom->epsilon;
  double zprd_slab = domain->zprd*slab_volfactor;
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
    qsum*dipole_r2 - qsum*qsum*zprd_slab*zprd_slab/12.0)/volume;
  const double qscale = qqrd2e * scale;

  if (eflag_global) energy += qscale * e_slabcorr;

  // per-atom energy

  if (eflag_atom) {
    double efact = qscale * MY_2PI/volume;
    for (int i = 0; i < nlocal; i++)
      eatom[i] += efact * eps[i]*q[i]*(x[i][2]*dipole_all - 0.5*(dipole_r2 +
        qsum*x[i][2]*x[i][2]) - qsum*zprd_slab*zprd_slab/12.0);
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
   memory usage of local arrays
------------------------------------------------------------------------- */

double PPPMDispDielectric::memory_usage()
{
  double bytes = PPPMDisp::memory_usage();
  bytes += nmax*3 * sizeof(double);
  bytes += nmax * sizeof(double);
  return bytes;
}
