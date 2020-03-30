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
   Contributing authors: William McDoniel (RWTH Aachen University)
------------------------------------------------------------------------- */

#include <mpi.h>
#include <cstdlib>
#include <cmath>
#include "pppm_disp_intel.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "modify.h"
#include "fft3d_wrap.h"
#include "gridcomm.h"
#include "math_const.h"
#include "math_special.h"
#include "memory.h"
#include "suffix.h"

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathSpecial;

#define MAXORDER   7
#define OFFSET 16384
#define SMALL 0.00001
#define LARGE 10000.0
#define EPS_HOC 1.0e-7

enum{GEOMETRIC,ARITHMETIC,SIXTHPOWER};
enum{REVERSE_RHO, REVERSE_RHO_G, REVERSE_RHO_A, REVERSE_RHO_NONE};
enum{FORWARD_IK, FORWARD_AD, FORWARD_IK_PERATOM, FORWARD_AD_PERATOM,
     FORWARD_IK_G, FORWARD_AD_G, FORWARD_IK_PERATOM_G, FORWARD_AD_PERATOM_G,
     FORWARD_IK_A, FORWARD_AD_A, FORWARD_IK_PERATOM_A, FORWARD_AD_PERATOM_A,
     FORWARD_IK_NONE, FORWARD_AD_NONE, FORWARD_IK_PERATOM_NONE,
     FORWARD_AD_PERATOM_NONE};

#ifdef FFT_SINGLE
#define ZEROF 0.0f
#define ONEF  1.0f
#else
#define ZEROF 0.0
#define ONEF  1.0
#endif

/* ---------------------------------------------------------------------- */

PPPMDispIntel::PPPMDispIntel(LAMMPS *lmp) : PPPMDisp(lmp)
{
  suffix_flag |= Suffix::INTEL;
  triclinic_support = 0;

  order = 7;
  order_6 = 7; //sets default stencil sizes to 7

  perthread_density = NULL;
  particle_ekx = particle_eky = particle_ekz = NULL;
  particle_ekx0 = particle_eky0 = particle_ekz0 = NULL;
  particle_ekx1 = particle_eky1 = particle_ekz1 = NULL;
  particle_ekx2 = particle_eky2 = particle_ekz2 = NULL;
  particle_ekx3 = particle_eky3 = particle_ekz3 = NULL;
  particle_ekx4 = particle_eky4 = particle_ekz4 = NULL;
  particle_ekx5 = particle_eky5 = particle_ekz5 = NULL;
  particle_ekx6 = particle_eky6 = particle_ekz6 = NULL;

  rho_lookup = drho_lookup = NULL;
  rho6_lookup = drho6_lookup = NULL;
  rho_points = 0;

  _use_table = _use_packing = _use_lrt = 0;
}

PPPMDispIntel::~PPPMDispIntel()
{
  memory->destroy(perthread_density);
  memory->destroy(particle_ekx);
  memory->destroy(particle_eky);
  memory->destroy(particle_ekz);

  memory->destroy(rho_lookup);
  memory->destroy(drho_lookup);
  memory->destroy(rho6_lookup);
  memory->destroy(drho6_lookup);
}

/* ----------------------------------------------------------------------
   called once before run
------------------------------------------------------------------------- */

void PPPMDispIntel::init()
{

  PPPMDisp::init();
  int ifix = modify->find_fix("package_intel");
  if (ifix < 0)
    error->all(FLERR,
               "The 'package intel' command is required for /intel styles");
  fix = static_cast<FixIntel *>(modify->fix[ifix]);

  #ifdef _LMP_INTEL_OFFLOAD
  _use_base = 0;
  if (fix->offload_balance() != 0.0) {
    _use_base = 1;
    return;
  }
  #endif

  fix->kspace_init_check();

  _use_lrt = fix->lrt();
  if (_use_lrt)
    error->all(FLERR,
               "LRT mode is currently not supported for pppm/disp/intel");


  // For vectorization, we need some padding in the end
  // The first thread computes on the global density
  if ((comm->nthreads > 1) && !_use_lrt) {
    memory->destroy(perthread_density);
    memory->create(perthread_density, comm->nthreads-1,
                   ngrid + INTEL_P3M_ALIGNED_MAXORDER,
                   "pppmdispintel:perthread_density");
  }

  _use_table = fix->pppm_table();
  if (_use_table) {
    rho_points = 5000;
    memory->destroy(rho_lookup);
    memory->create(rho_lookup, rho_points, INTEL_P3M_ALIGNED_MAXORDER,
                   "pppmdispintel:rho_lookup");
    memory->destroy(rho6_lookup);
    memory->create(rho6_lookup, rho_points, INTEL_P3M_ALIGNED_MAXORDER,
                   "pppmdispintel:rho6_lookup");

    if(differentiation_flag == 1) {
      memory->destroy(drho_lookup);
      memory->create(drho_lookup, rho_points, INTEL_P3M_ALIGNED_MAXORDER,
                     "pppmdispintel:drho_lookup");
      memory->destroy(drho6_lookup);
      memory->create(drho6_lookup, rho_points, INTEL_P3M_ALIGNED_MAXORDER,
                     "pppmdispintel:drho6_lookup");
    }
    precompute_rho();
  }
  if (order > INTEL_P3M_MAXORDER)
    error->all(FLERR,"PPPM order greater than supported by USER-INTEL\n");
}

/* ----------------------------------------------------------------------
   compute the PPPMDispIntel long-range force, energy, virial
------------------------------------------------------------------------- */

void PPPMDispIntel::compute(int eflag, int vflag)
{
  #ifdef _LMP_INTEL_OFFLOAD
  if (_use_base) {
    PPPMDisp::compute(eflag, vflag);
    return;
  }
  #endif
  int i;
  // convert atoms from box to lamda coords

  ev_init(eflag,vflag);

  if (evflag_atom && !peratom_allocate_flag) {
    allocate_peratom();
    if (function[0]) {
      cg_peratom->ghost_notify();
      cg_peratom->setup();
    }
    if (function[1] + function[2] + function[3]) {
      cg_peratom_6->ghost_notify();
      cg_peratom_6->setup();
    }
    peratom_allocate_flag = 1;
  }
  if (triclinic == 0) boxlo = domain->boxlo;
  else {
    boxlo = domain->boxlo_lamda;
    domain->x2lamda(atom->nlocal);
  }
  // extend size of per-atom arrays if necessary

  if (atom->nmax > nmax) {

    if (function[0]) memory->destroy(part2grid);
    if (function[1] + function[2] + function[3]) memory->destroy(part2grid_6);
    if (differentiation_flag == 1) {
      memory->destroy(particle_ekx);
      memory->destroy(particle_eky);
      memory->destroy(particle_ekz);
      if (function[2] == 1){
        memory->destroy(particle_ekx0);
        memory->destroy(particle_eky0);
        memory->destroy(particle_ekz0);
        memory->destroy(particle_ekx1);
        memory->destroy(particle_eky1);
        memory->destroy(particle_ekz1);
        memory->destroy(particle_ekx2);
        memory->destroy(particle_eky2);
        memory->destroy(particle_ekz2);
        memory->destroy(particle_ekx3);
        memory->destroy(particle_eky3);
        memory->destroy(particle_ekz3);
        memory->destroy(particle_ekx4);
        memory->destroy(particle_eky4);
        memory->destroy(particle_ekz4);
        memory->destroy(particle_ekx5);
        memory->destroy(particle_eky5);
        memory->destroy(particle_ekz5);
        memory->destroy(particle_ekx6);
        memory->destroy(particle_eky6);
        memory->destroy(particle_ekz6);
      }

    }
    nmax = atom->nmax;
    if (function[0]) memory->create(part2grid,nmax,3,"pppm/disp:part2grid");
    if (function[1] + function[2] + function[3])
      memory->create(part2grid_6,nmax,3,"pppm/disp:part2grid_6");
    if (differentiation_flag == 1) {
      memory->create(particle_ekx, nmax, "pppmdispintel:pekx");
      memory->create(particle_eky, nmax, "pppmdispintel:peky");
      memory->create(particle_ekz, nmax, "pppmdispintel:pekz");
      if (function[2] == 1){
        memory->create(particle_ekx0, nmax, "pppmdispintel:pekx0");
        memory->create(particle_eky0, nmax, "pppmdispintel:peky0");
        memory->create(particle_ekz0, nmax, "pppmdispintel:pekz0");
        memory->create(particle_ekx1, nmax, "pppmdispintel:pekx1");
        memory->create(particle_eky1, nmax, "pppmdispintel:peky1");
        memory->create(particle_ekz1, nmax, "pppmdispintel:pekz1");
        memory->create(particle_ekx2, nmax, "pppmdispintel:pekx2");
        memory->create(particle_eky2, nmax, "pppmdispintel:peky2");
        memory->create(particle_ekz2, nmax, "pppmdispintel:pekz2");
        memory->create(particle_ekx3, nmax, "pppmdispintel:pekx3");
        memory->create(particle_eky3, nmax, "pppmdispintel:peky3");
        memory->create(particle_ekz3, nmax, "pppmdispintel:pekz3");
        memory->create(particle_ekx4, nmax, "pppmdispintel:pekx4");
        memory->create(particle_eky4, nmax, "pppmdispintel:peky4");
        memory->create(particle_ekz4, nmax, "pppmdispintel:pekz4");
        memory->create(particle_ekx5, nmax, "pppmdispintel:pekx5");
        memory->create(particle_eky5, nmax, "pppmdispintel:peky5");
        memory->create(particle_ekz5, nmax, "pppmdispintel:pekz5");
        memory->create(particle_ekx6, nmax, "pppmdispintel:pekx6");
        memory->create(particle_eky6, nmax, "pppmdispintel:peky6");
        memory->create(particle_ekz6, nmax, "pppmdispintel:pekz6");
      }
    }
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

    //perform calculations for coulomb interactions only

    if (fix->precision() == FixIntel::PREC_MODE_MIXED) {
      particle_map<float,double>(delxinv, delyinv, delzinv, shift, part2grid,
                                 nupper, nlower, nxlo_out, nylo_out, nzlo_out,
                                 nxhi_out, nyhi_out, nzhi_out,
                                 fix->get_mixed_buffers());
      make_rho_c<float,double>(fix->get_mixed_buffers());
    } else if (fix->precision() == FixIntel::PREC_MODE_DOUBLE) {
      particle_map<double,double>(delxinv, delyinv, delzinv, shift, part2grid,
                                  nupper, nlower, nxlo_out, nylo_out,
                                  nzlo_out, nxhi_out, nyhi_out, nzhi_out,
                                  fix->get_double_buffers());
      make_rho_c<double,double>(fix->get_double_buffers());
    } else {
      particle_map<float,float>(delxinv, delyinv, delzinv, shift, part2grid,
                                nupper, nlower, nxlo_out, nylo_out, nzlo_out,
                                nxhi_out, nyhi_out, nzhi_out,
                                fix->get_single_buffers());
      make_rho_c<float,float>(fix->get_single_buffers());
    }

    cg->reverse_comm(this,REVERSE_RHO);

    brick2fft(nxlo_in, nylo_in, nzlo_in, nxhi_in, nyhi_in, nzhi_in,
              density_brick, density_fft, work1,remap);

    if (differentiation_flag == 1) {
      poisson_ad(work1, work2, density_fft, fft1, fft2,
                 nx_pppm, ny_pppm, nz_pppm, nfft,
                 nxlo_fft, nylo_fft, nzlo_fft, nxhi_fft, nyhi_fft, nzhi_fft,
                 nxlo_in, nylo_in, nzlo_in, nxhi_in, nyhi_in, nzhi_in,
                 energy_1, greensfn, virial_1, vg,vg2, u_brick, v0_brick,
                 v1_brick, v2_brick, v3_brick, v4_brick, v5_brick);

      cg->forward_comm(this,FORWARD_AD);

      if (fix->precision() == FixIntel::PREC_MODE_MIXED) {
        fieldforce_c_ad<float,double>(fix->get_mixed_buffers());
      } else if (fix->precision() == FixIntel::PREC_MODE_DOUBLE) {
        fieldforce_c_ad<double,double>(fix->get_double_buffers());
      } else {
        fieldforce_c_ad<float,float>(fix->get_single_buffers());
      }

      if (vflag_atom) cg_peratom->forward_comm(this, FORWARD_AD_PERATOM);

    } else {
      poisson_ik(work1, work2, density_fft, fft1, fft2,
                 nx_pppm, ny_pppm, nz_pppm, nfft,
                 nxlo_fft, nylo_fft, nzlo_fft, nxhi_fft, nyhi_fft, nzhi_fft,
                 nxlo_in, nylo_in, nzlo_in, nxhi_in, nyhi_in, nzhi_in,
                 energy_1, greensfn, fkx, fky, fkz,fkx2, fky2, fkz2,
                 vdx_brick, vdy_brick, vdz_brick, virial_1, vg,vg2,
                 u_brick, v0_brick, v1_brick, v2_brick, v3_brick, v4_brick,
                 v5_brick);

      cg->forward_comm(this, FORWARD_IK);

      if (fix->precision() == FixIntel::PREC_MODE_MIXED) {
        fieldforce_c_ik<float,double>(fix->get_mixed_buffers());
      } else if (fix->precision() == FixIntel::PREC_MODE_DOUBLE) {
        fieldforce_c_ik<double,double>(fix->get_double_buffers());
      } else {
        fieldforce_c_ik<float,float>(fix->get_single_buffers());
      }

      if (evflag_atom) cg_peratom->forward_comm(this, FORWARD_IK_PERATOM);
    }
    if (evflag_atom) fieldforce_c_peratom();
  }

  if (function[1]) {
    //perform calculations for geometric mixing

    if (fix->precision() == FixIntel::PREC_MODE_MIXED) {
      particle_map<float,double>(delxinv_6, delyinv_6, delzinv_6, shift_6,
                                 part2grid_6, nupper_6, nlower_6, nxlo_out_6,
                                 nylo_out_6, nzlo_out_6, nxhi_out_6,
                                 nyhi_out_6, nzhi_out_6,
                                 fix->get_mixed_buffers());
      make_rho_g<float,double>(fix->get_mixed_buffers());
    } else if (fix->precision() == FixIntel::PREC_MODE_DOUBLE) {
      particle_map<double,double>(delxinv_6, delyinv_6, delzinv_6, shift_6,
                                  part2grid_6, nupper_6, nlower_6, nxlo_out_6,
                                  nylo_out_6, nzlo_out_6, nxhi_out_6,
                                  nyhi_out_6, nzhi_out_6,
                                  fix->get_double_buffers());
      make_rho_g<double,double>(fix->get_double_buffers());
    } else {
      particle_map<float,float>(delxinv_6, delyinv_6, delzinv_6, shift_6,
                                part2grid_6, nupper_6, nlower_6, nxlo_out_6,
                                nylo_out_6, nzlo_out_6, nxhi_out_6,
                                nyhi_out_6, nzhi_out_6,
                                fix->get_single_buffers());
      make_rho_g<float,float>(fix->get_single_buffers());
    }


    cg_6->reverse_comm(this, REVERSE_RHO_G);

    brick2fft(nxlo_in_6, nylo_in_6, nzlo_in_6, nxhi_in_6, nyhi_in_6, nzhi_in_6,
              density_brick_g, density_fft_g, work1_6,remap_6);

    if (differentiation_flag == 1) {

      poisson_ad(work1_6, work2_6, density_fft_g, fft1_6, fft2_6,
                 nx_pppm_6, ny_pppm_6, nz_pppm_6, nfft_6,
                 nxlo_fft_6, nylo_fft_6, nzlo_fft_6, nxhi_fft_6,
                 nyhi_fft_6, nzhi_fft_6, nxlo_in_6, nylo_in_6, nzlo_in_6,
                 nxhi_in_6, nyhi_in_6, nzhi_in_6, energy_6, greensfn_6,
                 virial_6, vg_6, vg2_6, u_brick_g, v0_brick_g, v1_brick_g,
                 v2_brick_g, v3_brick_g, v4_brick_g, v5_brick_g);

      cg_6->forward_comm(this,FORWARD_AD_G);

    if (fix->precision() == FixIntel::PREC_MODE_MIXED) {
      fieldforce_g_ad<float,double>(fix->get_mixed_buffers());
    } else if (fix->precision() == FixIntel::PREC_MODE_DOUBLE) {
      fieldforce_g_ad<double,double>(fix->get_double_buffers());
    } else {
      fieldforce_g_ad<float,float>(fix->get_single_buffers());
    }

      if (vflag_atom) cg_peratom_6->forward_comm(this,FORWARD_AD_PERATOM_G);

    } else {
      poisson_ik(work1_6, work2_6, density_fft_g, fft1_6, fft2_6,
                 nx_pppm_6, ny_pppm_6, nz_pppm_6, nfft_6, nxlo_fft_6,
                 nylo_fft_6, nzlo_fft_6, nxhi_fft_6, nyhi_fft_6, nzhi_fft_6,
                 nxlo_in_6, nylo_in_6, nzlo_in_6, nxhi_in_6, nyhi_in_6,
                 nzhi_in_6, energy_6, greensfn_6, fkx_6, fky_6, fkz_6,
                 fkx2_6, fky2_6, fkz2_6, vdx_brick_g, vdy_brick_g,
                 vdz_brick_g, virial_6, vg_6, vg2_6, u_brick_g, v0_brick_g,
                 v1_brick_g, v2_brick_g, v3_brick_g, v4_brick_g, v5_brick_g);

      cg_6->forward_comm(this,FORWARD_IK_G);

    if (fix->precision() == FixIntel::PREC_MODE_MIXED) {
      fieldforce_g_ik<float,double>(fix->get_mixed_buffers());
    } else if (fix->precision() == FixIntel::PREC_MODE_DOUBLE) {
      fieldforce_g_ik<double,double>(fix->get_double_buffers());
    } else {
      fieldforce_g_ik<float,float>(fix->get_single_buffers());
    }


      if (evflag_atom) cg_peratom_6->forward_comm(this, FORWARD_IK_PERATOM_G);
    }
    if (evflag_atom) fieldforce_g_peratom();
  }

  if (function[2]) {
    //perform calculations for arithmetic mixing

    if (fix->precision() == FixIntel::PREC_MODE_MIXED) {
      particle_map<float,double>(delxinv_6, delyinv_6, delzinv_6, shift_6,
                                 part2grid_6, nupper_6, nlower_6,
                                 nxlo_out_6, nylo_out_6, nzlo_out_6,
                                 nxhi_out_6, nyhi_out_6, nzhi_out_6,
                                 fix->get_mixed_buffers());
      make_rho_a<float,double>(fix->get_mixed_buffers());
    } else if (fix->precision() == FixIntel::PREC_MODE_DOUBLE) {
      particle_map<double,double>(delxinv_6, delyinv_6, delzinv_6, shift_6,
                                  part2grid_6, nupper_6, nlower_6, nxlo_out_6,
                                  nylo_out_6, nzlo_out_6, nxhi_out_6,
                                  nyhi_out_6, nzhi_out_6,
                                  fix->get_double_buffers());
      make_rho_a<double,double>(fix->get_double_buffers());
    } else {
      particle_map<float,float>(delxinv_6, delyinv_6, delzinv_6, shift_6,
                                part2grid_6, nupper_6, nlower_6, nxlo_out_6,
                                nylo_out_6, nzlo_out_6, nxhi_out_6,
                                nyhi_out_6, nzhi_out_6,
                                fix->get_single_buffers());
      make_rho_a<float,float>(fix->get_single_buffers());
    }

    cg_6->reverse_comm(this, REVERSE_RHO_A);

    brick2fft_a();

    if ( differentiation_flag == 1) {

      poisson_ad(work1_6, work2_6, density_fft_a3, fft1_6, fft2_6,
                 nx_pppm_6, ny_pppm_6, nz_pppm_6, nfft_6, nxlo_fft_6,
                 nylo_fft_6, nzlo_fft_6, nxhi_fft_6, nyhi_fft_6, nzhi_fft_6,
                 nxlo_in_6, nylo_in_6, nzlo_in_6, nxhi_in_6, nyhi_in_6,
                 nzhi_in_6, energy_6, greensfn_6, virial_6, vg_6, vg2_6,
                 u_brick_a3, v0_brick_a3, v1_brick_a3, v2_brick_a3,
                 v3_brick_a3, v4_brick_a3, v5_brick_a3);
      poisson_2s_ad(density_fft_a0, density_fft_a6, u_brick_a0, v0_brick_a0,
                    v1_brick_a0, v2_brick_a0, v3_brick_a0, v4_brick_a0,
                    v5_brick_a0, u_brick_a6, v0_brick_a6, v1_brick_a6,
                    v2_brick_a6, v3_brick_a6, v4_brick_a6, v5_brick_a6);
      poisson_2s_ad(density_fft_a1, density_fft_a5, u_brick_a1, v0_brick_a1,
                    v1_brick_a1, v2_brick_a1, v3_brick_a1, v4_brick_a1,
                    v5_brick_a1, u_brick_a5, v0_brick_a5, v1_brick_a5,
                    v2_brick_a5, v3_brick_a5, v4_brick_a5, v5_brick_a5);
      poisson_2s_ad(density_fft_a2, density_fft_a4, u_brick_a2, v0_brick_a2,
                    v1_brick_a2, v2_brick_a2, v3_brick_a2, v4_brick_a2,
                    v5_brick_a2, u_brick_a4, v0_brick_a4, v1_brick_a4,
                    v2_brick_a4, v3_brick_a4, v4_brick_a4, v5_brick_a4);

      cg_6->forward_comm(this, FORWARD_AD_A);

    if (fix->precision() == FixIntel::PREC_MODE_MIXED) {
      fieldforce_a_ad<float,double>(fix->get_mixed_buffers());
    } else if (fix->precision() == FixIntel::PREC_MODE_DOUBLE) {
      fieldforce_a_ad<double,double>(fix->get_double_buffers());
    } else {
      fieldforce_a_ad<float,float>(fix->get_single_buffers());
    }

      if (evflag_atom) cg_peratom_6->forward_comm(this, FORWARD_AD_PERATOM_A);

    }  else {

      poisson_ik(work1_6, work2_6, density_fft_a3, fft1_6, fft2_6,
                 nx_pppm_6, ny_pppm_6, nz_pppm_6, nfft_6, nxlo_fft_6,
                 nylo_fft_6, nzlo_fft_6, nxhi_fft_6, nyhi_fft_6, nzhi_fft_6,
                 nxlo_in_6, nylo_in_6, nzlo_in_6, nxhi_in_6, nyhi_in_6,
                 nzhi_in_6, energy_6, greensfn_6, fkx_6, fky_6, fkz_6,fkx2_6,
                 fky2_6, fkz2_6, vdx_brick_a3, vdy_brick_a3, vdz_brick_a3,
                 virial_6, vg_6, vg2_6, u_brick_a3, v0_brick_a3, v1_brick_a3,
                 v2_brick_a3, v3_brick_a3, v4_brick_a3, v5_brick_a3);
      poisson_2s_ik(density_fft_a0, density_fft_a6, vdx_brick_a0,
                    vdy_brick_a0, vdz_brick_a0, vdx_brick_a6, vdy_brick_a6,
                    vdz_brick_a6, u_brick_a0, v0_brick_a0, v1_brick_a0,
                    v2_brick_a0, v3_brick_a0, v4_brick_a0, v5_brick_a0,
                    u_brick_a6, v0_brick_a6, v1_brick_a6, v2_brick_a6,
                    v3_brick_a6, v4_brick_a6, v5_brick_a6);
      poisson_2s_ik(density_fft_a1, density_fft_a5, vdx_brick_a1,
                    vdy_brick_a1, vdz_brick_a1, vdx_brick_a5, vdy_brick_a5,
                    vdz_brick_a5, u_brick_a1, v0_brick_a1, v1_brick_a1,
                    v2_brick_a1, v3_brick_a1, v4_brick_a1, v5_brick_a1,
                    u_brick_a5, v0_brick_a5, v1_brick_a5, v2_brick_a5,
                    v3_brick_a5, v4_brick_a5, v5_brick_a5);
      poisson_2s_ik(density_fft_a2, density_fft_a4, vdx_brick_a2,
                    vdy_brick_a2, vdz_brick_a2, vdx_brick_a4, vdy_brick_a4,
                    vdz_brick_a4, u_brick_a2, v0_brick_a2, v1_brick_a2,
                    v2_brick_a2, v3_brick_a2, v4_brick_a2, v5_brick_a2,
                    u_brick_a4, v0_brick_a4, v1_brick_a4, v2_brick_a4,
                    v3_brick_a4, v4_brick_a4, v5_brick_a4);

      cg_6->forward_comm(this, FORWARD_IK_A);

      if (fix->precision() == FixIntel::PREC_MODE_MIXED) {
        fieldforce_a_ik<float,double>(fix->get_mixed_buffers());
      } else if (fix->precision() == FixIntel::PREC_MODE_DOUBLE) {
        fieldforce_a_ik<double,double>(fix->get_double_buffers());
      } else {
        fieldforce_a_ik<float,float>(fix->get_single_buffers());
      }

      if (evflag_atom) cg_peratom_6->forward_comm(this, FORWARD_IK_PERATOM_A);
    }
    if (evflag_atom) fieldforce_a_peratom();
  }

  if (function[3]) {
    //perform calculations if no mixing rule applies

    if (fix->precision() == FixIntel::PREC_MODE_MIXED) {
      particle_map<float,double>(delxinv_6, delyinv_6, delzinv_6, shift_6,
                                 part2grid_6, nupper_6, nlower_6, nxlo_out_6,
                                 nylo_out_6, nzlo_out_6, nxhi_out_6,
                                 nyhi_out_6, nzhi_out_6,
                                 fix->get_mixed_buffers());
      make_rho_none<float,double>(fix->get_mixed_buffers());
    } else if (fix->precision() == FixIntel::PREC_MODE_DOUBLE) {
      particle_map<double,double>(delxinv_6, delyinv_6, delzinv_6, shift_6,
                                  part2grid_6, nupper_6, nlower_6, nxlo_out_6,
                                  nylo_out_6, nzlo_out_6, nxhi_out_6,
                                  nyhi_out_6, nzhi_out_6,
                                  fix->get_double_buffers());
      make_rho_none<double,double>(fix->get_double_buffers());
    } else {
      particle_map<float,float>(delxinv_6, delyinv_6, delzinv_6, shift_6,
                                part2grid_6, nupper_6, nlower_6, nxlo_out_6,
                                nylo_out_6, nzlo_out_6, nxhi_out_6,
                                nyhi_out_6, nzhi_out_6,
                                fix->get_single_buffers());
      make_rho_none<float,float>(fix->get_single_buffers());
    }

    cg_6->reverse_comm(this, REVERSE_RHO_NONE);

    brick2fft_none();

    if (differentiation_flag == 1) {

      int n = 0;
      for (int k = 0; k<nsplit_alloc/2; k++) {
        poisson_none_ad(n,n+1,density_fft_none[n],density_fft_none[n+1],
                        u_brick_none[n],u_brick_none[n+1],
                        v0_brick_none, v1_brick_none, v2_brick_none,
                        v3_brick_none, v4_brick_none, v5_brick_none);
        n += 2;
      }

      cg_6->forward_comm(this,FORWARD_AD_NONE);

    if (fix->precision() == FixIntel::PREC_MODE_MIXED) {
      fieldforce_none_ad<float,double>(fix->get_mixed_buffers());
    } else if (fix->precision() == FixIntel::PREC_MODE_DOUBLE) {
      fieldforce_none_ad<double,double>(fix->get_double_buffers());
    } else {
      fieldforce_none_ad<float,float>(fix->get_single_buffers());
    }

      if (vflag_atom) cg_peratom_6->forward_comm(this,FORWARD_AD_PERATOM_NONE);

    } else {
      int n = 0;
      for (int k = 0; k<nsplit_alloc/2; k++) {

        poisson_none_ik(n,n+1,density_fft_none[n], density_fft_none[n+1],
                        vdx_brick_none[n], vdy_brick_none[n],
                        vdz_brick_none[n], vdx_brick_none[n+1],
                        vdy_brick_none[n+1], vdz_brick_none[n+1],
                        u_brick_none, v0_brick_none, v1_brick_none,
                        v2_brick_none, v3_brick_none, v4_brick_none,
                        v5_brick_none);
        n += 2;
      }

      cg_6->forward_comm(this,FORWARD_IK_NONE);

    if (fix->precision() == FixIntel::PREC_MODE_MIXED) {
      fieldforce_none_ik<float,double>(fix->get_mixed_buffers());
    } else if (fix->precision() == FixIntel::PREC_MODE_DOUBLE) {
      fieldforce_none_ik<double,double>(fix->get_double_buffers());
    } else {
      fieldforce_none_ik<float,float>(fix->get_single_buffers());
    }

      if (evflag_atom)
        cg_peratom_6->forward_comm(this, FORWARD_IK_PERATOM_NONE);
    }
    if (evflag_atom) fieldforce_none_peratom();
  }

  // update qsum and qsqsum, if atom count has changed and energy needed

  if ((eflag_global || eflag_atom) && atom->natoms != natoms_original) {
    qsum_qsq();
    natoms_original = atom->natoms;
  }

  // sum energy across procs and add in volume-dependent term

  const double qscale = force->qqrd2e * scale;
  if (eflag_global) {
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
  }

  // sum virial across procs

  if (vflag_global) {
    double virial_all[6];
    MPI_Allreduce(virial_1,virial_all,6,MPI_DOUBLE,MPI_SUM,world);
    for (i = 0; i < 6; i++) virial[i] = 0.5*qscale*volume*virial_all[i];
    MPI_Allreduce(virial_6,virial_all,6,MPI_DOUBLE,MPI_SUM,world);
    for (i = 0; i < 6; i++) virial[i] += 0.5*volume*virial_all[i];
    if (function[1]+function[2]+function[3]){
      double a =  MY_PI*MY_PIS/(6*volume)*pow(g_ewald_6,3)*csumij;
      virial[0] -= a;
      virial[1] -= a;
      virial[2] -= a;
    }
  }

  if (eflag_atom) {
    if (function[0]) {
      double *q = atom->q;
      for (i = 0; i < atom->nlocal; i++) {
        eatom[i] -= qscale*g_ewald*q[i]*q[i]/MY_PIS + qscale*MY_PI2*q[i]*
          qsum / (g_ewald*g_ewald*volume); //coulomb self energy correction
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
      for (i = 0; i < atom->nlocal; i++) {
        tmp = atom->type[i];
        //dispersion self virial correction
        for (int n = 0; n < 3; n++) vatom[i][n] -= MY_PI*MY_PIS/(6*volume)*
                                      pow(g_ewald_6,3)*csumi[tmp];
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


/* ---------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   find center grid pt for each of my particles
   check that full stencil for the particle will fit in my 3d brick
   store central grid pt indices in part2grid array
------------------------------------------------------------------------- */

template<class flt_t, class acc_t>
void PPPMDispIntel::particle_map(double delx, double dely, double delz,
                                 double sft, int** p2g, int nup, int nlow,
                                 int nxlo, int nylo, int nzlo,
                                 int nxhi, int nyhi, int nzhi,
                                 IntelBuffers<flt_t,acc_t> * /*buffers*/)
{
  int nlocal = atom->nlocal;
  int nthr = comm->nthreads;

  if (!std::isfinite(boxlo[0]) || !std::isfinite(boxlo[1]) || !std::isfinite(boxlo[2]))
    error->one(FLERR,"Non-numeric box dimensions - simulation unstable");

  int flag = 0;

  #if defined(_OPENMP)
  #pragma omp parallel default(none) \
    shared(nlocal, nthr, delx, dely, delz, sft, p2g, nup, nlow, nxlo,\
           nylo, nzlo, nxhi, nyhi, nzhi) reduction(+:flag) if(!_use_lrt)
  #endif
  {
    double **x = atom->x;
    const flt_t lo0 = boxlo[0];
    const flt_t lo1 = boxlo[1];
    const flt_t lo2 = boxlo[2];
    const flt_t xi = delx;
    const flt_t yi = dely;
    const flt_t zi = delz;
    const flt_t fshift = sft;


    int iifrom, iito, tid;
    IP_PRE_omp_range_id_align(iifrom, iito, tid, nlocal, nthr, sizeof(ATOM_T));

    #if defined(LMP_SIMD_COMPILER)
    #pragma vector aligned
    #pragma simd reduction(+:flag)
    #endif
    for (int i = iifrom; i < iito; i++) {

    // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
    // current particle coord can be outside global and local box
    // add/subtract OFFSET to avoid int(-0.75) = 0 when want it to be -1

    int nx = static_cast<int> ((x[i][0]-lo0)*xi+fshift) - OFFSET;
    int ny = static_cast<int> ((x[i][1]-lo1)*yi+fshift) - OFFSET;
    int nz = static_cast<int> ((x[i][2]-lo2)*zi+fshift) - OFFSET;

    p2g[i][0] = nx;
    p2g[i][1] = ny;
    p2g[i][2] = nz;

    // check that entire stencil around nx,ny,nz will fit in my 3d brick

    if (nx+nlow < nxlo || nx+nup > nxhi ||
        ny+nlow < nylo || ny+nup > nyhi ||
        nz+nlow < nzlo || nz+nup > nzhi)
      flag = 1;
  }
  }

  if (flag) error->one(FLERR,"Out of range atoms - cannot compute PPPMDisp");
}

/* ----------------------------------------------------------------------
   create discretized "density" on section of global grid due to my particles
   density(x,y,z) = charge "density" at grid points of my 3d brick
   (nxlo:nxhi,nylo:nyhi,nzlo:nzhi) is extent of my brick (including ghosts)
   in global grid
------------------------------------------------------------------------- */

template<class flt_t, class acc_t, int use_table>
void PPPMDispIntel::make_rho_c(IntelBuffers<flt_t,acc_t> * /*buffers*/)
{
  // clear 3d density array

  FFT_SCALAR * _noalias global_density =
    &(density_brick[nzlo_out][nylo_out][nxlo_out]);

  // loop over my charges, add their contribution to nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt

  //double *q = atom->q;
  //double **x = atom->x;
  int nlocal = atom->nlocal;
  int nthr = comm->nthreads;

  #if defined(_OPENMP)
  #pragma omp parallel default(none) \
    shared(nthr, nlocal, global_density) if(!_use_lrt)
  #endif
  {
  double *q = atom->q;
  double **x = atom->x;

    const int nix = nxhi_out - nxlo_out + 1;
    const int niy = nyhi_out - nylo_out + 1;

    const flt_t lo0 = boxlo[0];
    const flt_t lo1 = boxlo[1];
    const flt_t lo2 = boxlo[2];
    const flt_t xi = delxinv;
    const flt_t yi = delyinv;
    const flt_t zi = delzinv;
    const flt_t fshift = shift;
    const flt_t fshiftone = shiftone;
    const flt_t fdelvolinv = delvolinv;

    int ifrom, ito, tid;
    IP_PRE_omp_range_id(ifrom, ito, tid, nlocal, nthr);
    FFT_SCALAR * _noalias my_density = tid == 0 ? global_density :
      perthread_density[tid - 1];
    // clear 3d density array
    memset(my_density, 0, ngrid * sizeof(FFT_SCALAR));

    for (int i = ifrom; i < ito; i++) {

      int nx = part2grid[i][0];
      int ny = part2grid[i][1];
      int nz = part2grid[i][2];

      int nysum = nlower + ny - nylo_out;
      int nxsum = nlower + nx - nxlo_out;
      int nzsum = (nlower + nz - nzlo_out)*nix*niy + nysum*nix + nxsum;

      FFT_SCALAR dx = nx+fshiftone - (x[i][0]-lo0)*xi;
      FFT_SCALAR dy = ny+fshiftone - (x[i][1]-lo1)*yi;
      FFT_SCALAR dz = nz+fshiftone - (x[i][2]-lo2)*zi;

      _alignvar(flt_t rho[3][INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};

      if (use_table) {
        dx = dx*half_rho_scale + half_rho_scale_plus;
        int idx = dx;
        dy = dy*half_rho_scale + half_rho_scale_plus;
        int idy = dy;
        dz = dz*half_rho_scale + half_rho_scale_plus;
        int idz = dz;
        #if defined(LMP_SIMD_COMPILER)
        #pragma simd
        #endif
        for(int k = 0; k < INTEL_P3M_ALIGNED_MAXORDER; k++) {
          rho[0][k] = rho_lookup[idx][k];
          rho[1][k] = rho_lookup[idy][k];
          rho[2][k] = rho_lookup[idz][k];
        }
      } else {
        #if defined(LMP_SIMD_COMPILER)
        #pragma simd
        #endif
        for (int k = nlower; k <= nupper; k++) {
          FFT_SCALAR r1,r2,r3;
          r1 = r2 = r3 = ZEROF;

          for (int l = order-1; l >= 0; l--) {
            r1 = rho_coeff[l][k] + r1*dx;
            r2 = rho_coeff[l][k] + r2*dy;
            r3 = rho_coeff[l][k] + r3*dz;
          }
          rho[0][k-nlower] = r1;
          rho[1][k-nlower] = r2;
          rho[2][k-nlower] = r3;
        }
      }

      FFT_SCALAR z0 = fdelvolinv * q[i];

      #if defined(LMP_SIMD_COMPILER)
      #pragma loop_count min(2), max(INTEL_P3M_ALIGNED_MAXORDER), avg(7)
      #endif
      for (int n = 0; n < order; n++) {
        int mz = n*nix*niy + nzsum;
        FFT_SCALAR y0 = z0*rho[2][n];
        #if defined(LMP_SIMD_COMPILER)
        #pragma loop_count min(2), max(INTEL_P3M_ALIGNED_MAXORDER), avg(7)
        #endif
        for (int m = 0; m < order; m++) {
          int mzy = m*nix + mz;
          FFT_SCALAR x0 = y0*rho[1][m];
          #if defined(LMP_SIMD_COMPILER)
          #pragma loop_count min(2), max(INTEL_P3M_ALIGNED_MAXORDER), avg(7)
          #pragma simd
          #endif
          for (int l = 0; l < order; l++) {
            int mzyx = l + mzy;
            my_density[mzyx] += x0*rho[0][l];
          }
        }
      }
    }
  }

  // reduce all the perthread_densities into global_density
  #if defined(_OPENMP)
  #pragma omp parallel default(none) \
    shared(nthr, global_density) if(!_use_lrt)
  #endif
  {
    int ifrom, ito, tid;
    IP_PRE_omp_range_id(ifrom, ito, tid, ngrid, nthr);

    #if defined(LMP_SIMD_COMPILER)
    #pragma simd
    #endif
    for (int i = ifrom; i < ito; i++) {
      for(int j = 1; j < nthr; j++) {
        global_density[i] += perthread_density[j-1][i];
      }
    }
  }
}

/* ----------------------------------------------------------------------
   create discretized "density" on section of global grid due to my particles
   density(x,y,z) = dispersion "density" at grid points of my 3d brick
   (nxlo:nxhi,nylo:nyhi,nzlo:nzhi) is extent of my brick (including ghosts)
   in global grid --- geometric mixing
------------------------------------------------------------------------- */

template<class flt_t, class acc_t, int use_table>
void PPPMDispIntel::make_rho_g(IntelBuffers<flt_t,acc_t> * /*buffers*/)
{
  // clear 3d density array

  FFT_SCALAR * _noalias global_density =
    &(density_brick_g[nzlo_out_6][nylo_out_6][nxlo_out_6]);

  // loop over my charges, add their contribution to nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt

  int nlocal = atom->nlocal;
  int nthr = comm->nthreads;

  #if defined(_OPENMP)
  #pragma omp parallel default(none) \
    shared(nthr, nlocal, global_density) if(!_use_lrt)
  #endif
  {
    int type;
    double **x = atom->x;

    const int nix = nxhi_out_6 - nxlo_out_6 + 1;
    const int niy = nyhi_out_6 - nylo_out_6 + 1;

    const flt_t lo0 = boxlo[0];
    const flt_t lo1 = boxlo[1];
    const flt_t lo2 = boxlo[2];
    const flt_t xi = delxinv_6;
    const flt_t yi = delyinv_6;
    const flt_t zi = delzinv_6;
    const flt_t fshift = shift_6;
    const flt_t fshiftone = shiftone_6;
    const flt_t fdelvolinv = delvolinv_6;

    int ifrom, ito, tid;
    IP_PRE_omp_range_id(ifrom, ito, tid, nlocal, nthr);
    FFT_SCALAR * _noalias my_density = tid == 0 ? global_density :
      perthread_density[tid - 1];

    // clear 3d density array
    memset(my_density, 0, ngrid_6 * sizeof(FFT_SCALAR));

    for (int i = ifrom; i < ito; i++) {

      int nx = part2grid_6[i][0];
      int ny = part2grid_6[i][1];
      int nz = part2grid_6[i][2];

      int nysum = nlower_6 + ny - nylo_out_6;
      int nxsum = nlower_6 + nx - nxlo_out_6;
      int nzsum = (nlower_6 + nz - nzlo_out_6)*nix*niy + nysum*nix + nxsum;

      FFT_SCALAR dx = nx+fshiftone - (x[i][0]-lo0)*xi;
      FFT_SCALAR dy = ny+fshiftone - (x[i][1]-lo1)*yi;
      FFT_SCALAR dz = nz+fshiftone - (x[i][2]-lo2)*zi;

      _alignvar(flt_t rho[3][INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};

      if (use_table) {
        dx = dx*half_rho_scale + half_rho_scale_plus;
        int idx = dx;
        dy = dy*half_rho_scale + half_rho_scale_plus;
        int idy = dy;
        dz = dz*half_rho_scale + half_rho_scale_plus;
        int idz = dz;
        #if defined(LMP_SIMD_COMPILER)
        #pragma simd
        #endif
        for(int k = 0; k < INTEL_P3M_ALIGNED_MAXORDER; k++) {
          rho[0][k] = rho6_lookup[idx][k];
          rho[1][k] = rho6_lookup[idy][k];
          rho[2][k] = rho6_lookup[idz][k];
        }
      } else {
        #if defined(LMP_SIMD_COMPILER)
        #pragma simd
        #endif
        for (int k = nlower_6; k <= nupper_6; k++) {
          FFT_SCALAR r1,r2,r3;
          r1 = r2 = r3 = ZEROF;

          for (int l = order_6-1; l >= 0; l--) {
            r1 = rho_coeff_6[l][k] + r1*dx;
            r2 = rho_coeff_6[l][k] + r2*dy;
            r3 = rho_coeff_6[l][k] + r3*dz;
          }
          rho[0][k-nlower_6] = r1;
          rho[1][k-nlower_6] = r2;
          rho[2][k-nlower_6] = r3;
        }
      }

      type = atom->type[i];
      FFT_SCALAR z0 = fdelvolinv * B[type];

      #if defined(LMP_SIMD_COMPILER)
      #pragma loop_count min(2), max(INTEL_P3M_ALIGNED_MAXORDER), avg(7)
      #endif
      for (int n = 0; n < order_6; n++) {
        int mz = n*nix*niy + nzsum;
        FFT_SCALAR y0 = z0*rho[2][n];
        #if defined(LMP_SIMD_COMPILER)
        #pragma loop_count min(2), max(INTEL_P3M_ALIGNED_MAXORDER), avg(7)
        #endif
        for (int m = 0; m < order_6; m++) {
          int mzy = m*nix + mz;
          FFT_SCALAR x0 = y0*rho[1][m];
          #if defined(LMP_SIMD_COMPILER)
          #pragma loop_count min(2), max(INTEL_P3M_ALIGNED_MAXORDER), avg(7)
          #pragma simd
          #endif
          for (int l = 0; l < order; l++) {
            int mzyx = l + mzy;
            my_density[mzyx] += x0*rho[0][l];
          }
        }
      }
    }
  }

  // reduce all the perthread_densities into global_density
  #if defined(_OPENMP)
  #pragma omp parallel default(none) \
    shared(nthr, global_density) if(!_use_lrt)
  #endif
  {
    int ifrom, ito, tid;
    IP_PRE_omp_range_id(ifrom, ito, tid, ngrid_6, nthr);

    #if defined(LMP_SIMD_COMPILER)
    #pragma simd
    #endif
    for (int i = ifrom; i < ito; i++) {
      for(int j = 1; j < nthr; j++) {
        global_density[i] += perthread_density[j-1][i];
      }
    }
  }

}

/* ----------------------------------------------------------------------
   create discretized "density" on section of global grid due to my particles
   density(x,y,z) = dispersion "density" at grid points of my 3d brick
   (nxlo:nxhi,nylo:nyhi,nzlo:nzhi) is extent of my brick (including ghosts)
   in global grid --- arithmetic mixing
------------------------------------------------------------------------- */

template<class flt_t, class acc_t, int use_table>
void PPPMDispIntel::make_rho_a(IntelBuffers<flt_t,acc_t> * /*buffers*/)
{
  // clear 3d density array

  memset(&(density_brick_a0[nzlo_out_6][nylo_out_6][nxlo_out_6]),0,
         ngrid_6*sizeof(FFT_SCALAR));
  memset(&(density_brick_a1[nzlo_out_6][nylo_out_6][nxlo_out_6]),0,
         ngrid_6*sizeof(FFT_SCALAR));
  memset(&(density_brick_a2[nzlo_out_6][nylo_out_6][nxlo_out_6]),0,
         ngrid_6*sizeof(FFT_SCALAR));
  memset(&(density_brick_a3[nzlo_out_6][nylo_out_6][nxlo_out_6]),0,
         ngrid_6*sizeof(FFT_SCALAR));
  memset(&(density_brick_a4[nzlo_out_6][nylo_out_6][nxlo_out_6]),0,
         ngrid_6*sizeof(FFT_SCALAR));
  memset(&(density_brick_a5[nzlo_out_6][nylo_out_6][nxlo_out_6]),0,
         ngrid_6*sizeof(FFT_SCALAR));
  memset(&(density_brick_a6[nzlo_out_6][nylo_out_6][nxlo_out_6]),0,
         ngrid_6*sizeof(FFT_SCALAR));

  // loop over my charges, add their contribution to nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt

  int nlocal = atom->nlocal;

    double **x = atom->x;

    const int nix = nxhi_out_6 - nxlo_out_6 + 1;
    const int niy = nyhi_out_6 - nylo_out_6 + 1;

    const flt_t lo0 = boxlo[0];
    const flt_t lo1 = boxlo[1];
    const flt_t lo2 = boxlo[2];
    const flt_t xi = delxinv_6;
    const flt_t yi = delyinv_6;
    const flt_t zi = delzinv_6;
    const flt_t fshift = shift_6;
    const flt_t fshiftone = shiftone_6;
    const flt_t fdelvolinv = delvolinv_6;

    for (int i = 0; i < nlocal; i++) {

      int nx = part2grid_6[i][0];
      int ny = part2grid_6[i][1];
      int nz = part2grid_6[i][2];

      int nxsum = nx + nlower_6;
      int nysum = ny + nlower_6;
      int nzsum = nz + nlower_6;

      FFT_SCALAR dx = nx+fshiftone - (x[i][0]-lo0)*xi;
      FFT_SCALAR dy = ny+fshiftone - (x[i][1]-lo1)*yi;
      FFT_SCALAR dz = nz+fshiftone - (x[i][2]-lo2)*zi;

      _alignvar(flt_t rho[3][INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};

      if (use_table) {
        dx = dx*half_rho_scale + half_rho_scale_plus;
        int idx = dx;
        dy = dy*half_rho_scale + half_rho_scale_plus;
        int idy = dy;
        dz = dz*half_rho_scale + half_rho_scale_plus;
        int idz = dz;
        #if defined(LMP_SIMD_COMPILER)
        #pragma simd
        #endif
        for(int k = 0; k < INTEL_P3M_ALIGNED_MAXORDER; k++) {
          rho[0][k] = rho6_lookup[idx][k];
          rho[1][k] = rho6_lookup[idy][k];
          rho[2][k] = rho6_lookup[idz][k];
        }
      } else {
        #if defined(LMP_SIMD_COMPILER)
        #pragma simd
        #endif
        for (int k = nlower_6; k <= nupper_6; k++) {
          FFT_SCALAR r1,r2,r3;
          r1 = r2 = r3 = ZEROF;

          for (int l = order_6-1; l >= 0; l--) {
            r1 = rho_coeff_6[l][k] + r1*dx;
            r2 = rho_coeff_6[l][k] + r2*dy;
            r3 = rho_coeff_6[l][k] + r3*dz;
          }
          rho[0][k-nlower_6] = r1;
          rho[1][k-nlower_6] = r2;
          rho[2][k-nlower_6] = r3;
        }
      }

      const int type = atom->type[i];
      FFT_SCALAR z0 = fdelvolinv;

      #if defined(LMP_SIMD_COMPILER)
      #pragma loop_count min(2), max(INTEL_P3M_ALIGNED_MAXORDER), avg(7)
      #endif
      for (int n = 0; n < order_6; n++) {
        int mz = n + nzsum;
        FFT_SCALAR y0 = z0*rho[2][n];
        #if defined(LMP_SIMD_COMPILER)
        #pragma loop_count min(2), max(INTEL_P3M_ALIGNED_MAXORDER), avg(7)
        #endif
        for (int m = 0; m < order_6; m++) {
          int my = m + nysum;
          FFT_SCALAR x0 = y0*rho[1][m];
          #if defined(LMP_SIMD_COMPILER)
          #pragma loop_count min(2), max(INTEL_P3M_ALIGNED_MAXORDER), avg(7)
          #pragma simd
          #endif
          for (int l = 0; l < order; l++) {
            int mx = l + nxsum;
            FFT_SCALAR w = x0*rho[0][l];
            density_brick_a0[mz][my][mx] += w*B[7*type];
            density_brick_a1[mz][my][mx] += w*B[7*type+1];
            density_brick_a2[mz][my][mx] += w*B[7*type+2];
            density_brick_a3[mz][my][mx] += w*B[7*type+3];
            density_brick_a4[mz][my][mx] += w*B[7*type+4];
            density_brick_a5[mz][my][mx] += w*B[7*type+5];
            density_brick_a6[mz][my][mx] += w*B[7*type+6];
          }
        }
      }
    }
}

/* ----------------------------------------------------------------------
   create discretized "density" on section of global grid due to my particles
   density(x,y,z) = dispersion "density" at grid points of my 3d brick
   (nxlo:nxhi,nylo:nyhi,nzlo:nzhi) is extent of my brick (including ghosts)
   in global grid --- case when mixing rules don't apply
------------------------------------------------------------------------- */

template<class flt_t, class acc_t, int use_table>
void PPPMDispIntel::make_rho_none(IntelBuffers<flt_t,acc_t> * /*buffers*/)
{

  FFT_SCALAR * _noalias global_density = &(density_brick_none[0][nzlo_out_6][nylo_out_6][nxlo_out_6]);

  // loop over my charges, add their contribution to nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt

  int nlocal = atom->nlocal;
  int nthr = comm->nthreads;

  #if defined(_OPENMP)
  #pragma omp parallel default(none) \
    shared(nthr, nlocal, global_density) if(!_use_lrt)
  #endif
  {
    int type;
    double **x = atom->x;

    const int nix = nxhi_out_6 - nxlo_out_6 + 1;
    const int niy = nyhi_out_6 - nylo_out_6 + 1;

    const flt_t lo0 = boxlo[0];
    const flt_t lo1 = boxlo[1];
    const flt_t lo2 = boxlo[2];
    const flt_t xi = delxinv_6;
    const flt_t yi = delyinv_6;
    const flt_t zi = delzinv_6;
    const flt_t fshift = shift_6;
    const flt_t fshiftone = shiftone_6;
    const flt_t fdelvolinv = delvolinv_6;

    int ifrom, ito, tid;
    IP_PRE_omp_range_id(ifrom, ito, tid, nlocal, nthr);
    FFT_SCALAR * _noalias my_density = tid == 0 ? global_density :
      perthread_density[tid - 1];
    // clear 3d density array
    memset(my_density, 0, ngrid_6 * sizeof(FFT_SCALAR));

    for (int i = ifrom; i < ito; i++) {

      int nx = part2grid_6[i][0];
      int ny = part2grid_6[i][1];
      int nz = part2grid_6[i][2];

      int nysum = nlower_6 + ny - nylo_out_6;
      int nxsum = nlower_6 + nx - nxlo_out_6;
      int nzsum = (nlower_6 + nz - nzlo_out_6)*nix*niy + nysum*nix + nxsum;

      FFT_SCALAR dx = nx+fshiftone - (x[i][0]-lo0)*xi;
      FFT_SCALAR dy = ny+fshiftone - (x[i][1]-lo1)*yi;
      FFT_SCALAR dz = nz+fshiftone - (x[i][2]-lo2)*zi;

      _alignvar(flt_t rho[3][INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};

      if (use_table) {
        dx = dx*half_rho_scale + half_rho_scale_plus;
        int idx = dx;
        dy = dy*half_rho_scale + half_rho_scale_plus;
        int idy = dy;
        dz = dz*half_rho_scale + half_rho_scale_plus;
        int idz = dz;
        #if defined(LMP_SIMD_COMPILER)
        #pragma simd
        #endif
        for(int k = 0; k < INTEL_P3M_ALIGNED_MAXORDER; k++) {
          rho[0][k] = rho6_lookup[idx][k];
          rho[1][k] = rho6_lookup[idy][k];
          rho[2][k] = rho6_lookup[idz][k];
        }
      } else {
        #if defined(LMP_SIMD_COMPILER)
        #pragma simd
        #endif
        for (int k = nlower_6; k <= nupper_6; k++) {
          FFT_SCALAR r1,r2,r3;
          r1 = r2 = r3 = ZEROF;

          for (int l = order_6-1; l >= 0; l--) {
            r1 = rho_coeff_6[l][k] + r1*dx;
            r2 = rho_coeff_6[l][k] + r2*dy;
            r3 = rho_coeff_6[l][k] + r3*dz;
          }
          rho[0][k-nlower_6] = r1;
          rho[1][k-nlower_6] = r2;
          rho[2][k-nlower_6] = r3;
        }
      }

      type = atom->type[i];
      FFT_SCALAR z0 = fdelvolinv;

      #if defined(LMP_SIMD_COMPILER)
      #pragma loop_count min(2), max(INTEL_P3M_ALIGNED_MAXORDER), avg(7)
      #endif
      for (int n = 0; n < order_6; n++) {
        int mz = n*nix*niy + nzsum;
        FFT_SCALAR y0 = z0*rho[2][n];
        #if defined(LMP_SIMD_COMPILER)
        #pragma loop_count min(2), max(INTEL_P3M_ALIGNED_MAXORDER), avg(7)
        #endif
        for (int m = 0; m < order_6; m++) {
          int mzy = m*nix + mz;
          FFT_SCALAR x0 = y0*rho[1][m];
          #if defined(LMP_SIMD_COMPILER)
          #pragma loop_count min(2), max(INTEL_P3M_ALIGNED_MAXORDER), avg(7)
          #pragma simd
          #endif
          for (int l = 0; l < order; l++) {
            int mzyx = l + mzy;
            FFT_SCALAR w0 = x0*rho[0][l];
            for(int k = 0; k < nsplit; k++)
              my_density[mzyx + k*ngrid_6] += x0*rho[0][l];
          }
        }
      }
    }
  }

  // reduce all the perthread_densities into global_density
  #if defined(_OPENMP)
  #pragma omp parallel default(none) \
    shared(nthr, global_density) if(!_use_lrt)
  #endif
  {
    int ifrom, ito, tid;
    IP_PRE_omp_range_id(ifrom, ito, tid, ngrid_6*nsplit, nthr);

    #if defined(LMP_SIMD_COMPILER)
    #pragma simd
    #endif
    for (int i = ifrom; i < ito; i++) {
      for(int j = 1; j < nthr; j++) {
        global_density[i] += perthread_density[j-1][i];
      }
    }
  }

}

/* ----------------------------------------------------------------------
   interpolate from grid to get electric field & force on my particles
   for ik scheme
------------------------------------------------------------------------- */

template<class flt_t, class acc_t, int use_table>
void PPPMDispIntel::fieldforce_c_ik(IntelBuffers<flt_t,acc_t> * /*buffers*/)
{

  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of E-field on particle

  //double *q = atom->q;
  //double **x = atom->x;
  //double **f = atom->f;

  int nlocal = atom->nlocal;
  int nthr = comm->nthreads;


  #if defined(_OPENMP)
  #pragma omp parallel default(none) \
    shared(nlocal, nthr) if(!_use_lrt)
  #endif
  {

    double *q = atom->q;
    double **x = atom->x;
    double **f = atom->f;

    const flt_t lo0 = boxlo[0];
    const flt_t lo1 = boxlo[1];
    const flt_t lo2 = boxlo[2];
    const flt_t xi = delxinv;
    const flt_t yi = delyinv;
    const flt_t zi = delzinv;
    const flt_t fshiftone = shiftone;
    const flt_t fqqrd2es = qqrd2e * scale;

    int ifrom, ito, tid;
    IP_PRE_omp_range_id(ifrom, ito, tid, nlocal, nthr);

    _alignvar(flt_t rho0[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
    _alignvar(flt_t rho1[INTEL_P3M_ALIGNED_MAXORDER] , 64)= {0};
    _alignvar(flt_t rho2[INTEL_P3M_ALIGNED_MAXORDER] , 64)= {0};

    for (int i = ifrom; i < ito; i++) {
      int nx = part2grid[i][0];
      int ny = part2grid[i][1];
      int nz = part2grid[i][2];

      int nxsum = nx + nlower;
      int nysum = ny + nlower;
      int nzsum = nz + nlower;;

      FFT_SCALAR dx = nx+fshiftone - (x[i][0]-lo0)*xi;
      FFT_SCALAR dy = ny+fshiftone - (x[i][1]-lo1)*yi;
      FFT_SCALAR dz = nz+fshiftone - (x[i][2]-lo2)*zi;

      if (use_table) {
        dx = dx*half_rho_scale + half_rho_scale_plus;
        int idx = dx;
        dy = dy*half_rho_scale + half_rho_scale_plus;
        int idy = dy;
        dz = dz*half_rho_scale + half_rho_scale_plus;
        int idz = dz;
        #if defined(LMP_SIMD_COMPILER)
        #pragma simd
        #endif
        for (int k = 0; k < INTEL_P3M_ALIGNED_MAXORDER; k++) {
          rho0[k] = rho_lookup[idx][k];
          rho1[k] = rho_lookup[idy][k];
          rho2[k] = rho_lookup[idz][k];
        }
      } else {
        #if defined(LMP_SIMD_COMPILER)
        #pragma simd
        #endif
        for (int k = nlower; k <= nupper; k++) {
          FFT_SCALAR r1 = rho_coeff[order-1][k];
          FFT_SCALAR r2 = rho_coeff[order-1][k];
          FFT_SCALAR r3 = rho_coeff[order-1][k];
          for (int l = order-2; l >= 0; l--) {
            r1 = rho_coeff[l][k] + r1*dx;
            r2 = rho_coeff[l][k] + r2*dy;
            r3 = rho_coeff[l][k] + r3*dz;
          }

          rho0[k-nlower] = r1;
          rho1[k-nlower] = r2;
          rho2[k-nlower] = r3;
        }
      }

      _alignvar(FFT_SCALAR ekx_arr[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
      _alignvar(FFT_SCALAR eky_arr[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
      _alignvar(FFT_SCALAR ekz_arr[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};

      #if defined(LMP_SIMD_COMPILER)
      #pragma loop_count min(2), max(INTEL_P3M_ALIGNED_MAXORDER), avg(7)
      #endif
      for (int n = 0; n < order; n++) {
        int mz = n+nzsum;
        FFT_SCALAR z0 = rho2[n];
        #if defined(LMP_SIMD_COMPILER)
        #pragma loop_count min(2), max(INTEL_P3M_ALIGNED_MAXORDER), avg(7)
        #endif
        for (int m = 0; m < order; m++) {
          int my = m+nysum;
          FFT_SCALAR y0 = z0*rho1[m];
          #if defined(LMP_SIMD_COMPILER)
          #pragma loop_count min(2), max(INTEL_P3M_ALIGNED_MAXORDER), avg(7)
          #pragma simd
          #endif
          for (int l = 0; l < order; l++) {
            int mx = l+nxsum;
            FFT_SCALAR x0 = y0*rho0[l];
              ekx_arr[l] -= x0*vdx_brick[mz][my][mx];
              eky_arr[l] -= x0*vdy_brick[mz][my][mx];
              ekz_arr[l] -= x0*vdz_brick[mz][my][mx];

          }
        }
      }

      FFT_SCALAR ekx, eky, ekz;
      ekx = eky = ekz = ZEROF;

      for (int l = 0; l < order; l++) {
        ekx += ekx_arr[l];
        eky += eky_arr[l];
        ekz += ekz_arr[l];
      }

      // convert E-field to force

      const flt_t qfactor = fqqrd2es * q[i];
      f[i][0] += qfactor*ekx;
      f[i][1] += qfactor*eky;
      if (slabflag != 2) f[i][2] += qfactor*ekz;
    }
  }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get electric field & force on my particles
   for ad scheme
------------------------------------------------------------------------- */

template<class flt_t, class acc_t, int use_table>
void PPPMDispIntel::fieldforce_c_ad(IntelBuffers<flt_t,acc_t> * /*buffers*/)
{

  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of E-field on particle

  //double *q = atom->q;
  //double **x = atom->x;
  //double **f = atom->f;

  int nlocal = atom->nlocal;
  int nthr = comm->nthreads;

  FFT_SCALAR * _noalias const particle_ekx = this->particle_ekx;
  FFT_SCALAR * _noalias const particle_eky = this->particle_eky;
  FFT_SCALAR * _noalias const particle_ekz = this->particle_ekz;

  #if defined(_OPENMP)
  #pragma omp parallel default(none) \
    shared(nlocal, nthr) if(!_use_lrt)
  #endif
  {

    double *prd;
    if (triclinic == 0) prd = domain->prd;
    else prd = domain->prd_lamda;

    double *q = atom->q;
    double **x = atom->x;
    double **f = atom->f;
    const flt_t ftwo_pi = MY_PI * 2.0;
    const flt_t ffour_pi = MY_PI * 4.0;

    const flt_t lo0 = boxlo[0];
    const flt_t lo1 = boxlo[1];
    const flt_t lo2 = boxlo[2];
    const flt_t xi = delxinv;
    const flt_t yi = delyinv;
    const flt_t zi = delzinv;
    const flt_t fshiftone = shiftone;
    const flt_t fqqrd2es = qqrd2e * scale;

    const double xprd = prd[0];
    const double yprd = prd[1];
    const double zprd = prd[2]*slab_volfactor;

    const flt_t hx_inv = nx_pppm/xprd;
    const flt_t hy_inv = ny_pppm/yprd;
    const flt_t hz_inv = nz_pppm/zprd;

    const flt_t fsf_coeff0 = sf_coeff[0];
    const flt_t fsf_coeff1 = sf_coeff[1];
    const flt_t fsf_coeff2 = sf_coeff[2];
    const flt_t fsf_coeff3 = sf_coeff[3];
    const flt_t fsf_coeff4 = sf_coeff[4];
    const flt_t fsf_coeff5 = sf_coeff[5];

    int ifrom, ito, tid;
    IP_PRE_omp_range_id(ifrom, ito, tid, nlocal, nthr);

    _alignvar(flt_t rho[3][INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
    _alignvar(flt_t drho[3][INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};

    for (int i = ifrom; i < ito; i++) {
      int nx = part2grid[i][0];
      int ny = part2grid[i][1];
      int nz = part2grid[i][2];
      FFT_SCALAR dx = nx+fshiftone - (x[i][0]-lo0)*xi;
      FFT_SCALAR dy = ny+fshiftone - (x[i][1]-lo1)*yi;
      FFT_SCALAR dz = nz+fshiftone - (x[i][2]-lo2)*zi;

      int nxsum = nx + nlower;
      int nysum = ny + nlower;
      int nzsum = nz + nlower;

      if (use_table) {
        dx = dx*half_rho_scale + half_rho_scale_plus;
        int idx = dx;
        dy = dy*half_rho_scale + half_rho_scale_plus;
        int idy = dy;
        dz = dz*half_rho_scale + half_rho_scale_plus;
        int idz = dz;

        #if defined(LMP_SIMD_COMPILER)
        #pragma simd
        #endif
        for(int k = 0; k < INTEL_P3M_ALIGNED_MAXORDER; k++) {
          rho[0][k] = rho_lookup[idx][k];
          rho[1][k] = rho_lookup[idy][k];
          rho[2][k] = rho_lookup[idz][k];
          drho[0][k] = drho_lookup[idx][k];
          drho[1][k] = drho_lookup[idy][k];
          drho[2][k] = drho_lookup[idz][k];
        }
      } else {
        #if defined(LMP_SIMD_COMPILER)
        #pragma simd
        #endif
        for (int k = nlower; k <= nupper; k++) {
          FFT_SCALAR r1,r2,r3,dr1,dr2,dr3;
          dr1 = dr2 = dr3 = ZEROF;

          r1 = rho_coeff[order-1][k];
          r2 = rho_coeff[order-1][k];
          r3 = rho_coeff[order-1][k];
          for (int l = order-2; l >= 0; l--) {
            r1 = rho_coeff[l][k] + r1 * dx;
            r2 = rho_coeff[l][k] + r2 * dy;
            r3 = rho_coeff[l][k] + r3 * dz;
            dr1 = drho_coeff[l][k] + dr1 * dx;
            dr2 = drho_coeff[l][k] + dr2 * dy;
            dr3 = drho_coeff[l][k] + dr3 * dz;
          }
          rho[0][k-nlower] = r1;
          rho[1][k-nlower] = r2;
          rho[2][k-nlower] = r3;
          drho[0][k-nlower] = dr1;
          drho[1][k-nlower] = dr2;
          drho[2][k-nlower] = dr3;
        }
      }
      _alignvar(FFT_SCALAR ekx[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
      _alignvar(FFT_SCALAR eky[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
      _alignvar(FFT_SCALAR ekz[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};

      particle_ekx[i] = particle_eky[i] = particle_ekz[i] = ZEROF;
      #if defined(LMP_SIMD_COMPILER)
      #pragma loop_count min(2), max(INTEL_P3M_ALIGNED_MAXORDER), avg(7)
      #endif
      for (int n = 0; n < order; n++) {
        int mz = n + nzsum;
        #if defined(LMP_SIMD_COMPILER)
        #pragma loop_count min(2), max(INTEL_P3M_ALIGNED_MAXORDER), avg(7)
        #endif
        for (int m = 0; m < order; m++) {
          int my = m + nysum;
          FFT_SCALAR ekx_p = rho[1][m] * rho[2][n];
          FFT_SCALAR eky_p = drho[1][m] * rho[2][n];
          FFT_SCALAR ekz_p = rho[1][m] * drho[2][n];
          #if defined(LMP_SIMD_COMPILER)
          #pragma loop_count min(2), max(INTEL_P3M_ALIGNED_MAXORDER), avg(7)
          #pragma simd
          #endif
          for (int l = 0; l < order; l++) {
            int mx = l + nxsum;
            ekx[l] += drho[0][l] * ekx_p * u_brick[mz][my][mx];
            eky[l] +=  rho[0][l] * eky_p * u_brick[mz][my][mx];
            ekz[l] +=  rho[0][l] * ekz_p * u_brick[mz][my][mx];
          }
        }
      }

      #if defined(LMP_SIMD_COMPILER)
      #pragma loop_count min(2), max(INTEL_P3M_ALIGNED_MAXORDER), avg(7)
      #endif
      for (int l = 0; l < order; l++){
        particle_ekx[i] += ekx[l];
        particle_eky[i] += eky[l];
        particle_ekz[i] += ekz[l];
      }
    }
    #if defined(LMP_SIMD_COMPILER)
    #pragma simd
    #endif
    for (int i = ifrom; i < ito; i++) {
      particle_ekx[i] *= hx_inv;
      particle_eky[i] *= hy_inv;
      particle_ekz[i] *= hz_inv;

      // convert E-field to force

      const flt_t qfactor = fqqrd2es * q[i];
      const flt_t twoqsq = (flt_t)2.0 * q[i] * q[i];

      const flt_t s1 = x[i][0] * hx_inv;
      const flt_t s2 = x[i][1] * hy_inv;
      const flt_t s3 = x[i][2] * hz_inv;
      flt_t sf = fsf_coeff0 * sin(ftwo_pi * s1);
      sf += fsf_coeff1 * sin(ffour_pi * s1);
      sf *= twoqsq;
      f[i][0] += qfactor * particle_ekx[i] - fqqrd2es * sf;

      sf = fsf_coeff2 * sin(ftwo_pi * s2);
      sf += fsf_coeff3 * sin(ffour_pi * s2);
      sf *= twoqsq;
      f[i][1] += qfactor * particle_eky[i] - fqqrd2es * sf;

      sf = fsf_coeff4 * sin(ftwo_pi * s3);
      sf += fsf_coeff5 * sin(ffour_pi * s3);
      sf *= twoqsq;

      if (slabflag != 2) f[i][2] += qfactor * particle_ekz[i] - fqqrd2es * sf;
    }
  }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get dispersion field & force on my particles
   for geometric mixing rule
------------------------------------------------------------------------- */

template<class flt_t, class acc_t, int use_table>
void PPPMDispIntel::fieldforce_g_ik(IntelBuffers<flt_t,acc_t> * /*buffers*/)
{

  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of dispersion field on particle

  int nlocal = atom->nlocal;
  int nthr = comm->nthreads;


  #if defined(_OPENMP)
  #pragma omp parallel default(none) \
    shared(nlocal, nthr) if(!_use_lrt)
  #endif
  {

    double lj;
    int type;
    double **x = atom->x;
    double **f = atom->f;

    const flt_t lo0 = boxlo[0];
    const flt_t lo1 = boxlo[1];
    const flt_t lo2 = boxlo[2];
    const flt_t xi = delxinv_6;
    const flt_t yi = delyinv_6;
    const flt_t zi = delzinv_6;
    const flt_t fshiftone = shiftone_6;

    int ifrom, ito, tid;
    IP_PRE_omp_range_id(ifrom, ito, tid, nlocal, nthr);

    _alignvar(flt_t rho0[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
    _alignvar(flt_t rho1[INTEL_P3M_ALIGNED_MAXORDER] , 64)= {0};
    _alignvar(flt_t rho2[INTEL_P3M_ALIGNED_MAXORDER] , 64)= {0};

    for (int i = ifrom; i < ito; i++) {
      int nx = part2grid_6[i][0];
      int ny = part2grid_6[i][1];
      int nz = part2grid_6[i][2];

      int nxsum = nx + nlower_6;
      int nysum = ny + nlower_6;
      int nzsum = nz + nlower_6;

      FFT_SCALAR dx = nx+fshiftone - (x[i][0]-lo0)*xi;
      FFT_SCALAR dy = ny+fshiftone - (x[i][1]-lo1)*yi;
      FFT_SCALAR dz = nz+fshiftone - (x[i][2]-lo2)*zi;

      if (use_table) {
        dx = dx*half_rho_scale + half_rho_scale_plus;
        int idx = dx;
        dy = dy*half_rho_scale + half_rho_scale_plus;
        int idy = dy;
        dz = dz*half_rho_scale + half_rho_scale_plus;
        int idz = dz;
        #if defined(LMP_SIMD_COMPILER)
        #pragma simd
        #endif
        for (int k = 0; k < INTEL_P3M_ALIGNED_MAXORDER; k++) {
          rho0[k] = rho6_lookup[idx][k];
          rho1[k] = rho6_lookup[idy][k];
          rho2[k] = rho6_lookup[idz][k];
        }
      } else {
        #if defined(LMP_SIMD_COMPILER)
        #pragma simd
        #endif
        for (int k = nlower_6; k <= nupper_6; k++) {
          FFT_SCALAR r1 = rho_coeff_6[order_6-1][k];
          FFT_SCALAR r2 = rho_coeff_6[order_6-1][k];
          FFT_SCALAR r3 = rho_coeff_6[order_6-1][k];
          for (int l = order_6-2; l >= 0; l--) {
            r1 = rho_coeff_6[l][k] + r1*dx;
            r2 = rho_coeff_6[l][k] + r2*dy;
            r3 = rho_coeff_6[l][k] + r3*dz;
          }

          rho0[k-nlower_6] = r1;
          rho1[k-nlower_6] = r2;
          rho2[k-nlower_6] = r3;
        }
      }

      _alignvar(FFT_SCALAR ekx_arr[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
      _alignvar(FFT_SCALAR eky_arr[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
      _alignvar(FFT_SCALAR ekz_arr[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};

      #if defined(LMP_SIMD_COMPILER)
      #pragma loop_count min(2), max(INTEL_P3M_ALIGNED_MAXORDER), avg(7)
      #endif
      for (int n = 0; n < order_6; n++) {
        int mz = n+nzsum;
        FFT_SCALAR z0 = rho2[n];
        #if defined(LMP_SIMD_COMPILER)
        #pragma loop_count min(2), max(INTEL_P3M_ALIGNED_MAXORDER), avg(7)
        #endif
        for (int m = 0; m < order_6; m++) {
          int my = m+nysum;
          FFT_SCALAR y0 = z0*rho1[m];
          #if defined(LMP_SIMD_COMPILER)
          #pragma loop_count min(2), max(INTEL_P3M_ALIGNED_MAXORDER), avg(7)
          #pragma simd
          #endif
          for (int l = 0; l < order; l++) {
            int mx = l+nxsum;
            FFT_SCALAR x0 = y0*rho0[l];
              ekx_arr[l] -= x0*vdx_brick_g[mz][my][mx];
              eky_arr[l] -= x0*vdy_brick_g[mz][my][mx];
              ekz_arr[l] -= x0*vdz_brick_g[mz][my][mx];

          }
        }
      }

      FFT_SCALAR ekx, eky, ekz;
      ekx = eky = ekz = ZEROF;

      for (int l = 0; l < order; l++) {
        ekx += ekx_arr[l];
        eky += eky_arr[l];
        ekz += ekz_arr[l];
      }

      // convert E-field to force

      type = atom->type[i];
      lj = B[type];
      f[i][0] += lj*ekx;
      f[i][1] += lj*eky;
      if (slabflag != 2) f[i][2] += lj*ekz;
    }
  }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get dispersion field & force on my particles
   for geometric mixing rule for ad scheme
------------------------------------------------------------------------- */

template<class flt_t, class acc_t, int use_table>
void PPPMDispIntel::fieldforce_g_ad(IntelBuffers<flt_t,acc_t> * /*buffers*/)
{

  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of dispersion field on particle

  int nlocal = atom->nlocal;
  int nthr = comm->nthreads;

  FFT_SCALAR * _noalias const particle_ekx = this->particle_ekx;
  FFT_SCALAR * _noalias const particle_eky = this->particle_eky;
  FFT_SCALAR * _noalias const particle_ekz = this->particle_ekz;

  #if defined(_OPENMP)
  #pragma omp parallel default(none) \
    shared(nlocal, nthr) if(!_use_lrt)
  #endif
  {

    double *prd;
    if (triclinic == 0) prd = domain->prd;
    else prd = domain->prd_lamda;

    double **x = atom->x;
    double **f = atom->f;
    const flt_t ftwo_pi = MY_PI * 2.0;
    const flt_t ffour_pi = MY_PI * 4.0;

    const flt_t lo0 = boxlo[0];
    const flt_t lo1 = boxlo[1];
    const flt_t lo2 = boxlo[2];
    const flt_t xi = delxinv_6;
    const flt_t yi = delyinv_6;
    const flt_t zi = delzinv_6;
    const flt_t fshiftone = shiftone_6;

    const double xprd = prd[0];
    const double yprd = prd[1];
    const double zprd = prd[2]*slab_volfactor;

    const flt_t hx_inv = nx_pppm_6/xprd;
    const flt_t hy_inv = ny_pppm_6/yprd;
    const flt_t hz_inv = nz_pppm_6/zprd;

    const flt_t fsf_coeff0 = sf_coeff_6[0];
    const flt_t fsf_coeff1 = sf_coeff_6[1];
    const flt_t fsf_coeff2 = sf_coeff_6[2];
    const flt_t fsf_coeff3 = sf_coeff_6[3];
    const flt_t fsf_coeff4 = sf_coeff_6[4];
    const flt_t fsf_coeff5 = sf_coeff_6[5];

    int ifrom, ito, tid;
    IP_PRE_omp_range_id(ifrom, ito, tid, nlocal, nthr);

    _alignvar(flt_t rho[3][INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
    _alignvar(flt_t drho[3][INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};

    for (int i = ifrom; i < ito; i++) {
      int nx = part2grid_6[i][0];
      int ny = part2grid_6[i][1];
      int nz = part2grid_6[i][2];
      FFT_SCALAR dx = nx+fshiftone - (x[i][0]-lo0)*xi;
      FFT_SCALAR dy = ny+fshiftone - (x[i][1]-lo1)*yi;
      FFT_SCALAR dz = nz+fshiftone - (x[i][2]-lo2)*zi;

      int nxsum = nx + nlower_6;
      int nysum = ny + nlower_6;
      int nzsum = nz + nlower_6;

      if (use_table) {
        dx = dx*half_rho_scale + half_rho_scale_plus;
        int idx = dx;
        dy = dy*half_rho_scale + half_rho_scale_plus;
        int idy = dy;
        dz = dz*half_rho_scale + half_rho_scale_plus;
        int idz = dz;

        #if defined(LMP_SIMD_COMPILER)
        #pragma simd
        #endif
        for(int k = 0; k < INTEL_P3M_ALIGNED_MAXORDER; k++) {
          rho[0][k] = rho6_lookup[idx][k];
          rho[1][k] = rho6_lookup[idy][k];
          rho[2][k] = rho6_lookup[idz][k];
          drho[0][k] = drho6_lookup[idx][k];
          drho[1][k] = drho6_lookup[idy][k];
          drho[2][k] = drho6_lookup[idz][k];
        }
      } else {
        #if defined(LMP_SIMD_COMPILER)
        #pragma simd
        #endif
        for (int k = nlower_6; k <= nupper_6; k++) {
          FFT_SCALAR r1,r2,r3,dr1,dr2,dr3;
          dr1 = dr2 = dr3 = ZEROF;

          r1 = rho_coeff_6[order_6-1][k];
          r2 = rho_coeff_6[order_6-1][k];
          r3 = rho_coeff_6[order_6-1][k];
          for (int l = order_6-2; l >= 0; l--) {
            r1 = rho_coeff_6[l][k] + r1 * dx;
            r2 = rho_coeff_6[l][k] + r2 * dy;
            r3 = rho_coeff_6[l][k] + r3 * dz;
            dr1 = drho_coeff_6[l][k] + dr1 * dx;
            dr2 = drho_coeff_6[l][k] + dr2 * dy;
            dr3 = drho_coeff_6[l][k] + dr3 * dz;
          }
          rho[0][k-nlower_6] = r1;
          rho[1][k-nlower_6] = r2;
          rho[2][k-nlower_6] = r3;
          drho[0][k-nlower_6] = dr1;
          drho[1][k-nlower_6] = dr2;
          drho[2][k-nlower_6] = dr3;
        }
      }
      _alignvar(FFT_SCALAR ekx[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
      _alignvar(FFT_SCALAR eky[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
      _alignvar(FFT_SCALAR ekz[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};

      particle_ekx[i] = particle_eky[i] = particle_ekz[i] = ZEROF;
      #if defined(LMP_SIMD_COMPILER)
      #pragma loop_count min(2), max(INTEL_P3M_ALIGNED_MAXORDER), avg(7)
      #endif
      for (int n = 0; n < order_6; n++) {
        int mz = n + nzsum;
        #if defined(LMP_SIMD_COMPILER)
        #pragma loop_count min(2), max(INTEL_P3M_ALIGNED_MAXORDER), avg(7)
        #endif
        for (int m = 0; m < order_6; m++) {
          int my = m + nysum;
          FFT_SCALAR ekx_p = rho[1][m] * rho[2][n];
          FFT_SCALAR eky_p = drho[1][m] * rho[2][n];
          FFT_SCALAR ekz_p = rho[1][m] * drho[2][n];
          #if defined(LMP_SIMD_COMPILER)
          #pragma loop_count min(2), max(INTEL_P3M_ALIGNED_MAXORDER), avg(7)
          #pragma simd
          #endif
          for (int l = 0; l < order; l++) {
            int mx = l + nxsum;
            ekx[l] += drho[0][l] * ekx_p * u_brick_g[mz][my][mx];
            eky[l] +=  rho[0][l] * eky_p * u_brick_g[mz][my][mx];
            ekz[l] +=  rho[0][l] * ekz_p * u_brick_g[mz][my][mx];
          }
        }
      }

      #if defined(LMP_SIMD_COMPILER)
      #pragma loop_count min(2), max(INTEL_P3M_ALIGNED_MAXORDER), avg(7)
      #endif
      for (int l = 0; l < order; l++){
        particle_ekx[i] += ekx[l];
        particle_eky[i] += eky[l];
        particle_ekz[i] += ekz[l];
      }
    }
    #if defined(LMP_SIMD_COMPILER)
    #pragma simd
    #endif
    for (int i = ifrom; i < ito; i++) {
      particle_ekx[i] *= hx_inv;
      particle_eky[i] *= hy_inv;
      particle_ekz[i] *= hz_inv;

      // convert E-field to force

      const int type = atom->type[i];
      const flt_t lj = B[type];
      const flt_t twoljsq = 2.*lj*lj;

      const flt_t s1 = x[i][0] * hx_inv;
      const flt_t s2 = x[i][1] * hy_inv;
      const flt_t s3 = x[i][2] * hz_inv;
      flt_t sf = fsf_coeff0 * sin(ftwo_pi * s1);
      sf += fsf_coeff1 * sin(ffour_pi * s1);
      sf *= twoljsq;
      f[i][0] += lj * particle_ekx[i] - sf;

      sf = fsf_coeff2 * sin(ftwo_pi * s2);
      sf += fsf_coeff3 * sin(ffour_pi * s2);
      sf *= twoljsq;
      f[i][1] += lj * particle_eky[i] - sf;

      sf = fsf_coeff4 * sin(ftwo_pi * s3);
      sf += fsf_coeff5 * sin(ffour_pi * s3);
      sf *= twoljsq;

      if (slabflag != 2) f[i][2] += lj * particle_ekz[i] -  sf;
    }
  }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get dispersion field & force on my particles
   for arithmetic mixing rule and ik scheme
------------------------------------------------------------------------- */

template<class flt_t, class acc_t, int use_table>
void PPPMDispIntel::fieldforce_a_ik(IntelBuffers<flt_t,acc_t> * /*buffers*/)
{

  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of dispersion field on particle

  int nlocal = atom->nlocal;
  int nthr = comm->nthreads;


  #if defined(_OPENMP)
  #pragma omp parallel default(none) \
    shared(nlocal, nthr) if(!_use_lrt)
  #endif
  {
    double **x = atom->x;
    double **f = atom->f;

    const flt_t lo0 = boxlo[0];
    const flt_t lo1 = boxlo[1];
    const flt_t lo2 = boxlo[2];
    const flt_t xi = delxinv_6;
    const flt_t yi = delyinv_6;
    const flt_t zi = delzinv_6;
    const flt_t fshiftone = shiftone_6;

    int ifrom, ito, tid;
    IP_PRE_omp_range_id(ifrom, ito, tid, nlocal, nthr);

    _alignvar(flt_t rho0[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
    _alignvar(flt_t rho1[INTEL_P3M_ALIGNED_MAXORDER] , 64)= {0};
    _alignvar(flt_t rho2[INTEL_P3M_ALIGNED_MAXORDER] , 64)= {0};

    for (int i = ifrom; i < ito; i++) {
      int nx = part2grid_6[i][0];
      int ny = part2grid_6[i][1];
      int nz = part2grid_6[i][2];

      int nxsum = nx + nlower_6;
      int nysum = ny + nlower_6;
      int nzsum = nz + nlower_6;

      FFT_SCALAR dx = nx+fshiftone - (x[i][0]-lo0)*xi;
      FFT_SCALAR dy = ny+fshiftone - (x[i][1]-lo1)*yi;
      FFT_SCALAR dz = nz+fshiftone - (x[i][2]-lo2)*zi;

      if (use_table) {
        dx = dx*half_rho_scale + half_rho_scale_plus;
        int idx = dx;
        dy = dy*half_rho_scale + half_rho_scale_plus;
        int idy = dy;
        dz = dz*half_rho_scale + half_rho_scale_plus;
        int idz = dz;
        #if defined(LMP_SIMD_COMPILER)
        #pragma simd
        #endif
        for (int k = 0; k < INTEL_P3M_ALIGNED_MAXORDER; k++) {
          rho0[k] = rho6_lookup[idx][k];
          rho1[k] = rho6_lookup[idy][k];
          rho2[k] = rho6_lookup[idz][k];
        }
      } else {
        #if defined(LMP_SIMD_COMPILER)
        #pragma simd
        #endif
        for (int k = nlower_6; k <= nupper_6; k++) {
          FFT_SCALAR r1 = rho_coeff_6[order_6-1][k];
          FFT_SCALAR r2 = rho_coeff_6[order_6-1][k];
          FFT_SCALAR r3 = rho_coeff_6[order_6-1][k];
          for (int l = order_6-2; l >= 0; l--) {
            r1 = rho_coeff_6[l][k] + r1*dx;
            r2 = rho_coeff_6[l][k] + r2*dy;
            r3 = rho_coeff_6[l][k] + r3*dz;
          }

          rho0[k-nlower_6] = r1;
          rho1[k-nlower_6] = r2;
          rho2[k-nlower_6] = r3;
        }
      }

      _alignvar(FFT_SCALAR ekx0_arr[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
      _alignvar(FFT_SCALAR eky0_arr[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
      _alignvar(FFT_SCALAR ekz0_arr[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
      _alignvar(FFT_SCALAR ekx1_arr[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
      _alignvar(FFT_SCALAR eky1_arr[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
      _alignvar(FFT_SCALAR ekz1_arr[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
      _alignvar(FFT_SCALAR ekx2_arr[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
      _alignvar(FFT_SCALAR eky2_arr[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
      _alignvar(FFT_SCALAR ekz2_arr[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
      _alignvar(FFT_SCALAR ekx3_arr[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
      _alignvar(FFT_SCALAR eky3_arr[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
      _alignvar(FFT_SCALAR ekz3_arr[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
      _alignvar(FFT_SCALAR ekx4_arr[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
      _alignvar(FFT_SCALAR eky4_arr[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
      _alignvar(FFT_SCALAR ekz4_arr[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
      _alignvar(FFT_SCALAR ekx5_arr[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
      _alignvar(FFT_SCALAR eky5_arr[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
      _alignvar(FFT_SCALAR ekz5_arr[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
      _alignvar(FFT_SCALAR ekx6_arr[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
      _alignvar(FFT_SCALAR eky6_arr[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
      _alignvar(FFT_SCALAR ekz6_arr[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};


      #if defined(LMP_SIMD_COMPILER)
      #pragma loop_count min(2), max(INTEL_P3M_ALIGNED_MAXORDER), avg(7)
      #endif
      for (int n = 0; n < order_6; n++) {
        int mz = n+nzsum;
        FFT_SCALAR z0 = rho2[n];
        #if defined(LMP_SIMD_COMPILER)
        #pragma loop_count min(2), max(INTEL_P3M_ALIGNED_MAXORDER), avg(7)
        #endif
        for (int m = 0; m < order_6; m++) {
          int my = m+nysum;
          FFT_SCALAR y0 = z0*rho1[m];
          #if defined(LMP_SIMD_COMPILER)
          #pragma loop_count min(2), max(INTEL_P3M_ALIGNED_MAXORDER), avg(7)
          #pragma simd
          #endif
          for (int l = 0; l < order; l++) {
            int mx = l+nxsum;
            FFT_SCALAR x0 = y0*rho0[l];
              ekx0_arr[l] -= x0*vdx_brick_a0[mz][my][mx];
              eky0_arr[l] -= x0*vdy_brick_a0[mz][my][mx];
              ekz0_arr[l] -= x0*vdz_brick_a0[mz][my][mx];
              ekx1_arr[l] -= x0*vdx_brick_a1[mz][my][mx];
              eky1_arr[l] -= x0*vdy_brick_a1[mz][my][mx];
              ekz1_arr[l] -= x0*vdz_brick_a1[mz][my][mx];
              ekx2_arr[l] -= x0*vdx_brick_a2[mz][my][mx];
              eky2_arr[l] -= x0*vdy_brick_a2[mz][my][mx];
              ekz2_arr[l] -= x0*vdz_brick_a2[mz][my][mx];
              ekx3_arr[l] -= x0*vdx_brick_a3[mz][my][mx];
              eky3_arr[l] -= x0*vdy_brick_a3[mz][my][mx];
              ekz3_arr[l] -= x0*vdz_brick_a3[mz][my][mx];
              ekx4_arr[l] -= x0*vdx_brick_a4[mz][my][mx];
              eky4_arr[l] -= x0*vdy_brick_a4[mz][my][mx];
              ekz4_arr[l] -= x0*vdz_brick_a4[mz][my][mx];
              ekx5_arr[l] -= x0*vdx_brick_a5[mz][my][mx];
              eky5_arr[l] -= x0*vdy_brick_a5[mz][my][mx];
              ekz5_arr[l] -= x0*vdz_brick_a5[mz][my][mx];
              ekx6_arr[l] -= x0*vdx_brick_a6[mz][my][mx];
              eky6_arr[l] -= x0*vdy_brick_a6[mz][my][mx];
              ekz6_arr[l] -= x0*vdz_brick_a6[mz][my][mx];
          }
        }
      }

      FFT_SCALAR ekx0, eky0, ekz0, ekx1, eky1, ekz1, ekx2, eky2, ekz2;
      FFT_SCALAR ekx3, eky3, ekz3, ekx4, eky4, ekz4, ekx5, eky5, ekz5;
      FFT_SCALAR ekx6, eky6, ekz6;
      ekx0 = eky0 = ekz0 = ZEROF;
      ekx1 = eky1 = ekz1 = ZEROF;
      ekx2 = eky2 = ekz2 = ZEROF;
      ekx3 = eky3 = ekz3 = ZEROF;
      ekx4 = eky4 = ekz4 = ZEROF;
      ekx5 = eky5 = ekz5 = ZEROF;
      ekx6 = eky6 = ekz6 = ZEROF;

      for (int l = 0; l < order; l++) {
        ekx0 += ekx0_arr[l];
        eky0 += eky0_arr[l];
        ekz0 += ekz0_arr[l];
        ekx1 += ekx1_arr[l];
        eky1 += eky1_arr[l];
        ekz1 += ekz1_arr[l];
        ekx2 += ekx2_arr[l];
        eky2 += eky2_arr[l];
        ekz2 += ekz2_arr[l];
        ekx3 += ekx3_arr[l];
        eky3 += eky3_arr[l];
        ekz3 += ekz3_arr[l];
        ekx4 += ekx4_arr[l];
        eky4 += eky4_arr[l];
        ekz4 += ekz4_arr[l];
        ekx5 += ekx5_arr[l];
        eky5 += eky5_arr[l];
        ekz5 += ekz5_arr[l];
        ekx6 += ekx6_arr[l];
        eky6 += eky6_arr[l];
        ekz6 += ekz6_arr[l];
      }

      // convert D-field to force

      const int type = atom->type[i];
      const FFT_SCALAR lj0 = B[7*type+6];
      const FFT_SCALAR lj1 = B[7*type+5];
      const FFT_SCALAR lj2 = B[7*type+4];
      const FFT_SCALAR lj3 = B[7*type+3];
      const FFT_SCALAR lj4 = B[7*type+2];
      const FFT_SCALAR lj5 = B[7*type+1];
      const FFT_SCALAR lj6 = B[7*type];

      f[i][0] += lj0*ekx0 + lj1*ekx1 + lj2*ekx2 + lj3*ekx3 +
        lj4*ekx4 + lj5*ekx5 + lj6*ekx6;
      f[i][1] += lj0*eky0 + lj1*eky1 + lj2*eky2 + lj3*eky3 +
        lj4*eky4 + lj5*eky5 + lj6*eky6;
      if (slabflag != 2) f[i][2] += lj0*ekz0 + lj1*ekz1 + lj2*ekz2 +
                           lj3*ekz3 + lj4*ekz4 + lj5*ekz5 + lj6*ekz6;
    }
  }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get dispersion field & force on my particles
   for arithmetic mixing rule for the ad scheme
------------------------------------------------------------------------- */

template<class flt_t, class acc_t, int use_table>
void PPPMDispIntel::fieldforce_a_ad(IntelBuffers<flt_t,acc_t> * /*buffers*/)
{

  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of dispersion field on particle

  int nlocal = atom->nlocal;
  int nthr = comm->nthreads;

  FFT_SCALAR * _noalias const particle_ekx0 = this->particle_ekx0;
  FFT_SCALAR * _noalias const particle_eky0 = this->particle_eky0;
  FFT_SCALAR * _noalias const particle_ekz0 = this->particle_ekz0;
  FFT_SCALAR * _noalias const particle_ekx1 = this->particle_ekx1;
  FFT_SCALAR * _noalias const particle_eky1 = this->particle_eky1;
  FFT_SCALAR * _noalias const particle_ekz1 = this->particle_ekz1;
  FFT_SCALAR * _noalias const particle_ekx2 = this->particle_ekx2;
  FFT_SCALAR * _noalias const particle_eky2 = this->particle_eky2;
  FFT_SCALAR * _noalias const particle_ekz2 = this->particle_ekz2;
  FFT_SCALAR * _noalias const particle_ekx3 = this->particle_ekx3;
  FFT_SCALAR * _noalias const particle_eky3 = this->particle_eky3;
  FFT_SCALAR * _noalias const particle_ekz3 = this->particle_ekz3;
  FFT_SCALAR * _noalias const particle_ekx4 = this->particle_ekx4;
  FFT_SCALAR * _noalias const particle_eky4 = this->particle_eky4;
  FFT_SCALAR * _noalias const particle_ekz4 = this->particle_ekz4;
  FFT_SCALAR * _noalias const particle_ekx5 = this->particle_ekx5;
  FFT_SCALAR * _noalias const particle_eky5 = this->particle_eky5;
  FFT_SCALAR * _noalias const particle_ekz5 = this->particle_ekz5;
  FFT_SCALAR * _noalias const particle_ekx6 = this->particle_ekx6;
  FFT_SCALAR * _noalias const particle_eky6 = this->particle_eky6;
  FFT_SCALAR * _noalias const particle_ekz6 = this->particle_ekz6;

  #if defined(_OPENMP)
  #pragma omp parallel default(none) \
    shared(nlocal, nthr) if(!_use_lrt)
  #endif
  {

    double *prd;
    if (triclinic == 0) prd = domain->prd;
    else prd = domain->prd_lamda;

    double **x = atom->x;
    double **f = atom->f;
    const flt_t ftwo_pi = MY_PI * 2.0;
    const flt_t ffour_pi = MY_PI * 4.0;

    const flt_t lo0 = boxlo[0];
    const flt_t lo1 = boxlo[1];
    const flt_t lo2 = boxlo[2];
    const flt_t xi = delxinv_6;
    const flt_t yi = delyinv_6;
    const flt_t zi = delzinv_6;
    const flt_t fshiftone = shiftone_6;

    const double xprd = prd[0];
    const double yprd = prd[1];
    const double zprd = prd[2]*slab_volfactor;

    const flt_t hx_inv = nx_pppm_6/xprd;
    const flt_t hy_inv = ny_pppm_6/yprd;
    const flt_t hz_inv = nz_pppm_6/zprd;

    const flt_t fsf_coeff0 = sf_coeff_6[0];
    const flt_t fsf_coeff1 = sf_coeff_6[1];
    const flt_t fsf_coeff2 = sf_coeff_6[2];
    const flt_t fsf_coeff3 = sf_coeff_6[3];
    const flt_t fsf_coeff4 = sf_coeff_6[4];
    const flt_t fsf_coeff5 = sf_coeff_6[5];

    int ifrom, ito, tid;
    IP_PRE_omp_range_id(ifrom, ito, tid, nlocal, nthr);

    _alignvar(flt_t rho[3][INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
    _alignvar(flt_t drho[3][INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};

    for (int i = ifrom; i < ito; i++) {
      int nx = part2grid_6[i][0];
      int ny = part2grid_6[i][1];
      int nz = part2grid_6[i][2];
      FFT_SCALAR dx = nx+fshiftone - (x[i][0]-lo0)*xi;
      FFT_SCALAR dy = ny+fshiftone - (x[i][1]-lo1)*yi;
      FFT_SCALAR dz = nz+fshiftone - (x[i][2]-lo2)*zi;

      int nxsum = nx + nlower_6;
      int nysum = ny + nlower_6;
      int nzsum = nz + nlower_6;

      if (use_table) {
        dx = dx*half_rho_scale + half_rho_scale_plus;
        int idx = dx;
        dy = dy*half_rho_scale + half_rho_scale_plus;
        int idy = dy;
        dz = dz*half_rho_scale + half_rho_scale_plus;
        int idz = dz;

        #if defined(LMP_SIMD_COMPILER)
        #pragma simd
        #endif
        for(int k = 0; k < INTEL_P3M_ALIGNED_MAXORDER; k++) {
          rho[0][k] = rho6_lookup[idx][k];
          rho[1][k] = rho6_lookup[idy][k];
          rho[2][k] = rho6_lookup[idz][k];
          drho[0][k] = drho6_lookup[idx][k];
          drho[1][k] = drho6_lookup[idy][k];
          drho[2][k] = drho6_lookup[idz][k];
        }
      } else {
        #if defined(LMP_SIMD_COMPILER)
        #pragma simd
        #endif
        for (int k = nlower_6; k <= nupper_6; k++) {
          FFT_SCALAR r1,r2,r3,dr1,dr2,dr3;
          dr1 = dr2 = dr3 = ZEROF;

          r1 = rho_coeff_6[order_6-1][k];
          r2 = rho_coeff_6[order_6-1][k];
          r3 = rho_coeff_6[order_6-1][k];
          for (int l = order_6-2; l >= 0; l--) {
            r1 = rho_coeff_6[l][k] + r1 * dx;
            r2 = rho_coeff_6[l][k] + r2 * dy;
            r3 = rho_coeff_6[l][k] + r3 * dz;
            dr1 = drho_coeff_6[l][k] + dr1 * dx;
            dr2 = drho_coeff_6[l][k] + dr2 * dy;
            dr3 = drho_coeff_6[l][k] + dr3 * dz;
          }
          rho[0][k-nlower_6] = r1;
          rho[1][k-nlower_6] = r2;
          rho[2][k-nlower_6] = r3;
          drho[0][k-nlower_6] = dr1;
          drho[1][k-nlower_6] = dr2;
          drho[2][k-nlower_6] = dr3;
        }
      }
      _alignvar(FFT_SCALAR ekx0[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
      _alignvar(FFT_SCALAR eky0[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
      _alignvar(FFT_SCALAR ekz0[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
      _alignvar(FFT_SCALAR ekx1[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
      _alignvar(FFT_SCALAR eky1[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
      _alignvar(FFT_SCALAR ekz1[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
      _alignvar(FFT_SCALAR ekx2[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
      _alignvar(FFT_SCALAR eky2[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
      _alignvar(FFT_SCALAR ekz2[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
      _alignvar(FFT_SCALAR ekx3[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
      _alignvar(FFT_SCALAR eky3[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
      _alignvar(FFT_SCALAR ekz3[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
      _alignvar(FFT_SCALAR ekx4[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
      _alignvar(FFT_SCALAR eky4[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
      _alignvar(FFT_SCALAR ekz4[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
      _alignvar(FFT_SCALAR ekx5[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
      _alignvar(FFT_SCALAR eky5[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
      _alignvar(FFT_SCALAR ekz5[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
      _alignvar(FFT_SCALAR ekx6[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
      _alignvar(FFT_SCALAR eky6[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
      _alignvar(FFT_SCALAR ekz6[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};

      particle_ekx0[i] = particle_eky0[i] = particle_ekz0[i] = ZEROF;
      particle_ekx1[i] = particle_eky1[i] = particle_ekz1[i] = ZEROF;
      particle_ekx2[i] = particle_eky2[i] = particle_ekz2[i] = ZEROF;
      particle_ekx3[i] = particle_eky3[i] = particle_ekz3[i] = ZEROF;
      particle_ekx4[i] = particle_eky4[i] = particle_ekz4[i] = ZEROF;
      particle_ekx5[i] = particle_eky5[i] = particle_ekz5[i] = ZEROF;
      particle_ekx6[i] = particle_eky6[i] = particle_ekz6[i] = ZEROF;

      #if defined(LMP_SIMD_COMPILER)
      #pragma loop_count min(2), max(INTEL_P3M_ALIGNED_MAXORDER), avg(7)
      #endif
      for (int n = 0; n < order_6; n++) {
        int mz = n + nzsum;
        #if defined(LMP_SIMD_COMPILER)
        #pragma loop_count min(2), max(INTEL_P3M_ALIGNED_MAXORDER), avg(7)
        #endif
        for (int m = 0; m < order_6; m++) {
          int my = m + nysum;
          FFT_SCALAR ekx_p = rho[1][m] * rho[2][n];
          FFT_SCALAR eky_p = drho[1][m] * rho[2][n];
          FFT_SCALAR ekz_p = rho[1][m] * drho[2][n];
          #if defined(LMP_SIMD_COMPILER)
          #pragma loop_count min(2), max(INTEL_P3M_ALIGNED_MAXORDER), avg(7)
          #pragma simd
          #endif
          for (int l = 0; l < order; l++) {
            int mx = l + nxsum;
            FFT_SCALAR x0 = drho[0][l] * ekx_p;
            FFT_SCALAR y0 = rho[0][l] * eky_p;
            FFT_SCALAR z0 = rho[0][l] * ekz_p;

            ekx0[l] +=  x0 * u_brick_a0[mz][my][mx];
            eky0[l] +=  y0 * u_brick_a0[mz][my][mx];
            ekz0[l] +=  z0 * u_brick_a0[mz][my][mx];
            ekx1[l] +=  x0 * u_brick_a1[mz][my][mx];
            eky1[l] +=  y0 * u_brick_a1[mz][my][mx];
            ekz1[l] +=  z0 * u_brick_a1[mz][my][mx];
            ekx2[l] +=  x0 * u_brick_a2[mz][my][mx];
            eky2[l] +=  y0 * u_brick_a2[mz][my][mx];
            ekz2[l] +=  z0 * u_brick_a2[mz][my][mx];
            ekx3[l] +=  x0 * u_brick_a3[mz][my][mx];
            eky3[l] +=  y0 * u_brick_a3[mz][my][mx];
            ekz3[l] +=  z0 * u_brick_a3[mz][my][mx];
            ekx4[l] +=  x0 * u_brick_a4[mz][my][mx];
            eky4[l] +=  y0 * u_brick_a4[mz][my][mx];
            ekz4[l] +=  z0 * u_brick_a4[mz][my][mx];
            ekx5[l] +=  x0 * u_brick_a5[mz][my][mx];
            eky5[l] +=  y0 * u_brick_a5[mz][my][mx];
            ekz5[l] +=  z0 * u_brick_a5[mz][my][mx];
            ekx6[l] +=  x0 * u_brick_a6[mz][my][mx];
            eky6[l] +=  y0 * u_brick_a6[mz][my][mx];
            ekz6[l] +=  z0 * u_brick_a6[mz][my][mx];
          }
        }
      }

      #if defined(LMP_SIMD_COMPILER)
      #pragma loop_count min(2), max(INTEL_P3M_ALIGNED_MAXORDER), avg(7)
      #endif
      for (int l = 0; l < order; l++){
        particle_ekx0[i] += ekx0[l];
        particle_eky0[i] += eky0[l];
        particle_ekz0[i] += ekz0[l];
        particle_ekx1[i] += ekx1[l];
        particle_eky1[i] += eky1[l];
        particle_ekz1[i] += ekz1[l];
        particle_ekx2[i] += ekx2[l];
        particle_eky2[i] += eky2[l];
        particle_ekz2[i] += ekz2[l];
        particle_ekx3[i] += ekx3[l];
        particle_eky3[i] += eky3[l];
        particle_ekz3[i] += ekz3[l];
        particle_ekx4[i] += ekx4[l];
        particle_eky4[i] += eky4[l];
        particle_ekz4[i] += ekz4[l];
        particle_ekx5[i] += ekx5[l];
        particle_eky5[i] += eky5[l];
        particle_ekz5[i] += ekz5[l];
        particle_ekx6[i] += ekx6[l];
        particle_eky6[i] += eky6[l];
        particle_ekz6[i] += ekz6[l];
      }
    }
    #if defined(LMP_SIMD_COMPILER)
    #pragma simd
    #endif
    for (int i = ifrom; i < ito; i++) {
      particle_ekx0[i] *= hx_inv;
      particle_eky0[i] *= hy_inv;
      particle_ekz0[i] *= hz_inv;
      particle_ekx1[i] *= hx_inv;
      particle_eky1[i] *= hy_inv;
      particle_ekz1[i] *= hz_inv;
      particle_ekx2[i] *= hx_inv;
      particle_eky2[i] *= hy_inv;
      particle_ekz2[i] *= hz_inv;
      particle_ekx3[i] *= hx_inv;
      particle_eky3[i] *= hy_inv;
      particle_ekz3[i] *= hz_inv;
      particle_ekx4[i] *= hx_inv;
      particle_eky4[i] *= hy_inv;
      particle_ekz4[i] *= hz_inv;
      particle_ekx5[i] *= hx_inv;
      particle_eky5[i] *= hy_inv;
      particle_ekz5[i] *= hz_inv;
      particle_ekx6[i] *= hx_inv;
      particle_eky6[i] *= hy_inv;
      particle_ekz6[i] *= hz_inv;

      // convert D-field to force

      const int type = atom->type[i];
      const FFT_SCALAR lj0 = B[7*type+6];
      const FFT_SCALAR lj1 = B[7*type+5];
      const FFT_SCALAR lj2 = B[7*type+4];
      const FFT_SCALAR lj3 = B[7*type+3];
      const FFT_SCALAR lj4 = B[7*type+2];
      const FFT_SCALAR lj5 = B[7*type+1];
      const FFT_SCALAR lj6 = B[7*type];

      const flt_t s1 = x[i][0] * hx_inv;
      const flt_t s2 = x[i][1] * hy_inv;
      const flt_t s3 = x[i][2] * hz_inv;
      flt_t sf = fsf_coeff0 * sin(ftwo_pi * s1);
      sf += fsf_coeff1 * sin(ffour_pi * s1);
      sf *= 4*lj0*lj6 + 4*lj1*lj5 + 4*lj2*lj4 + 2*lj3*lj3;
      f[i][0] += lj0*particle_ekx0[i] + lj1*particle_ekx1[i] +
        lj2*particle_ekx2[i] + lj3*particle_ekx3[i] + lj4*particle_ekx4[i] +
        lj5*particle_ekx5[i] + lj6*particle_ekx6[i] - sf;

      sf = fsf_coeff2 * sin(ftwo_pi * s2);
      sf += fsf_coeff3 * sin(ffour_pi * s2);
      sf *= 4*lj0*lj6 + 4*lj1*lj5 + 4*lj2*lj4 + 2*lj3*lj3;
      f[i][1] += lj0*particle_eky0[i] + lj1*particle_eky1[i] +
        lj2*particle_eky2[i] + lj3*particle_eky3[i] + lj4*particle_eky4[i] +
        lj5*particle_eky5[i] + lj6*particle_eky6[i] - sf;

      sf = fsf_coeff4 * sin(ftwo_pi * s3);
      sf += fsf_coeff5 * sin(ffour_pi * s3);
      sf *= 4*lj0*lj6 + 4*lj1*lj5 + 4*lj2*lj4 + 2*lj3*lj3;
      if (slabflag != 2)
      f[i][2] += lj0*particle_ekz0[i] + lj1*particle_ekz1[i] +
        lj2*particle_ekz2[i] + lj3*particle_ekz3[i] + lj4*particle_ekz4[i] +
        lj5*particle_ekz5[i] + lj6*particle_ekz6[i] - sf;
    }
  }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get dispersion field & force on my particles
   for no mixing rule and ik scheme
------------------------------------------------------------------------- */

template<class flt_t, class acc_t, int use_table>
void PPPMDispIntel::fieldforce_none_ik(IntelBuffers<flt_t,acc_t> * /*buffers*/)
{

  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of dispersion field on particle

  int nlocal = atom->nlocal;
  int nthr = comm->nthreads;


  #if defined(_OPENMP)
  #pragma omp parallel default(none) \
    shared(nlocal, nthr) if(!_use_lrt)
  #endif
  {

    double lj;
    int type;
    double **x = atom->x;
    double **f = atom->f;

    const flt_t lo0 = boxlo[0];
    const flt_t lo1 = boxlo[1];
    const flt_t lo2 = boxlo[2];
    const flt_t xi = delxinv_6;
    const flt_t yi = delyinv_6;
    const flt_t zi = delzinv_6;
    const flt_t fshiftone = shiftone_6;

    int ifrom, ito, tid;
    IP_PRE_omp_range_id(ifrom, ito, tid, nlocal, nthr);

    _alignvar(flt_t rho0[INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
    _alignvar(flt_t rho1[INTEL_P3M_ALIGNED_MAXORDER] , 64)= {0};
    _alignvar(flt_t rho2[INTEL_P3M_ALIGNED_MAXORDER] , 64)= {0};

    for (int i = ifrom; i < ito; i++) {
      int nx = part2grid_6[i][0];
      int ny = part2grid_6[i][1];
      int nz = part2grid_6[i][2];

      int nxsum = nx + nlower_6;
      int nysum = ny + nlower_6;
      int nzsum = nz + nlower_6;

      FFT_SCALAR dx = nx+fshiftone - (x[i][0]-lo0)*xi;
      FFT_SCALAR dy = ny+fshiftone - (x[i][1]-lo1)*yi;
      FFT_SCALAR dz = nz+fshiftone - (x[i][2]-lo2)*zi;

      if (use_table) {
        dx = dx*half_rho_scale + half_rho_scale_plus;
        int idx = dx;
        dy = dy*half_rho_scale + half_rho_scale_plus;
        int idy = dy;
        dz = dz*half_rho_scale + half_rho_scale_plus;
        int idz = dz;
        #if defined(LMP_SIMD_COMPILER)
        #pragma simd
        #endif
        for (int k = 0; k < INTEL_P3M_ALIGNED_MAXORDER; k++) {
          rho0[k] = rho6_lookup[idx][k];
          rho1[k] = rho6_lookup[idy][k];
          rho2[k] = rho6_lookup[idz][k];
        }
      } else {
        #if defined(LMP_SIMD_COMPILER)
        #pragma simd
        #endif
        for (int k = nlower_6; k <= nupper_6; k++) {
          FFT_SCALAR r1 = rho_coeff_6[order_6-1][k];
          FFT_SCALAR r2 = rho_coeff_6[order_6-1][k];
          FFT_SCALAR r3 = rho_coeff_6[order_6-1][k];
          for (int l = order_6-2; l >= 0; l--) {
            r1 = rho_coeff_6[l][k] + r1*dx;
            r2 = rho_coeff_6[l][k] + r2*dy;
            r3 = rho_coeff_6[l][k] + r3*dz;
          }

          rho0[k-nlower_6] = r1;
          rho1[k-nlower_6] = r2;
          rho2[k-nlower_6] = r3;
        }
      }


      _alignvar(FFT_SCALAR ekx_arr[nsplit*INTEL_P3M_ALIGNED_MAXORDER],64);
      _alignvar(FFT_SCALAR eky_arr[nsplit*INTEL_P3M_ALIGNED_MAXORDER],64);
      _alignvar(FFT_SCALAR ekz_arr[nsplit*INTEL_P3M_ALIGNED_MAXORDER],64);

      for (int k = 0; k < nsplit*INTEL_P3M_ALIGNED_MAXORDER; k++) {
        ekx_arr[k] = eky_arr[k] = ekz_arr[k] = ZEROF;
      }

      for (int k = 0; k < nsplit; k++) {
        #if defined(LMP_SIMD_COMPILER)
        #pragma loop_count min(2), max(INTEL_P3M_ALIGNED_MAXORDER), avg(7)
        #endif
        for (int n = 0; n < order_6; n++) {
          int mz = n+nzsum;
          FFT_SCALAR z0 = rho2[n];
          #if defined(LMP_SIMD_COMPILER)
          #pragma loop_count min(2), max(INTEL_P3M_ALIGNED_MAXORDER), avg(7)
          #endif
          for (int m = 0; m < order_6; m++) {
            int my = m+nysum;
            FFT_SCALAR y0 = z0*rho1[m];
            #if defined(LMP_SIMD_COMPILER)
            #pragma loop_count min(2), max(INTEL_P3M_ALIGNED_MAXORDER), avg(7)
            #pragma simd
            #endif
            for (int l = 0; l < order; l++) {
              int mx = l+nxsum;
              FFT_SCALAR x0 = y0*rho0[l];
              ekx_arr[k*INTEL_P3M_ALIGNED_MAXORDER + l] -=
                x0*vdx_brick_none[k][mz][my][mx];
              eky_arr[k*INTEL_P3M_ALIGNED_MAXORDER + l] -=
                x0*vdy_brick_none[k][mz][my][mx];
              ekz_arr[k*INTEL_P3M_ALIGNED_MAXORDER + l] -=
                x0*vdz_brick_none[k][mz][my][mx];
            }
          }
        }
      }

      _alignvar(FFT_SCALAR ekx[nsplit], 64);
      _alignvar(FFT_SCALAR eky[nsplit], 64);
      _alignvar(FFT_SCALAR ekz[nsplit], 64);
      for (int k = 0; k < nsplit; k++) {
        ekx[k] = eky[k] = ekz[k] = ZEROF;
      }

      for (int l = 0; l < order; l++) {
        for (int k = 0; k < nsplit; k++) {
          ekx[k] += ekx_arr[k*INTEL_P3M_ALIGNED_MAXORDER + l];
          eky[k] += eky_arr[k*INTEL_P3M_ALIGNED_MAXORDER + l];
          ekz[k] += ekz_arr[k*INTEL_P3M_ALIGNED_MAXORDER + l];
        }
      }

      // convert E-field to force

      type = atom->type[i];
      for (int k = 0; k < nsplit; k++) {
        lj = B[nsplit*type + k];
        f[i][0] += lj*ekx[k];
        f[i][1] += lj*eky[k];
        if (slabflag != 2) f[i][2] += lj*ekz[k];
      }
    }
  }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get dispersion field & force on my particles
   for no mixing rule for the ad scheme
------------------------------------------------------------------------- */

template<class flt_t, class acc_t, int use_table>
void PPPMDispIntel::fieldforce_none_ad(IntelBuffers<flt_t,acc_t> * /*buffers*/)
{
  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of dispersion field on particle

  int nlocal = atom->nlocal;
  int nthr = comm->nthreads;

   #if defined(_OPENMP)
   #pragma omp parallel default(none)           \
     shared(nlocal, nthr) if(!_use_lrt)
   #endif
  {

    double *prd;
    if (triclinic == 0) prd = domain->prd;
    else prd = domain->prd_lamda;

    double **x = atom->x;
    double **f = atom->f;
    const flt_t ftwo_pi = MY_PI * 2.0;
    const flt_t ffour_pi = MY_PI * 4.0;

    const flt_t lo0 = boxlo[0];
    const flt_t lo1 = boxlo[1];
    const flt_t lo2 = boxlo[2];
    const flt_t xi = delxinv_6;
    const flt_t yi = delyinv_6;
    const flt_t zi = delzinv_6;
    const flt_t fshiftone = shiftone_6;

    const double xprd = prd[0];
    const double yprd = prd[1];
    const double zprd = prd[2]*slab_volfactor;

    const flt_t hx_inv = nx_pppm_6/xprd;
    const flt_t hy_inv = ny_pppm_6/yprd;
    const flt_t hz_inv = nz_pppm_6/zprd;

    const flt_t fsf_coeff0 = sf_coeff_6[0];
    const flt_t fsf_coeff1 = sf_coeff_6[1];
    const flt_t fsf_coeff2 = sf_coeff_6[2];
    const flt_t fsf_coeff3 = sf_coeff_6[3];
    const flt_t fsf_coeff4 = sf_coeff_6[4];
    const flt_t fsf_coeff5 = sf_coeff_6[5];

    int ifrom, ito, tid;
    IP_PRE_omp_range_id(ifrom, ito, tid, nlocal, nthr);

    _alignvar(flt_t rho[3][INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
    _alignvar(flt_t drho[3][INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};

    for (int i = ifrom; i < ito; i++) {
      int nx = part2grid_6[i][0];
      int ny = part2grid_6[i][1];
      int nz = part2grid_6[i][2];
      FFT_SCALAR dx = nx+fshiftone - (x[i][0]-lo0)*xi;
      FFT_SCALAR dy = ny+fshiftone - (x[i][1]-lo1)*yi;
      FFT_SCALAR dz = nz+fshiftone - (x[i][2]-lo2)*zi;

      int nxsum = nx + nlower_6;
      int nysum = ny + nlower_6;
      int nzsum = nz + nlower_6;

      if (use_table) {
        dx = dx*half_rho_scale + half_rho_scale_plus;
        int idx = dx;
        dy = dy*half_rho_scale + half_rho_scale_plus;
        int idy = dy;
        dz = dz*half_rho_scale + half_rho_scale_plus;
        int idz = dz;

        #if defined(LMP_SIMD_COMPILER)
        #pragma simd
        #endif
        for(int k = 0; k < INTEL_P3M_ALIGNED_MAXORDER; k++) {
          rho[0][k] = rho6_lookup[idx][k];
          rho[1][k] = rho6_lookup[idy][k];
          rho[2][k] = rho6_lookup[idz][k];
          drho[0][k] = drho6_lookup[idx][k];
          drho[1][k] = drho6_lookup[idy][k];
          drho[2][k] = drho6_lookup[idz][k];
        }
      } else {
        #if defined(LMP_SIMD_COMPILER)
        #pragma simd
        #endif
        for (int k = nlower_6; k <= nupper_6; k++) {
          FFT_SCALAR r1,r2,r3,dr1,dr2,dr3;
          dr1 = dr2 = dr3 = ZEROF;

          r1 = rho_coeff_6[order_6-1][k];
          r2 = rho_coeff_6[order_6-1][k];
          r3 = rho_coeff_6[order_6-1][k];
          for (int l = order_6-2; l >= 0; l--) {
            r1 = rho_coeff_6[l][k] + r1 * dx;
            r2 = rho_coeff_6[l][k] + r2 * dy;
            r3 = rho_coeff_6[l][k] + r3 * dz;
            dr1 = drho_coeff_6[l][k] + dr1 * dx;
            dr2 = drho_coeff_6[l][k] + dr2 * dy;
            dr3 = drho_coeff_6[l][k] + dr3 * dz;
          }
          rho[0][k-nlower_6] = r1;
          rho[1][k-nlower_6] = r2;
          rho[2][k-nlower_6] = r3;
          drho[0][k-nlower_6] = dr1;
          drho[1][k-nlower_6] = dr2;
          drho[2][k-nlower_6] = dr3;
        }
      }
      _alignvar(FFT_SCALAR ekx[nsplit*INTEL_P3M_ALIGNED_MAXORDER], 64);
      _alignvar(FFT_SCALAR eky[nsplit*INTEL_P3M_ALIGNED_MAXORDER], 64);
      _alignvar(FFT_SCALAR ekz[nsplit*INTEL_P3M_ALIGNED_MAXORDER], 64);

      for (int k = 0; k < nsplit*INTEL_P3M_ALIGNED_MAXORDER; k++) {
        ekx[k]=eky[k]=ekz[k]=ZEROF;
      }

      for (int k = 0; k < nsplit; k++) {
        particle_ekx[i] = particle_eky[i] = particle_ekz[i] = ZEROF;
        #if defined(LMP_SIMD_COMPILER)
        #pragma loop_count min(2), max(INTEL_P3M_ALIGNED_MAXORDER), avg(7)
        #endif
        for (int n = 0; n < order_6; n++) {
          int mz = n + nzsum;
          #if defined(LMP_SIMD_COMPILER)
          #pragma loop_count min(2), max(INTEL_P3M_ALIGNED_MAXORDER), avg(7)
          #endif
          for (int m = 0; m < order_6; m++) {
            int my = m + nysum;
            FFT_SCALAR ekx_p = rho[1][m] * rho[2][n];
            FFT_SCALAR eky_p = drho[1][m] * rho[2][n];
            FFT_SCALAR ekz_p = rho[1][m] * drho[2][n];
            #if defined(LMP_SIMD_COMPILER)
            #pragma loop_count min(2), max(INTEL_P3M_ALIGNED_MAXORDER), avg(7)
            #pragma simd
            #endif
            for (int l = 0; l < order; l++) {
              int mx = l + nxsum;
              ekx[k*INTEL_P3M_ALIGNED_MAXORDER+l] += drho[0][l] * ekx_p *
                u_brick_none[k][mz][my][mx];
              eky[k*INTEL_P3M_ALIGNED_MAXORDER+l] +=  rho[0][l] * eky_p *
                u_brick_none[k][mz][my][mx];
              ekz[k*INTEL_P3M_ALIGNED_MAXORDER+l] +=  rho[0][l] * ekz_p *
                u_brick_none[k][mz][my][mx];
            }
          }
        }
      }

      _alignvar(FFT_SCALAR ekx_tot[nsplit], 64);
      _alignvar(FFT_SCALAR eky_tot[nsplit], 64);
      _alignvar(FFT_SCALAR ekz_tot[nsplit], 64);
      for (int k = 0; k < nsplit; k++) {
        ekx_tot[k] = eky_tot[k] = ekz_tot[k] = ZEROF;
      }

      #if defined(LMP_SIMD_COMPILER)
      #pragma loop_count min(2), max(INTEL_P3M_ALIGNED_MAXORDER), avg(7)
      #endif
      for (int l = 0; l < order; l++){
        for (int k = 0; k < nsplit; k++) {
          ekx_tot[k] += ekx[k*INTEL_P3M_ALIGNED_MAXORDER+l];
          eky_tot[k] += eky[k*INTEL_P3M_ALIGNED_MAXORDER+l];
          ekz_tot[k] += ekz[k*INTEL_P3M_ALIGNED_MAXORDER+l];
        }
      }

      for (int k = 0; k < nsplit; k++) {
        ekx_tot[k] *= hx_inv;
        eky_tot[k] *= hy_inv;
        ekz_tot[k] *= hz_inv;
      }
      // convert D-field to force

      const int type = atom->type[i];

      const flt_t s1 = x[i][0] * hx_inv;
      const flt_t s2 = x[i][1] * hy_inv;
      const flt_t s3 = x[i][2] * hz_inv;
      flt_t sf1 = fsf_coeff0 * sin(ftwo_pi * s1);
      sf1 += fsf_coeff1 * sin(ffour_pi * s1);

      flt_t sf2 = fsf_coeff2 * sin(ftwo_pi * s2);
      sf2 += fsf_coeff3 * sin(ffour_pi * s2);

      flt_t sf3 = fsf_coeff4 * sin(ftwo_pi * s3);
      sf3 += fsf_coeff5 * sin(ffour_pi * s3);
      for (int k = 0; k < nsplit; k++) {
        const flt_t lj = B[nsplit*type + k];
        const flt_t twoljsq = lj*lj * B[k] * 2;
        flt_t sf = sf1*twoljsq;
        f[i][0] += lj * ekx_tot[k] - sf;
        sf = sf2*twoljsq;
        f[i][1] += lj * eky_tot[k] - sf;
        sf = sf3*twoljsq;
        if (slabflag != 2) f[i][2] += lj * ekz_tot[k] -  sf;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   precompute rho coefficients as a lookup table to save time in make_rho
   and fieldforce.  Instead of doing this polynomial for every atom 6 times
   per time step, precompute it for some number of points.
------------------------------------------------------------------------- */

void PPPMDispIntel::precompute_rho()
{

  half_rho_scale = (rho_points - 1.)/2.;
  half_rho_scale_plus = half_rho_scale + 0.5;

  for (int i = 0; i < rho_points; i++) {
    FFT_SCALAR dx = -1. + 1./half_rho_scale * (FFT_SCALAR)i;
    #if defined(LMP_SIMD_COMPILER)
    #pragma simd
    #endif
    for (int k=nlower; k<=nupper;k++){
      FFT_SCALAR r1 = ZEROF;
      for(int l=order-1; l>=0; l--){
        r1 = rho_coeff[l][k] + r1*dx;
      }
      rho_lookup[i][k-nlower] = r1;
    }
    for (int k = nupper-nlower+1; k < INTEL_P3M_ALIGNED_MAXORDER; k++) {
      rho_lookup[i][k] = 0;
    }
    if (differentiation_flag == 1) {
      #if defined(LMP_SIMD_COMPILER)
      #pragma simd
      #endif
      for(int k=nlower; k<=nupper;k++){
        FFT_SCALAR r1 = ZEROF;
        for(int l=order-2; l>=0; l--){
          r1 = drho_coeff[l][k] + r1*dx;
        }
        drho_lookup[i][k-nlower] = r1;
      }
      for (int k = nupper-nlower+1; k < INTEL_P3M_ALIGNED_MAXORDER; k++) {
        drho_lookup[i][k] = 0;
      }
    }
  }
  for (int i = 0; i < rho_points; i++) {
    FFT_SCALAR dx = -1. + 1./half_rho_scale * (FFT_SCALAR)i;
    #if defined(LMP_SIMD_COMPILER)
    #pragma simd
    #endif
    for (int k=nlower_6; k<=nupper_6;k++){
      FFT_SCALAR r1 = ZEROF;
      for(int l=order_6-1; l>=0; l--){
        r1 = rho_coeff_6[l][k] + r1*dx;
      }
      rho6_lookup[i][k-nlower_6] = r1;
    }
    for (int k = nupper_6-nlower_6+1; k < INTEL_P3M_ALIGNED_MAXORDER; k++) {
      rho6_lookup[i][k] = 0;
    }
    if (differentiation_flag == 1) {
      #if defined(LMP_SIMD_COMPILER)
      #pragma simd
      #endif
      for(int k=nlower_6; k<=nupper_6;k++){
        FFT_SCALAR r1 = ZEROF;
        for(int l=order_6-2; l>=0; l--){
          r1 = drho_coeff_6[l][k] + r1*dx;
        }
        drho6_lookup[i][k-nlower_6] = r1;
      }
      for (int k = nupper_6-nlower_6+1; k < INTEL_P3M_ALIGNED_MAXORDER; k++) {
        drho6_lookup[i][k] = 0;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   Returns 0 if Intel optimizations for PPPM ignored due to offload
------------------------------------------------------------------------- */

#ifdef _LMP_INTEL_OFFLOAD
int PPPMDispIntel::use_base() {
  return _use_base;
}
#endif
