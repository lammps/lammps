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
   Contributing authors: Rodrigo Canales (RWTH Aachen University)
                         W. Michael Brown (Intel)
------------------------------------------------------------------------- */

#include <mpi.h>
#include <stdlib.h>
#include <math.h>
#include "pppm_intel.h"
#include "atom.h"
#include "error.h"
#include "gridcomm.h"
#include "math_const.h"
#include "math_special.h"
#include "memory.h"
#include "suffix.h"

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathSpecial;

#define MAXORDER 7
#define OFFSET 16384
#define LARGE 10000.0
#define SMALL 0.00001
#define EPS_HOC 1.0e-7

enum{REVERSE_RHO};
enum{FORWARD_IK,FORWARD_AD,FORWARD_IK_PERATOM,FORWARD_AD_PERATOM};

#ifdef FFT_SINGLE
#define ZEROF 0.0f
#define ONEF  1.0f
#else
#define ZEROF 0.0
#define ONEF  1.0
#endif

/* ---------------------------------------------------------------------- */

PPPMIntel::PPPMIntel(LAMMPS *lmp, int narg, char **arg) : PPPM(lmp, narg, arg)
{
  suffix_flag |= Suffix::INTEL;
}

PPPMIntel::~PPPMIntel()
{
}

/* ----------------------------------------------------------------------
   called once before run
------------------------------------------------------------------------- */

void PPPMIntel::init()
{
  PPPM::init();
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

  if (order > INTEL_P3M_MAXORDER)
    error->all(FLERR,"PPPM order greater than supported by USER-INTEL\n");

  /*
  if (fix->precision() == FixIntel::PREC_MODE_MIXED)
    pack_force_const(force_const_single, fix->get_mixed_buffers());
  else if (fix->precision() == FixIntel::PREC_MODE_DOUBLE)
    pack_force_const(force_const_double, fix->get_double_buffers());
  else
    pack_force_const(force_const_single, fix->get_single_buffers());
  */
}

/* ----------------------------------------------------------------------
   compute the PPPMIntel long-range force, energy, virial
------------------------------------------------------------------------- */

void PPPMIntel::compute(int eflag, int vflag)
{
  #ifdef _LMP_INTEL_OFFLOAD
  if (_use_base) {
    PPPM::compute(eflag, vflag);
    return;
  }
  #endif
  compute_first(eflag,vflag);
  compute_second(eflag,vflag);
}

/* ---------------------------------------------------------------------- */

void PPPMIntel::compute_first(int eflag, int vflag)
{
  int i,j;

  // set energy/virial flags
  // invoke allocate_peratom() if needed for first time

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = evflag_atom = eflag_global = vflag_global =
         eflag_atom = vflag_atom = 0;

  if (evflag_atom && !peratom_allocate_flag) {
    allocate_peratom();
    cg_peratom->ghost_notify();
    cg_peratom->setup();
  }

  // if atom count has changed, update qsum and qsqsum

  if (atom->natoms != natoms_original) {
    qsum_qsq();
    natoms_original = atom->natoms;
  }

  // return if there are no charges

  if (qsqsum == 0.0) return;

  // convert atoms from box to lamda coords

  if (triclinic == 0) boxlo = domain->boxlo;
  else {
    boxlo = domain->boxlo_lamda;
    domain->x2lamda(atom->nlocal);
  }

  // extend size of per-atom arrays if necessary

  if (atom->nmax > nmax) {
    memory->destroy(part2grid);
    nmax = atom->nmax;
    memory->create(part2grid,nmax,3,"pppm:part2grid");
  }

  // find grid points for all my particles
  // map my particle charge onto my local 3d density grid

  if (fix->precision() == FixIntel::PREC_MODE_MIXED) {
    particle_map<float,double>(fix->get_mixed_buffers());
    make_rho<float,double>(fix->get_mixed_buffers());
  } else if (fix->precision() == FixIntel::PREC_MODE_DOUBLE) {
    particle_map<double,double>(fix->get_double_buffers());
    make_rho<double,double>(fix->get_double_buffers());
  } else {
    particle_map<float,float>(fix->get_single_buffers());
    make_rho<float,float>(fix->get_single_buffers());
  }

  // all procs communicate density values from their ghost cells
  //   to fully sum contribution in their 3d bricks
  // remap from 3d decomposition to FFT decomposition

  cg->reverse_comm(this,REVERSE_RHO);
  brick2fft();

  // compute potential gradient on my FFT grid and
  //   portion of e_long on this proc's FFT grid
  // return gradients (electric fields) in 3d brick decomposition
  // also performs per-atom calculations via poisson_peratom()

  poisson();

  // all procs communicate E-field values
  // to fill ghost cells surrounding their 3d bricks

  if (differentiation_flag == 1) cg->forward_comm(this,FORWARD_AD);
  else cg->forward_comm(this,FORWARD_IK);

  // extra per-atom energy/virial communication

  if (evflag_atom) {
    if (differentiation_flag == 1 && vflag_atom)
      cg_peratom->forward_comm(this,FORWARD_AD_PERATOM);
    else if (differentiation_flag == 0)
      cg_peratom->forward_comm(this,FORWARD_IK_PERATOM);
  }
}

/* ---------------------------------------------------------------------- */

void PPPMIntel::compute_second(int eflag, int vflag)
{
  int i,j;

  // calculate the force on my particles

  if (differentiation_flag == 1) {
    if (fix->precision() == FixIntel::PREC_MODE_MIXED)
      fieldforce_ad<float,double>(fix->get_mixed_buffers());
    else if (fix->precision() == FixIntel::PREC_MODE_DOUBLE)
      fieldforce_ad<double,double>(fix->get_double_buffers());
    else
      fieldforce_ad<float,float>(fix->get_single_buffers());
  } else {
    if (fix->precision() == FixIntel::PREC_MODE_MIXED)
      fieldforce_ik<float,double>(fix->get_mixed_buffers());
    else if (fix->precision() == FixIntel::PREC_MODE_DOUBLE)
      fieldforce_ik<double,double>(fix->get_double_buffers());
    else
      fieldforce_ik<float,float>(fix->get_single_buffers());
  }

  // extra per-atom energy/virial communication

  if (evflag_atom) fieldforce_peratom();

  // sum global energy across procs and add in volume-dependent term

  const double qscale = qqrd2e * scale;

  if (eflag_global) {
    double energy_all;
    MPI_Allreduce(&energy,&energy_all,1,MPI_DOUBLE,MPI_SUM,world);
    energy = energy_all;

    energy *= 0.5*volume;
    energy -= g_ewald*qsqsum/MY_PIS +
      MY_PI2*qsum*qsum / (g_ewald*g_ewald*volume);
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
    int nlocal = atom->nlocal;
    int ntotal = nlocal;
    if (tip4pflag) ntotal += atom->nghost;

    if (eflag_atom) {
      for (i = 0; i < nlocal; i++) {
        eatom[i] *= 0.5;
        eatom[i] -= g_ewald*q[i]*q[i]/MY_PIS + MY_PI2*q[i]*qsum /
          (g_ewald*g_ewald*volume);
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
   find center grid pt for each of my particles
   check that full stencil for the particle will fit in my 3d brick
   store central grid pt indices in part2grid array
------------------------------------------------------------------------- */

template<class flt_t, class acc_t>
void PPPMIntel::particle_map(IntelBuffers<flt_t,acc_t> *buffers)
{
  int nx,ny,nz;

  ATOM_T * _noalias const x = buffers->get_x(0);
  int nlocal = atom->nlocal;

  int flag = 0;

  if (!ISFINITE(boxlo[0]) || !ISFINITE(boxlo[1]) || !ISFINITE(boxlo[2]))
    error->one(FLERR,"Non-numeric box dimensions - simulation unstable");

  const flt_t lo0 = boxlo[0];
  const flt_t lo1 = boxlo[1];
  const flt_t lo2 = boxlo[2];
  const flt_t xi = delxinv;
  const flt_t yi = delyinv;
  const flt_t zi = delzinv;
  const flt_t fshift = shift;

  #if defined(LMP_SIMD_COMPILER)
  #pragma vector aligned
  #pragma simd
  #endif
  for (int i = 0; i < nlocal; i++) {

    // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
    // current particle coord can be outside global and local box
    // add/subtract OFFSET to avoid int(-0.75) = 0 when want it to be -1

    nx = static_cast<int> ((x[i].x-lo0)*xi+fshift) - OFFSET;
    ny = static_cast<int> ((x[i].y-lo1)*yi+fshift) - OFFSET;
    nz = static_cast<int> ((x[i].z-lo2)*zi+fshift) - OFFSET;

    part2grid[i][0] = nx;
    part2grid[i][1] = ny;
    part2grid[i][2] = nz;

    // check that entire stencil around nx,ny,nz will fit in my 3d brick

    if (nx+nlower < nxlo_out || nx+nupper > nxhi_out ||
        ny+nlower < nylo_out || ny+nupper > nyhi_out ||
        nz+nlower < nzlo_out || nz+nupper > nzhi_out)
      flag = 1;
  }

  if (flag) error->one(FLERR,"Out of range atoms - cannot compute PPPM");
}


/* ----------------------------------------------------------------------
   create discretized "density" on section of global grid due to my particles
   density(x,y,z) = charge "density" at grid points of my 3d brick
   (nxlo:nxhi,nylo:nyhi,nzlo:nzhi) is extent of my brick (including ghosts)
   in global grid
------------------------------------------------------------------------- */

template<class flt_t, class acc_t>
void PPPMIntel::make_rho(IntelBuffers<flt_t,acc_t> *buffers)
{
  // clear 3d density array

  memset(&(density_brick[nzlo_out][nylo_out][nxlo_out]),0,
         ngrid*sizeof(FFT_SCALAR));

  // loop over my charges, add their contribution to nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt

  ATOM_T * _noalias const x = buffers->get_x(0);
  flt_t * _noalias const q = buffers->get_q(0);
  int nlocal = atom->nlocal;

  const flt_t lo0 = boxlo[0];
  const flt_t lo1 = boxlo[1];
  const flt_t lo2 = boxlo[2];
  const flt_t xi = delxinv;
  const flt_t yi = delyinv;
  const flt_t zi = delzinv;
  const flt_t fshift = shift;
  const flt_t fshiftone = shiftone;
  const flt_t fdelvolinv = delvolinv;

  for (int i = 0; i < nlocal; i++) {

    int nx = part2grid[i][0];
    int ny = part2grid[i][1];
    int nz = part2grid[i][2];
    FFT_SCALAR dx = nx+fshiftone - (x[i].x-lo0)*xi;
    FFT_SCALAR dy = ny+fshiftone - (x[i].y-lo1)*yi;
    FFT_SCALAR dz = nz+fshiftone - (x[i].z-lo2)*zi;

    flt_t rho[3][INTEL_P3M_MAXORDER];

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

    FFT_SCALAR z0 = fdelvolinv * q[i];
    for (int n = nlower; n <= nupper; n++) {
      int mz = n+nz;
      FFT_SCALAR y0 = z0*rho[2][n-nlower];
      for (int m = nlower; m <= nupper; m++) {
        int my = m+ny;
        FFT_SCALAR x0 = y0*rho[1][m-nlower];
        for (int l = nlower; l <= nupper; l++) {
          int mx = l+nx;
          density_brick[mz][my][mx] += x0*rho[0][l-nlower];
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get electric field & force on my particles for ik
------------------------------------------------------------------------- */

template<class flt_t, class acc_t>
void PPPMIntel::fieldforce_ik(IntelBuffers<flt_t,acc_t> *buffers)
{
  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of E-field on particle

  ATOM_T * _noalias const x = buffers->get_x(0);
  flt_t * _noalias const q = buffers->get_q(0);
  FORCE_T * _noalias const f = buffers->get_f();
  int nlocal = atom->nlocal;

  const flt_t lo0 = boxlo[0];
  const flt_t lo1 = boxlo[1];
  const flt_t lo2 = boxlo[2];
  const flt_t xi = delxinv;
  const flt_t yi = delyinv;
  const flt_t zi = delzinv;
  const flt_t fshiftone = shiftone;
  const flt_t fqqrd2es = qqrd2e * scale;

  #if defined(LMP_SIMD_COMPILER)
  #pragma vector aligned nontemporal
  #pragma simd
  #endif
  for (int i = 0; i < nlocal; i++) {
    int nx = part2grid[i][0];
    int ny = part2grid[i][1];
    int nz = part2grid[i][2];
    FFT_SCALAR dx = nx+fshiftone - (x[i].x-lo0)*xi;
    FFT_SCALAR dy = ny+fshiftone - (x[i].y-lo1)*yi;
    FFT_SCALAR dz = nz+fshiftone - (x[i].z-lo2)*zi;

    flt_t rho[3][INTEL_P3M_MAXORDER];

    for (int k = nlower; k <= nupper; k++) {
      FFT_SCALAR r1 = rho_coeff[order-1][k];
      FFT_SCALAR r2 = rho_coeff[order-1][k];
      FFT_SCALAR r3 = rho_coeff[order-1][k];
      for (int l = order-2; l >= 0; l--) {
        r1 = rho_coeff[l][k] + r1*dx;
        r2 = rho_coeff[l][k] + r2*dy;
        r3 = rho_coeff[l][k] + r3*dz;
      }
      rho[0][k-nlower] = r1;
      rho[1][k-nlower] = r2;
      rho[2][k-nlower] = r3;
    }

    FFT_SCALAR ekx, eky, ekz;
    ekx = eky = ekz = ZEROF;
    for (int n = nlower; n <= nupper; n++) {
      int mz = n+nz;
      FFT_SCALAR z0 = rho[2][n-nlower];
      for (int m = nlower; m <= nupper; m++) {
        int my = m+ny;
        FFT_SCALAR y0 = z0*rho[1][m-nlower];
        for (int l = nlower; l <= nupper; l++) {
          int mx = l+nx;
          FFT_SCALAR x0 = y0*rho[0][l-nlower];
          ekx -= x0*vdx_brick[mz][my][mx];
          eky -= x0*vdy_brick[mz][my][mx];
          ekz -= x0*vdz_brick[mz][my][mx];
        }
      }
    }

    // convert E-field to force

    const flt_t qfactor = fqqrd2es * q[i];
    f[i].x += qfactor*ekx;
    f[i].y += qfactor*eky;
    if (slabflag != 2) f[i].z += qfactor*ekz;
  }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get electric field & force on my particles for ad
------------------------------------------------------------------------- */

template<class flt_t, class acc_t>
void PPPMIntel::fieldforce_ad(IntelBuffers<flt_t,acc_t> *buffers)
{
  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of E-field on particle

  ATOM_T * _noalias const x = buffers->get_x(0);
  const flt_t * _noalias const q = buffers->get_q(0);
  FORCE_T * _noalias const f = buffers->get_f();
  int nlocal = atom->nlocal;

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

  const double *prd = domain->prd;
  const double xprd = prd[0];
  const double yprd = prd[1];
  const double zprd = prd[2];

  const flt_t hx_inv = nx_pppm/xprd;
  const flt_t hy_inv = ny_pppm/yprd;
  const flt_t hz_inv = nz_pppm/zprd;

  const flt_t fsf_coeff0 = sf_coeff[0];
  const flt_t fsf_coeff1 = sf_coeff[1];
  const flt_t fsf_coeff2 = sf_coeff[2];
  const flt_t fsf_coeff3 = sf_coeff[3];
  const flt_t fsf_coeff4 = sf_coeff[4];
  const flt_t fsf_coeff5 = sf_coeff[5];

  #if defined(LMP_SIMD_COMPILER)
  #pragma vector aligned nontemporal
  #pragma simd
  #endif
  for (int i = 0; i < nlocal; i++) {
    int nx = part2grid[i][0];
    int ny = part2grid[i][1];
    int nz = part2grid[i][2];
    FFT_SCALAR dx = nx+fshiftone - (x[i].x-lo0)*xi;
    FFT_SCALAR dy = ny+fshiftone - (x[i].y-lo1)*yi;
    FFT_SCALAR dz = nz+fshiftone - (x[i].z-lo2)*zi;

    flt_t rho[3][INTEL_P3M_MAXORDER];
    flt_t drho[3][INTEL_P3M_MAXORDER];

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

    FFT_SCALAR ekx, eky, ekz;
    ekx = eky = ekz = ZEROF;
    for (int n = nlower; n <= nupper; n++) {
      int mz = n+nz;
      for (int m = nlower; m <= nupper; m++) {
        int my = m+ny;
        FFT_SCALAR ekx_p = rho[1][m-nlower] * rho[2][n-nlower];
        FFT_SCALAR eky_p = drho[1][m-nlower] * rho[2][n-nlower];
        FFT_SCALAR ekz_p = rho[1][m-nlower] * drho[2][n-nlower];
        for (int l = nlower; l <= nupper; l++) {
          int mx = l+nx;
          ekx += drho[0][l-nlower] * ekx_p * u_brick[mz][my][mx];
          eky += rho[0][l-nlower] * eky_p * u_brick[mz][my][mx];
          ekz += rho[0][l-nlower] * ekz_p * u_brick[mz][my][mx];
        }
      }
    }
    ekx *= hx_inv;
    eky *= hy_inv;
    ekz *= hz_inv;

    // convert E-field to force

    const flt_t qfactor = fqqrd2es * q[i];
    const flt_t twoqsq = (flt_t)2.0 * q[i] * q[i];

    const flt_t s1 = x[i].x * hx_inv;
    const flt_t s2 = x[i].y * hy_inv;
    const flt_t s3 = x[i].z * hz_inv;
    flt_t sf = fsf_coeff0 * sin(ftwo_pi * s1);
    sf += fsf_coeff1 * sin(ffour_pi * s1);
    sf *= twoqsq;
    f[i].x += qfactor * ekx - fqqrd2es * sf;

    sf = fsf_coeff2 * sin(ftwo_pi * s2);
    sf += fsf_coeff3 * sin(ffour_pi * s2);
    sf *= twoqsq;
    f[i].y += qfactor * eky - fqqrd2es * sf;

    sf = fsf_coeff4 * sin(ftwo_pi * s3);
    sf += fsf_coeff5 * sin(ffour_pi * s3);
    sf *= twoqsq;

    if (slabflag != 2) f[i].z += qfactor * ekz - fqqrd2es * sf;
  }
}

/* ----------------------------------------------------------------------
  Pack data into intel package buffers if using LRT mode
------------------------------------------------------------------------- */

void PPPMIntel::pack_buffers()
{
  fix->start_watch(TIME_PACK);
  #if defined(_OPENMP)
  #pragma omp parallel default(none)
  #endif
  {
    int ifrom, ito, tid;
    IP_PRE_omp_range_id_align(ifrom, ito, tid, atom->nlocal+atom->nghost,
                              comm->nthreads, 
                              sizeof(IntelBuffers<float,double>::atom_t));
    if (fix->precision() == FixIntel::PREC_MODE_MIXED)
      fix->get_mixed_buffers()->thr_pack(ifrom,ito,1);
    else if (fix->precision() == FixIntel::PREC_MODE_DOUBLE)
      fix->get_double_buffers()->thr_pack(ifrom,ito,1);
    else
      fix->get_single_buffers()->thr_pack(ifrom,ito,1);
  }
  fix->stop_watch(TIME_PACK);
}

/* ----------------------------------------------------------------------
   Returns 0 if Intel optimizations for PPPM ignored due to offload
------------------------------------------------------------------------- */

#ifdef _LMP_INTEL_OFFLOAD
int PPPMIntel::use_base() {
  return _use_base;
}
#endif
