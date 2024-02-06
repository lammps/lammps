// clang-format off
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
   Contributing authors: William McDoniel (RWTH Aachen University)
                         Rodrigo Canales (RWTH Aachen University)
                         Markus Hoehnerbach (RWTH Aachen University)
                         W. Michael Brown (Intel)
------------------------------------------------------------------------- */

#include "pppm_intel.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "grid3d.h"
#include "math_const.h"
#include "math_special.h"
#include "memory.h"
#include "modify.h"
#include "suffix.h"

#include <cmath>

#include "omp_compat.h"

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathSpecial;

static constexpr int OFFSET = 16384;
static constexpr FFT_SCALAR ZEROF = 0.0;

enum { REVERSE_RHO };
enum { FORWARD_IK, FORWARD_AD, FORWARD_IK_PERATOM, FORWARD_AD_PERATOM };

/* ---------------------------------------------------------------------- */

PPPMIntel::PPPMIntel(LAMMPS *lmp) : PPPM(lmp)
{
  suffix_flag |= Suffix::INTEL;

  order = 7; //sets default stencil size to 7

  perthread_density = nullptr;
  particle_ekx = particle_eky = particle_ekz = nullptr;

  rho_lookup = drho_lookup = nullptr;
  rho_points = 0;

  _use_table = _use_lrt = 0;
}

PPPMIntel::~PPPMIntel()
{
  memory->destroy(perthread_density);
  memory->destroy(particle_ekx);
  memory->destroy(particle_eky);
  memory->destroy(particle_ekz);

  memory->destroy(rho_lookup);
  memory->destroy(drho_lookup);
}

/* ----------------------------------------------------------------------
   called once before run
------------------------------------------------------------------------- */

void PPPMIntel::init()
{
  PPPM::init();
  fix = static_cast<FixIntel *>(modify->get_fix_by_id("package_intel"));
  if (!fix) error->all(FLERR, "The 'package intel' command is required for /intel styles");

  #ifdef _LMP_INTEL_OFFLOAD
  _use_base = 0;
  if (fix->offload_balance() != 0.0) {
    _use_base = 1;
    return;
  }
  #endif

  fix->kspace_init_check();

  _use_lrt = fix->lrt();

  // For vectorization, we need some padding in the end
  // The first thread computes on the global density
  if ((comm->nthreads > 1) && !_use_lrt) {
    memory->destroy(perthread_density);
    memory->create(perthread_density, comm->nthreads-1,
                   ngrid + INTEL_P3M_ALIGNED_MAXORDER,
                   "pppmintel:perthread_density");
  }

  _use_table = fix->pppm_table();
  if (_use_table) {
    rho_points = 5000;
    memory->destroy(rho_lookup);
    memory->create(rho_lookup, rho_points, INTEL_P3M_ALIGNED_MAXORDER,
                   "pppmintel:rho_lookup");
    if (differentiation_flag == 1) {
      memory->destroy(drho_lookup);
      memory->create(drho_lookup, rho_points, INTEL_P3M_ALIGNED_MAXORDER,
                     "pppmintel:drho_lookup");
    }
    precompute_rho();
  }

  if (order > INTEL_P3M_MAXORDER)
    error->all(FLERR,"PPPM order greater than supported by INTEL\n");

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
  // set energy/virial flags
  // invoke allocate_peratom() if needed for first time

  ev_init(eflag,vflag);

  if (evflag_atom && !peratom_allocate_flag) allocate_peratom();

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
    if (differentiation_flag == 1) {
      memory->destroy(particle_ekx);
      memory->destroy(particle_eky);
      memory->destroy(particle_ekz);
    }
    nmax = atom->nmax;
    memory->create(part2grid,nmax,3,"pppm:part2grid");
    if (differentiation_flag == 1) {
      memory->create(particle_ekx, nmax, "pppmintel:pekx");
      memory->create(particle_eky, nmax, "pppmintel:peky");
      memory->create(particle_ekz, nmax, "pppmintel:pekz");
    }
  }

  // find grid points for all my particles
  // map my particle charge onto my local 3d density grid
  // optimized versions can only be used for orthogonal boxes

  if (triclinic) {
    PPPM::particle_map();
    PPPM::make_rho();
  } else {

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
  }

  // all procs communicate density values from their ghost cells
  //   to fully sum contribution in their 3d bricks
  // remap from 3d decomposition to FFT decomposition

  gc->reverse_comm(Grid3d::KSPACE,this,REVERSE_RHO,1,sizeof(FFT_SCALAR),
                   gc_buf1,gc_buf2,MPI_FFT_SCALAR);
  brick2fft();

  // compute potential gradient on my FFT grid and
  //   portion of e_long on this proc's FFT grid
  // return gradients (electric fields) in 3d brick decomposition
  // also performs per-atom calculations via poisson_peratom()

  if (differentiation_flag == 1) poisson_ad();
  else poisson_ik();

  // all procs communicate E-field values
  // to fill ghost cells surrounding their 3d bricks

  if (differentiation_flag == 1)
    gc->forward_comm(Grid3d::KSPACE,this,FORWARD_AD,1,sizeof(FFT_SCALAR),
                     gc_buf1,gc_buf2,MPI_FFT_SCALAR);
  else
    gc->forward_comm(Grid3d::KSPACE,this,FORWARD_IK,3,sizeof(FFT_SCALAR),
                     gc_buf1,gc_buf2,MPI_FFT_SCALAR);

  // extra per-atom energy/virial communication

  if (evflag_atom) {
    if (differentiation_flag == 1 && vflag_atom)
      gc->forward_comm(Grid3d::KSPACE,this,FORWARD_AD_PERATOM,6,sizeof(FFT_SCALAR),
                       gc_buf1,gc_buf2,MPI_FFT_SCALAR);
    else if (differentiation_flag == 0)
      gc->forward_comm(Grid3d::KSPACE,this,FORWARD_IK_PERATOM,7,sizeof(FFT_SCALAR),
                       gc_buf1,gc_buf2,MPI_FFT_SCALAR);
  }
}

/* ---------------------------------------------------------------------- */

void PPPMIntel::compute_second(int /*eflag*/, int /*vflag*/)
{
  int i,j;

  // calculate the force on my particles
  // optimized versions can only be used for orthogonal boxes

  if (triclinic) {
    PPPM::fieldforce();
  } else {
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
  ATOM_T * _noalias const x = buffers->get_x(0);
  int nlocal = atom->nlocal;
  int nthr;
  if (_use_lrt)
    nthr = 1;
  else
    nthr = comm->nthreads;

  int flag = 0;

  if (!std::isfinite(boxlo[0]) || !std::isfinite(boxlo[1]) || !std::isfinite(boxlo[2]))
    error->one(FLERR,"Non-numeric box dimensions - simulation unstable");

  #if defined(_OPENMP)
  #pragma omp parallel LMP_DEFAULT_NONE \
    shared(nlocal, nthr) reduction(+:flag) if (!_use_lrt)
  #endif
  {
    const flt_t lo0 = boxlo[0];
    const flt_t lo1 = boxlo[1];
    const flt_t lo2 = boxlo[2];
    const flt_t xi = delxinv;
    const flt_t yi = delyinv;
    const flt_t zi = delzinv;
    const flt_t fshift = shift;

    int iifrom, iito, tid;
    IP_PRE_omp_range_id_align(iifrom, iito, tid, nlocal, nthr, sizeof(ATOM_T));

    #if defined(LMP_SIMD_COMPILER)
#if defined(USE_OMP_SIMD)
    #pragma omp simd reduction(+:flag)
#else
    #pragma simd reduction(+:flag)
#endif
    #pragma vector aligned
    #endif
    for (int i = iifrom; i < iito; i++) {

      // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
      // current particle coord can be outside global and local box
      // add/subtract OFFSET to avoid int(-0.75) = 0 when want it to be -1
      int nx = static_cast<int> ((x[i].x-lo0)*xi+fshift) - OFFSET;
      int ny = static_cast<int> ((x[i].y-lo1)*yi+fshift) - OFFSET;
      int nz = static_cast<int> ((x[i].z-lo2)*zi+fshift) - OFFSET;

      part2grid[i][0] = nx;
      part2grid[i][1] = ny;
      part2grid[i][2] = nz;

      // check that entire stencil around nx,ny,nz will fit in my 3d brick

      if (nx+nlower < nxlo_out || nx+nupper > nxhi_out ||
          ny+nlower < nylo_out || ny+nupper > nyhi_out ||
          nz+nlower < nzlo_out || nz+nupper > nzhi_out)
        flag = 1;
    }
  }

  if (flag) error->one(FLERR,"Out of range atoms - cannot compute PPPM");
}


/* ----------------------------------------------------------------------
   create discretized "density" on section of global grid due to my particles
   density(x,y,z) = charge "density" at grid points of my 3d brick
   (nxlo:nxhi,nylo:nyhi,nzlo:nzhi) is extent of my brick (including ghosts)
   in global grid
------------------------------------------------------------------------- */

template<class flt_t, class acc_t, int use_table>
void PPPMIntel::make_rho(IntelBuffers<flt_t,acc_t> *buffers)
{
  FFT_SCALAR * _noalias global_density =
    &(density_brick[nzlo_out][nylo_out][nxlo_out]);

  // loop over my charges, add their contribution to nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt

  ATOM_T * _noalias const x = buffers->get_x(0);
  flt_t * _noalias const q = buffers->get_q(0);
  int nlocal = atom->nlocal;
  int nthr;
  if (_use_lrt)
    nthr = 1;
  else
    nthr = comm->nthreads;

  #if defined(_OPENMP)
  #pragma omp parallel LMP_DEFAULT_NONE \
    shared(nthr, nlocal, global_density) if (!_use_lrt)
  #endif
  {
    const int nix = nxhi_out - nxlo_out + 1;
    const int niy = nyhi_out - nylo_out + 1;

    const flt_t lo0 = boxlo[0];
    const flt_t lo1 = boxlo[1];
    const flt_t lo2 = boxlo[2];
    const flt_t xi = delxinv;
    const flt_t yi = delyinv;
    const flt_t zi = delzinv;
    const flt_t fshiftone = shiftone;
    const flt_t fdelvolinv = delvolinv;

    int ifrom, ito, tid;
    IP_PRE_omp_range_id(ifrom, ito, tid, nlocal, nthr);
    FFT_SCALAR * _noalias my_density = tid == 0 ?
      global_density : perthread_density[tid - 1];
    // clear 3d density array
    memset(my_density, 0, ngrid * sizeof(FFT_SCALAR));

    for (int i = ifrom; i < ito; i++) {

      int nx = part2grid[i][0];
      int ny = part2grid[i][1];
      int nz = part2grid[i][2];

      int nysum = nlower + ny - nylo_out;
      int nxsum = nlower + nx - nxlo_out;
      int nzsum = (nlower + nz - nzlo_out)*nix*niy + nysum*nix + nxsum;

      FFT_SCALAR dx = nx+fshiftone - (x[i].x-lo0)*xi;
      FFT_SCALAR dy = ny+fshiftone - (x[i].y-lo1)*yi;
      FFT_SCALAR dz = nz+fshiftone - (x[i].z-lo2)*zi;

      _alignvar(flt_t rho[3][INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};

      if (use_table) {
        dx = dx*half_rho_scale + half_rho_scale_plus;
        int idx = dx;
        dy = dy*half_rho_scale + half_rho_scale_plus;
        int idy = dy;
        dz = dz*half_rho_scale + half_rho_scale_plus;
        int idz = dz;
        #if defined(LMP_SIMD_COMPILER)
#if defined(USE_OMP_SIMD)
        #pragma omp simd
#else
        #pragma simd
#endif
        #endif
        for (int k = 0; k < INTEL_P3M_ALIGNED_MAXORDER; k++) {
          rho[0][k] = rho_lookup[idx][k];
          rho[1][k] = rho_lookup[idy][k];
          rho[2][k] = rho_lookup[idz][k];
        }
      } else {
        #if defined(LMP_SIMD_COMPILER)
#if defined(USE_OMP_SIMD)
        #pragma omp simd
#else
        #pragma simd
#endif
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
#if defined(USE_OMP_SIMD)
          #pragma omp simd
#else
          #pragma simd
#endif
          #endif
          for (int l = 0; l < INTEL_P3M_ALIGNED_MAXORDER; l++) {
            int mzyx = l + mzy;
            my_density[mzyx] += x0*rho[0][l];
          }
        }
      }
    }
  }

  // reduce all the perthread_densities into global_density
  if (nthr > 1) {
    #if defined(_OPENMP)
    #pragma omp parallel LMP_DEFAULT_NONE \
      shared(nthr, global_density) if (!_use_lrt)
    #endif
    {
      int ifrom, ito, tid;
      IP_PRE_omp_range_id(ifrom, ito, tid, ngrid, nthr);

      #if defined(LMP_SIMD_COMPILER)
#if defined(USE_OMP_SIMD)
      #pragma omp simd
#else
      #pragma simd
#endif
      #endif
      for (int i = ifrom; i < ito; i++) {
        for (int j = 1; j < nthr; j++) {
          global_density[i] += perthread_density[j-1][i];
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get electric field & force on my particles for ik
------------------------------------------------------------------------- */

template<class flt_t, class acc_t, int use_table>
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
  int nthr;
  if (_use_lrt)
    nthr = 1;
  else
    nthr = comm->nthreads;

  if (fix->need_zero(0)) {
    int zl = nlocal;
    if (force->newton_pair) zl += atom->nghost;
    memset(f, 0, zl * sizeof(FORCE_T));
  }

  #if defined(_OPENMP)
  #pragma omp parallel LMP_DEFAULT_NONE \
    shared(nlocal, nthr) if (!_use_lrt)
  #endif
  {
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

    _alignvar(flt_t rho0[2 * INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
    _alignvar(flt_t rho1[INTEL_P3M_ALIGNED_MAXORDER] , 64)= {0};
    _alignvar(flt_t rho2[INTEL_P3M_ALIGNED_MAXORDER] , 64)= {0};

    for (int i = ifrom; i < ito; i++) {
      int nx = part2grid[i][0];
      int ny = part2grid[i][1];
      int nz = part2grid[i][2];

      int nxsum = nx + nlower;
      int nysum = ny + nlower;
      int nzsum = nz + nlower;

      FFT_SCALAR dx = nx+fshiftone - (x[i].x-lo0)*xi;
      FFT_SCALAR dy = ny+fshiftone - (x[i].y-lo1)*yi;
      FFT_SCALAR dz = nz+fshiftone - (x[i].z-lo2)*zi;

      if (use_table) {
        dx = dx*half_rho_scale + half_rho_scale_plus;
        int idx = dx;
        dy = dy*half_rho_scale + half_rho_scale_plus;
        int idy = dy;
        dz = dz*half_rho_scale + half_rho_scale_plus;
        int idz = dz;
        #if defined(LMP_SIMD_COMPILER)
#if defined(USE_OMP_SIMD)
        #pragma omp simd
#else
        #pragma simd
#endif
        #endif
        for (int k = 0; k < INTEL_P3M_ALIGNED_MAXORDER; k++) {
          rho0[k] = rho_lookup[idx][k];
          rho1[k] = rho_lookup[idy][k];
          rho2[k] = rho_lookup[idz][k];
        }
      } else {
        #if defined(LMP_SIMD_COMPILER)
#if defined(USE_OMP_SIMD)
        #pragma omp simd
#else
        #pragma simd
#endif
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
#if defined(USE_OMP_SIMD)
          #pragma omp simd
#else
          #pragma simd
#endif
          #endif
          for (int l = 0; l < INTEL_P3M_ALIGNED_MAXORDER; l++) {
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
      f[i].x += qfactor*ekx;
      f[i].y += qfactor*eky;
      if (slabflag != 2) f[i].z += qfactor*ekz;
    }
  }
}

/* ----------------------------------------------------------------------
   interpolate from grid to get electric field & force on my particles for ad
------------------------------------------------------------------------- */

template<class flt_t, class acc_t, int use_table>
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
  int nthr;
  if (_use_lrt)
    nthr = 1;
  else
    nthr = comm->nthreads;

  FFT_SCALAR * _noalias const particle_ekx = this->particle_ekx;
  FFT_SCALAR * _noalias const particle_eky = this->particle_eky;
  FFT_SCALAR * _noalias const particle_ekz = this->particle_ekz;

  if (fix->need_zero(0)) {
    int zl = nlocal;
    if (force->newton_pair) zl += atom->nghost;
    memset(f, 0, zl * sizeof(FORCE_T));
  }

  #if defined(_OPENMP)
  #pragma omp parallel LMP_DEFAULT_NONE \
    shared(nlocal, nthr) if (!_use_lrt)
  #endif
  {
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

    int ifrom, ito, tid;
    IP_PRE_omp_range_id(ifrom, ito, tid, nlocal, nthr);

    _alignvar(flt_t rho[3][INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};
    _alignvar(flt_t drho[3][INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};

    for (int i = ifrom; i < ito; i++) {
      int nx = part2grid[i][0];
      int ny = part2grid[i][1];
      int nz = part2grid[i][2];
      FFT_SCALAR dx = nx+fshiftone - (x[i].x-lo0)*xi;
      FFT_SCALAR dy = ny+fshiftone - (x[i].y-lo1)*yi;
      FFT_SCALAR dz = nz+fshiftone - (x[i].z-lo2)*zi;

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
#if defined(USE_OMP_SIMD)
        #pragma omp simd
#else
        #pragma simd
#endif
        #endif
        for (int k = 0; k < INTEL_P3M_ALIGNED_MAXORDER; k++) {
          rho[0][k] = rho_lookup[idx][k];
          rho[1][k] = rho_lookup[idy][k];
          rho[2][k] = rho_lookup[idz][k];
          drho[0][k] = drho_lookup[idx][k];
          drho[1][k] = drho_lookup[idy][k];
          drho[2][k] = drho_lookup[idz][k];
        }
      } else {
        #if defined(LMP_SIMD_COMPILER)
#if defined(USE_OMP_SIMD)
        #pragma omp simd
#else
        #pragma simd
#endif
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
#if defined(USE_OMP_SIMD)
          #pragma omp simd
#else
          #pragma simd
#endif
          #endif
          for (int l = 0; l < INTEL_P3M_ALIGNED_MAXORDER; l++) {
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
      for (int l = 0; l < order; l++) {
        particle_ekx[i] += ekx[l];
        particle_eky[i] += eky[l];
        particle_ekz[i] += ekz[l];
      }
    }

    #if defined(LMP_SIMD_COMPILER)
#if defined(USE_OMP_SIMD)
    #pragma omp simd
#else
    #pragma simd
#endif
    #endif
    for (int i = ifrom; i < ito; i++) {
      i = IP_PRE_dword_index(i);
      particle_ekx[i] *= hx_inv;
      particle_eky[i] *= hy_inv;
      particle_ekz[i] *= hz_inv;

      // convert E-field to force

      const flt_t qfactor = fqqrd2es * q[i];
      const flt_t twoqsq = (flt_t)2.0 * q[i] * q[i];

      const flt_t s1 = x[i].x * hx_inv;
      const flt_t s2 = x[i].y * hy_inv;
      const flt_t s3 = x[i].z * hz_inv;
      flt_t sf = fsf_coeff0 * std::sin(ftwo_pi * s1);
      sf += fsf_coeff1 * std::sin(ffour_pi * s1);
      sf *= twoqsq;
      f[i].x += qfactor * particle_ekx[i] - fqqrd2es * sf;

      sf = fsf_coeff2 * std::sin(ftwo_pi * s2);
      sf += fsf_coeff3 * std::sin(ffour_pi * s2);
      sf *= twoqsq;
      f[i].y += qfactor * particle_eky[i] - fqqrd2es * sf;

      sf = fsf_coeff4 * std::sin(ftwo_pi * s3);
      sf += fsf_coeff5 * std::sin(ffour_pi * s3);
      sf *= twoqsq;

      if (slabflag != 2) f[i].z += qfactor * particle_ekz[i] - fqqrd2es * sf;
    }
  }
}

/* ----------------------------------------------------------------------
   precompute rho coefficients as a lookup table to save time in make_rho
   and fieldforce.  Instead of doing this polynomial for every atom 6 times
   per time step, precompute it for some number of points.
------------------------------------------------------------------------- */

void PPPMIntel::precompute_rho()
{

  half_rho_scale = (rho_points - 1.)/2.;
  half_rho_scale_plus = half_rho_scale + 0.5;

  for (int i = 0; i < rho_points; i++) {
    FFT_SCALAR dx = -1. + 1./half_rho_scale * (FFT_SCALAR)i;
    #if defined(LMP_SIMD_COMPILER)
#if defined(USE_OMP_SIMD)
    #pragma omp simd
#else
    #pragma simd
#endif
    #endif
    for (int k=nlower; k<=nupper;k++) {
      FFT_SCALAR r1 = ZEROF;
      for (int l=order-1; l>=0; l--) {
        r1 = rho_coeff[l][k] + r1*dx;
      }
      rho_lookup[i][k-nlower] = r1;
    }
    for (int k = nupper-nlower+1; k < INTEL_P3M_ALIGNED_MAXORDER; k++) {
      rho_lookup[i][k] = 0;
    }
    if (differentiation_flag == 1) {
      #if defined(LMP_SIMD_COMPILER)
#if defined(USE_OMP_SIMD)
      #pragma omp simd
#else
      #pragma simd
#endif
      #endif
      for (int k=nlower; k<=nupper;k++) {
        FFT_SCALAR r1 = ZEROF;
        for (int l=order-2; l>=0; l--) {
          r1 = drho_coeff[l][k] + r1*dx;
        }
        drho_lookup[i][k-nlower] = r1;
      }
      for (int k = nupper-nlower+1; k < INTEL_P3M_ALIGNED_MAXORDER; k++) {
        drho_lookup[i][k] = 0;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local arrays
------------------------------------------------------------------------- */

double PPPMIntel::memory_usage()
{
  double bytes = PPPM::memory_usage();
  if ((comm->nthreads > 1) && !_use_lrt) {
    bytes += (double)(comm->nthreads - 1) * (ngrid + INTEL_P3M_ALIGNED_MAXORDER) *
      sizeof(FFT_SCALAR);
  }
  if (differentiation_flag == 1) {
    bytes += (double)3 * nmax * sizeof(FFT_SCALAR);
  }
  if (_use_table) {
    bytes += (double)rho_points * INTEL_P3M_ALIGNED_MAXORDER * sizeof(FFT_SCALAR);
    if (differentiation_flag == 1) {
      bytes += (double)rho_points * INTEL_P3M_ALIGNED_MAXORDER * sizeof(FFT_SCALAR);
    }
  }
  return bytes;
}

/* ----------------------------------------------------------------------
  Pack data into intel package buffers if using LRT mode
------------------------------------------------------------------------- */

void PPPMIntel::pack_buffers()
{
  fix->start_watch(TIME_PACK);
  int packthreads;
  if (comm->nthreads > INTEL_HTHREADS) packthreads = comm->nthreads;
  else packthreads = 1;
  #if defined(_OPENMP)
  #pragma omp parallel if (packthreads > 1)
  #endif
  {
    int ifrom, ito, tid;
    IP_PRE_omp_range_id_align(ifrom, ito, tid, atom->nlocal+atom->nghost,
                              packthreads,
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
   Allocate density_brick with extra padding for vector writes
------------------------------------------------------------------------- */

void PPPMIntel::allocate()
{
  PPPM::allocate();
  memory->destroy3d_offset(density_brick,nzlo_out,nylo_out,nxlo_out);
  create3d_offset(density_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                  nxlo_out,nxhi_out,"pppm:density_brick");

  if (differentiation_flag == 1) {
    memory->destroy3d_offset(u_brick,nzlo_out,nylo_out,nxlo_out);
    create3d_offset(u_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                    nxlo_out,nxhi_out,"pppm:u_brick");
  } else {
    memory->destroy3d_offset(vdx_brick,nzlo_out,nylo_out,nxlo_out);
    memory->destroy3d_offset(vdy_brick,nzlo_out,nylo_out,nxlo_out);
    memory->destroy3d_offset(vdz_brick,nzlo_out,nylo_out,nxlo_out);
    create3d_offset(vdx_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                    nxlo_out,nxhi_out,"pppm:vdx_brick");
    create3d_offset(vdy_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                    nxlo_out,nxhi_out,"pppm:vdy_brick");
    create3d_offset(vdz_brick,nzlo_out,nzhi_out,nylo_out,nyhi_out,
                    nxlo_out,nxhi_out,"pppm:vdz_brick");
  }
}

/* ----------------------------------------------------------------------
   Create 3D-offset allocation with extra padding for vector writes
------------------------------------------------------------------------- */

FFT_SCALAR *** PPPMIntel::create3d_offset(FFT_SCALAR ***&array, int n1lo,
                                          int n1hi, int n2lo, int n2hi,
                                          int n3lo, int n3hi,
                                          const char *name)
{
  int n1 = n1hi - n1lo + 1;
  int n2 = n2hi - n2lo + 1;
  int n3 = n3hi - n3lo + 1;

  bigint nbytes = ((bigint) sizeof(FFT_SCALAR)) * n1*n2*n3 +
    INTEL_P3M_ALIGNED_MAXORDER*2;
  auto data = (FFT_SCALAR *) memory->smalloc(nbytes,name);
  nbytes = ((bigint) sizeof(FFT_SCALAR *)) * n1*n2;
  auto plane = (FFT_SCALAR **) memory->smalloc(nbytes,name);
  nbytes = ((bigint) sizeof(FFT_SCALAR **)) * n1;
  array = (FFT_SCALAR ***) memory->smalloc(nbytes,name);

  bigint m;
  bigint n = 0;
  for (int i = 0; i < n1; i++) {
    m = ((bigint) i) * n2;
    array[i] = &plane[m];
    for (int j = 0; j < n2; j++) {
      plane[m+j] = &data[n];
      n += n3;
    }
  }

  m = ((bigint) n1) * n2;
  for (bigint i = 0; i < m; i++) array[0][i] -= n3lo;
  for (int i = 0; i < n1; i++) array[i] -= n2lo;
  array -= n1lo;
  return array;
}

/* ----------------------------------------------------------------------
   Returns 0 if Intel optimizations for PPPM ignored due to offload
------------------------------------------------------------------------- */

#ifdef _LMP_INTEL_OFFLOAD
int PPPMIntel::use_base() {
  return _use_base;
}
#endif

/* ----------------------------------------------------------------------
   allows usage in derived classes (pppm/electrode/intel)
------------------------------------------------------------------------- */
template void PPPMIntel::particle_map<float,double>(IntelBuffers<float,double> *buffers);
template void PPPMIntel::particle_map<double,double>(IntelBuffers<double,double> *buffers);
template void PPPMIntel::particle_map<float,float>(IntelBuffers<float,float> *buffers);
template void PPPMIntel::make_rho<float,double,0>(IntelBuffers<float,double> *buffers);
template void PPPMIntel::make_rho<double,double,0>(IntelBuffers<double,double> *buffers);
template void PPPMIntel::make_rho<float,float,0>(IntelBuffers<float,float> *buffers);
template void PPPMIntel::make_rho<float,double,1>(IntelBuffers<float,double> *buffers);
template void PPPMIntel::make_rho<double,double,1>(IntelBuffers<double,double> *buffers);
template void PPPMIntel::make_rho<float,float,1>(IntelBuffers<float,float> *buffers);
