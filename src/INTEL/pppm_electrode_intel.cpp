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

#include "pppm_electrode_intel.h"

#include "angle.h"
#include "atom.h"
#include "bond.h"
#include "citeme.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fft3d_wrap.h"
#include "force.h"
#include "grid3d.h"
#include "math_const.h"
#include "math_special.h"
#include "memory.h"
#include "neighbor.h"
#include "omp_compat.h"
#include "pair.h"
#include "pppm_intel.h"
#include "remap_wrap.h"
#include "slab_dipole.h"
#include "update.h"
#include "wire_dipole.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace std;

#define MAXORDER 7
#define OFFSET 16384
#define LARGE 10000.0
#define SMALL 0.00001
#define EPS_HOC 1.0e-7

enum { REVERSE_RHO };
enum { FORWARD_IK, FORWARD_AD, FORWARD_IK_PERATOM, FORWARD_AD_PERATOM };
enum : bool { ELECTRODE = true, ELECTROLYTE = false };

#ifdef FFT_SINGLE
#define ZEROF 0.0f
#define ONEF 1.0f
#else
#define ZEROF 0.0
#define ONEF 1.0
#endif

static const char cite_pppm_electrode[] =
    "kspace_style pppm/electrode command:\n\n"
    "@article{Ahrens2021,\n"
    "author = {Ahrens-Iwers, Ludwig J.V. and Mei{\\ss}ner, Robert H.},\n"
    "doi = {10.1063/5.0063381},\n"
    "title = {{Constant potential simulations on a mesh}},\n"
    "journal = {Journal of Chemical Physics},\n"
    "year = {2021}\n"
    "volume = {155},\n"
    "pages = {104104},\n"
    "}\n";

PPPMElectrodeIntel::PPPMElectrodeIntel(LAMMPS *lmp) :
    PPPMIntel(lmp), ElectrodeKSpace(), electrolyte_density_brick(nullptr),
    electrolyte_density_fft(nullptr), boundcorr(nullptr)
{
  if (lmp->citeme) lmp->citeme->add(cite_pppm_electrode);

  group_group_enable = 0;
  electrolyte_density_brick = nullptr;
  electrolyte_density_fft = nullptr;
  compute_vector_called = false;
  last_source_grpbit = 1 << 0;    // default to "all"
  last_invert_source = false;
}

PPPMElectrodeIntel::~PPPMElectrodeIntel()
{
  memory->destroy3d_offset(electrolyte_density_brick, nzlo_out, nylo_out, nxlo_out);
  memory->destroy(electrolyte_density_fft);
  if ((differentiation_flag != 1) && !peratom_allocate_flag)
    memory->destroy3d_offset(u_brick, nzlo_out, nylo_out, nxlo_out);
}

void PPPMElectrodeIntel::init()
{
  PPPMIntel::init();

  // PPPM/electrode/intel - specific checks
  if (slabflag == 3)
    error->all(FLERR, "Cannot (yet) use PPPM/electrode/intel with 'kspace_modify slab ew2d'");

  triclinic_check();
  triclinic = domain->triclinic;
  if (triclinic) error->all(FLERR, "Cannot (yet) use PPPM/electrode with triclinic box ");
  if (domain->dimension == 2) error->all(FLERR, "Cannot use PPPM/electrode with 2d simulation");

  if (!atom->q_flag) error->all(FLERR, "KSpace style requires atom attribute q");

  if (slabflag == 0 && wireflag == 0 && domain->nonperiodic > 0)
    error->all(FLERR, "Cannot use non-periodic boundaries with PPPM/electrode");
  if (slabflag) {
    if (domain->xperiodic != 1 || domain->yperiodic != 1 || domain->boundary[2][0] != 1 ||
        domain->boundary[2][1] != 1)
      error->all(FLERR, "Incorrect boundaries with slab PPPM/electrode");
  } else if (wireflag) {
    if (domain->zperiodic != 1 || domain->boundary[0][0] != 1 || domain->boundary[0][1] != 1 ||
        domain->boundary[1][0] != 1 || domain->boundary[1][1] != 1)
      error->all(FLERR, "Incorrect boundaries with wire PPPM/electrode");
  }
  compute_step = -1;
}

/* ----------------------------------------------------------------------
   adjust PPPM coeffs, called initially and whenever volume has changed
------------------------------------------------------------------------- */

void PPPMElectrodeIntel::setup()
{
  // PPPM/electrode/intel - specific checks

  if (slabflag == 0 && wireflag == 0 && domain->nonperiodic > 0)
    error->all(FLERR, "Cannot use non-periodic boundaries with PPPM/electrode");
  if (slabflag) {
    if (domain->xperiodic != 1 || domain->yperiodic != 1 || domain->boundary[2][0] != 1 ||
        domain->boundary[2][1] != 1)
      error->all(FLERR, "Incorrect boundaries with slab PPPM/electrode");
  } else if (wireflag) {
    if (domain->zperiodic != 1 || domain->boundary[0][0] != 1 || domain->boundary[0][1] != 1 ||
        domain->boundary[1][0] != 1 || domain->boundary[1][1] != 1)
      error->all(FLERR, "Incorrect boundaries with wire PPPM/electrode");
  }

  double *prd = domain->prd;

  // volume-dependent factors
  // adjust z dimension for 2d slab PPPM
  // z dimension for 3d PPPM is zprd since slab_volfactor = 1.0

  double xprd = prd[0];
  double yprd = prd[1];
  double zprd = prd[2];
  double xprd_wire = xprd * wire_volfactor;
  double yprd_wire = yprd * wire_volfactor;
  double zprd_slab = zprd * slab_volfactor;
  volume = xprd_wire * yprd_wire * zprd_slab;

  // wire_volfactor hacks to reuse compute_gf code
  prd[0] *= wire_volfactor;
  prd[1] *= wire_volfactor;
  PPPMIntel::setup();
  prd[0] /= wire_volfactor;
  prd[1] /= wire_volfactor;

}

void PPPMElectrodeIntel::compute(int eflag, int vflag)
{
#ifdef _LMP_INTEL_OFFLOAD
  if (_use_base) {
    error->all(FLERR, "Cannot use pppm/electrode/intel with offload");
    // PPPM::compute(eflag, vflag);
    // would work if the above line referred to PPPMElectrode
    // but the required multiple inheritances would be insane
  }
#endif

  ev_init(eflag, vflag);

  if (evflag_atom && !peratom_allocate_flag) allocate_peratom();

  // update qsum and qsqsum
  qsum_qsq();

  // return if there are no charges
  if (qsqsum == 0.0) return;

  // convert atoms from box to lamda coords
  if (triclinic == 0)
    boxlo = domain->boxlo;
  else {
    boxlo = domain->boxlo_lamda;
    domain->x2lamda(atom->nlocal);
  }

  start_compute();

  // find grid points for all my particles
  // map my particle charge onto my local 3d density grid
  // optimized versions can only be used for orthogonal boxes

  if (compute_vector_called) {
    // electrolyte_density_brick is filled, so we can
    // grab only electrode atoms
    switch (fix->precision()) {
      case FixIntel::PREC_MODE_MIXED:
        make_rho_in_brick<float, double>(fix->get_mixed_buffers(), last_source_grpbit,
                                         density_brick, !last_invert_source);
        break;
      case FixIntel::PREC_MODE_DOUBLE:
        make_rho_in_brick<double, double>(fix->get_double_buffers(), last_source_grpbit,
                                          density_brick, !last_invert_source);
        break;
      default:
        make_rho_in_brick<float, float>(fix->get_single_buffers(), last_source_grpbit,
                                        density_brick, !last_invert_source);
    }
    gc->reverse_comm(Grid3d::KSPACE, this, REVERSE_RHO, 1, sizeof(FFT_SCALAR), gc_buf1, gc_buf2,
                     MPI_FFT_SCALAR);
    for (int nz = nzlo_out; nz <= nzhi_out; nz++)
      for (int ny = nylo_out; ny <= nyhi_out; ny++)
        for (int nx = nxlo_out; nx <= nxhi_out; nx++) {
          density_brick[nz][ny][nx] += electrolyte_density_brick[nz][ny][nx];
        }
  } else {
    switch (fix->precision()) {
      case FixIntel::PREC_MODE_MIXED:
        PPPMIntel::make_rho<float, double>(fix->get_mixed_buffers());
        break;
      case FixIntel::PREC_MODE_DOUBLE:
        PPPMIntel::make_rho<double, double>(fix->get_double_buffers());
        break;
      default:
        PPPMIntel::make_rho<float, float>(fix->get_single_buffers());
    }
    // all procs communicate density values from their ghost cells
    //   to fully sum contribution in their 3d bricks
    // remap from 3d decomposition to FFT decomposition

    gc->reverse_comm(Grid3d::KSPACE, this, REVERSE_RHO, 1, sizeof(FFT_SCALAR), gc_buf1, gc_buf2,
                     MPI_FFT_SCALAR);
  }

  brick2fft();

  // compute potential gradient on my FFT grid and
  //   portion of e_long on this proc's FFT grid
  // return gradients (electric fields) in 3d brick decomposition
  // also performs per-atom calculations via poisson_peratom()

  if (differentiation_flag == 1)
    poisson_ad();
  else
    poisson_ik();

  // all procs communicate E-field values
  // to fill ghost cells surrounding their 3d bricks

  if (differentiation_flag == 1)
    gc->forward_comm(Grid3d::KSPACE, this, FORWARD_AD, 1, sizeof(FFT_SCALAR), gc_buf1, gc_buf2,
                     MPI_FFT_SCALAR);
  else
    gc->forward_comm(Grid3d::KSPACE, this, FORWARD_IK, 3, sizeof(FFT_SCALAR), gc_buf1, gc_buf2,
                     MPI_FFT_SCALAR);

  // extra per-atom energy/virial communication

  if (evflag_atom) {
    if (differentiation_flag == 1 && vflag_atom)
      gc->forward_comm(Grid3d::KSPACE, this, FORWARD_AD_PERATOM, 6, sizeof(FFT_SCALAR), gc_buf1,
                       gc_buf2, MPI_FFT_SCALAR);
    else if (differentiation_flag == 0)
      gc->forward_comm(Grid3d::KSPACE, this, FORWARD_IK_PERATOM, 7, sizeof(FFT_SCALAR), gc_buf1,
                       gc_buf2, MPI_FFT_SCALAR);
  }
  int tempslabflag = slabflag;
  slabflag = 0;    // bypass compute_second's slabcorr()
  PPPMIntel::compute_second(eflag, vflag);
  slabflag = tempslabflag;
  boundcorr->compute_corr(qsum,  eflag_atom, eflag_global, energy, eatom);
  compute_vector_called = false;
}

void PPPMElectrodeIntel::start_compute()
{
  if (compute_step < update->ntimestep) {
    if (compute_step == -1) setup();
    boxlo = domain->boxlo;
    // extend size of per-atom arrays if necessary
    if (atom->nmax > nmax) {
      memory->destroy(part2grid);
      if (differentiation_flag == 1) {
        memory->destroy(particle_ekx);
        memory->destroy(particle_eky);
        memory->destroy(particle_ekz);
      }
      nmax = atom->nmax;
      memory->create(part2grid, nmax, 3, "pppm:part2grid");
      if (differentiation_flag == 1) {
        memory->create(particle_ekx, nmax, "pppmintel:pekx");
        memory->create(particle_eky, nmax, "pppmintel:peky");
        memory->create(particle_ekz, nmax, "pppmintel:pekz");
      }
    }
    switch (fix->precision()) {
      case FixIntel::PREC_MODE_MIXED:
        PPPMIntel::particle_map<float, double>(fix->get_mixed_buffers());
        break;
      case FixIntel::PREC_MODE_DOUBLE:
        PPPMIntel::particle_map<double, double>(fix->get_double_buffers());
        break;
      default:
        PPPMIntel::particle_map<float, float>(fix->get_single_buffers());
    }
    compute_step = update->ntimestep;
  }
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */
void PPPMElectrodeIntel::compute_vector(double *vec, int sensor_grpbit, int source_grpbit,
                                        bool invert_source)
{
  start_compute();

  last_source_grpbit = source_grpbit;
  last_invert_source = invert_source;

  // temporarily store and switch pointers so we can use brick2fft() for
  // electrolyte density (without writing an additional function)
  FFT_SCALAR ***density_brick_real = density_brick;
  FFT_SCALAR *density_fft_real = density_fft;
  if (neighbor->ago != 0) pack_buffers(); // since midstep positions may be outdated
  switch (fix->precision()) {
    case FixIntel::PREC_MODE_MIXED:
      make_rho_in_brick<float, double>(fix->get_mixed_buffers(), source_grpbit,
                                       electrolyte_density_brick, invert_source);
      break;
    case FixIntel::PREC_MODE_DOUBLE:
      make_rho_in_brick<double, double>(fix->get_double_buffers(), source_grpbit,
                                        electrolyte_density_brick, invert_source);
      break;
    default:
      make_rho_in_brick<float, float>(fix->get_single_buffers(), source_grpbit,
                                      electrolyte_density_brick, invert_source);
  }
  density_brick = electrolyte_density_brick;
  density_fft = electrolyte_density_fft;
  gc->reverse_comm(Grid3d::KSPACE, this, REVERSE_RHO, 1, sizeof(FFT_SCALAR), gc_buf1, gc_buf2,
                   MPI_FFT_SCALAR);
  brick2fft();
  // switch back pointers
  density_brick = density_brick_real;
  density_fft = density_fft_real;

  // transform electrolyte charge density (r -> k) (complex conjugate)
  for (int i = 0, n = 0; i < nfft; i++) {
    work1[n++] = electrolyte_density_fft[i];
    work1[n++] = ZEROF;
  }
  fft1->compute(work1, work1, -1);

  // k->r FFT of Green's * electrolyte density = u_brick
  for (int i = 0, n = 0; i < nfft; i++) {
    work2[n] = work1[n] * greensfn[i];
    n++;
    work2[n] = work1[n] * greensfn[i];
    n++;
  }
  fft2->compute(work2, work2, 1);

  for (int k = nzlo_in, n = 0; k <= nzhi_in; k++)
    for (int j = nylo_in; j <= nyhi_in; j++)
      for (int i = nxlo_in; i <= nxhi_in; i++) {
        u_brick[k][j][i] = work2[n];
        n += 2;
      }

  gc->forward_comm(Grid3d::KSPACE, this, FORWARD_AD, 1, sizeof(FFT_SCALAR), gc_buf1, gc_buf2,
                   MPI_FFT_SCALAR);

  switch (fix->precision()) {
    case FixIntel::PREC_MODE_MIXED:
      project_psi<float, double>(fix->get_mixed_buffers(), vec, sensor_grpbit);
      break;
    case FixIntel::PREC_MODE_DOUBLE:
      project_psi<double, double>(fix->get_double_buffers(), vec, sensor_grpbit);
      break;
    default:
      project_psi<float, float>(fix->get_single_buffers(), vec, sensor_grpbit);
  }
  compute_vector_called = true;
}

// project u_brick with weight matrix
template <class flt_t, class acc_t, int use_table>
void PPPMElectrodeIntel::project_psi(IntelBuffers<flt_t, acc_t> *buffers, double *vec,
                                     int sensor_grpbit)
{
  ATOM_T *_noalias const x = buffers->get_x(0);

  int const nlocal = atom->nlocal;
  int nthr;

  if (_use_lrt)
    nthr = 1;
  else
    nthr = comm->nthreads;

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(nthr, vec, sensor_grpbit) if (!_use_lrt)
#endif
  {
    int *mask = atom->mask;
    const flt_t scaleinv = 1.0 / (nx_pppm * ny_pppm * nz_pppm);

    const flt_t lo0 = boxlo[0];
    const flt_t lo1 = boxlo[1];
    const flt_t lo2 = boxlo[2];
    const flt_t xi = delxinv;
    const flt_t yi = delyinv;
    const flt_t zi = delzinv;
    const flt_t fshiftone = shiftone;

    int ifrom, ito, tid;
    IP_PRE_omp_range_id(ifrom, ito, tid, nlocal, nthr);

    for (int i = ifrom; i < ito; i++) {
      if (!(mask[i] & sensor_grpbit)) continue;

      double v = 0.;
      // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
      // (dx,dy,dz) = distance to "lower left" grid pt
      // (mx,my,mz) = global coords of moving stencil pt

      int nx = part2grid[i][0];
      int ny = part2grid[i][1];
      int nz = part2grid[i][2];

      FFT_SCALAR dx = nx + fshiftone - (x[i].x - lo0) * xi;
      FFT_SCALAR dy = ny + fshiftone - (x[i].y - lo1) * yi;
      FFT_SCALAR dz = nz + fshiftone - (x[i].z - lo2) * zi;

      _alignvar(flt_t rho[3][INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};

      if (use_table) {
        dx = dx * half_rho_scale + half_rho_scale_plus;
        int idx = dx;
        dy = dy * half_rho_scale + half_rho_scale_plus;
        int idy = dy;
        dz = dz * half_rho_scale + half_rho_scale_plus;
        int idz = dz;
#if defined(LMP_SIMD_COMPILER)
#pragma simd
#endif
        for (int k = 0; k < INTEL_P3M_ALIGNED_MAXORDER; k++) {
          rho[0][k] = rho_lookup[idx][k];
          rho[1][k] = rho_lookup[idy][k];
          rho[2][k] = rho_lookup[idz][k];
        }
      } else {
#if defined(LMP_SIMD_COMPILER)
#pragma simd
#endif
        for (int k = nlower; k <= nupper; k++) {
          FFT_SCALAR r1, r2, r3;
          r1 = r2 = r3 = ZEROF;

          for (int l = order - 1; l >= 0; l--) {
            r1 = rho_coeff[l][k] + r1 * dx;
            r2 = rho_coeff[l][k] + r2 * dy;
            r3 = rho_coeff[l][k] + r3 * dz;
          }
          rho[0][k - nlower] = r1;
          rho[1][k - nlower] = r2;
          rho[2][k - nlower] = r3;
        }
      }

#if defined(LMP_SIMD_COMPILER)
#pragma loop_count min(2), max(INTEL_P3M_ALIGNED_MAXORDER), avg(7)
#endif
      for (int n = 0; n < order; n++) {
        int miz = nlower + n + nz;
        FFT_SCALAR z0 = rho[2][n];
#if defined(LMP_SIMD_COMPILER)
#pragma loop_count min(2), max(INTEL_P3M_ALIGNED_MAXORDER), avg(7)
#endif
        for (int m = 0; m < order; m++) {
          int miy = nlower + m + ny;
          FFT_SCALAR y0 = z0 * rho[1][m];
#if defined(LMP_SIMD_COMPILER)
#pragma loop_count min(2), max(INTEL_P3M_ALIGNED_MAXORDER), avg(7)
#pragma simd
#endif
          for (int l = 0; l < order; l++) {
            int mix = nlower + l + nx;
            v += y0 * rho[0][l] * u_brick[miz][miy][mix];
          }
        }
      }
      vec[i] += v * scaleinv;
    }
  }
}
/* ----------------------------------------------------------------------
  ---------------------------------------------------------------------  */
void PPPMElectrodeIntel::compute_matrix(bigint *imat, double **matrix, bool timer_flag)
{
  // TODO replace compute with required setup
  compute(1, 0);

  // fft green's function k -> r
  double *greens_real;
  memory->create(greens_real, nz_pppm * ny_pppm * nx_pppm, "pppm/electrode:greens_real");
  memset(greens_real, 0, nz_pppm * ny_pppm * nx_pppm * sizeof(double));
  for (int i = 0, n = 0; i < nfft; i++) {
    work2[n++] = greensfn[i];
    work2[n++] = ZEROF;
  }
  fft2->compute(work2, work2, -1);
  for (int k = nzlo_in, n = 0; k <= nzhi_in; k++)
    for (int j = nylo_in; j <= nyhi_in; j++)
      for (int i = nxlo_in; i <= nxhi_in; i++) {
        greens_real[ny_pppm * nx_pppm * k + nx_pppm * j + i] = work2[n];
        n += 2;
      }
  MPI_Allreduce(MPI_IN_PLACE, greens_real, nz_pppm * ny_pppm * nx_pppm, MPI_DOUBLE, MPI_SUM, world);
  int const nlocal = atom->nlocal;
  int nmat = std::count_if(&imat[0], &imat[nlocal], [](int x) {
    return x >= 0;
  });
  MPI_Allreduce(MPI_IN_PLACE, &nmat, 1, MPI_INT, MPI_SUM, world);
  double **x_ele;
  memory->create(x_ele, nmat, 3, "pppm/electrode:x_ele");
  memset(&(x_ele[0][0]), 0, nmat * 3 * sizeof(double));
  double **x = atom->x;
  for (int i = 0; i < nlocal; i++) {
    int ipos = imat[i];
    if (ipos < 0) continue;
    for (int dim = 0; dim < 3; dim++) x_ele[ipos][dim] = x[i][dim];
  }
  MPI_Allreduce(MPI_IN_PLACE, &(x_ele[0][0]), nmat * 3, MPI_DOUBLE, MPI_SUM, world);

  if (conp_one_step)
    one_step_multiplication(imat, greens_real, x_ele, matrix, nmat, timer_flag);
  else
    two_step_multiplication(imat, greens_real, x_ele, matrix, nmat, timer_flag);
  memory->destroy(greens_real);
  memory->destroy(x_ele);
}

/* ----------------------------------------------------------------------*/

void PPPMElectrodeIntel::one_step_multiplication(bigint *imat, double *greens_real, double **x_ele,
                                                 double **matrix, int const nmat, bool timer_flag)
{
  // map green's function in real space from mesh to particle positions
  // with matrix multiplication 'W^T G W' in one steps. Uses less memory than
  // two_step_multiplication
  //
  int const nlocal = atom->nlocal;
  double **x = atom->x;
  MPI_Barrier(world);
  double step1_time = MPI_Wtime();

  // precalculate rho_1d for local electrode
  std::vector<int> j_list;
  for (int j = 0; j < nlocal; j++) {
    int jpos = imat[j];
    if (jpos < 0) continue;
    j_list.push_back(j);
  }
  int const nj_local = j_list.size();

  FFT_SCALAR ***rho1d_j;
  memory->create(rho1d_j, nj_local, 3, order, "pppm/electrode:rho1d_j");

  _alignvar(FFT_SCALAR rho[3][INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};

  for (int jlist_pos = 0; jlist_pos < nj_local; jlist_pos++) {
    int j = j_list[jlist_pos];
    int njx = part2grid[j][0];
    int njy = part2grid[j][1];
    int njz = part2grid[j][2];
    FFT_SCALAR const djx = njx + shiftone - (x[j][0] - boxlo[0]) * delxinv;
    FFT_SCALAR const djy = njy + shiftone - (x[j][1] - boxlo[1]) * delyinv;
    FFT_SCALAR const djz = njz + shiftone - (x[j][2] - boxlo[2]) * delzinv;
    if (_use_table) {
      int idx = (int) (djx * half_rho_scale + half_rho_scale_plus);
      int idy = (int) (djy * half_rho_scale + half_rho_scale_plus);
      int idz = (int) (djz * half_rho_scale + half_rho_scale_plus);
#if defined(LMP_SIMD_COMPILER)
#pragma simd
#endif
      for (int k = 0; k < INTEL_P3M_ALIGNED_MAXORDER; k++) {
        rho[0][k] = rho_lookup[idx][k];
        rho[1][k] = rho_lookup[idy][k];
        rho[2][k] = rho_lookup[idz][k];
      }
    } else {
#if defined(LMP_SIMD_COMPILER)
#pragma simd
#endif
      for (int k = nlower; k <= nupper; k++) {
        FFT_SCALAR r1, r2, r3;
        r1 = r2 = r3 = ZEROF;

        for (int l = order - 1; l >= 0; l--) {
          r1 = rho_coeff[l][k] + r1 * djx;
          r2 = rho_coeff[l][k] + r2 * djy;
          r3 = rho_coeff[l][k] + r3 * djz;
        }
        rho[0][k - nlower] = r1;
        rho[1][k - nlower] = r2;
        rho[2][k - nlower] = r3;
      }
    }
    for (int dim = 0; dim < 3; dim++) {
      for (int oi = 0; oi < order; oi++) { rho1d_j[jlist_pos][dim][oi] = rho[dim][oi]; }
    }
  }

  // nested loops over weights of electrode atoms i and j
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  int const order2 = INTEL_P3M_ALIGNED_MAXORDER * INTEL_P3M_ALIGNED_MAXORDER;
  int const order6 = order2 * order2 * order2;
  _alignvar(double amesh[order6], 64) = {0};
  for (int ipos = 0; ipos < nmat; ipos++) {
    double *_noalias xi_ele = x_ele[ipos];
    // new calculation for nx, ny, nz because part2grid available for nlocal,
    // only
    int nix = static_cast<int>((xi_ele[0] - boxlo[0]) * delxinv + shift) - OFFSET;
    int niy = static_cast<int>((xi_ele[1] - boxlo[1]) * delyinv + shift) - OFFSET;
    int niz = static_cast<int>((xi_ele[2] - boxlo[2]) * delzinv + shift) - OFFSET;
    FFT_SCALAR dix = nix + shiftone - (xi_ele[0] - boxlo[0]) * delxinv;
    FFT_SCALAR diy = niy + shiftone - (xi_ele[1] - boxlo[1]) * delyinv;
    FFT_SCALAR diz = niz + shiftone - (xi_ele[2] - boxlo[2]) * delzinv;
    if (_use_table) {
      int idx = (int) (dix * half_rho_scale + half_rho_scale_plus);
      int idy = (int) (diy * half_rho_scale + half_rho_scale_plus);
      int idz = (int) (diz * half_rho_scale + half_rho_scale_plus);
#if defined(LMP_SIMD_COMPILER)
#pragma simd
#endif
      for (int k = 0; k < INTEL_P3M_ALIGNED_MAXORDER; k++) {
        rho[0][k] = rho_lookup[idx][k];
        rho[1][k] = rho_lookup[idy][k];
        rho[2][k] = rho_lookup[idz][k];
      }
    } else {
#if defined(LMP_SIMD_COMPILER)
#pragma simd
#endif
      for (int k = nlower; k <= nupper; k++) {
        FFT_SCALAR r1, r2, r3;
        r1 = r2 = r3 = ZEROF;

        for (int l = order - 1; l >= 0; l--) {
          r1 = rho_coeff[l][k] + r1 * dix;
          r2 = rho_coeff[l][k] + r2 * diy;
          r3 = rho_coeff[l][k] + r3 * diz;
        }
        rho[0][k - nlower] = r1;
        rho[1][k - nlower] = r2;
        rho[2][k - nlower] = r3;
      }
    }
    int njx = -1;
    int njy = -1;
    int njz = -1;    // force initial build_amesh
    for (int jlist_pos = 0; jlist_pos < nj_local; jlist_pos++) {
      int j = j_list[jlist_pos];
      int jpos = imat[j];
      if ((ipos < jpos) == !((ipos - jpos) % 2)) continue;
      double aij = 0.;
      if (njx != part2grid[j][0] || njy != part2grid[j][1] || njz != part2grid[j][2]) {
        njx = part2grid[j][0];
        njy = part2grid[j][1];
        njz = part2grid[j][2];
        build_amesh(njx - nix, njy - niy, njz - niz, amesh, greens_real);
      }
      int ind_amesh = 0;
      for (int ni = 0; ni < order; ni++) {
        FFT_SCALAR const iz0 = rho[2][ni];
        for (int nj = 0; nj < order; nj++) {
          FFT_SCALAR const jz0 = rho1d_j[jlist_pos][2][nj];
          for (int mi = 0; mi < order; mi++) {
            FFT_SCALAR const iy0 = iz0 * rho[1][mi];
            for (int mj = 0; mj < order; mj++) {
              FFT_SCALAR const jy0 = jz0 * rho1d_j[jlist_pos][1][mj];
              for (int li = 0; li < order; li++) {
                FFT_SCALAR const ix0 = iy0 * rho[0][li];
                double aij_xscan = 0.;
                for (int lj = 0; lj < order; lj++) {
                  aij_xscan += amesh[ind_amesh] * rho1d_j[jlist_pos][0][lj];
                  ind_amesh++;
                }
                aij += (double) ix0 * jy0 * aij_xscan;
              }
            }
          }
        }
      }
      matrix[ipos][jpos] += aij / volume;
      if (ipos != jpos) matrix[jpos][ipos] += aij / volume;
    }
  }
  MPI_Barrier(world);
  memory->destroy(rho1d_j);
  if (timer_flag && (comm->me == 0))
    utils::logmesg(lmp, "Single step time: {:.4g} s\n", MPI_Wtime() - step1_time);
}

/* ----------------------------------------------------------------------*/

void PPPMElectrodeIntel::build_amesh(const int dx,    // = njx - nix
                                     const int dy,    // = njy - niy
                                     const int dz,    // = njz - niz
                                     double *amesh, double *const greens_real)
{
  auto fmod = [](int x, int n) {    // fast unsigned mod
    int r = abs(x);
    while (r >= n) r -= n;
    return r;
  };
  int ind_amesh = 0;

  for (int iz = 0; iz < order; iz++)
    for (int jz = 0; jz < order; jz++) {
      int const mz = fmod(dz + jz - iz, nz_pppm) * nx_pppm * ny_pppm;
      for (int iy = 0; iy < order; iy++)
        for (int jy = 0; jy < order; jy++) {
          int const my = fmod(dy + jy - iy, ny_pppm) * nx_pppm;
          for (int ix = 0; ix < order; ix++)
            for (int jx = 0; jx < order; jx++) {
              int const mx = fmod(dx + jx - ix, nx_pppm);
              amesh[ind_amesh] = greens_real[mz + my + mx];
              ind_amesh++;
            }
        }
    }
}
/* ----------------------------------------------------------------------*/

void PPPMElectrodeIntel::two_step_multiplication(bigint *imat, double *greens_real, double **x_ele,
                                                 double **matrix, int const nmat, bool timer_flag)
{
  // map green's function in real space from mesh to particle positions
  // with matrix multiplication 'W^T G W' in two steps. gw is result of
  // first multiplication.
  int const nlocal = atom->nlocal;
  MPI_Barrier(world);
  double step1_time = MPI_Wtime();
  int nx_ele = nxhi_out - nxlo_out + 1;    // nx_pppm + order + 1;
  int ny_ele = nyhi_out - nylo_out + 1;    // ny_pppm + order + 1;
  int nz_ele = nzhi_out - nzlo_out + 1;    // nz_pppm + order + 1;
  int nxyz = nx_ele * ny_ele * nz_ele;
  vector<vector<double>> gw(nmat, vector<double>(nxyz, 0.));

  _alignvar(FFT_SCALAR rho[3][INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};

  // loops over weights of electrode atoms and weights of complete grid
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  for (int ipos = 0; ipos < nmat; ipos++) {
    double *_noalias xi_ele = x_ele[ipos];
    // new calculation for nx, ny, nz because part2grid available for
    // nlocal, only
    int nix = static_cast<int>((xi_ele[0] - boxlo[0]) * delxinv + shift) - OFFSET;
    int niy = static_cast<int>((xi_ele[1] - boxlo[1]) * delyinv + shift) - OFFSET;
    int niz = static_cast<int>((xi_ele[2] - boxlo[2]) * delzinv + shift) - OFFSET;
    FFT_SCALAR dx = nix + shiftone - (xi_ele[0] - boxlo[0]) * delxinv;
    FFT_SCALAR dy = niy + shiftone - (xi_ele[1] - boxlo[1]) * delyinv;
    FFT_SCALAR dz = niz + shiftone - (xi_ele[2] - boxlo[2]) * delzinv;
    if (_use_table) {
      dx = dx * half_rho_scale + half_rho_scale_plus;
      int idx = dx;
      dy = dy * half_rho_scale + half_rho_scale_plus;
      int idy = dy;
      dz = dz * half_rho_scale + half_rho_scale_plus;
      int idz = dz;
#if defined(LMP_SIMD_COMPILER)
#pragma simd
#endif
      for (int k = 0; k < INTEL_P3M_ALIGNED_MAXORDER; k++) {
        rho[0][k] = rho_lookup[idx][k];
        rho[1][k] = rho_lookup[idy][k];
        rho[2][k] = rho_lookup[idz][k];
      }
    } else {
#if defined(LMP_SIMD_COMPILER)
#pragma simd
#endif
      for (int k = nlower; k <= nupper; k++) {
        FFT_SCALAR r1, r2, r3;
        r1 = r2 = r3 = ZEROF;

        for (int l = order - 1; l >= 0; l--) {
          r1 = rho_coeff[l][k] + r1 * dx;
          r2 = rho_coeff[l][k] + r2 * dy;
          r3 = rho_coeff[l][k] + r3 * dz;
        }
        rho[0][k - nlower] = r1;
        rho[1][k - nlower] = r2;
        rho[2][k - nlower] = r3;
      }
    }
    for (int ni = nlower; ni <= nupper; ni++) {
      double iz0 = rho[2][ni - nlower];
      int miz = ni + niz;
      for (int mi = nlower; mi <= nupper; mi++) {
        double iy0 = iz0 * rho[1][mi - nlower];
        int miy = mi + niy;
        for (int li = nlower; li <= nupper; li++) {
          int mix = li + nix;
          double ix0 = iy0 * rho[0][li - nlower];
          for (int mjz = nzlo_out; mjz <= nzhi_out; mjz++) {
            int mz = abs(mjz - miz) % nz_pppm;
            for (int mjy = nylo_out; mjy <= nyhi_out; mjy++) {
              int my = abs(mjy - miy) % ny_pppm;
              for (int mjx = nxlo_out; mjx <= nxhi_out; mjx++) {
                int mx = abs(mjx - mix) % nx_pppm;
                gw[ipos][nx_ele * ny_ele * (mjz - nzlo_out) + nx_ele * (mjy - nylo_out) +
                         (mjx - nxlo_out)] +=
                    ix0 * greens_real[mz * nx_pppm * ny_pppm + my * nx_pppm + mx];
              }
            }
          }
        }
      }
    }
  }
  MPI_Barrier(world);
  if (timer_flag && (comm->me == 0))
    utils::logmesg(lmp, "step 1 time: {:.4g} s\n", MPI_Wtime() - step1_time);

  // nested loop over electrode atoms i and j and stencil of i
  double step2_time = MPI_Wtime();
  double **x = atom->x;
  for (int i = 0; i < nlocal; i++) {
    int ipos = imat[i];
    if (ipos < 0) continue;
    int nix = part2grid[i][0];
    int niy = part2grid[i][1];
    int niz = part2grid[i][2];
    FFT_SCALAR dix = nix + shiftone - (x[i][0] - boxlo[0]) * delxinv;
    FFT_SCALAR diy = niy + shiftone - (x[i][1] - boxlo[1]) * delyinv;
    FFT_SCALAR diz = niz + shiftone - (x[i][2] - boxlo[2]) * delzinv;
    if (_use_table) {
      dix = dix * half_rho_scale + half_rho_scale_plus;
      int idx = dix;
      diy = diy * half_rho_scale + half_rho_scale_plus;
      int idy = diy;
      diz = diz * half_rho_scale + half_rho_scale_plus;
      int idz = diz;
#if defined(LMP_SIMD_COMPILER)
#pragma simd
#endif
      for (int k = 0; k < INTEL_P3M_ALIGNED_MAXORDER; k++) {
        rho[0][k] = (double) rho_lookup[idx][k];
        rho[1][k] = (double) rho_lookup[idy][k];
        rho[2][k] = (double) rho_lookup[idz][k];
      }
    } else {
#if defined(LMP_SIMD_COMPILER)
#pragma simd
#endif
      for (int k = nlower; k <= nupper; k++) {
        FFT_SCALAR r1, r2, r3;
        r1 = r2 = r3 = ZEROF;

        for (int l = order - 1; l >= 0; l--) {
          r1 = rho_coeff[l][k] + r1 * dix;
          r2 = rho_coeff[l][k] + r2 * diy;
          r3 = rho_coeff[l][k] + r3 * diz;
        }
        rho[0][k - nlower] = (double) r1;
        rho[1][k - nlower] = (double) r2;
        rho[2][k - nlower] = (double) r3;
      }
    }
    for (int jpos = 0; jpos < nmat; jpos++) {
      double aij = 0.;
      for (int ni = nlower; ni <= nupper; ni++) {
        double iz0 = rho[2][ni - nlower];
        int miz = ni + niz;
        for (int mi = nlower; mi <= nupper; mi++) {
          double iy0 = iz0 * rho[1][mi - nlower];
          int miy = mi + niy;
          for (int li = nlower; li <= nupper; li++) {
            int mix = li + nix;
            double ix0 = iy0 * rho[0][li - nlower];
            int miz0 = miz - nzlo_out;
            int miy0 = miy - nylo_out;
            int mix0 = mix - nxlo_out;
            aij += ix0 * gw[jpos][nx_ele * ny_ele * miz0 + nx_ele * miy0 + mix0];
          }
        }
      }
      matrix[ipos][jpos] += aij / volume;
    }
  }
  MPI_Barrier(world);
  if (timer_flag && (comm->me == 0))
    utils::logmesg(lmp, "step 2 time: {:.4g} s\n", MPI_Wtime() - step2_time);
}

template <class flt_t, class acc_t, int use_table>
void PPPMElectrodeIntel::make_rho_in_brick(IntelBuffers<flt_t, acc_t> *buffers, int source_grpbit,
                                           FFT_SCALAR ***scratch_brick, bool invert_source)
{

  FFT_SCALAR *_noalias global_density = &(scratch_brick[nzlo_out][nylo_out][nxlo_out]);

  ATOM_T *_noalias const x = buffers->get_x(0);
  flt_t *_noalias const q = buffers->get_q(0);
  int nlocal = atom->nlocal;
  int nthr;
  if (_use_lrt)
    nthr = 1;
  else
    nthr = comm->nthreads;

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(nthr, nlocal, global_density, source_grpbit, \
                                                 invert_source) if (!_use_lrt)
#endif
  {
    int *mask = atom->mask;
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
    FFT_SCALAR *_noalias my_density = tid == 0 ? global_density : perthread_density[tid - 1];
    // clear 3d density array
    memset(my_density, 0, ngrid * sizeof(FFT_SCALAR));

    for (int i = ifrom; i < ito; i++) {
      bool const i_in_source = !!(mask[i] & source_grpbit) ^ invert_source;
      if (!i_in_source) continue;

      int nx = part2grid[i][0];
      int ny = part2grid[i][1];
      int nz = part2grid[i][2];

      int nysum = nlower + ny - nylo_out;
      int nxsum = nlower + nx - nxlo_out;
      int nzsum = (nlower + nz - nzlo_out) * nix * niy + nysum * nix + nxsum;

      FFT_SCALAR dx = nx + fshiftone - (x[i].x - lo0) * xi;
      FFT_SCALAR dy = ny + fshiftone - (x[i].y - lo1) * yi;
      FFT_SCALAR dz = nz + fshiftone - (x[i].z - lo2) * zi;

      _alignvar(flt_t rho[3][INTEL_P3M_ALIGNED_MAXORDER], 64) = {0};

      if (use_table) {
        dx = dx * half_rho_scale + half_rho_scale_plus;
        int idx = dx;
        dy = dy * half_rho_scale + half_rho_scale_plus;
        int idy = dy;
        dz = dz * half_rho_scale + half_rho_scale_plus;
        int idz = dz;
#if defined(LMP_SIMD_COMPILER)
#pragma simd
#endif
        for (int k = 0; k < INTEL_P3M_ALIGNED_MAXORDER; k++) {
          rho[0][k] = rho_lookup[idx][k];
          rho[1][k] = rho_lookup[idy][k];
          rho[2][k] = rho_lookup[idz][k];
        }
      } else {
#if defined(LMP_SIMD_COMPILER)
#pragma simd
#endif
        for (int k = nlower; k <= nupper; k++) {
          FFT_SCALAR r1, r2, r3;
          r1 = r2 = r3 = ZEROF;

          for (int l = order - 1; l >= 0; l--) {
            r1 = rho_coeff[l][k] + r1 * dx;
            r2 = rho_coeff[l][k] + r2 * dy;
            r3 = rho_coeff[l][k] + r3 * dz;
          }
          rho[0][k - nlower] = r1;
          rho[1][k - nlower] = r2;
          rho[2][k - nlower] = r3;
        }
      }

      FFT_SCALAR z0 = fdelvolinv * q[i];

#if defined(LMP_SIMD_COMPILER)
#pragma loop_count min(2), max(INTEL_P3M_ALIGNED_MAXORDER), avg(7)
#endif
      for (int n = 0; n < order; n++) {
        int mz = n * nix * niy + nzsum;
        FFT_SCALAR y0 = z0 * rho[2][n];
#if defined(LMP_SIMD_COMPILER)
#pragma loop_count min(2), max(INTEL_P3M_ALIGNED_MAXORDER), avg(7)
#endif
        for (int m = 0; m < order; m++) {
          int mzy = m * nix + mz;
          FFT_SCALAR x0 = y0 * rho[1][m];
#if defined(LMP_SIMD_COMPILER)
#pragma simd
#endif
          for (int l = 0; l < INTEL_P3M_ALIGNED_MAXORDER; l++) {
            int mzyx = l + mzy;
            my_density[mzyx] += x0 * rho[0][l];
          }
        }
      }
    }
  }

  // reduce all the perthread_densities into global_density
  if (nthr > 1) {
#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(nthr, global_density) if (!_use_lrt)
#endif
    {
      int ifrom, ito, tid;
      IP_PRE_omp_range_id(ifrom, ito, tid, ngrid, nthr);

#if defined(LMP_SIMD_COMPILER)
#pragma simd
#endif
      for (int i = ifrom; i < ito; i++) {
        for (int j = 1; j < nthr; j++) { global_density[i] += perthread_density[j - 1][i]; }
      }
    }
  }
}

void PPPMElectrodeIntel::compute_group_group(int /*groupbit_A*/, int /*groupbit_B*/,
                                             int /*AA_flag*/)
{
  error->all(FLERR, "group group interaction not implemented in pppm/electrode yet");
}

void PPPMElectrodeIntel::compute_matrix_corr(bigint *imat, double **matrix)
{
  boundcorr->matrix_corr(imat, matrix);
}

void PPPMElectrodeIntel::compute_vector_corr(double *vec, int sensor_grpbit, int source_grpbit,
                                             bool invert_source)
{
  boundcorr->vector_corr(vec, sensor_grpbit, source_grpbit, invert_source);
}

void PPPMElectrodeIntel::allocate()
{
  if (slabflag == 1) {
    // EW3Dc dipole correction
    boundcorr = new SlabDipole(lmp);
  } else if (wireflag == 1) {
    // EW3Dc wire correction
    boundcorr = new WireDipole(lmp);
  } else {
    // base BoundaryCorrection -- used for ffield
    boundcorr = new BoundaryCorrection(lmp);
  }

  PPPM::allocate();
  /* ----------------------------------------------------------------------
     Allocate density_brick with extra padding for vector writes
  ------------------------------------------------------------------------- */
  PPPMIntel::create3d_offset(electrolyte_density_brick, nzlo_out, nzhi_out, nylo_out, nyhi_out,
                             nxlo_out, nxhi_out, "pppm/electrode:electrolyte_density_brick");
  memory->create(electrolyte_density_fft, nfft_both, "pppm/electrode:electrolyte_density_fft");
  memory->destroy3d_offset(density_brick, nzlo_out, nylo_out, nxlo_out);
  PPPMIntel::create3d_offset(density_brick, nzlo_out, nzhi_out, nylo_out, nyhi_out, nxlo_out,
                             nxhi_out, "pppm:density_brick");

  if ((differentiation_flag == 1) || u_brick != nullptr) {
    memory->destroy3d_offset(u_brick, nzlo_out, nylo_out, nxlo_out);
  }
  if (differentiation_flag != 1) {
    memory->destroy3d_offset(vdx_brick, nzlo_out, nylo_out, nxlo_out);
    memory->destroy3d_offset(vdy_brick, nzlo_out, nylo_out, nxlo_out);
    memory->destroy3d_offset(vdz_brick, nzlo_out, nylo_out, nxlo_out);
    PPPMIntel::create3d_offset(vdx_brick, nzlo_out, nzhi_out, nylo_out, nyhi_out, nxlo_out,
                               nxhi_out, "pppm:vdx_brick");
    PPPMIntel::create3d_offset(vdy_brick, nzlo_out, nzhi_out, nylo_out, nyhi_out, nxlo_out,
                               nxhi_out, "pppm:vdy_brick");
    PPPMIntel::create3d_offset(vdz_brick, nzlo_out, nzhi_out, nylo_out, nyhi_out, nxlo_out,
                               nxhi_out, "pppm:vdz_brick");
  }
  PPPMIntel::create3d_offset(u_brick, nzlo_out, nzhi_out, nylo_out, nyhi_out, nxlo_out, nxhi_out,
                             "pppm:u_brick");    // use u_brick in project_psi
}

void PPPMElectrodeIntel::allocate_peratom()
{
  // duplicated to avoid reallocating u_brick
  peratom_allocate_flag = 1;

  memory->create3d_offset(v0_brick, nzlo_out, nzhi_out, nylo_out, nyhi_out, nxlo_out, nxhi_out,
                          "pppm:v0_brick");

  memory->create3d_offset(v1_brick, nzlo_out, nzhi_out, nylo_out, nyhi_out, nxlo_out, nxhi_out,
                          "pppm:v1_brick");
  memory->create3d_offset(v2_brick, nzlo_out, nzhi_out, nylo_out, nyhi_out, nxlo_out, nxhi_out,
                          "pppm:v2_brick");
  memory->create3d_offset(v3_brick, nzlo_out, nzhi_out, nylo_out, nyhi_out, nxlo_out, nxhi_out,
                          "pppm:v3_brick");
  memory->create3d_offset(v4_brick, nzlo_out, nzhi_out, nylo_out, nyhi_out, nxlo_out, nxhi_out,
                          "pppm:v4_brick");
  memory->create3d_offset(v5_brick, nzlo_out, nzhi_out, nylo_out, nyhi_out, nxlo_out, nxhi_out,
                          "pppm:v5_brick");

  // use same GC ghost grid object for peratom grid communication
  // but need to reallocate a larger gc_buf1 and gc_buf2

  if (differentiation_flag)
    npergrid = 6;
  else
    npergrid = 7;

  memory->destroy(gc_buf1);
  memory->destroy(gc_buf2);
  memory->create(gc_buf1, npergrid * ngc_buf1, "pppm:gc_buf1");
  memory->create(gc_buf2, npergrid * ngc_buf2, "pppm:gc_buf2");
}

void PPPMElectrodeIntel::deallocate()
{
  if (boundcorr != nullptr) delete boundcorr;
  // duplicated to always deallocate u_brick
  memory->destroy3d_offset(density_brick, nzlo_out, nylo_out, nxlo_out);
  memory->destroy3d_offset(u_brick, nzlo_out, nylo_out, nxlo_out);

  if (differentiation_flag == 1) {
    memory->destroy(sf_precoeff1);
    memory->destroy(sf_precoeff2);
    memory->destroy(sf_precoeff3);
    memory->destroy(sf_precoeff4);
    memory->destroy(sf_precoeff5);
    memory->destroy(sf_precoeff6);
  } else {
    memory->destroy3d_offset(vdx_brick, nzlo_out, nylo_out, nxlo_out);
    memory->destroy3d_offset(vdy_brick, nzlo_out, nylo_out, nxlo_out);
    memory->destroy3d_offset(vdz_brick, nzlo_out, nylo_out, nxlo_out);
  }

  memory->destroy(density_fft);
  memory->destroy(greensfn);
  memory->destroy(work1);
  memory->destroy(work2);
  memory->destroy(vg);

  if (triclinic == 0) {
    memory->destroy1d_offset(fkx, nxlo_fft);
    memory->destroy1d_offset(fky, nylo_fft);
    memory->destroy1d_offset(fkz, nzlo_fft);
  } else {
    memory->destroy(fkx);
    memory->destroy(fky);
    memory->destroy(fkz);
  }

  memory->destroy(gf_b);
  if (stagger_flag) gf_b = nullptr;
  memory->destroy2d_offset(rho1d, -order_allocated / 2);
  memory->destroy2d_offset(drho1d, -order_allocated / 2);
  memory->destroy2d_offset(rho_coeff, (1 - order_allocated) / 2);
  memory->destroy2d_offset(drho_coeff, (1 - order_allocated) / 2);

  delete fft1;
  delete fft2;
  delete remap;
  delete gc;
  memory->destroy(gc_buf1);
  memory->destroy(gc_buf2);
}

/* ----------------------------------------------------------------------
  Pack charge data into intel package buffers after updates
------------------------------------------------------------------------- */

void PPPMElectrodeIntel::pack_buffers_q()
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
      fix->get_mixed_buffers()->thr_pack_q(ifrom,ito);
    else if (fix->precision() == FixIntel::PREC_MODE_DOUBLE)
      fix->get_double_buffers()->thr_pack_q(ifrom,ito);
    else
      fix->get_single_buffers()->thr_pack_q(ifrom,ito);
  }
  fix->stop_watch(TIME_PACK);
}
