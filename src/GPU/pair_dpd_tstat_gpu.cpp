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
   Contributing author: Trung Dac Nguyen (ORNL)
------------------------------------------------------------------------- */

#include "pair_dpd_tstat_gpu.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "gpu_extra.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "suffix.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;

// External functions from cuda library for atom decomposition

int dpd_tstat_gpu_init(const int ntypes, double **cutsq, double **host_a0, double **host_gamma,
                       double **host_sigma, double **host_cut, double *special_lj, const int inum,
                       const int nall, const int max_nbors, const int maxspecial,
                       const double cell_size, int &gpu_mode, FILE *screen);
void dpd_tstat_gpu_clear();
int **dpd_tstat_gpu_compute_n(const int ago, const int inum_full, const int nall, double **host_x,
                              int *host_type, double *sublo, double *subhi, tagint *tag,
                              int **nspecial, tagint **special, const bool eflag, const bool vflag,
                              const bool eatom, const bool vatom, int &host_start, int **ilist,
                              int **jnum, const double cpu_time, bool &success, double **host_v,
                              const double dtinvsqrt, const int seed, const int timestep,
                              double *boxlo, double *prd);
void dpd_tstat_gpu_compute(const int ago, const int inum_full, const int nall, double **host_x,
                           int *host_type, int *ilist, int *numj, int **firstneigh,
                           const bool eflag, const bool vflag, const bool eatom, const bool vatom,
                           int &host_start, const double cpu_time, bool &success, tagint *tag,
                           double **host_v, const double dtinvsqrt, const int seed,
                           const int timestep, const int nlocal, double *boxlo, double *prd);
void dpd_tstat_gpu_update_coeff(int ntypes, double **host_a0, double **host_gamma,
                                double **host_sigma, double **host_cut);
double dpd_tstat_gpu_bytes();

#define EPSILON 1.0e-10

//#define _USE_UNIFORM_SARU_LCG
//#define _USE_UNIFORM_SARU_TEA8
//#define _USE_GAUSSIAN_SARU_LCG

#if !defined(_USE_UNIFORM_SARU_LCG) && !defined(_USE_UNIFORM_SARU_TEA8) && \
    !defined(_USE_GAUSSIAN_SARU_LCG)
#define _USE_UNIFORM_SARU_LCG
#endif

// References:
// 1. Y. Afshar, F. Schmid, A. Pishevar, S. Worley, Comput. Phys. Comm. 184 (2013), 1119–1128.
// 2. C. L. Phillips, J. A. Anderson, S. C. Glotzer, Comput. Phys. Comm. 230 (2011), 7191-7201.
// PRNG period = 3666320093*2^32 ~ 2^64 ~ 10^19

#define LCGA 0x4beb5d59    // Full period 32 bit LCG
#define LCGC 0x2600e1f7
#define oWeylPeriod 0xda879add    // Prime period 3666320093
#define oWeylOffset 0x8009d14b
#define TWO_N32 0.232830643653869628906250e-9f /* 2^-32 */

// specifically implemented for steps = 1; high = 1.0; low = -1.0
// returns uniformly distributed random numbers u in [-1.0;1.0]
// using the inherent LCG, then multiply u with sqrt(3) to "match"
// with a normal random distribution.
// Afshar et al. mutlplies u in [-0.5;0.5] with sqrt(12)
// Curly brackets to make variables local to the scope.
#ifdef _USE_UNIFORM_SARU_LCG
#define numtyp double
#define SQRT3 (numtyp) 1.7320508075688772935274463
#define saru(seed1, seed2, seed, timestep, randnum)                                \
  {                                                                                \
    unsigned int seed3 = seed + timestep;                                          \
    seed3 ^= (seed1 << 7) ^ (seed2 >> 6);                                          \
    seed2 += (seed1 >> 4) ^ (seed3 >> 15);                                         \
    seed1 ^= (seed2 << 9) + (seed3 << 8);                                          \
    seed3 ^= 0xA5366B4D * ((seed2 >> 11) ^ (seed1 << 1));                          \
    seed2 += 0x72BE1579 * ((seed1 << 4) ^ (seed3 >> 16));                          \
    seed1 ^= 0x3F38A6ED * ((seed3 >> 5) ^ (((signed int) seed2) >> 22));           \
    seed2 += seed1 * seed3;                                                        \
    seed1 += seed3 ^ (seed2 >> 2);                                                 \
    seed2 ^= ((signed int) seed2) >> 17;                                           \
    unsigned int state = 0x79dedea3 * (seed1 ^ (((signed int) seed1) >> 14));      \
    unsigned int wstate = (state + seed2) ^ (((signed int) state) >> 8);           \
    state = state + (wstate * (wstate ^ 0xdddf97f5));                              \
    wstate = 0xABCB96F7 + (wstate >> 1);                                           \
    state = LCGA * state + LCGC;                                                   \
    wstate = wstate + oWeylOffset + ((((signed int) wstate) >> 31) & oWeylPeriod); \
    unsigned int v = (state ^ (state >> 26)) + wstate;                             \
    unsigned int s = (signed int) ((v ^ (v >> 20)) * 0x6957f5a7);                  \
    randnum = SQRT3 * (s * TWO_N32 * (numtyp) 2.0 - (numtyp) 1.0);                 \
  }
#endif

// specifically implemented for steps = 1; high = 1.0; low = -1.0
// returns uniformly distributed random numbers u in [-1.0;1.0] using TEA8
// then multiply u with sqrt(3) to "match" with a normal random distribution
// Afshar et al. mutlplies u in [-0.5;0.5] with sqrt(12)
#ifdef _USE_UNIFORM_SARU_TEA8
#define numtyp double
#define SQRT3 (numtyp) 1.7320508075688772935274463
#define k0 0xA341316C
#define k1 0xC8013EA4
#define k2 0xAD90777D
#define k3 0x7E95761E
#define delta 0x9e3779b9
#define rounds 8
#define saru(seed1, seed2, seed, timestep, randnum)                           \
  {                                                                           \
    unsigned int seed3 = seed + timestep;                                     \
    seed3 ^= (seed1 << 7) ^ (seed2 >> 6);                                     \
    seed2 += (seed1 >> 4) ^ (seed3 >> 15);                                    \
    seed1 ^= (seed2 << 9) + (seed3 << 8);                                     \
    seed3 ^= 0xA5366B4D * ((seed2 >> 11) ^ (seed1 << 1));                     \
    seed2 += 0x72BE1579 * ((seed1 << 4) ^ (seed3 >> 16));                     \
    seed1 ^= 0x3F38A6ED * ((seed3 >> 5) ^ (((signed int) seed2) >> 22));      \
    seed2 += seed1 * seed3;                                                   \
    seed1 += seed3 ^ (seed2 >> 2);                                            \
    seed2 ^= ((signed int) seed2) >> 17;                                      \
    unsigned int state = 0x79dedea3 * (seed1 ^ (((signed int) seed1) >> 14)); \
    unsigned int wstate = (state + seed2) ^ (((signed int) state) >> 8);      \
    state = state + (wstate * (wstate ^ 0xdddf97f5));                         \
    wstate = 0xABCB96F7 + (wstate >> 1);                                      \
    unsigned int sum = 0;                                                     \
    for (int i = 0; i < rounds; i++) {                                        \
      sum += delta;                                                           \
      state += ((wstate << 4) + k0) ^ (wstate + sum) ^ ((wstate >> 5) + k1);  \
      wstate += ((state << 4) + k2) ^ (state + sum) ^ ((state >> 5) + k3);    \
    }                                                                         \
    unsigned int v = (state ^ (state >> 26)) + wstate;                        \
    unsigned int s = (signed int) ((v ^ (v >> 20)) * 0x6957f5a7);             \
    randnum = SQRT3 * (s * TWO_N32 * (numtyp) 2.0 - (numtyp) 1.0);            \
  }
#endif

// specifically implemented for steps = 1; high = 1.0; low = -1.0
// returns two uniformly distributed random numbers r1 and r2 in [-1.0;1.0],
// and uses the polar method (Marsaglia's) to transform to a normal random value
// This is used to compared with CPU DPD using RandMars::gaussian()
#ifdef _USE_GAUSSIAN_SARU_LCG
#define numtyp double
#define saru(seed1, seed2, seed, timestep, randnum)                                  \
  {                                                                                  \
    unsigned int seed3 = seed + timestep;                                            \
    seed3 ^= (seed1 << 7) ^ (seed2 >> 6);                                            \
    seed2 += (seed1 >> 4) ^ (seed3 >> 15);                                           \
    seed1 ^= (seed2 << 9) + (seed3 << 8);                                            \
    seed3 ^= 0xA5366B4D * ((seed2 >> 11) ^ (seed1 << 1));                            \
    seed2 += 0x72BE1579 * ((seed1 << 4) ^ (seed3 >> 16));                            \
    seed1 ^= 0x3F38A6ED * ((seed3 >> 5) ^ (((signed int) seed2) >> 22));             \
    seed2 += seed1 * seed3;                                                          \
    seed1 += seed3 ^ (seed2 >> 2);                                                   \
    seed2 ^= ((signed int) seed2) >> 17;                                             \
    unsigned int state = 0x12345678;                                                 \
    unsigned int wstate = 12345678;                                                  \
    state = 0x79dedea3 * (seed1 ^ (((signed int) seed1) >> 14));                     \
    wstate = (state + seed2) ^ (((signed int) state) >> 8);                          \
    state = state + (wstate * (wstate ^ 0xdddf97f5));                                \
    wstate = 0xABCB96F7 + (wstate >> 1);                                             \
    unsigned int v, s;                                                               \
    numtyp r1, r2, rsq;                                                              \
    while (1) {                                                                      \
      state = LCGA * state + LCGC;                                                   \
      wstate = wstate + oWeylOffset + ((((signed int) wstate) >> 31) & oWeylPeriod); \
      v = (state ^ (state >> 26)) + wstate;                                          \
      s = (signed int) ((v ^ (v >> 20)) * 0x6957f5a7);                               \
      r1 = s * TWO_N32 * (numtyp) 2.0 - (numtyp) 1.0;                                \
      state = LCGA * state + LCGC;                                                   \
      wstate = wstate + oWeylOffset + ((((signed int) wstate) >> 31) & oWeylPeriod); \
      v = (state ^ (state >> 26)) + wstate;                                          \
      s = (signed int) ((v ^ (v >> 20)) * 0x6957f5a7);                               \
      r2 = s * TWO_N32 * (numtyp) 2.0 - (numtyp) 1.0;                                \
      rsq = r1 * r1 + r2 * r2;                                                       \
      if (rsq < (numtyp) 1.0) break;                                                 \
    }                                                                                \
    numtyp fac = sqrt((numtyp) -2.0 * log(rsq) / rsq);                               \
    randnum = r2 * fac;                                                              \
  }
#endif

/* ---------------------------------------------------------------------- */

PairDPDTstatGPU::PairDPDTstatGPU(LAMMPS *lmp) : PairDPDTstat(lmp), gpu_mode(GPU_FORCE)
{
  respa_enable = 0;
  reinitflag = 0;
  cpu_time = 0.0;
  suffix_flag |= Suffix::GPU;
  GPU_EXTRA::gpu_ready(lmp->modify, lmp->error);
}

/* ----------------------------------------------------------------------
   free all arrays
------------------------------------------------------------------------- */

PairDPDTstatGPU::~PairDPDTstatGPU()
{
  dpd_tstat_gpu_clear();
}

/* ---------------------------------------------------------------------- */

void PairDPDTstatGPU::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  // adjust sigma if target T is changing

  if (t_start != t_stop) {
    double delta = update->ntimestep - update->beginstep;
    if (delta != 0.0) delta /= update->endstep - update->beginstep;
    temperature = t_start + delta * (t_stop - t_start);
    double boltz = force->boltz;
    for (int i = 1; i <= atom->ntypes; i++)
      for (int j = i; j <= atom->ntypes; j++)
        sigma[i][j] = sigma[j][i] = sqrt(2.0 * boltz * temperature * gamma[i][j]);

    dpd_tstat_gpu_update_coeff(atom->ntypes + 1, a0, gamma, sigma, cut);
  }

  int nall = atom->nlocal + atom->nghost;
  int inum, host_start;

  double dtinvsqrt = 1.0 / sqrt(update->dt);

  bool success = true;
  int *ilist, *numneigh, **firstneigh;
  if (gpu_mode != GPU_FORCE) {
    double sublo[3], subhi[3];
    if (domain->triclinic == 0) {
      sublo[0] = domain->sublo[0];
      sublo[1] = domain->sublo[1];
      sublo[2] = domain->sublo[2];
      subhi[0] = domain->subhi[0];
      subhi[1] = domain->subhi[1];
      subhi[2] = domain->subhi[2];
    } else {
      domain->bbox(domain->sublo_lamda, domain->subhi_lamda, sublo, subhi);
    }
    inum = atom->nlocal;
    firstneigh = dpd_tstat_gpu_compute_n(
        neighbor->ago, inum, nall, atom->x, atom->type, sublo, subhi, atom->tag, atom->nspecial,
        atom->special, eflag, vflag, eflag_atom, vflag_atom, host_start, &ilist, &numneigh,
        cpu_time, success, atom->v, dtinvsqrt, seed, update->ntimestep, domain->boxlo, domain->prd);
  } else {
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;
    dpd_tstat_gpu_compute(neighbor->ago, inum, nall, atom->x, atom->type, ilist, numneigh,
                          firstneigh, eflag, vflag, eflag_atom, vflag_atom, host_start, cpu_time,
                          success, atom->tag, atom->v, dtinvsqrt, seed, update->ntimestep,
                          atom->nlocal, domain->boxlo, domain->prd);
  }
  if (!success) error->one(FLERR, "Insufficient memory on accelerator");

  if (atom->molecular != Atom::ATOMIC && neighbor->ago == 0)
    neighbor->build_topology();
  if (host_start < inum) {
    cpu_time = platform::walltime();
    cpu_compute(host_start, inum, eflag, vflag, ilist, numneigh, firstneigh);
    cpu_time = platform::walltime() - cpu_time;
  }
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairDPDTstatGPU::init_style()
{

  // Repeat cutsq calculation because done after call to init_style
  double maxcut = -1.0;
  double mcut;
  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = i; j <= atom->ntypes; j++) {
      if (setflag[i][j] != 0 || (setflag[i][i] != 0 && setflag[j][j] != 0)) {
        mcut = init_one(i, j);
        mcut *= mcut;
        if (mcut > maxcut) maxcut = mcut;
        cutsq[i][j] = cutsq[j][i] = mcut;
      } else
        cutsq[i][j] = cutsq[j][i] = 0.0;
    }
  }
  double cell_size = sqrt(maxcut) + neighbor->skin;

  int maxspecial = 0;
  if (atom->molecular != Atom::ATOMIC) maxspecial = atom->maxspecial;
  int mnf = 5e-2 * neighbor->oneatom;
  int success = dpd_tstat_gpu_init(atom->ntypes + 1, cutsq, a0, gamma, sigma, cut,
                                   force->special_lj, atom->nlocal, atom->nlocal + atom->nghost,
                                   mnf, maxspecial, cell_size, gpu_mode, screen);
  GPU_EXTRA::check_flag(success, error, world);

  if (gpu_mode == GPU_FORCE) neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ---------------------------------------------------------------------- */

double PairDPDTstatGPU::memory_usage()
{
  double bytes = Pair::memory_usage();
  return bytes + dpd_tstat_gpu_bytes();
}

/* ---------------------------------------------------------------------- */

void PairDPDTstatGPU::cpu_compute(int start, int inum, int /* eflag */, int /* vflag */, int *ilist,
                                  int *numneigh, int **firstneigh)
{
  int i, j, ii, jj, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, fpair;
  double vxtmp, vytmp, vztmp, delvx, delvy, delvz;
  double rsq, r, rinv, dot, wd, randnum, factor_dpd, factor_sqrt;
  int *jlist;
  tagint itag, jtag;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *type = atom->type;
  tagint *tag = atom->tag;
  double *special_lj = force->special_lj;
  double dtinvsqrt = 1.0 / sqrt(update->dt);
  int timestep = (int) update->ntimestep;

  // loop over neighbors of my atoms

  for (ii = start; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    vxtmp = v[i][0];
    vytmp = v[i][1];
    vztmp = v[i][2];
    itype = type[i];
    itag = tag[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_dpd = special_lj[sbmask(j)];
      factor_sqrt = special_sqrt[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];
      jtag = tag[j];

      if (rsq < cutsq[itype][jtype]) {
        r = sqrt(rsq);
        if (r < EPSILON) continue;    // r can be 0.0 in DPD systems
        rinv = 1.0 / r;
        delvx = vxtmp - v[j][0];
        delvy = vytmp - v[j][1];
        delvz = vztmp - v[j][2];
        dot = delx * delvx + dely * delvy + delz * delvz;
        wd = 1.0 - r / cut[itype][jtype];

        unsigned int tag1 = itag, tag2 = jtag;
        if (tag1 > tag2) {
          tag1 = jtag;
          tag2 = itag;
        }

        randnum = 0.0;
        saru(tag1, tag2, seed, timestep, randnum);

        // conservative force = a0 * wd
        // drag force = -gamma * wd^2 * (delx dot delv) / r
        // random force = sigma * wd * rnd * dtinvsqrt;

        fpair = -factor_dpd * gamma[itype][jtype] * wd * wd * dot * rinv;
        fpair += factor_sqrt * sigma[itype][jtype] * wd * randnum * dtinvsqrt;
        fpair *= rinv;

        f[i][0] += delx * fpair;
        f[i][1] += dely * fpair;
        f[i][2] += delz * fpair;

        if (evflag) ev_tally_full(i, 0.0, 0.0, fpair, delx, dely, delz);
      }
    }
  }
}
