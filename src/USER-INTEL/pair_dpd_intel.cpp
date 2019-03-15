/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: W. Michael Brown (Intel)
                        Shun Xu (Computer Network Information Center, CAS)
------------------------------------------------------------------------- */

#include <cmath>
#include "pair_dpd_intel.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "suffix.h"
using namespace LAMMPS_NS;

#define LMP_MKL_RNG VSL_BRNG_MT19937
#define FC_PACKED1_T typename ForceConst<flt_t>::fc_packed1
#define IEPSILON 1.0e10

/* ---------------------------------------------------------------------- */

PairDPDIntel::PairDPDIntel(LAMMPS *lmp) :
  PairDPD(lmp)
{
  suffix_flag |= Suffix::INTEL;
  respa_enable = 0;
  random_thread = NULL;
  _nrandom_thread = 0;
}

/* ---------------------------------------------------------------------- */

PairDPDIntel::~PairDPDIntel()
{
  #if defined(_OPENMP)
  if (_nrandom_thread) {
    #ifdef LMP_USE_MKL_RNG
    for (int i = 0; i < _nrandom_thread; i++)
      vslDeleteStream(&random_thread[i]);
    #else
    for (int i = 1; i < _nrandom_thread; i++)
      delete random_thread[i];
    #endif
  }
  #endif
  delete []random_thread;
}

/* ---------------------------------------------------------------------- */

void PairDPDIntel::compute(int eflag, int vflag)
{
  if (fix->precision() == FixIntel::PREC_MODE_MIXED)
    compute<float,double>(eflag, vflag, fix->get_mixed_buffers(),
                          force_const_single);
  else if (fix->precision() == FixIntel::PREC_MODE_DOUBLE)
    compute<double,double>(eflag, vflag, fix->get_double_buffers(),
                           force_const_double);
  else
    compute<float,float>(eflag, vflag, fix->get_single_buffers(),
                         force_const_single);

  fix->balance_stamp();
  vflag_fdotr = 0;
}

template <class flt_t, class acc_t>
void PairDPDIntel::compute(int eflag, int vflag,
                           IntelBuffers<flt_t,acc_t> *buffers,
                           const ForceConst<flt_t> &fc)
{
  if (eflag || vflag) {
    ev_setup(eflag, vflag);
  } else evflag = vflag_fdotr = 0;

  const int inum = list->inum;
  const int nthreads = comm->nthreads;
  const int host_start = fix->host_start_pair();
  const int offload_end = fix->offload_end_pair();
  const int ago = neighbor->ago;

  if (ago != 0 && fix->separate_buffers() == 0) {
    fix->start_watch(TIME_PACK);

    int packthreads;
    if (nthreads > INTEL_HTHREADS) packthreads = nthreads;
    else packthreads = 1;
    #if defined(_OPENMP)
    #pragma omp parallel if(packthreads > 1)
    #endif
    {
      int ifrom, ito, tid;
      IP_PRE_omp_range_id_align(ifrom, ito, tid, atom->nlocal + atom->nghost,
                                packthreads, sizeof(ATOM_T));
      buffers->thr_pack(ifrom,ito,ago);
    }
    fix->stop_watch(TIME_PACK);
  }

  int ovflag = 0;
  if (vflag_fdotr) ovflag = 2;
  else if (vflag) ovflag = 1;
  if (_onetype) {
    if (eflag) {
      if (force->newton_pair) {
        eval<1,1,1>(1, ovflag, buffers, fc, 0, offload_end);
        eval<1,1,1>(0, ovflag, buffers, fc, host_start, inum);
      } else {
        eval<1,1,0>(1, ovflag, buffers, fc, 0, offload_end);
        eval<1,1,0>(0, ovflag, buffers, fc, host_start, inum);
      }
    } else {
      if (force->newton_pair) {
        eval<1,0,1>(1, ovflag, buffers, fc, 0, offload_end);
        eval<1,0,1>(0, ovflag, buffers, fc, host_start, inum);
      } else {
        eval<1,0,0>(1, ovflag, buffers, fc, 0, offload_end);
        eval<1,0,0>(0, ovflag, buffers, fc, host_start, inum);
      }
    }
  } else {
    if (eflag) {
      if (force->newton_pair) {
        eval<0,1,1>(1, ovflag, buffers, fc, 0, offload_end);
        eval<0,1,1>(0, ovflag, buffers, fc, host_start, inum);
      } else {
        eval<0,1,0>(1, ovflag, buffers, fc, 0, offload_end);
        eval<0,1,0>(0, ovflag, buffers, fc, host_start, inum);
      }
    } else {
      if (force->newton_pair) {
        eval<0,0,1>(1, ovflag, buffers, fc, 0, offload_end);
        eval<0,0,1>(0, ovflag, buffers, fc, host_start, inum);
      } else {
        eval<0,0,0>(1, ovflag, buffers, fc, 0, offload_end);
        eval<0,0,0>(0, ovflag, buffers, fc, host_start, inum);
      }
    }
  }
}

template <int ONETYPE, int EFLAG, int NEWTON_PAIR, class flt_t, class acc_t>
void PairDPDIntel::eval(const int offload, const int vflag,
                        IntelBuffers<flt_t,acc_t> *buffers,
                        const ForceConst<flt_t> &fc,
                        const int astart, const int aend)
{
  const int inum = aend - astart;
  if (inum == 0) return;
  int nlocal, nall, minlocal;
  fix->get_buffern(offload, nlocal, nall, minlocal);

  const int ago = neighbor->ago;
  IP_PRE_pack_separate_buffers(fix, buffers, ago, offload, nlocal, nall);

  ATOM_T * _noalias const x = buffers->get_x(offload);
  typedef struct { double x, y, z; } lmp_vt;
  lmp_vt *v = (lmp_vt *)atom->v[0];
  const flt_t dtinvsqrt = 1.0/sqrt(update->dt);

  const int * _noalias const ilist = list->ilist;
  const int * _noalias const numneigh = list->numneigh;
  const int ** _noalias const firstneigh = (const int **)list->firstneigh;
  const FC_PACKED1_T * _noalias const param = fc.param[0];
  const flt_t * _noalias const special_lj = fc.special_lj;
  int * _noalias const rngi_thread = fc.rngi;
  const int rng_size = buffers->get_max_nbors();

  const int ntypes = atom->ntypes + 1;
  const int eatom = this->eflag_atom;

  // Determine how much data to transfer
  int x_size, q_size, f_stride, ev_size, separate_flag;
  IP_PRE_get_transfern(ago, NEWTON_PAIR, EFLAG, vflag,
                       buffers, offload, fix, separate_flag,
                       x_size, q_size, ev_size, f_stride);

  int tc;
  FORCE_T * _noalias f_start;
  acc_t * _noalias ev_global;
  IP_PRE_get_buffers(offload, buffers, fix, tc, f_start, ev_global);
  const int nthreads = tc;
  int *overflow = fix->get_off_overflow_flag();
  {
    #if defined(__MIC__) && defined(_LMP_INTEL_OFFLOAD)
    *timer_compute = MIC_Wtime();
    #endif

    IP_PRE_repack_for_offload(NEWTON_PAIR, separate_flag, nlocal, nall,
                              f_stride, x, 0);

    acc_t oevdwl, ov0, ov1, ov2, ov3, ov4, ov5;
    if (EFLAG || vflag)
      oevdwl = ov0 = ov1 = ov2 = ov3 = ov4 = ov5 = (acc_t)0;
    if (NEWTON_PAIR == 0 && inum != nlocal)
      memset(f_start, 0, f_stride * sizeof(FORCE_T));

    // loop over neighbors of my atoms
    #if defined(_OPENMP)
    #pragma omp parallel reduction(+:oevdwl,ov0,ov1,ov2,ov3,ov4,ov5)
    #endif
    {
      int iifrom, iip, iito, tid;
      IP_PRE_omp_stride_id(iifrom, iip, iito, tid, inum, nthreads);
      iifrom += astart;
      iito += astart;

      #ifdef LMP_USE_MKL_RNG
      VSLStreamStatePtr *my_random = &(random_thread[tid]);
      #else
      RanMars *my_random = random_thread[tid];
      #endif
      flt_t *my_rand_buffer = fc.rand_buffer_thread[tid];
      int rngi = rngi_thread[tid];

      int foff;
      if (NEWTON_PAIR) foff = tid * f_stride - minlocal;
      else foff = -minlocal;
      FORCE_T * _noalias const f = f_start + foff;
      if (NEWTON_PAIR) memset(f + minlocal, 0, f_stride * sizeof(FORCE_T));

      flt_t icut, a0, gamma, sigma;
      if (ONETYPE) {
        icut = param[_onetype].icut;
        a0 = param[_onetype].a0;
        gamma = param[_onetype].gamma;
        sigma = param[_onetype].sigma;
      }
      for (int ii = iifrom; ii < iito; ii += iip) {
        const int i = ilist[ii];
        int itype, ptr_off;
        const FC_PACKED1_T * _noalias parami;
        if (!ONETYPE) {
          itype = x[i].w;
          ptr_off = itype * ntypes;
          parami = param + ptr_off;
        }

        const int * _noalias const jlist = firstneigh[i];
        int jnum = numneigh[i];
        IP_PRE_neighbor_pad(jnum, offload);

        acc_t fxtmp, fytmp, fztmp, fwtmp;
        acc_t sevdwl, sv0, sv1, sv2, sv3, sv4, sv5;

        const flt_t xtmp = x[i].x;
        const flt_t ytmp = x[i].y;
        const flt_t ztmp = x[i].z;
        const flt_t vxtmp = v[i].x;
        const flt_t vytmp = v[i].y;
        const flt_t vztmp = v[i].z;
        fxtmp = fytmp = fztmp = (acc_t)0;
        if (EFLAG) fwtmp = sevdwl = (acc_t)0;
        if (NEWTON_PAIR == 0)
          if (vflag==1) sv0 = sv1 = sv2 = sv3 = sv4 = sv5 = (acc_t)0;

        if (rngi + jnum > rng_size) {
          #ifdef LMP_USE_MKL_RNG
          if (sizeof(flt_t) == sizeof(float))
            vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, *my_random, rngi,
                          (float*)my_rand_buffer, (float)0.0, (float)1.0 );
          else
            vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, *my_random, rngi,
                          (double*)my_rand_buffer, 0.0, 1.0 );
          #else
          for (int jj = 0; jj < rngi; jj++)
            my_rand_buffer[jj] = my_random->gaussian();
          #endif
          rngi = 0;
        }

        #if defined(LMP_SIMD_COMPILER)
        #pragma vector aligned
        #pragma simd reduction(+:fxtmp, fytmp, fztmp, fwtmp, sevdwl, \
                                 sv0, sv1, sv2, sv3, sv4, sv5)
        #endif
        for (int jj = 0; jj < jnum; jj++) {
          flt_t forcelj, evdwl;
          forcelj = evdwl = (flt_t)0.0;

          int j, jtype, sbindex;
          if (!ONETYPE) {
            sbindex = jlist[jj] >> SBBITS & 3;
            j = jlist[jj] & NEIGHMASK;
          } else
            j = jlist[jj];

          const flt_t delx = xtmp - x[j].x;
          const flt_t dely = ytmp - x[j].y;
          const flt_t delz = ztmp - x[j].z;
          if (!ONETYPE) {
            jtype = x[j].w;
            icut = parami[jtype].icut;
          }
          const flt_t rsq = delx * delx + dely * dely + delz * delz;
          const flt_t rinv = (flt_t)1.0/sqrt(rsq);

          if (rinv > icut) {
            flt_t factor_dpd;
            if (!ONETYPE) factor_dpd = special_lj[sbindex];

            flt_t delvx = vxtmp - v[j].x;
            flt_t delvy = vytmp - v[j].y;
            flt_t delvz = vztmp - v[j].z;
            flt_t dot = delx*delvx + dely*delvy + delz*delvz;
            flt_t randnum = my_rand_buffer[jj];

            flt_t iwd = rinv - icut;
            if (rinv > (flt_t)IEPSILON) iwd = (flt_t)0.0;

            if (!ONETYPE) {
              a0 = parami[jtype].a0;
              gamma = parami[jtype].gamma;
              sigma = parami[jtype].sigma;
            }
            flt_t fpair = a0 - iwd * gamma * dot + sigma * randnum * dtinvsqrt;
            if (!ONETYPE) fpair *= factor_dpd;
            fpair *= iwd;

            const flt_t fpx = fpair * delx;
            fxtmp += fpx;
            if (NEWTON_PAIR) f[j].x -= fpx;
            const flt_t fpy = fpair * dely;
            fytmp += fpy;
            if (NEWTON_PAIR) f[j].y -= fpy;
            const flt_t fpz = fpair * delz;
            fztmp += fpz;
            if (NEWTON_PAIR) f[j].z -= fpz;

            if (EFLAG) {
              flt_t cut = (flt_t)1.0/icut;
              flt_t r = (flt_t)1.0/rinv;
              evdwl = (flt_t)0.5 * a0 * (cut - (flt_t)2.0*r + rsq * icut);
              if (!ONETYPE) evdwl *= factor_dpd;
              sevdwl += evdwl;
              if (eatom) {
                fwtmp += (flt_t)0.5 * evdwl;
                if (NEWTON_PAIR)
                  f[j].w += (flt_t)0.5 * evdwl;
              }
            }

            if (NEWTON_PAIR == 0)
              IP_PRE_ev_tally_nborv(vflag, delx, dely, delz, fpx, fpy, fpz);
          } // if rsq
        } // for jj
        if (NEWTON_PAIR) {
          f[i].x += fxtmp;
          f[i].y += fytmp;
          f[i].z += fztmp;
        } else {
          f[i].x = fxtmp;
          f[i].y = fytmp;
          f[i].z = fztmp;
        }

        IP_PRE_ev_tally_atom(NEWTON_PAIR, EFLAG, vflag, f, fwtmp);
        rngi += jnum;
      } // for ii

      IP_PRE_fdotr_reduce_omp(NEWTON_PAIR, nall, minlocal, nthreads, f_start,
                              f_stride, x, offload, vflag, ov0, ov1, ov2, ov3,
                              ov4, ov5);
      rngi_thread[tid] = rngi;
    } // end omp

    IP_PRE_fdotr_reduce(NEWTON_PAIR, nall, nthreads, f_stride, vflag,
                        ov0, ov1, ov2, ov3, ov4, ov5);

    if (EFLAG || vflag) {
      if (NEWTON_PAIR == 0) {
        oevdwl *= (acc_t)0.5;
        ov0 *= (acc_t)0.5;
        ov1 *= (acc_t)0.5;
        ov2 *= (acc_t)0.5;
        ov3 *= (acc_t)0.5;
        ov4 *= (acc_t)0.5;
        ov5 *= (acc_t)0.5;
      }
      ev_global[0] = oevdwl;
      ev_global[1] = (acc_t)0.0;
      ev_global[2] = ov0;
      ev_global[3] = ov1;
      ev_global[4] = ov2;
      ev_global[5] = ov3;
      ev_global[6] = ov4;
      ev_global[7] = ov5;
    }
    #if defined(__MIC__) && defined(_LMP_INTEL_OFFLOAD)
    *timer_compute = MIC_Wtime() - *timer_compute;
    #endif
  } // end offload

  if (offload)
    fix->stop_watch(TIME_OFFLOAD_LATENCY);
  else
    fix->stop_watch(TIME_HOST_PAIR);

  if (EFLAG || vflag)
    fix->add_result_array(f_start, ev_global, offload, eatom, 0, vflag);
  else
    fix->add_result_array(f_start, 0, offload);
}

/* ----------------------------------------------------------------------
   global settings
   ------------------------------------------------------------------------- */

void PairDPDIntel::settings(int narg, char **arg) {
  #if defined(_OPENMP)
  if (_nrandom_thread) {
    #ifdef LMP_USE_MKL_RNG
    for (int i = 0; i < _nrandom_thread; i++)
      vslDeleteStream(&random_thread[i]);
    #else
    for (int i = 1; i < _nrandom_thread; i++)
      delete random_thread[i];
    #endif
  }
  delete []random_thread;
  #endif
  PairDPD::settings(narg,arg);
  _nrandom_thread = comm->nthreads;

  #ifdef LMP_USE_MKL_RNG

  random_thread=new VSLStreamStatePtr[comm->nthreads];
  #if defined(_OPENMP)
  #pragma omp parallel
  {
    int tid = omp_get_thread_num();
    vslNewStream(&random_thread[tid], LMP_MKL_RNG,
                 seed + comm->me + comm->nprocs * tid );
  }
  #endif

  #else

  random_thread =new RanMars*[comm->nthreads];
  random_thread[0] = random;
  #if defined(_OPENMP)
  #pragma omp parallel
  {
    int tid = omp_get_thread_num();
    if (tid > 0)
      random_thread[tid] = new RanMars(lmp, seed+comm->me+comm->nprocs*tid);
  }
  #endif

  #endif
}

/* ---------------------------------------------------------------------- */

void PairDPDIntel::init_style()
{
  PairDPD::init_style();
  if (force->newton_pair == 0) {
    neighbor->requests[neighbor->nrequest-1]->half = 0;
    neighbor->requests[neighbor->nrequest-1]->full = 1;
  }
  neighbor->requests[neighbor->nrequest-1]->intel = 1;

  int ifix = modify->find_fix("package_intel");
  if (ifix < 0)
    error->all(FLERR,
               "The 'package intel' command is required for /intel styles");
  fix = static_cast<FixIntel *>(modify->fix[ifix]);

  fix->pair_init_check();
  #ifdef _LMP_INTEL_OFFLOAD
  if (fix->offload_balance() != 0.0)
    error->all(FLERR,
          "Offload for dpd/intel is not yet available. Set balance to 0.");
  #endif

  if (fix->precision() == FixIntel::PREC_MODE_MIXED)
    pack_force_const(force_const_single, fix->get_mixed_buffers());
  else if (fix->precision() == FixIntel::PREC_MODE_DOUBLE)
    pack_force_const(force_const_double, fix->get_double_buffers());
  else
    pack_force_const(force_const_single, fix->get_single_buffers());
}

/* ---------------------------------------------------------------------- */

template <class flt_t, class acc_t>
void PairDPDIntel::pack_force_const(ForceConst<flt_t> &fc,
                                    IntelBuffers<flt_t,acc_t> *buffers)
{
  _onetype = 0;

  int tp1 = atom->ntypes + 1;
  fc.set_ntypes(tp1,comm->nthreads,buffers->get_max_nbors(),memory,_cop);

  // Repeat cutsq calculation because done after call to init_style
  int mytypes = 0;
  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = i; j <= atom->ntypes; j++) {
      if (setflag[i][j] != 0 || (setflag[i][i] != 0 && setflag[j][j] != 0)) {
        mytypes++;
        _onetype = i * tp1 + j;
      }
      double cut = init_one(i,j);
      cutsq[i][j] = cutsq[j][i] = cut*cut;
      double icut = 1.0 / cut;
      fc.param[i][j].icut = fc.param[j][i].icut = icut;
    }
  }
  if (mytypes > 1 || atom->molecular) _onetype = 0;

  for (int i = 0; i < 4; i++) {
    fc.special_lj[i] = force->special_lj[i];
    fc.special_lj[0] = 1.0;
  }

  for (int i = 1; i < tp1; i++) {
    for (int j = 1; j < tp1; j++) {
      fc.param[i][j].a0 = a0[i][j];
      fc.param[i][j].gamma = gamma[i][j];
      fc.param[i][j].sigma = sigma[i][j];
    }
  }
}

/* ---------------------------------------------------------------------- */

template <class flt_t>
void PairDPDIntel::ForceConst<flt_t>::set_ntypes(const int ntypes,
                                                 const int nthreads,
                                                 const int max_nbors,
                                                 Memory *memory,
                                                 const int cop) {
  if (ntypes != _ntypes) {
    if (_ntypes > 0) {
      _memory->destroy(param);
      _memory->destroy(rand_buffer_thread);
      _memory->destroy(rngi);
    }
    if (ntypes > 0) {
      _cop = cop;
      memory->create(param,ntypes,ntypes,"fc.param");
      memory->create(rand_buffer_thread, nthreads, max_nbors,
                     "fc.rand_buffer_thread");
      memory->create(rngi,nthreads,"fc.param");
      for (int i = 0; i < nthreads; i++) rngi[i] = max_nbors;
    }
  }
  _ntypes = ntypes;
  _memory = memory;
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
   ------------------------------------------------------------------------- */

void PairDPDIntel::read_restart_settings(FILE *fp)
{
  #if defined(_OPENMP)
  if (_nrandom_thread) {
    #ifdef LMP_USE_MKL_RNG
    for (int i = 0; i < _nrandom_thread; i++)
      vslDeleteStream(&random_thread[i]);
    #else
    for (int i = 1; i < _nrandom_thread; i++)
      delete random_thread[i];
    #endif
  }
  delete []random_thread;
  #endif
  PairDPD::read_restart_settings(fp);
  _nrandom_thread = comm->nthreads;

  #ifdef LMP_USE_MKL_RNG

  random_thread=new VSLStreamStatePtr[comm->nthreads];
  #if defined(_OPENMP)
  #pragma omp parallel
  {
    int tid = omp_get_thread_num();
    vslNewStream(&random_thread[tid], LMP_MKL_RNG,
                 seed + comm->me + comm->nprocs * tid );
  }
  #endif

  #else

  random_thread =new RanMars*[comm->nthreads];
  random_thread[0] = random;
  #if defined(_OPENMP)
  #pragma omp parallel
  {
    int tid = omp_get_thread_num();
    if (tid > 0)
      random_thread[tid] = new RanMars(lmp, seed+comm->me+comm->nprocs*tid);
  }
  #endif

  #endif
}
