/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: W. Michael Brown (Intel)
------------------------------------------------------------------------- */

#include "math.h"
#include "pair_lj_cut_intel.h"
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

#define FC_PACKED1_T typename ForceConst<flt_t>::fc_packed1
#define FC_PACKED2_T typename ForceConst<flt_t>::fc_packed2

/* ---------------------------------------------------------------------- */

PairLJCutIntel::PairLJCutIntel(LAMMPS *lmp) :
  PairLJCut(lmp)
{
  suffix_flag |= Suffix::INTEL;
  respa_enable = 0;
  cut_respa = NULL;
}

/* ---------------------------------------------------------------------- */

void PairLJCutIntel::compute(int eflag, int vflag)
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
void PairLJCutIntel::compute(int eflag, int vflag,
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
    if (ago != 0) {
      #if defined(_OPENMP)
      #pragma omp parallel default(none) shared(eflag,vflag,buffers,fc)
      #endif
      {
        int ifrom, ito, tid;
	IP_PRE_omp_range_id_align(ifrom, ito, tid, atom->nlocal + atom->nghost,
				  nthreads, sizeof(ATOM_T));
	buffers->thr_pack(ifrom,ito,ago);
      }
    }
    fix->stop_watch(TIME_PACK);
  }

  if (evflag || vflag_fdotr) {
    int ovflag = 0;
    if (vflag_fdotr) ovflag = 2;
    else if (vflag) ovflag = 1;
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
    if (force->newton_pair) {
      eval<0,0,1>(1, 0, buffers, fc, 0, offload_end);
      eval<0,0,1>(0, 0, buffers, fc, host_start, inum);
    } else {
      eval<0,0,0>(1, 0, buffers, fc, 0, offload_end);
      eval<0,0,0>(0, 0, buffers, fc, host_start, inum);
    }
  }
}

template <int EVFLAG, int EFLAG, int NEWTON_PAIR, class flt_t, class acc_t>
void PairLJCutIntel::eval(const int offload, const int vflag,
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

  ATOM_T * restrict const x = buffers->get_x(offload);

  const int * restrict const numneigh = list->numneigh;
  const int * restrict const cnumneigh = buffers->cnumneigh(list);
  const int * restrict const firstneigh = buffers->firstneigh(list);
  const flt_t * restrict const special_lj = fc.special_lj;
  const FC_PACKED1_T * restrict const ljc12o = fc.ljc12o[0];
  const FC_PACKED2_T * restrict const lj34 = fc.lj34[0];

  const int ntypes = atom->ntypes + 1;
  const int eatom = this->eflag_atom;

  // Determine how much data to transfer
  int x_size, q_size, f_stride, ev_size, separate_flag;
  IP_PRE_get_transfern(ago, NEWTON_PAIR, EVFLAG, EFLAG, vflag,
		       buffers, offload, fix, separate_flag,
		       x_size, q_size, ev_size, f_stride);

  int tc;
  FORCE_T * restrict f_start;
  acc_t * restrict ev_global;
  IP_PRE_get_buffers(offload, buffers, fix, tc, f_start, ev_global);
  const int nthreads = tc;
  int *overflow = fix->get_off_overflow_flag();
  {
    #ifdef __MIC__
    *timer_compute = MIC_Wtime();
    #endif

    IP_PRE_repack_for_offload(NEWTON_PAIR, separate_flag, nlocal, nall, 
			      f_stride, x, 0);

    acc_t oevdwl, ov0, ov1, ov2, ov3, ov4, ov5;
    if (EVFLAG) {
      oevdwl = (acc_t)0;
      if (vflag) ov0 = ov1 = ov2 = ov3 = ov4 = ov5 = (acc_t)0;
    }

    // loop over neighbors of my atoms
    #if defined(_OPENMP)
    #pragma omp parallel default(none) \
      shared(f_start,f_stride,nlocal,nall,minlocal) \
      reduction(+:oevdwl,ov0,ov1,ov2,ov3,ov4,ov5)
    #endif
    {
      int iifrom, iito, tid;
      IP_PRE_omp_range_id(iifrom, iito, tid, inum, nthreads);
      iifrom += astart;
      iito += astart;

      FORCE_T * restrict const f = f_start - minlocal + (tid * f_stride);
      memset(f + minlocal, 0, f_stride * sizeof(FORCE_T));

      for (int i = iifrom; i < iito; ++i) {
        const int itype = x[i].w;

        const int ptr_off = itype * ntypes;
        const FC_PACKED1_T * restrict const ljc12oi = ljc12o + ptr_off;
        const FC_PACKED2_T * restrict const lj34i = lj34 + ptr_off;

        const int * restrict const jlist = firstneigh + cnumneigh[i];
        const int jnum = numneigh[i];

        acc_t fxtmp, fytmp, fztmp, fwtmp;
        acc_t sevdwl, sv0, sv1, sv2, sv3, sv4, sv5;

        const flt_t xtmp = x[i].x;
        const flt_t ytmp = x[i].y;
        const flt_t ztmp = x[i].z;
        fxtmp = fytmp = fztmp = (acc_t)0;
        if (EVFLAG) {
          if (EFLAG) fwtmp = sevdwl = (acc_t)0;
          if (vflag==1) sv0 = sv1 = sv2 = sv3 = sv4 = sv5 = (acc_t)0;
        }

        #pragma vector aligned
	#pragma simd reduction(+:fxtmp, fytmp, fztmp, fwtmp, sevdwl, \
	                       sv0, sv1, sv2, sv3, sv4, sv5)
        for (int jj = 0; jj < jnum; jj++) {
          flt_t forcelj, evdwl;
          forcelj = evdwl = (flt_t)0.0;

          const int sbindex = jlist[jj] >> SBBITS & 3;
          const int j = jlist[jj] & NEIGHMASK;
          const flt_t delx = xtmp - x[j].x;
          const flt_t dely = ytmp - x[j].y;
          const flt_t delz = ztmp - x[j].z;
          const int jtype = x[j].w;
          const flt_t rsq = delx * delx + dely * dely + delz * delz;

          #ifdef __MIC__
          if (rsq < ljc12oi[jtype].cutsq) {
	  #endif
            flt_t factor_lj = special_lj[sbindex];
            flt_t r2inv = 1.0 / rsq;
            flt_t r6inv = r2inv * r2inv * r2inv;
            #ifndef __MIC__
	    if (rsq > ljc12oi[jtype].cutsq) r6inv = (flt_t)0.0;
	    #endif
            forcelj = r6inv * (ljc12oi[jtype].lj1 * r6inv - ljc12oi[jtype].lj2);
            flt_t fpair = factor_lj * forcelj * r2inv;

            fxtmp += delx * fpair;
            fytmp += dely * fpair;
            fztmp += delz * fpair;
            if (NEWTON_PAIR || j < nlocal) {
              f[j].x -= delx * fpair;
              f[j].y -= dely * fpair;
              f[j].z -= delz * fpair;
            }

            if (EVFLAG) {
              flt_t ev_pre = (flt_t)0;
              if (NEWTON_PAIR || i<nlocal)
                ev_pre += (flt_t)0.5;
              if (NEWTON_PAIR || j<nlocal)
                ev_pre += (flt_t)0.5;

              if (EFLAG) {
                evdwl = r6inv * (lj34i[jtype].lj3 * r6inv-lj34i[jtype].lj4) -
                    ljc12oi[jtype].offset;
                evdwl *= factor_lj;
                sevdwl += ev_pre*evdwl;
                if (eatom) {
                  if (NEWTON_PAIR || i < nlocal)
                    fwtmp += 0.5 * evdwl;
                  if (NEWTON_PAIR || j < nlocal)
                    f[j].w += 0.5 * evdwl;
                }
              }

	      IP_PRE_ev_tally_nbor(vflag, ev_pre, fpair,
				   delx, dely, delz);
            }
          #ifdef __MIC__
          } // if rsq
          #endif
        } // for jj
        f[i].x += fxtmp;
        f[i].y += fytmp;
        f[i].z += fztmp;
        
	IP_PRE_ev_tally_atom(EVFLAG, EFLAG, vflag, f, fwtmp);
      } // for ii

      #if defined(_OPENMP)
      #pragma omp barrier
      #endif
      IP_PRE_fdotr_acc_force(NEWTON_PAIR, EVFLAG,  EFLAG, vflag, eatom, nall,
			     nlocal, minlocal, nthreads, f_start, f_stride, 
			     x);
    } // end omp
    if (EVFLAG) {
      if (EFLAG) {
        ev_global[0] = oevdwl;
	ev_global[1] = (acc_t)0.0;
      }
      if (vflag) {
        ev_global[2] = ov0;
        ev_global[3] = ov1;
        ev_global[4] = ov2;
        ev_global[5] = ov3;
        ev_global[6] = ov4;
        ev_global[7] = ov5;
      }
    }
    #ifdef __MIC__
    *timer_compute = MIC_Wtime() - *timer_compute;
    #endif
  } // end offload

  if (offload)
    fix->stop_watch(TIME_OFFLOAD_LATENCY);
  else
    fix->stop_watch(TIME_HOST_PAIR);

  if (EVFLAG)
    fix->add_result_array(f_start, ev_global, offload, eatom);
  else
    fix->add_result_array(f_start, 0, offload);
}

/* ---------------------------------------------------------------------- */

void PairLJCutIntel::init_style()
{
  PairLJCut::init_style();
  neighbor->requests[neighbor->nrequest-1]->intel = 1;

  int ifix = modify->find_fix("package_intel");
  if (ifix < 0)
    error->all(FLERR,
               "The 'package intel' command is required for /intel styles");
  fix = static_cast<FixIntel *>(modify->fix[ifix]);

  #ifdef _LMP_INTEL_OFFLOAD
  if (fix->offload_balance() != 0.0)
    error->all(FLERR,
          "Offload for lj/cut/intel is not yet available. Set balance to 0.");
  #endif
  if (fix->precision() == FixIntel::PREC_MODE_MIXED) {
    fix->get_mixed_buffers()->free_all_nbor_buffers();
    pack_force_const(force_const_single, fix->get_mixed_buffers());
  } else if (fix->precision() == FixIntel::PREC_MODE_DOUBLE) {
    fix->get_double_buffers()->free_all_nbor_buffers();
    pack_force_const(force_const_double, fix->get_double_buffers());
  } else {
    fix->get_single_buffers()->free_all_nbor_buffers();
    pack_force_const(force_const_single, fix->get_single_buffers());
  }
}

/* ---------------------------------------------------------------------- */

template <class flt_t, class acc_t>
void PairLJCutIntel::pack_force_const(ForceConst<flt_t> &fc,
                                      IntelBuffers<flt_t,acc_t> *buffers)
{
  int tp1 = atom->ntypes + 1;
  fc.set_ntypes(tp1,memory,_cop);
  buffers->set_ntypes(tp1);
  flt_t **cutneighsq = buffers->get_cutneighsq();

  // Repeat cutsq calculation because done after call to init_style
  double cut, cutneigh;
  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = i; j <= atom->ntypes; j++) {
      if (setflag[i][j] != 0 || (setflag[i][i] != 0 && setflag[j][j] != 0)) {
        cut = init_one(i,j);
        cutneigh = cut + neighbor->skin;
        cutsq[i][j] = cutsq[j][i] = cut*cut;
        cutneighsq[i][j] = cutneighsq[j][i] = cutneigh * cutneigh;
      }
    }
  }

  for (int i = 0; i < 4; i++) {
    fc.special_lj[i] = force->special_lj[i];
    fc.special_lj[0] = 1.0;
  }

  for (int i = 0; i < tp1; i++) {
    for (int j = 0; j < tp1; j++) {
      fc.ljc12o[i][j].lj1 = lj1[i][j];
      fc.ljc12o[i][j].lj2 = lj2[i][j];
      fc.lj34[i][j].lj3 = lj3[i][j];
      fc.lj34[i][j].lj4 = lj4[i][j];
      fc.ljc12o[i][j].cutsq = cutsq[i][j];
      fc.ljc12o[i][j].offset = offset[i][j];
    }
  }
}

/* ---------------------------------------------------------------------- */

template <class flt_t>
void PairLJCutIntel::ForceConst<flt_t>::set_ntypes(const int ntypes,
                                                   Memory *memory,
						   const int cop) {
  if (ntypes != _ntypes) {
    if (_ntypes > 0) {
      fc_packed1 *oljc12o = ljc12o[0];
      fc_packed2 *olj34 = lj34[0];

      _memory->destroy(oljc12o);
      _memory->destroy(olj34);
    }
    if (ntypes > 0) {
      _cop = cop;
      memory->create(ljc12o,ntypes,ntypes,"fc.c12o");
      memory->create(lj34,ntypes,ntypes,"fc.lj34");
    }
  }
  _ntypes = ntypes;
  _memory = memory;
}
