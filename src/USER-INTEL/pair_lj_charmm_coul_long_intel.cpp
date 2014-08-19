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
#include "pair_lj_charmm_coul_long_intel.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "group.h"
#include "kspace.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "suffix.h"
using namespace LAMMPS_NS;

#define LJ_T typename IntelBuffers<flt_t,flt_t>::vec4_t
#define TABLE_T typename ForceConst<flt_t>::table_t

/* ---------------------------------------------------------------------- */

PairLJCharmmCoulLongIntel::PairLJCharmmCoulLongIntel(LAMMPS *lmp) :
  PairLJCharmmCoulLong(lmp)
{
  suffix_flag |= Suffix::INTEL;
  respa_enable = 0;
  cut_respa = NULL;
}

/* ---------------------------------------------------------------------- */

PairLJCharmmCoulLongIntel::~PairLJCharmmCoulLongIntel()
{
}

/* ---------------------------------------------------------------------- */

void PairLJCharmmCoulLongIntel::compute(int eflag, int vflag)
{
  if (fix->precision()==FixIntel::PREC_MODE_MIXED)
    compute<float,double>(eflag, vflag, fix->get_mixed_buffers(), 
                          force_const_single);
  else if (fix->precision()==FixIntel::PREC_MODE_DOUBLE)
    compute<double,double>(eflag, vflag, fix->get_double_buffers(),
                           force_const_double);
  else
    compute<float,float>(eflag, vflag, fix->get_single_buffers(),
                         force_const_single);

  fix->balance_stamp();
  vflag_fdotr = 0;
}

template <class flt_t, class acc_t>
void PairLJCharmmCoulLongIntel::compute(int eflag, int vflag,
					IntelBuffers<flt_t,acc_t> *buffers,
					const ForceConst<flt_t> &fc)
{
  if (eflag || vflag) {
    ev_setup(eflag,vflag);
  } else evflag = vflag_fdotr = 0;

  const int inum = list->inum;
  const int nthreads = comm->nthreads;
  const int host_start = fix->host_start_pair();
  const int offload_end = fix->offload_end_pair();
  const int ago = neighbor->ago;
  
  if (ago != 0 && fix->separate_buffers() == 0) {
    fix->start_watch(TIME_PACK);
    #if defined(_OPENMP)
    #pragma omp parallel default(none) shared(eflag,vflag,buffers,fc)
    #endif
    {
      int ifrom, ito, tid;
      IP_PRE_omp_range_id_align(ifrom, ito, tid, atom->nlocal+atom->nghost, 
			      nthreads, sizeof(ATOM_T));
      buffers->thr_pack(ifrom,ito,ago);
    }
    fix->stop_watch(TIME_PACK);
  }
  
  // -------------------- Regular version
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

/* ---------------------------------------------------------------------- */

template <int EVFLAG, int EFLAG, int NEWTON_PAIR, class flt_t, class acc_t>
void PairLJCharmmCoulLongIntel::eval(const int offload, const int vflag,
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
  flt_t * restrict const q = buffers->get_q(offload);

  const int * restrict const numneigh = list->numneigh;
  const int * restrict const cnumneigh = buffers->cnumneigh(list);
  const int * restrict const firstneigh = buffers->firstneigh(list);

  const flt_t * restrict const special_coul = fc.special_coul;
  const flt_t * restrict const special_lj = fc.special_lj;
  const flt_t qqrd2e = force->qqrd2e;
  const flt_t inv_denom_lj = (flt_t)1.0/denom_lj;

  const flt_t * restrict const cutsq = fc.cutsq[0];
  const LJ_T * restrict const lj = fc.lj[0];
  const TABLE_T * restrict const table = fc.table;
  const flt_t * restrict const etable = fc.etable;
  const flt_t * restrict const detable = fc.detable;
  const flt_t * restrict const ctable = fc.ctable;
  const flt_t * restrict const dctable = fc.dctable;
  const flt_t cut_ljsq = fc.cut_ljsq;
  const flt_t cut_lj_innersq = fc.cut_lj_innersq;
  const flt_t cut_coulsq = fc.cut_coulsq;
  const flt_t g_ewald = fc.g_ewald;
  const flt_t tabinnersq = fc.tabinnersq;

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
  #ifdef _LMP_INTEL_OFFLOAD
  int *overflow = fix->get_off_overflow_flag();
  double *timer_compute = fix->off_watch_pair();
  // Redeclare as local variables for offload
  const int ncoultablebits = this->ncoultablebits;
  const int ncoulmask = this->ncoulmask;
  const int ncoulshiftbits = this->ncoulshiftbits;
  #ifdef INTEL_ALLOW_TABLE
  #define ITABLE_IN in(table,etable,detable:length(0) alloc_if(0) free_if(0)) \
                    in(ctable,dctable:length(0) alloc_if(0) free_if(0)) \
                    in(ncoultablebits,tabinnersq,ncoulmask,ncoulshiftbits)
  #else
  #define ITABLE_IN
  #endif

  if (offload) fix->start_watch(TIME_OFFLOAD_LATENCY);
  #pragma offload target(mic:_cop) if(offload) \
    in(special_lj,special_coul:length(0) alloc_if(0) free_if(0)) \
    in(cutsq,lj:length(0) alloc_if(0) free_if(0)) \
    in(firstneigh:length(0) alloc_if(0) free_if(0)) \
    in(cnumneigh:length(0) alloc_if(0) free_if(0)) \
    in(numneigh:length(0) alloc_if(0) free_if(0)) \
    in(x:length(x_size) alloc_if(0) free_if(0)) \
    in(q:length(q_size) alloc_if(0) free_if(0)) \
    in(overflow:length(0) alloc_if(0) free_if(0)) \
    in(nthreads,qqrd2e,g_ewald,inum,nall,ntypes,cut_coulsq,vflag,eatom) \
    in(f_stride,separate_flag,offload) \
    in(astart,cut_ljsq,cut_lj_innersq,nlocal,inv_denom_lj,minlocal) \
    out(f_start:length(f_stride) alloc_if(0) free_if(0)) \
    out(ev_global:length(ev_size) alloc_if(0) free_if(0)) \
    out(timer_compute:length(1) alloc_if(0) free_if(0)) \
    ITABLE_IN signal(f_start)
  #endif
  {
    #ifdef __MIC__
    *timer_compute = MIC_Wtime();
    #endif

    IP_PRE_repack_for_offload(NEWTON_PAIR, separate_flag, nlocal, nall, 
			      f_stride, x, q);

    acc_t oevdwl, oecoul, ov0, ov1, ov2, ov3, ov4, ov5;
    if (EVFLAG) {
      oevdwl = oecoul = (acc_t)0;
      if (vflag) ov0 = ov1 = ov2 = ov3 = ov4 = ov5 = (acc_t)0;
    }

    // loop over neighbors of my atoms
    #if defined(_OPENMP)
    #pragma omp parallel default(none) \
      shared(f_start,f_stride,nlocal,nall,minlocal) \
      reduction(+:oevdwl,oecoul,ov0,ov1,ov2,ov3,ov4,ov5)
    #endif
    {
      int iifrom, iito, tid;
      IP_PRE_omp_range_id(iifrom, iito, tid, inum, nthreads);
      iifrom += astart;
      iito += astart;

      FORCE_T * restrict const f = f_start - minlocal + (tid * f_stride);
      memset(f + minlocal, 0, f_stride * sizeof(FORCE_T));
      flt_t cutboth = cut_coulsq;

      for (int i = iifrom; i < iito; ++i) {
	//        const int i = ilist[ii];
        const int itype = x[i].w;

        const int ptr_off = itype * ntypes;
        const flt_t * restrict const cutsqi = cutsq + ptr_off;
        const LJ_T * restrict const lji = lj + ptr_off;

        const int   * restrict const jlist = firstneigh + cnumneigh[i];
        const int jnum = numneigh[i];

        acc_t fxtmp,fytmp,fztmp,fwtmp;
	acc_t sevdwl, secoul, sv0, sv1, sv2, sv3, sv4, sv5;

        const flt_t xtmp = x[i].x;
        const flt_t ytmp = x[i].y;
        const flt_t ztmp = x[i].z;
        const flt_t qtmp = q[i];
        fxtmp = fytmp = fztmp = (acc_t)0;
        if (EVFLAG) {
	  if (EFLAG) fwtmp = sevdwl = secoul = (acc_t)0;
	  if (vflag==1) sv0 = sv1 = sv2 = sv3 = sv4 = sv5 = (acc_t)0;
	}

	#pragma vector aligned
	#pragma simd reduction(+:fxtmp, fytmp, fztmp, fwtmp, sevdwl, secoul, \
	                       sv0, sv1, sv2, sv3, sv4, sv5)
        for (int jj = 0; jj < jnum; jj++) {
          flt_t forcecoul, forcelj, evdwl, ecoul;
          forcecoul = forcelj = evdwl = ecoul = (flt_t)0.0;

          const int sbindex = jlist[jj] >> SBBITS & 3;
          const int j = jlist[jj] & NEIGHMASK;

          const flt_t delx = xtmp - x[j].x;
          const flt_t dely = ytmp - x[j].y;
          const flt_t delz = ztmp - x[j].z;
          const int jtype = x[j].w;
          const flt_t rsq = delx * delx + dely * dely + delz * delz;
          const flt_t r2inv = (flt_t)1.0 / rsq;

	  #ifdef __MIC__
	  if (rsq < cut_coulsq) {
          #endif
            #ifdef INTEL_ALLOW_TABLE
            if (!ncoultablebits || rsq <= tabinnersq) {
            #endif
              const flt_t A1 =  0.254829592;
              const flt_t A2 = -0.284496736;
              const flt_t A3 =  1.421413741;
              const flt_t A4 = -1.453152027;
              const flt_t A5 =  1.061405429;
              const flt_t EWALD_F = 1.12837917;
              const flt_t INV_EWALD_P = 1.0 / 0.3275911;

              const flt_t r = sqrt(rsq);
              const flt_t grij = g_ewald * r;
              const flt_t expm2 = exp(-grij * grij);
              const flt_t t = INV_EWALD_P / (INV_EWALD_P + grij);
              const flt_t erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
              const flt_t prefactor = qqrd2e * qtmp * q[j] / r;
              forcecoul = prefactor * (erfc + EWALD_F * grij * expm2);
              if (EFLAG) ecoul = prefactor * erfc;
              if (sbindex) {
                const flt_t adjust = ((flt_t)1.0 - special_coul[sbindex])*
                    prefactor;
                forcecoul -= adjust;
                if (EFLAG) ecoul -= adjust;
              }
            #ifdef INTEL_ALLOW_TABLE
            } else {
              float rsq_lookup = rsq;
              const int itable = (__intel_castf32_u32(rsq_lookup) &
                  ncoulmask) >> ncoulshiftbits;
              const flt_t fraction = (rsq_lookup - table[itable].r) *
                  table[itable].dr;

              const flt_t tablet = table[itable].f +
                  fraction * table[itable].df;
              forcecoul = qtmp * q[j] * tablet;
              if (EFLAG) ecoul = qtmp * q[j] * (etable[itable] +
                  fraction * detable[itable]);
              if (sbindex) {
                const flt_t table2 = ctable[itable] +
                    fraction * dctable[itable];
                const flt_t prefactor = qtmp * q[j] * table2;
                const flt_t adjust = ((flt_t)1.0 - special_coul[sbindex]) *
                    prefactor;
                forcecoul -= adjust;
                if (EFLAG) ecoul -= adjust;
              }
            }
            #endif
	  #ifdef __MIC__
	  }
	  #endif

	  #ifdef __MIC__
	  if (rsq < cut_ljsq) {
	  #endif
            flt_t r6inv = r2inv * r2inv * r2inv;
            forcelj = r6inv * (lji[jtype].x * r6inv - lji[jtype].y);
            if (EFLAG) evdwl = r6inv*(lji[jtype].z * r6inv - lji[jtype].w);

	    #ifdef __MIC__
	    if (rsq > cut_lj_innersq) {
	    #endif
              const flt_t drsq = cut_ljsq - rsq;
              const flt_t cut2 = (rsq - cut_lj_innersq) * drsq;
              const flt_t switch1 = drsq * (drsq * drsq + (flt_t)3.0 * cut2) *
                  inv_denom_lj;
              const flt_t switch2 = (flt_t)12.0 * rsq * cut2 * inv_denom_lj;
              if (EFLAG) {
		#ifndef __MIC__
		if (rsq > cut_lj_innersq) {
		#endif
                  forcelj = forcelj * switch1 + evdwl * switch2;
                  evdwl *= switch1;
		#ifndef __MIC__
		}
		#endif
              } else {
                const flt_t philj = r6inv * (lji[jtype].z*r6inv -
                    lji[jtype].w);
		#ifndef __MIC__
		if (rsq > cut_lj_innersq)
		#endif
                  forcelj =  forcelj * switch1 + philj * switch2;
              }
	    #ifdef __MIC__
	    }
	    #endif

            if (sbindex) {
              const flt_t factor_lj = special_lj[sbindex];
              forcelj *= factor_lj;
              if (EFLAG) evdwl *= factor_lj;
            }
	  #ifdef __MIC__
	  }
	  #else
	  if (rsq > cut_coulsq) { forcecoul = (flt_t)0.0; ecoul = (flt_t)0.0; }
	  if (rsq > cut_ljsq) { forcelj = (flt_t)0.0; evdwl = (flt_t)0.0; }
	  #endif

	  #ifdef __MIC__
	  if (rsq < cut_coulsq) {
	  #endif
            const flt_t fpair = (forcecoul + forcelj) * r2inv;
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
              if (NEWTON_PAIR || i < nlocal)
                ev_pre += (flt_t)0.5;
              if (NEWTON_PAIR || j < nlocal)
                ev_pre += (flt_t)0.5;

              if (EFLAG) {
                sevdwl += ev_pre * evdwl;
                secoul += ev_pre * ecoul;
                if (eatom) {
                  if (NEWTON_PAIR || i < nlocal)
                    fwtmp += (flt_t)0.5 * evdwl + (flt_t)0.5 * ecoul;
                  if (NEWTON_PAIR || j < nlocal) 
                    f[j].w += (flt_t)0.5 * evdwl + (flt_t)0.5 * ecoul;
                }
              }

	      IP_PRE_ev_tally_nbor(vflag, ev_pre, fpair,
				   delx, dely, delz);
            }
	  #ifdef __MIC__
	  }
	  #endif
        } // for jj
        f[i].x += fxtmp;
        f[i].y += fytmp;
        f[i].z += fztmp;

	IP_PRE_ev_tally_atomq(EVFLAG, EFLAG, vflag, f, fwtmp);
      } // for ii

      #if defined(_OPENMP)
      #pragma omp barrier
      #endif
      IP_PRE_fdotr_acc_force(NEWTON_PAIR, EVFLAG,  EFLAG, vflag, eatom, nall,
			     nlocal, minlocal, nthreads, f_start, f_stride, 
			     x);
    } // end of omp parallel region
    if (EVFLAG) {
      if (EFLAG) {
        ev_global[0] = oevdwl;
        ev_global[1] = oecoul;
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
  } // end of offload region

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

void PairLJCharmmCoulLongIntel::init_style()
{
  PairLJCharmmCoulLong::init_style();
  neighbor->requests[neighbor->nrequest-1]->intel = 1;

  int ifix = modify->find_fix("package_intel");
  if (ifix < 0)
    error->all(FLERR,
               "The 'package intel' command is required for /intel styles");
  fix = static_cast<FixIntel *>(modify->fix[ifix]);
  
  #ifdef _LMP_INTEL_OFFLOAD
  fix->set_offload_affinity();
  _cop = fix->coprocessor_number();
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

template <class flt_t, class acc_t>
void PairLJCharmmCoulLongIntel::pack_force_const(ForceConst<flt_t> &fc,
                                          IntelBuffers<flt_t,acc_t> *buffers)
{
  int tp1 = atom->ntypes + 1;
  int ntable = 1;
  if (ncoultablebits)
    for (int i = 0; i < ncoultablebits; i++) ntable *= 2;

  fc.set_ntypes(tp1, ntable, memory, _cop);
  buffers->set_ntypes(tp1);
  flt_t **cutneighsq = buffers->get_cutneighsq();

  // Repeat cutsq calculation because done after call to init_style
  double cut, cutneigh;
  if (cut_lj > cut_coul)
    error->all(FLERR,
	 "Intel varient of lj/charmm/coul/long expects lj cutoff<=coulombic");
  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = i; j <= atom->ntypes; j++) {
      if (setflag[i][j] != 0 || (setflag[i][i] != 0 && setflag[j][j] != 0)) {
        cut = init_one(i, j);
        cutneigh = cut + neighbor->skin;
        cutsq[i][j] = cutsq[j][i] = cut*cut;
        cutneighsq[i][j] = cutneighsq[j][i] = cutneigh * cutneigh;
      }
    }
  }

  cut_lj_innersq = cut_lj_inner * cut_lj_inner;
  cut_ljsq = cut_lj * cut_lj;
  cut_coulsq = cut_coul * cut_coul;
  cut_bothsq = MAX(cut_ljsq, cut_coulsq);

  fc.g_ewald = force->kspace->g_ewald;
  fc.tabinnersq = tabinnersq;
  fc.cut_coulsq = cut_coulsq;
  fc.cut_ljsq = cut_ljsq;
  fc.cut_lj_innersq = cut_lj_innersq;

  for (int i = 0; i < 4; i++) {
    fc.special_lj[i] = force->special_lj[i];
    fc.special_coul[i] = force->special_coul[i];
    fc.special_coul[0] = 1.0;
    fc.special_lj[0] = 1.0;
  }

  for (int i = 0; i < tp1; i++) {
    for (int j = 0; j < tp1; j++) {
      fc.lj[i][j].x = lj1[i][j];
      fc.lj[i][j].y = lj2[i][j];
      fc.lj[i][j].z = lj3[i][j];
      fc.lj[i][j].w = lj4[i][j];
      fc.cutsq[i][j] = cutsq[i][j];
    }
  }

  if (ncoultablebits) {
    for (int i = 0; i < ntable; i++) {
      fc.table[i].r = rtable[i];
      fc.table[i].dr = drtable[i];
      fc.table[i].f = ftable[i];
      fc.table[i].df = dftable[i];
      fc.etable[i] = etable[i];
      fc.detable[i] = detable[i];
      fc.ctable[i] = ctable[i];
      fc.dctable[i] = dctable[i];
    }
  }

  #ifdef _LMP_INTEL_OFFLOAD
  if (_cop < 0) return;
  flt_t * special_lj = fc.special_lj;
  flt_t * special_coul = fc.special_coul;
  flt_t * cutsq = fc.cutsq[0];
  LJ_T * lj = fc.lj[0];
  TABLE_T * table = fc.table;
  flt_t * etable = fc.etable;
  flt_t * detable = fc.detable;
  flt_t * ctable = fc.ctable;
  flt_t * dctable = fc.dctable;
  flt_t * ocutneighsq = cutneighsq[0];
  int tp1sq = tp1 * tp1;
  #pragma offload_transfer target(mic:_cop) \
    in(special_lj, special_coul: length(4) alloc_if(0) free_if(0)) \
    in(cutsq,lj: length(tp1sq) alloc_if(0) free_if(0)) \
    in(table: length(ntable) alloc_if(0) free_if(0)) \
    in(etable,detable,ctable,dctable: length(ntable) alloc_if(0) free_if(0)) \
    in(ocutneighsq: length(tp1sq) alloc_if(0) free_if(0))
  #endif
}

/* ---------------------------------------------------------------------- */

template <class flt_t>
void PairLJCharmmCoulLongIntel::ForceConst<flt_t>::set_ntypes(const int ntypes,
                                                              const int ntable,
                                                              Memory *memory,
							      const int cop) {
  if ( (ntypes != _ntypes || ntable != _ntable) ) {
    if (_ntypes > 0) {
      #ifdef _LMP_INTEL_OFFLOAD
      flt_t * ospecial_lj = special_lj;
      flt_t * ospecial_coul = special_coul;
      flt_t * ocutsq = cutsq[0];
      typename IntelBuffers<flt_t,flt_t>::vec4_t * olj = lj[0];
      table_t * otable = table;
      flt_t * oetable = etable;
      flt_t * odetable = detable;
      flt_t * octable = ctable;
      flt_t * odctable = dctable;
      if (ospecial_lj != NULL && ocutsq != NULL && olj != NULL &&
          otable != NULL && oetable != NULL && odetable != NULL &&
          octable != NULL && odctable != NULL && ospecial_coul != NULL &&
	  cop >= 0) {
        #pragma offload_transfer target(mic:cop) \
          nocopy(ospecial_lj, ospecial_coul: alloc_if(0) free_if(1)) \
	  nocopy(ocutsq, olj: alloc_if(0) free_if(1)) \
	  nocopy(otable: alloc_if(0) free_if(1)) \
	  nocopy(oetable, odetable, octable, odctable: alloc_if(0) free_if(1))
      }
      #endif

      _memory->destroy(cutsq);
      _memory->destroy(lj);
      _memory->destroy(table);
      _memory->destroy(etable);
      _memory->destroy(detable);
      _memory->destroy(ctable);
      _memory->destroy(dctable);
    }
    if (ntypes > 0) {
      _cop = cop;
      memory->create(cutsq,ntypes,ntypes,"fc.cutsq");
      memory->create(lj,ntypes,ntypes,"fc.lj");
      memory->create(table,ntable,"pair:fc.table");
      memory->create(etable,ntable,"pair:fc.etable");
      memory->create(detable,ntable,"pair:fc.detable");
      memory->create(ctable,ntable,"pair:fc.ctable");
      memory->create(dctable,ntable,"pair:fc.dctable");

      #ifdef _LMP_INTEL_OFFLOAD
      flt_t * ospecial_lj = special_lj;
      flt_t * ospecial_coul = special_coul;
      flt_t * ocutsq = cutsq[0];
      typename IntelBuffers<flt_t,flt_t>::vec4_t * olj = lj[0];
      table_t * otable = table;
      flt_t * oetable = etable;
      flt_t * odetable = detable;
      flt_t * octable = ctable;
      flt_t * odctable = dctable;
      int tp1sq = ntypes*ntypes;
      if (ospecial_lj != NULL && ocutsq != NULL && olj != NULL &&
          otable !=NULL && oetable != NULL && odetable != NULL &&
          octable != NULL && odctable != NULL && ospecial_coul != NULL &&
	  cop >= 0) {
        #pragma offload_transfer target(mic:cop) \
          nocopy(ospecial_lj: length(4) alloc_if(1) free_if(0)) \
          nocopy(ospecial_coul: length(4) alloc_if(1) free_if(0)) \
          nocopy(ocutsq,olj: length(tp1sq) alloc_if(1) free_if(0)) \
          nocopy(otable: length(ntable) alloc_if(1) free_if(0)) \
          nocopy(oetable,odetable: length(ntable) alloc_if(1) free_if(0)) \
          nocopy(octable,odctable: length(ntable) alloc_if(1) free_if(0))
      }
      #endif
    }
  }
  _ntypes=ntypes;
  _ntable=ntable;
  _memory=memory;
}
