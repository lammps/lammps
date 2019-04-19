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
   Contributing author: Markus Höhnerbach (RWTH)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_tersoff_intel.h"
#include "atom.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

// Currently Intel compiler is required for this pair style.
// For convenience, base class routines are called if not using Intel compiler.
#ifndef __INTEL_COMPILER
using namespace LAMMPS_NS;

PairTersoffIntel::PairTersoffIntel(LAMMPS *lmp) : PairTersoff(lmp)
{
}

void PairTersoffIntel::compute(int eflag, int vflag)
{
  PairTersoff::compute(eflag, vflag);
}

void PairTersoffIntel::init_style()
{
  if (comm->me == 0) {
    error->warning(FLERR, "Tersoff/intel currently requires intel compiler. "
                   "Using MANYBODY version.");
  }
  PairTersoff::init_style();
}

#else

#ifdef _LMP_INTEL_OFFLOAD
#pragma offload_attribute(push,target(mic))
#endif

#include "intel_intrinsics.h"
#include "math_const.h"

#ifdef _LMP_INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif

#include "group.h"
#include "kspace.h"
#include "modify.h"
#include "suffix.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairTersoffIntel::PairTersoffIntel(LAMMPS *lmp) : PairTersoff(lmp)
{
  suffix_flag |= Suffix::INTEL;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

// Dispatch the requested precision
void PairTersoffIntel::compute(int eflag, int vflag)
{
  if (fix->precision()==FixIntel::PREC_MODE_MIXED) {
    compute<float,double>(eflag, vflag, fix->get_mixed_buffers(),
                          force_const_single);
  } else if (fix->precision()==FixIntel::PREC_MODE_DOUBLE) {
    compute<double,double>(eflag, vflag, fix->get_double_buffers(),
                           force_const_double);
  } else {
    compute<float,float>(eflag, vflag, fix->get_single_buffers(),
                         force_const_single);
  }
  fix->balance_stamp();
  vflag_fdotr = 0;
}

// Dispatch the extent of computation:
//  do we need to calculate energy/virial
template <class flt_t, class acc_t>
void PairTersoffIntel::compute(int eflag, int vflag,
                                     IntelBuffers<flt_t,acc_t> *buffers,
                                     const ForceConst<flt_t> &fc)
{
  ev_init(eflag,vflag);
  if (vflag_atom)
    error->all(FLERR,"USER-INTEL package does not support per-atom stress");

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
  if (eflag) {
    eval<1>(1, ovflag, buffers, fc, 0, offload_end);
    eval<1>(0, ovflag, buffers, fc, host_start, inum);
  } else {
    eval<0>(1, ovflag, buffers, fc, 0, offload_end);
    eval<0>(0, ovflag, buffers, fc, host_start, inum);
  }
}

#ifdef _LMP_INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

// The complete Tersoff computation kernel is encapsulated here
//  everything is static, the class just serves as a unit of organization
template<class flt_t, class acc_t, lmp_intel::CalculationMode mic, bool pack_i>
struct IntelKernelTersoff : public lmp_intel::vector_routines<flt_t, acc_t, mic> {
  // instantiate the vector library and import the types
  typedef typename lmp_intel::vector_routines<flt_t, acc_t, mic> v;
  typedef typename v::fvec fvec;
  typedef typename v::ivec ivec;
  typedef typename v::bvec bvec;
  typedef typename v::avec avec;
  typedef typename v::iarr iarr;
  typedef typename v::farr farr;
  typedef typename v::aarr aarr;
  typedef typename PairTersoffIntel::ForceConst<flt_t>::c_inner_t c_inner_t;
  typedef typename PairTersoffIntel::ForceConst<flt_t>::c_outer_t c_outer_t;

  // for descriptions of these methods, please have a look at the original code
  // what's done in here is that they are inlined and vectorized
  // attractive() also provides an option to compute zeta as well
  static fvec zeta_vector(
      const c_inner_t * param,
      ivec xjw, bvec mask,
      fvec vrij, fvec rsq2,
      fvec vdijx, fvec vdijy, fvec vdijz,
      fvec dikx, fvec diky, fvec dikz
  );
  static void force_zeta_vector(
      const c_outer_t * param,
      ivec xjw,
      bvec mask,
      fvec vrijsq, fvec vzeta_ij,
      fvec *vfpair, fvec *vprefactor, int EVDWL, fvec *vevdwl,
      bvec vmask_repulsive
  );
  template<bool ZETA>
  static void attractive_vector(
      const c_inner_t * param,
      ivec xjw,
      bvec mask,
      fvec vprefactor,
      fvec vrijsq, fvec rsq2,
      fvec vdijx, fvec vdijy, fvec vdijz,
      fvec dikx, fvec diky, fvec dikz,
      fvec *fix, fvec *fiy, fvec *fiz,
      fvec *fjx, fvec *fjy, fvec *fjz,
      fvec *fkx, fvec *fky, fvec *fkz,
      fvec *zeta
  );

  // perform the actual computation
  template<bool EFLAG>
  static void kernel(
      int iito, int iifrom, int eatom, int vflag,
      const int * _noalias const numneigh,
      const int * _noalias const numneighhalf,
      const int * _noalias const cnumneigh,
      const int * _noalias const firstneigh, int ntypes,
      typename IntelBuffers<flt_t,acc_t>::atom_t * _noalias const x,
      const int * _noalias const ilist,
      const c_inner_t * _noalias const c_inner,
      const c_outer_t * _noalias const c_outer,
      typename IntelBuffers<flt_t,acc_t>::vec3_acc_t * _noalias const f,
      acc_t *evdwl
  );

  // perform one step of calculation, pass in i-j pairs of atoms (is, js)
  template<int EFLAG>
  static void kernel_step(
      int eatom, int vflag,
      const int * _noalias const numneigh,
      iarr vcnumneigh,
      const int * _noalias const firstneigh,
      int ntypes,
      typename IntelBuffers<flt_t,acc_t>::atom_t * _noalias const x,
      const c_inner_t * _noalias const c_inner,
      const c_outer_t * _noalias const c_outer,
      typename IntelBuffers<flt_t,acc_t>::vec3_acc_t * _noalias const f,
      avec *vsevdwl, int compress_idx, iarr is, iarr js, bvec vmask_repulsive
  );

  // perform one step of calculation, as opposed to the previous method now
  //  with fixed i and a number of js
  template<int EFLAG>
  static void kernel_step_const_i(
    int eatom, int vflag,
    const int * _noalias const numneigh, const int * _noalias const cnumneigh,
    const int * _noalias const firstneigh, int ntypes,
    typename IntelBuffers<flt_t,acc_t>::atom_t * _noalias const x,
    const c_inner_t * _noalias const c_inner,
    const c_outer_t * _noalias const c_outer,
    typename IntelBuffers<flt_t,acc_t>::vec3_acc_t * _noalias const f,
    avec *vsevdwl, int compress_idx, int ii, int i, iarr js,
    bvec vmask_repulsive
  );
};

#ifdef _LMP_INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif

/* ---------------------------------------------------------------------- */

// Dispatch to correct kernel instatiation and perform all the work neccesary
//  for offloading. In this routine we enter the Phi.
// This method is nearly identical to what happens in the other /intel styles
template <int EFLAG, class flt_t, class acc_t>
void PairTersoffIntel::eval(const int offload, const int vflag,
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
  tagint * _noalias tag = this->atom->tag;
  flt_t * _noalias const q = buffers->get_q(offload);

  const int * _noalias const ilist = list->ilist;
  const int * _noalias const numneigh = list->numneigh;
  const int * _noalias const firstneigh = list->firstneigh[ilist[0]];

  int *nhalf, *cnum;
  buffers->get_list_data3(list, nhalf, cnum);
  const int * _noalias const numneighhalf = nhalf;
  const int * _noalias const cnumneigh = cnum;

  typedef typename ForceConst<flt_t>::c_inner_t c_inner_t;
  typedef typename ForceConst<flt_t>::c_outer_t c_outer_t;
  typedef typename ForceConst<flt_t>::c_cutoff_t c_cutoff_t;
  const c_outer_t * _noalias const c_outer = fc.c_outer[0];
  const c_inner_t * _noalias const c_inner = fc.c_inner[0][0];
  const c_cutoff_t * _noalias const c_inner_cutoff = fc.c_cutoff_inner[0][0];

  const int ntypes = atom->ntypes + 1;
  const int eatom = this->eflag_atom;

  // Determine how much data to transfer
  int x_size, q_size, f_stride, ev_size, separate_flag;
  IP_PRE_get_transfern(ago, 1, EFLAG, vflag,
                       buffers, offload, fix, separate_flag,
                       x_size, q_size, ev_size, f_stride);

  int tc;
  FORCE_T * _noalias f_start;
  acc_t * _noalias ev_global;
  IP_PRE_get_buffers(offload, buffers, fix, tc, f_start, ev_global);

  const int nthreads = tc;

  #ifdef _LMP_INTEL_OFFLOAD
  int *overflow = fix->get_off_overflow_flag();
  double *timer_compute = fix->off_watch_pair();
  if (offload) fix->start_watch(TIME_OFFLOAD_LATENCY);
  #pragma offload target(mic:_cop) if(offload) \
    in(c_inner, c_outer :length(0) alloc_if(0) free_if(0)) \
    in(c_inner_cutoff :length(0) alloc_if(0) free_if(0)) \
    in(firstneigh:length(0) alloc_if(0) free_if(0)) \
    in(cnumneigh:length(0) alloc_if(0) free_if(0)) \
    in(numneigh:length(0) alloc_if(0) free_if(0)) \
    in(numneighhalf:length(0) alloc_if(0) free_if(0)) \
    in(x:length(x_size) alloc_if(0) free_if(0)) \
    in(ilist:length(0) alloc_if(0) free_if(0)) \
    in(overflow:length(0) alloc_if(0) free_if(0)) \
    in(astart,nthreads,inum,nall,ntypes,vflag,eatom) \
    in(f_stride,nlocal,minlocal,separate_flag,offload) \
    out(f_start:length(f_stride) alloc_if(0) free_if(0)) \
    out(ev_global:length(ev_size) alloc_if(0) free_if(0)) \
    out(timer_compute:length(1) alloc_if(0) free_if(0)) \
    signal(f_start)
  #endif
  {
    #ifdef _LMP_INTEL_OFFLOAD
    #ifdef __MIC__
    *timer_compute = MIC_Wtime();
    #endif
    #endif

    IP_PRE_repack_for_offload(1, separate_flag, nlocal, nall,
                              f_stride, x, 0);

    acc_t oevdwl, ov0, ov1, ov2, ov3, ov4, ov5;
    if (EFLAG || vflag)
      oevdwl = ov0 = ov1 = ov2 = ov3 = ov4 = ov5 = (acc_t)0;

    // loop over neighbors of my atoms
    #if defined(_OPENMP)
    #pragma omp parallel reduction(+:oevdwl,ov0,ov1,ov2,ov3,ov4,ov5)
    #endif
    {
      int iifrom, iito, tid;
      IP_PRE_omp_range_id(iifrom, iito, tid, inum, nthreads);
      iifrom += astart;
      iito += astart;

      FORCE_T * _noalias const f = f_start - minlocal + (tid * f_stride);
      memset(f + minlocal, 0, f_stride * sizeof(FORCE_T));

      {
        acc_t sevdwl;
        sevdwl = 0.;
        #define ARGS iito, iifrom, eatom, vflag, numneigh, numneighhalf, \
                     cnumneigh, firstneigh, ntypes, x, ilist, c_inner, \
                     c_outer, f, &sevdwl
        // Pick the variable i algorithm under specific conditions
        // do use scalar algorithm with very short vectors
        int VL = lmp_intel::vector_routines<flt_t,acc_t,lmp_intel::mode>::VL;
        bool pack_i = VL >= 8 &&
          lmp_intel::vector_traits<lmp_intel::mode>::support_integer_and_gather_ops;
        bool use_scalar = VL < 4;
        if (use_scalar) {
          IntelKernelTersoff<flt_t,acc_t,lmp_intel::NONE,false>::kernel<EFLAG>(ARGS);
        } else if (pack_i) {
          IntelKernelTersoff<flt_t,acc_t,lmp_intel::mode,true >::kernel<EFLAG>(ARGS);
        } else {
          IntelKernelTersoff<flt_t,acc_t,lmp_intel::mode,false>::kernel<EFLAG>(ARGS);
        }
        if (EFLAG) oevdwl += sevdwl;
      }

      IP_PRE_fdotr_reduce_omp(1, nall, minlocal, nthreads, f_start,
                              f_stride, x, offload, vflag, ov0, ov1, ov2, ov3,
                              ov4, ov5);
    } // end of omp parallel region

    IP_PRE_fdotr_reduce(1, nall, nthreads, f_stride, vflag,
                        ov0, ov1, ov2, ov3, ov4, ov5);

    if (EFLAG || vflag) {
      ev_global[0] = oevdwl;
      ev_global[1] = 0.0;
      ev_global[2] = ov0;
      ev_global[3] = ov1;
      ev_global[4] = ov2;
      ev_global[5] = ov3;
      ev_global[6] = ov4;
      ev_global[7] = ov5;
    }

    #ifdef _LMP_INTEL_OFFLOAD
    #ifdef __MIC__
    *timer_compute = MIC_Wtime() - *timer_compute;
    #endif
    #endif
  } // end of offload region

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
   init specific to this pair style
------------------------------------------------------------------------- */

// As in any other /intel pair style
void PairTersoffIntel::init_style()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style Tersoff requires atom IDs");
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style Tersoff requires newton pair on");

  // need a full neighbor list

  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->intel = 1;

  int ifix = modify->find_fix("package_intel");
  if (ifix < 0)
    error->all(FLERR,
               "The 'package intel' command is required for /intel styles");
  fix = static_cast<FixIntel *>(modify->fix[ifix]);

  fix->pair_init_check();
  fix->three_body_neighbor(1);
  #ifdef _LMP_INTEL_OFFLOAD
  _cop = fix->coprocessor_number();
  #endif
  if (fix->precision() == FixIntel::PREC_MODE_MIXED) {
    pack_force_const(force_const_single, fix->get_mixed_buffers());
    fix->get_mixed_buffers()->need_tag(1);
  } else if (fix->precision() == FixIntel::PREC_MODE_DOUBLE) {
    fix->get_double_buffers()->need_tag(1);
    pack_force_const(force_const_double, fix->get_double_buffers());
  } else {
    pack_force_const(force_const_single, fix->get_single_buffers());
    fix->get_single_buffers()->need_tag(1);
  }
}

// As in any other /intel pair style
template <class flt_t, class acc_t>
void PairTersoffIntel::pack_force_const(ForceConst<flt_t> &fc,
                                          IntelBuffers<flt_t,acc_t> *buffers)
{
  int tp1 = atom->ntypes + 1;

  fc.set_ntypes(tp1, memory, _cop);

  // Repeat cutsq calculation because done after call to init_style
  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = i; j <= atom->ntypes; j++) {
      double cut;
      if (setflag[i][j] != 0 || (setflag[i][i] != 0 && setflag[j][j] != 0))
        cut = init_one(i, j);
      else
        cut = 0.0;
      cutsq[i][j] = cutsq[j][i] = cut*cut;
    }
  }

  for (int i = 1; i < tp1; i++) {
    for (int j = 1; j < tp1; j++) {
      fc.c_inner_loop[i][j][0].d2 = 1.0;
      fc.c_inner_loop[i][0][j].d2 = 1.0;
      fc.c_inner_loop[0][i][j].d2 = 1.0;
      for (int k = 1; k < tp1; k++) {
        Param * param = &params[elem2param[map[i]][map[j]][map[k]]];
        fc.c_cutoff_inner[i][k][j].cutsq = static_cast<flt_t>(param->cutsq);
        fc.c_inner_loop[i][j][k].lam3 = static_cast<flt_t>(param->lam3);
        fc.c_inner_loop[i][j][k].bigr = static_cast<flt_t>(param->bigr);
        fc.c_inner_loop[i][j][k].bigd = static_cast<flt_t>(param->bigd);
        fc.c_inner_loop[i][j][k].c2 = static_cast<flt_t>(param->c * param->c);
        fc.c_inner_loop[i][j][k].d2 = static_cast<flt_t>(param->d * param->d);
        fc.c_inner_loop[i][j][k].h = static_cast<flt_t>(param->h);
        fc.c_inner_loop[i][j][k].gamma = static_cast<flt_t>(param->gamma);
        fc.c_inner_loop[i][j][k].powermint = static_cast<flt_t>(param->powermint);

        fc.c_inner[i][j][k].cutsq = static_cast<flt_t>(param->cutsq);
        fc.c_inner[i][j][k].lam3 = static_cast<flt_t>(param->lam3);
        fc.c_inner[i][j][k].bigr = static_cast<flt_t>(param->bigr);
        fc.c_inner[i][j][k].bigd = static_cast<flt_t>(param->bigd);
        fc.c_inner[i][j][k].c2 = static_cast<flt_t>(param->c * param->c);
        fc.c_inner[i][j][k].d2 = static_cast<flt_t>(param->d * param->d);
        fc.c_inner[i][j][k].h = static_cast<flt_t>(param->h);
        fc.c_inner[i][j][k].gamma = static_cast<flt_t>(param->gamma);
        fc.c_inner[i][j][k].powermint = static_cast<flt_t>(param->powermint);

      }
      Param * param = &params[elem2param[map[i]][map[j]][map[j]]];
      fc.c_cutoff_outer[i][j].cutsq = static_cast<flt_t>(param->cutsq);
      fc.c_first_loop[i][j].bigr = static_cast<flt_t>(param->bigr);
      fc.c_first_loop[i][j].bigd = static_cast<flt_t>(param->bigd);
      fc.c_first_loop[i][j].lam1 = static_cast<flt_t>(param->lam1);
      fc.c_first_loop[i][j].biga = static_cast<flt_t>(param->biga);
      fc.c_second_loop[i][j].lam2 = static_cast<flt_t>(param->lam2);
      fc.c_second_loop[i][j].beta = static_cast<flt_t>(param->beta);
      fc.c_second_loop[i][j].bigb = static_cast<flt_t>(param->bigb);
      fc.c_second_loop[i][j].powern = static_cast<flt_t>(param->powern);
      fc.c_second_loop[i][j].c1 = static_cast<flt_t>(param->c1);
      fc.c_second_loop[i][j].c2 = static_cast<flt_t>(param->c2);
      fc.c_second_loop[i][j].c3 = static_cast<flt_t>(param->c3);
      fc.c_second_loop[i][j].c4 = static_cast<flt_t>(param->c4);

      fc.c_outer[i][j].cutsq = static_cast<flt_t>(param->cutsq);
      fc.c_outer[i][j].bigr = static_cast<flt_t>(param->bigr);
      fc.c_outer[i][j].bigd = static_cast<flt_t>(param->bigd);
      fc.c_outer[i][j].lam1 = static_cast<flt_t>(param->lam1);
      fc.c_outer[i][j].biga = static_cast<flt_t>(param->biga);
      fc.c_outer[i][j].lam2 = static_cast<flt_t>(param->lam2);
      fc.c_outer[i][j].beta = static_cast<flt_t>(param->beta);
      fc.c_outer[i][j].bigb = static_cast<flt_t>(param->bigb);
      fc.c_outer[i][j].powern = static_cast<flt_t>(param->powern);
      fc.c_outer[i][j].c1 = static_cast<flt_t>(param->c1);
      fc.c_outer[i][j].c2 = static_cast<flt_t>(param->c2);
      fc.c_outer[i][j].c3 = static_cast<flt_t>(param->c3);
      fc.c_outer[i][j].c4 = static_cast<flt_t>(param->c4);

    }
    fc.c_outer[i][0].cutsq = 0.;
  }

  #ifdef _LMP_INTEL_OFFLOAD
  if (_cop < 0) return;
  typename ForceConst<flt_t>::c_first_loop_t * c_first_loop = fc.c_first_loop[0];
  typename ForceConst<flt_t>::c_cutoff_t * c_cutoff_outer = fc.c_cutoff_outer[0];
  typename ForceConst<flt_t>::c_outer_t * c_outer = fc.c_outer[0];
  typename ForceConst<flt_t>::c_second_loop_t * c_second_loop = fc.c_second_loop[0];
  typename ForceConst<flt_t>::c_inner_loop_t * c_inner_loop = fc.c_inner_loop[0][0];
  typename ForceConst<flt_t>::c_cutoff_t * c_cutoff_inner = fc.c_cutoff_inner[0][0];
  typename ForceConst<flt_t>::c_inner_t * c_inner = fc.c_inner[0][0];
  size_t VL = 512 / 8 / sizeof(flt_t);
  int ntypes = tp1;
  int ntypes_pad = ntypes + VL - ntypes % VL;
  int tp1sq = tp1 * tp1;
  int tp1cb = tp1 * tp1 * tp1;
  int tp1cb_pad = tp1 * tp1 * ntypes_pad;
  #pragma offload_transfer target(mic:_cop) \
    in(c_first_loop, c_second_loop, c_cutoff_outer, c_outer : length(tp1sq) alloc_if(0) free_if(0)) \
    in(c_inner : length(tp1cb) alloc_if(0) free_if(0)) \
    in(c_cutoff_inner : length(tp1cb_pad) alloc_if(0) free_if(0))
  #endif
}

/* ---------------------------------------------------------------------- */

// As in any other /intel pair style
template <class flt_t>
void PairTersoffIntel::ForceConst<flt_t>::set_ntypes(const int ntypes,
                                                           Memory *memory,
                                                           const int cop) {
  if ( (ntypes != _ntypes) ) {
    if (_ntypes > 0) {
      #ifdef _LMP_INTEL_OFFLOAD
      c_first_loop_t * oc_first_loop = c_first_loop[0];
      c_second_loop_t * oc_second_loop = c_second_loop[0];
      c_inner_loop_t * oc_inner_loop = c_inner_loop[0][0];
      c_cutoff_t * oc_cutoff_inner = c_cutoff_inner[0][0];
      c_cutoff_t * oc_cutoff_outer = c_cutoff_outer[0];
      c_inner_t * oc_inner = c_inner[0][0];
      c_outer_t * oc_outer = c_outer[0];
      if (c_first_loop != NULL && c_second_loop != NULL &&
          c_inner_loop != NULL &&  _cop >= 0) {

        #pragma offload_transfer target(mic:cop) \
          nocopy(oc_first_loop, oc_second_loop, oc_inner_loop: alloc_if(0) free_if(1)) \
          nocopy(oc_cutoff_outer, oc_cutoff_inner: alloc_if(0) free_if(1)) \
          nocopy(oc_inner, oc_outer: alloc_if(0) free_if(0))
      }
      #endif
      _memory->destroy(c_first_loop);
      _memory->destroy(c_second_loop);
      _memory->destroy(c_inner_loop);
      _memory->destroy(c_cutoff_outer);
      _memory->destroy(c_cutoff_inner);
      _memory->destroy(c_inner);
      _memory->destroy(c_outer);
    }
    if (ntypes > 0) {
      _cop = cop;
      size_t VL = 512 / 8 / sizeof(flt_t);
      int ntypes_pad = ntypes + VL - ntypes % VL;
      memory->create(c_first_loop,ntypes,ntypes,"fc.c_first_loop");
      memory->create(c_second_loop,ntypes,ntypes,"fc.c_second_loop");
      memory->create(c_cutoff_outer,ntypes,ntypes,"fc.c_cutoff_outer");
      memory->create(c_inner_loop,ntypes,ntypes,ntypes,"fc.c_inner_loop");
      memory->create(c_cutoff_inner,ntypes,ntypes,ntypes_pad,"fc.c_cutoff_inner");
      memory->create(c_inner,ntypes,ntypes,ntypes,"fc.c_inner");
      memory->create(c_outer,ntypes,ntypes,"fc.c_outer");
      #ifdef _LMP_INTEL_OFFLOAD
      c_first_loop_t * oc_first_loop = c_first_loop[0];
      c_second_loop_t * oc_second_loop = c_second_loop[0];
      c_cutoff_t * oc_cutoff_outer = c_cutoff_outer[0];
      c_inner_loop_t * oc_inner_loop = c_inner_loop[0][0];
      c_cutoff_t * oc_cutoff_inner = c_cutoff_inner[0][0];
      c_inner_t * oc_inner = c_inner[0][0];
      c_outer_t * oc_outer = c_outer[0];
      int tp1sq = ntypes * ntypes;
      int tp1cb = ntypes * ntypes * ntypes;
      int tp1cb_pad = ntypes * ntypes * ntypes_pad;
      if (oc_first_loop != NULL && oc_second_loop != NULL &&
          oc_inner_loop != NULL && cop >= 0) {
        #pragma offload_transfer target(mic:cop) \
          nocopy(oc_first_loop: length(tp1sq) alloc_if(1) free_if(0)) \
          nocopy(oc_second_loop: length(tp1sq) alloc_if(1) free_if(0)) \
          nocopy(oc_cutoff_outer: length(tp1sq) alloc_if(1) free_if(0)) \
          nocopy(oc_outer: length(tp1sq) alloc_if(1) free_if(0)) \
          nocopy(oc_inner_loop: length(tp1cb) alloc_if(1) free_if(0)) \
          nocopy(oc_inner: length(tp1cb) alloc_if(1) free_if(0)) \
          nocopy(oc_cutoff_inner: length(tp1cb_pad) alloc_if(1) free_if(0))
      }
      #endif
    }
  }
  _ntypes=ntypes;
  _memory=memory;
}

#ifdef _LMP_INTEL_OFFLOAD
#pragma offload_attribute(push,target(mic))
#endif

// The factor up to which we do caching
static const int N_CACHE = 8;

template<class flt_t, class acc_t, lmp_intel::CalculationMode mic, bool pack_i>
template<int EFLAG>
void IntelKernelTersoff<flt_t, acc_t, mic, pack_i>::kernel_step(
    int eatom, int vflag,
    const int * _noalias const numneigh, iarr vcnumneigh,
    const int * _noalias const firstneigh, int ntypes,
    typename IntelBuffers<flt_t,acc_t>::atom_t * _noalias const x,
    const typename PairTersoffIntel::ForceConst<flt_t>::c_inner_t * _noalias const c_inner,
    const typename PairTersoffIntel::ForceConst<flt_t>::c_outer_t * _noalias const c_outer,
    typename IntelBuffers<flt_t,acc_t>::vec3_acc_t * _noalias const f,
    avec *vsevdwl,
    int compress_idx,
    iarr is,
    iarr js,
    bvec vmask_repulsive
) {
  ivec v_i4floats((int) (4 * sizeof(typename v::fscal)));
  ivec v_i1(1);
  fvec v_2(0.);
  fvec v_0_5(0.5);
  ivec v_i0(0);
  ivec v_i_ntypes(ntypes);
  ivec v_i_NEIGHMASK(NEIGHMASK);

  farr fx, fy, fz, fw;
  int cache_idx = 0;
  fvec vfkx_cache[N_CACHE];
  fvec vfky_cache[N_CACHE];
  fvec vfkz_cache[N_CACHE];
  ivec vks_cache[N_CACHE];
  bvec vmask_cache[N_CACHE];
  ivec vkks_final_cache;
  bvec vmask_final_cache;
  iarr ts;
  // compute all the stuff we know from i and j
  // TDO: We could extract this from the driver routine
  ivec vis = v::int_mullo(v_i4floats, v::int_load_vl(is));
  ivec vjs = v::int_mullo(v_i4floats, v::int_load_vl(js));
  ivec vcnumneigh_i = v::int_load_vl(vcnumneigh);
  bvec vmask = v::mask_enable_lower(compress_idx);
  fvec vx_i = v::zero(), vy_i = v::zero(), vz_i = v::zero();
  ivec vw_i = v_i0;
  v::gather_x(vis, vmask, x, &vx_i, &vy_i, &vz_i, &vw_i);
  fvec vx_j = v::zero(), vy_j = v::zero(), vz_j = v::zero();
  ivec vw_j = v_i0;
  v::gather_x(vjs, vmask, x, &vx_j, &vy_j, &vz_j, &vw_j);
  fvec vdx_ij = vx_j - vx_i, vdy_ij = vy_j - vy_i, vdz_ij = vz_j - vz_i;
  fvec vrijsq = vdx_ij * vdx_ij + vdy_ij *  vdy_ij + vdz_ij * vdz_ij;
  fvec vrij = sqrt(vrijsq);
  ivec vis_orig = v::int_load_vl(is);
  ivec vnumneigh_i = v::int_gather<4>(v_i0, vmask, vis_orig, numneigh);
  ivec vc_idx_ij = v::int_mullo(v_i4floats, vw_j + v::int_mullo(v_i_ntypes, vw_i));

  fvec vzeta = v::zero();
  fvec vfxtmp = v::zero(), vfytmp = v::zero(), vfztmp = v::zero();
  fvec vfjxtmp = v::zero(), vfjytmp = v::zero(), vfjztmp = v::zero();
  // This piece of code faciliates the traversal of the k loop assuming
  //  nothing about i. As such, it uses masking to avoid superfluous loads
  //  and fast-forwards each lane until work is available.
  // This is useful because we can not make assumptions as to where in the
  //  neighbor list the atoms within the cutoff might be.
  // We also implement the caching in here, i.e. collect force contributions
  //  due to zeta.
  // This means that you will see four loops:
  // 1. the loop that does zeta calculation and caches the force contributions
  // 2. the loop that processes the remaining zeta calculations
  // 3. the loop that updates the force based on the cached force contributions
  // 4. the loop that computes force contributions for the remainder
  {
    ivec vkks = v_i0;
    bvec vactive_mask = vmask;
    bvec veff_old_mask(0);
    ivec vks, vw_k;
    fvec vx_k, vy_k, vz_k, vcutsq;
    while (! v::mask_testz(vactive_mask) && cache_idx < N_CACHE) {
      bvec vnew_mask = vactive_mask & ~ veff_old_mask;
      vks = v::int_mullo(v_i4floats, v_i_NEIGHMASK &
          v::int_gather<4>(vks, vactive_mask, vkks + vcnumneigh_i, firstneigh));
      v::gather_x(vks, vnew_mask, x, &vx_k, &vy_k, &vz_k, &vw_k);
      fvec vdx_ik = (vx_k - vx_i);
      fvec vdy_ik = (vy_k - vy_i);
      fvec vdz_ik = (vz_k - vz_i);
      fvec vrsq = vdx_ik * vdx_ik + vdy_ik *  vdy_ik + vdz_ik * vdz_ik;
      ivec vc_idx = v::int_mullo(v_i4floats, vw_k) + v::int_mullo(v_i_ntypes, vc_idx_ij);
      vcutsq = v::gather<4>(vcutsq, vnew_mask, vc_idx, c_inner);
      bvec vcutoff_mask = v::cmplt(vrsq, vcutsq);
      bvec vsame_mask = v::int_cmpneq(vjs, vks);
      bvec veff_mask = vcutoff_mask & vsame_mask & vactive_mask;
      if (v::mask_testz(~(veff_mask | ~vactive_mask))) {
        fvec vzeta_contrib;
        fvec vfix, vfiy, vfiz;
        fvec vfjx, vfjy, vfjz;
        fvec vfkx, vfky, vfkz;

        attractive_vector<true>(c_inner,vc_idx,veff_mask,fvec(1.),
            vrij,vrsq,vdx_ij,vdy_ij,vdz_ij,vdx_ik,vdy_ik,vdz_ik,
            &vfix,&vfiy,&vfiz,
            &vfjx,&vfjy,&vfjz,
            &vfkx,&vfky,&vfkz,
            &vzeta_contrib);
        vfxtmp = v::mask_add(vfxtmp, veff_mask, vfxtmp, vfix);
        vfytmp = v::mask_add(vfytmp, veff_mask, vfytmp, vfiy);
        vfztmp = v::mask_add(vfztmp, veff_mask, vfztmp, vfiz);
        vfjxtmp = v::mask_add(vfjxtmp, veff_mask, vfjxtmp, vfjx);
        vfjytmp = v::mask_add(vfjytmp, veff_mask, vfjytmp, vfjy);
        vfjztmp = v::mask_add(vfjztmp, veff_mask, vfjztmp, vfjz);

        vfkx_cache[cache_idx] = vfkx;
        vfky_cache[cache_idx] = vfky;
        vfkz_cache[cache_idx] = vfkz;
        vks_cache[cache_idx] = vks;
        vmask_cache[cache_idx] = veff_mask;
        cache_idx += 1;

        vzeta = v::mask_add(vzeta, veff_mask, vzeta, vzeta_contrib);
        vkks = vkks + v_i1;
        veff_old_mask = bvec(0);
      } else {
        vkks = v::int_mask_add(vkks, ~veff_mask, vkks, v_i1);
        veff_old_mask = veff_mask;
      }
      vactive_mask &= v::int_cmplt(vkks, vnumneigh_i);
    }
    vkks_final_cache = vkks;
    vmask_final_cache = vactive_mask;
    while (! v::mask_testz(vactive_mask)) {
      bvec vnew_mask = vactive_mask & ~ veff_old_mask;
      vks = v::int_mullo(v_i4floats, v_i_NEIGHMASK &
          v::int_gather<4>(vks, vactive_mask, vkks + vcnumneigh_i, firstneigh));
      v::gather_x(vks, vnew_mask, x, &vx_k, &vy_k, &vz_k, &vw_k);
      fvec vdx_ik = (vx_k - vx_i);
      fvec vdy_ik = (vy_k - vy_i);
      fvec vdz_ik = (vz_k - vz_i);
      fvec vrsq = vdx_ik * vdx_ik + vdy_ik *  vdy_ik + vdz_ik * vdz_ik;
      ivec vc_idx = v::int_mullo(v_i4floats, vw_k) + v::int_mullo(v_i_ntypes, vc_idx_ij);
      vcutsq = v::gather<4>(vcutsq, vnew_mask, vc_idx, c_inner);
      bvec vcutoff_mask = v::cmplt(vrsq, vcutsq);
      bvec vsame_mask = v::int_cmpneq(vjs, vks);
      bvec veff_mask = vcutoff_mask & vsame_mask & vactive_mask;
      if (v::mask_testz(~(veff_mask | ~vactive_mask))) {
        fvec vzeta_contrib;
        vzeta_contrib = zeta_vector(c_inner,vc_idx,veff_mask,vrij,vrsq,vdx_ij,vdy_ij,vdz_ij,vdx_ik,vdy_ik,vdz_ik);
        vzeta = v::mask_add(vzeta, veff_mask, vzeta, vzeta_contrib);
        vkks = vkks + v_i1;
        veff_old_mask = bvec(0);
      } else {
        vkks = v::int_mask_add(vkks, ~veff_mask, vkks, v_i1);
        veff_old_mask = veff_mask;
      }
      vactive_mask &= v::int_cmplt(vkks, vnumneigh_i);
    }
  }
  fvec vfpair, vevdwl, vprefactor, vfwtmp, vfjwtmp;
  force_zeta_vector(c_outer, vc_idx_ij, vmask, vrij, vzeta, &vfpair, &vprefactor, EFLAG, &vevdwl, vmask_repulsive);
  vfxtmp = vfxtmp * vprefactor + vdx_ij * vfpair;
  vfytmp = vfytmp * vprefactor + vdy_ij * vfpair;
  vfztmp = vfztmp * vprefactor + vdz_ij * vfpair;
  vfjxtmp = vfjxtmp * vprefactor - vdx_ij * vfpair;
  vfjytmp = vfjytmp * vprefactor - vdy_ij * vfpair;
  vfjztmp = vfjztmp * vprefactor - vdz_ij * vfpair;

  if (EFLAG) {
    *vsevdwl = v::acc_mask_add(*vsevdwl, vmask, *vsevdwl, vevdwl);
    if (eatom) {
      v::store(fw, (v_0_5 * vevdwl));
    }
  }
  {
    while (cache_idx-- > 0) {
      fvec vfkx = vprefactor * vfkx_cache[cache_idx];
      fvec vfky = vprefactor * vfky_cache[cache_idx];
      fvec vfkz = vprefactor * vfkz_cache[cache_idx];
      ivec vks = vks_cache[cache_idx];
      bvec veff_mask = vmask_cache[cache_idx];
      v::store(fx, vfkx);
      v::store(fy, vfky);
      v::store(fz, vfkz);
      v::int_store(ts, vks);
      for (int t = 0; t < v::VL; t++) {
        if (v::mask_test_at(veff_mask, t)) {
          int t_ = ts[t] / (4 * sizeof(typename v::fscal));
          f[t_].x += fx[t];
          f[t_].y += fy[t];
          f[t_].z += fz[t];
        }
      }
    }
    ivec vkks = vkks_final_cache;
    bvec vactive_mask = vmask_final_cache;
    bvec veff_old_mask(0);
    ivec vks, vw_k;
    fvec vx_k, vy_k, vz_k, vcutsq;
    while (! v::mask_testz(vactive_mask)) {
      bvec vnew_mask = vactive_mask & ~ veff_old_mask;
      vks = v::int_mullo(v_i4floats, v_i_NEIGHMASK &
          v::int_gather<4>(vks, vactive_mask, vkks + vcnumneigh_i, firstneigh));
      v::gather_x(vks, vnew_mask, x, &vx_k, &vy_k, &vz_k, &vw_k);
      fvec vdx_ik = vx_k - vx_i;
      fvec vdy_ik = vy_k - vy_i;
      fvec vdz_ik = vz_k - vz_i;
      fvec vrsq = vdx_ik * vdx_ik + vdy_ik *  vdy_ik + vdz_ik * vdz_ik;
      ivec vc_idx = v::int_mullo(v_i4floats, vw_k) + v::int_mullo(v_i_ntypes, vc_idx_ij);
      vcutsq = v::gather<4>(vcutsq, vnew_mask, vc_idx, c_inner);
      bvec vcutoff_mask = v::cmplt(vrsq, vcutsq);
      bvec vsame_mask = v::int_cmpneq(vjs, vks);
      bvec veff_mask = vcutoff_mask & vsame_mask & vactive_mask;
      if (v::mask_testz(~(veff_mask | ~vactive_mask))) {
        fvec vfix, vfiy, vfiz;
        fvec vfjx, vfjy, vfjz;
        fvec vfkx, vfky, vfkz;

        attractive_vector<false>(c_inner,vc_idx,veff_mask,vprefactor,
            vrij,vrsq,vdx_ij,vdy_ij,vdz_ij,vdx_ik,vdy_ik,vdz_ik,
            &vfix,&vfiy,&vfiz,
            &vfjx,&vfjy,&vfjz,
            &vfkx,&vfky,&vfkz,
            0);
        vfxtmp = v::mask_add(vfxtmp, veff_mask, vfxtmp, vfix);
        vfytmp = v::mask_add(vfytmp, veff_mask, vfytmp, vfiy);
        vfztmp = v::mask_add(vfztmp, veff_mask, vfztmp, vfiz);
        vfjxtmp = v::mask_add(vfjxtmp, veff_mask, vfjxtmp, vfjx);
        vfjytmp = v::mask_add(vfjytmp, veff_mask, vfjytmp, vfjy);
        vfjztmp = v::mask_add(vfjztmp, veff_mask, vfjztmp, vfjz);
        v::store(fx, vfkx);
        v::store(fy, vfky);
        v::store(fz, vfkz);
        v::int_store(ts, vks);
        for (int t = 0; t < v::VL; t++) {
          if (v::mask_test_at(veff_mask, t)) {
            int t_ = ts[t] / (4 * sizeof(typename v::fscal));
            f[t_].x += fx[t];
            f[t_].y += fy[t];
            f[t_].z += fz[t];
          }
        }
        vkks = vkks + v_i1;
        veff_old_mask = bvec(0);
      } else {
        vkks = v::int_mask_add(vkks, ~veff_mask, vkks, v_i1);
        veff_old_mask = veff_mask;
      }
      vactive_mask &= v::int_cmplt(vkks, vnumneigh_i);
    } // while (vactive_mask != 0)
  } // section k
  // We can not make any assumptions regarding conflicts.
  // So we sequentialize this.
  // TDO: Once AVX-512 is around check out VPCONFLICT
  v::store(fx, vfjxtmp);
  v::store(fy, vfjytmp);
  v::store(fz, vfjztmp);
  for (int t = 0; t < compress_idx; t++) {
    int t_ = js[t];
    f[t_].x += fx[t];
    f[t_].y += fy[t];
    f[t_].z += fz[t];
    if (EFLAG && eatom) {
      f[t_].w += fw[t];
    }
  }
  v::store(fx, vfxtmp);
  v::store(fy, vfytmp);
  v::store(fz, vfztmp);
  for (int t = 0; t < compress_idx; t++) {
    int t_ = is[t];
    f[t_].x += fx[t];
    f[t_].y += fy[t];
    f[t_].z += fz[t];
    if (EFLAG && eatom) {
      f[t_].w += fw[t];
    }
  }
}

// Specialized kernel step for fixed i, means that we don't have to use the
//  convoluted iteration scheme above, as the loop variables are uniform.
template<class flt_t, class acc_t, lmp_intel::CalculationMode mic, bool pack_i>
template<int EFLAG>
void IntelKernelTersoff<flt_t,acc_t,mic, pack_i>::kernel_step_const_i(
    int eatom, int vflag,
    const int * _noalias const numneigh, const int * _noalias const cnumneigh,
    const int * _noalias const firstneigh, int ntypes,
    typename IntelBuffers<flt_t,acc_t>::atom_t * _noalias const x,
    const typename PairTersoffIntel::ForceConst<flt_t>::c_inner_t * _noalias const c_inner,
    const typename PairTersoffIntel::ForceConst<flt_t>::c_outer_t * _noalias const c_outer,
    typename IntelBuffers<flt_t,acc_t>::vec3_acc_t * _noalias const f,
    avec *vsevdwl,
    int compress_idx,
    int ii,
    int i,
    iarr js,
    bvec vmask_repulsive
) {
  typedef typename v::fvec fvec;
  typedef typename v::ivec ivec;
  typedef typename v::bvec bvec;
  typedef typename v::farr farr;
  typedef typename v::iarr iarr;
  typedef typename v::avec avec;
  typedef typename v::aarr aarr;

  ivec v_i4floats((int) (4 * sizeof(typename v::fscal)));
  ivec v_i1(1), v_i0(0), v_i_ntypes(ntypes), v_i_NEIGHMASK(NEIGHMASK);
  fvec v_0_5(0.5);

  int cache_idx = 0;
  fvec vfkx_cache[N_CACHE];
  fvec vfky_cache[N_CACHE];
  fvec vfkz_cache[N_CACHE];
  int k_cache[N_CACHE];
  bvec vmask_cache[N_CACHE];
  int kk_final_cache;

  aarr fx, fy, fz, fw;
  iarr ts;

  bvec vmask = v::mask_enable_lower(compress_idx);
  fvec vx_i(x[i].x), vy_i(x[i].y), vz_i(x[i].z);
  int w_i = x[i].w;

  ivec vjs = v::int_mullo(v_i4floats, v::int_load_vl(js));
  fvec vx_j = v::zero(), vy_j = v::zero(), vz_j = v::zero();
  ivec vw_j = v_i0;
  v::gather_x(vjs, vmask, x, &vx_j, &vy_j, &vz_j, &vw_j);

  fvec vdx_ij = vx_j - vx_i, vdy_ij = vy_j - vy_i, vdz_ij = vz_j - vz_i;
  fvec vrijsq = vdx_ij * vdx_ij + vdy_ij *  vdy_ij + vdz_ij * vdz_ij;
  fvec vrij = sqrt(vrijsq);

  int cnumneigh_i = cnumneigh[ii];
  int numneigh_i = numneigh[i];
  ivec vc_idx_j = v::int_mullo(v_i4floats, vw_j);
  ivec vc_idx_j_ntypes = v::int_mullo(v_i_ntypes, vc_idx_j);

  avec vzeta = v::acc_zero();
  avec vfxtmp = v::acc_zero(), vfytmp = v::acc_zero(), vfztmp = v::acc_zero();
  avec vfjxtmp = v::acc_zero(), vfjytmp = v::acc_zero(), vfjztmp = v::acc_zero();

  // Same structure as kernel_step, just simpler as the loops all iterate over
  //  the same k
  int kk = 0;
  for (; kk < numneigh_i && cache_idx < N_CACHE; kk++) {
    int k = firstneigh[kk + cnumneigh_i] & NEIGHMASK;
    fvec vx_k(x[k].x);
    fvec vy_k(x[k].y);
    fvec vz_k(x[k].z);
    int w_k = x[k].w;
    fvec vdx_ik = vx_k - vx_i;
    fvec vdy_ik = vy_k - vy_i;
    fvec vdz_ik = vz_k - vz_i;
    fvec vrsq = vdx_ik * vdx_ik + vdy_ik * vdy_ik + vdz_ik * vdz_ik;
    fvec vcutsq = v::gather<4>(v::zero(), vmask, vc_idx_j_ntypes, &c_inner[ntypes * ntypes * w_i + w_k]);
    bvec vcutoff_mask = v::cmplt(vrsq, vcutsq);
    bvec vsame_mask = v::int_cmpneq(vjs, ivec(static_cast<int>(4 * sizeof(typename v::fscal) * k)));
    bvec veff_mask = vcutoff_mask & vsame_mask & vmask;
    if (! v::mask_testz(veff_mask)) {
      fvec vzeta_contrib;
      fvec vfix, vfiy, vfiz;
      fvec vfjx, vfjy, vfjz;
      fvec vfkx, vfky, vfkz;

      attractive_vector<true>(&c_inner[ntypes * ntypes * w_i + w_k],vc_idx_j_ntypes,veff_mask,fvec(1.),
          vrij,vrsq,vdx_ij,vdy_ij,vdz_ij,vdx_ik,vdy_ik,vdz_ik,
          &vfix,&vfiy,&vfiz,
          &vfjx,&vfjy,&vfjz,
          &vfkx,&vfky,&vfkz,
          &vzeta_contrib);
      vfxtmp  = v::acc_mask_add(vfxtmp, veff_mask, vfxtmp, vfix);
      vfytmp  = v::acc_mask_add(vfytmp, veff_mask, vfytmp, vfiy);
      vfztmp  = v::acc_mask_add(vfztmp, veff_mask, vfztmp, vfiz);
      vfjxtmp = v::acc_mask_add(vfjxtmp, veff_mask, vfjxtmp, vfjx);
      vfjytmp = v::acc_mask_add(vfjytmp, veff_mask, vfjytmp, vfjy);
      vfjztmp = v::acc_mask_add(vfjztmp, veff_mask, vfjztmp, vfjz);

      vfkx_cache[cache_idx] = v::mask_add(v::zero(), veff_mask, vfkx, v::zero());
      vfky_cache[cache_idx] = v::mask_add(v::zero(), veff_mask, vfky, v::zero());
      vfkz_cache[cache_idx] = v::mask_add(v::zero(), veff_mask, vfkz, v::zero());
      vmask_cache[cache_idx] = veff_mask;
      k_cache[cache_idx] = k;
      cache_idx += 1;

      vzeta = v::acc_mask_add(vzeta, veff_mask, vzeta, vzeta_contrib);
    }
  }
  kk_final_cache = kk;
  for (; kk < numneigh_i; kk++) {
    int k = firstneigh[kk + cnumneigh_i] & NEIGHMASK;
    fvec vx_k(x[k].x);
    fvec vy_k(x[k].y);
    fvec vz_k(x[k].z);
    int w_k = x[k].w;
    fvec vdx_ik = vx_k - vx_i;
    fvec vdy_ik = vy_k - vy_i;
    fvec vdz_ik = vz_k - vz_i;
    fvec vrsq = vdx_ik * vdx_ik + vdy_ik * vdy_ik + vdz_ik * vdz_ik;
    fvec vcutsq = v::gather<4>(v::zero(), vmask, vc_idx_j_ntypes, &c_inner[ntypes * ntypes * w_i + w_k]);
    bvec vcutoff_mask = v::cmplt(vrsq, vcutsq);
    bvec vsame_mask = v::int_cmpneq(vjs, ivec(static_cast<int>(4 * sizeof(typename v::fscal) * k)));
    bvec veff_mask = vcutoff_mask & vsame_mask & vmask;
    if (! v::mask_testz(veff_mask)) {
      fvec vzeta_contrib = zeta_vector(&c_inner[ntypes * ntypes * w_i + w_k], vc_idx_j_ntypes, veff_mask, vrij, vrsq,
          vdx_ij,vdy_ij,vdz_ij,vdx_ik,vdy_ik,vdz_ik);
      vzeta = v::acc_mask_add(vzeta, veff_mask, vzeta, vzeta_contrib);
    }
  }
  fvec vfpair, vevdwl, vprefactor, vfwtmp;
  force_zeta_vector(&c_outer[ntypes * w_i], vc_idx_j, vmask, vrij, vzeta, &vfpair, &vprefactor, EFLAG, &vevdwl, vmask_repulsive);
  avec vaprefactor(vprefactor);
  vfxtmp  = vfxtmp  * vaprefactor + avec(vdx_ij * vfpair);
  vfytmp  = vfytmp  * vaprefactor + avec(vdy_ij * vfpair);
  vfztmp  = vfztmp  * vaprefactor + avec(vdz_ij * vfpair);
  vfjxtmp = vfjxtmp * vaprefactor - avec(vdx_ij * vfpair);
  vfjytmp = vfjytmp * vaprefactor - avec(vdy_ij * vfpair);
  vfjztmp = vfjztmp * vaprefactor - avec(vdz_ij * vfpair);

  if (EFLAG) {
    *vsevdwl = v::acc_mask_add(*vsevdwl, vmask, *vsevdwl, vevdwl);
    if (eatom) {
      vfwtmp = v_0_5 * vevdwl;
      v::store(fw, vfwtmp);
    }
  }
  while (cache_idx-- > 0) {
    fvec vfkx = vprefactor * vfkx_cache[cache_idx];
    fvec vfky = vprefactor * vfky_cache[cache_idx];
    fvec vfkz = vprefactor * vfkz_cache[cache_idx];
    int k = k_cache[cache_idx];
    bvec veff_mask = vmask_cache[cache_idx];
    f[k].x += v::reduce_add(v::mask_add(v::zero(), veff_mask, vfkx, v::zero()));
    f[k].y += v::reduce_add(v::mask_add(v::zero(), veff_mask, vfky, v::zero()));
    f[k].z += v::reduce_add(v::mask_add(v::zero(), veff_mask, vfkz, v::zero()));
  }
  for (int kk = kk_final_cache; kk < numneigh_i; kk++) {
    int k = firstneigh[kk + cnumneigh_i] & NEIGHMASK;
    fvec vx_k(x[k].x);
    fvec vy_k(x[k].y);
    fvec vz_k(x[k].z);
    int w_k = x[k].w;
    fvec vdx_ik = vx_k - vx_i;
    fvec vdy_ik = vy_k - vy_i;
    fvec vdz_ik = vz_k - vz_i;
    fvec vrsq = vdx_ik * vdx_ik + vdy_ik * vdy_ik + vdz_ik * vdz_ik;
    fvec vcutsq = v::gather<4>(v::zero(), vmask, vc_idx_j_ntypes, &c_inner[ntypes * ntypes * w_i + w_k].cutsq);
    bvec vcutoff_mask = v::cmplt(vrsq, vcutsq);
    bvec vsame_mask = v::int_cmpneq(vjs, ivec(static_cast<int>(4 * sizeof(typename v::fscal) * k)));
    bvec veff_mask = vcutoff_mask & vsame_mask & vmask;
    if (! v::mask_testz(veff_mask)) {
       fvec vfix, vfiy, vfiz;
       fvec vfjx, vfjy, vfjz;
       fvec vfkx, vfky, vfkz;

       attractive_vector<false>(&c_inner[ntypes * ntypes * w_i + w_k],vc_idx_j_ntypes,veff_mask,vprefactor,
           vrij,vrsq,vdx_ij,vdy_ij,vdz_ij,vdx_ik,vdy_ik,vdz_ik,
           &vfix,&vfiy,&vfiz,
           &vfjx,&vfjy,&vfjz,
           &vfkx,&vfky,&vfkz,
           0);
       vfxtmp  = v::acc_mask_add(vfxtmp, veff_mask, vfxtmp, vfix);
       vfytmp  = v::acc_mask_add(vfytmp, veff_mask, vfytmp, vfiy);
       vfztmp  = v::acc_mask_add(vfztmp, veff_mask, vfztmp, vfiz);
       vfjxtmp = v::acc_mask_add(vfjxtmp, veff_mask, vfjxtmp, vfjx);
       vfjytmp = v::acc_mask_add(vfjytmp, veff_mask, vfjytmp, vfjy);
       vfjztmp = v::acc_mask_add(vfjztmp, veff_mask, vfjztmp, vfjz);
       f[k].x += v::reduce_add(v::mask_add(v::zero(), veff_mask, vfkx, v::zero()));
       f[k].y += v::reduce_add(v::mask_add(v::zero(), veff_mask, vfky, v::zero()));
       f[k].z += v::reduce_add(v::mask_add(v::zero(), veff_mask, vfkz, v::zero()));
    }
  }
  // TDO: This could be a scatter
  v::acc_store(fx, vfjxtmp);
  v::acc_store(fy, vfjytmp);
  v::acc_store(fz, vfjztmp);
  for (int t = 0; t < compress_idx; t++) {
    int t_ = js[t];
    f[t_].x += fx[t];
    f[t_].y += fy[t];
    f[t_].z += fz[t];
    if (EFLAG && eatom) {
      f[t_].w += fw[t];
    }
  }
  f[i].x += v::acc_reduce_add(v::acc_mask_add(v::acc_zero(), vmask, vfxtmp, v::zero()));
  f[i].y += v::acc_reduce_add(v::acc_mask_add(v::acc_zero(), vmask, vfytmp, v::zero()));
  f[i].z += v::acc_reduce_add(v::acc_mask_add(v::acc_zero(), vmask, vfztmp, v::zero()));
  if (EFLAG && eatom) {
    f[i].z += v::acc_reduce_add(v::acc_mask_add(v::acc_zero(), vmask, vfwtmp, v::zero()));
  }
}

template<class flt_t, class acc_t, lmp_intel::CalculationMode mic, bool pack_i>
template<bool EFLAG>
void IntelKernelTersoff<flt_t,acc_t,mic, pack_i>::kernel(
    int iito, int iifrom, int eatom, int vflag,
    const int * _noalias const numneigh,
    const int * _noalias const numneighhalf,
    const int * _noalias const cnumneigh,
    const int * _noalias const firstneigh, int ntypes,
    typename IntelBuffers<flt_t,acc_t>::atom_t * _noalias const x,
    const int * _noalias const ilist,
    const c_inner_t * _noalias const c_inner,
    const c_outer_t * _noalias const c_outer,
    typename IntelBuffers<flt_t,acc_t>::vec3_acc_t * _noalias const f,
    acc_t *evdwl
) {
  int compress_idx = 0;
  int ii, jj;
  iarr is, js, vc;
  avec vsevdwl = v::acc_zero();
  ivec v_i4floats(static_cast<int>(sizeof(typename v::fscal) * 4));
  ivec vj, v_NEIGHMASK(NEIGHMASK);
  bvec vmask_repulsive(0);
  iarr repulsive_flag = {0};
  // If you want to get the very most out of this, please uncomment.
  // Consider getting a coffee or doing something else.
  // Also good for heating.
  //#pragma forceinline recursive
  for (ii = iifrom; ii < iito; ii++) {
    // Right now this loop is scalar, to allow for the compiler to do
    //  its prefetching magic.
    int i = ilist[ii];
    int w_i = x[i].w;
    flt_t x_i = x[i].x;
    flt_t y_i = x[i].y;
    flt_t z_i = x[i].z;
    int jlist_off_i = cnumneigh[ii];
    int jnum = numneigh[i];
    for (jj = 0; jj < jnum; jj++) {
      int j = firstneigh[jlist_off_i + jj] & NEIGHMASK;
      int w_j = x[j].w;
      flt_t dx_ij = x[j].x - x_i;
      flt_t dy_ij = x[j].y - y_i;
      flt_t dz_ij = x[j].z - z_i;
      flt_t rsq = dx_ij*dx_ij + dy_ij*dy_ij + dz_ij*dz_ij;
      flt_t cutsq = c_outer[w_i * ntypes + w_j].cutsq;
      if (rsq < cutsq) {
        is[compress_idx] = ii;
        js[compress_idx] = j;
        vc[compress_idx] = jlist_off_i;
        if (jj < numneighhalf[ii])
          repulsive_flag[compress_idx] = 1;
        compress_idx += 1;
      }
      if (pack_i) {
        if (compress_idx == v::VL) {
          vmask_repulsive = v::int_cmpneq(v::int_load_vl(repulsive_flag), ivec(0));
          kernel_step<EFLAG>(
              eatom, vflag,
              numneigh, vc, firstneigh, ntypes,
              x, c_inner, c_outer, f,
              &vsevdwl, compress_idx,
              is, js, vmask_repulsive
          );
          compress_idx = 0;
          v::int_clear_arr(repulsive_flag);
        }
      } else {
        if (compress_idx == v::VL || (compress_idx > 0 && jj == jnum-1)) {
          vmask_repulsive = v::int_cmpneq(v::int_load_vl(repulsive_flag), ivec(0));
          kernel_step_const_i<EFLAG>(
              eatom, vflag,
              numneigh, cnumneigh, firstneigh, ntypes,
              x, c_inner, c_outer, f,
              &vsevdwl, compress_idx,
              ii, i, js, vmask_repulsive
          );
          compress_idx = 0;
          v::int_clear_arr(repulsive_flag);
        }
      }
    }
  }
  if (compress_idx > 0) {
        vmask_repulsive = v::int_cmpneq(v::int_load_vl(repulsive_flag), ivec(0));
        IntelKernelTersoff::kernel_step<EFLAG>(
            eatom, vflag,
            numneigh, vc, firstneigh, ntypes,
            x, c_inner, c_outer, f,
            &vsevdwl, compress_idx,
            is, js, vmask_repulsive
        );
  }
  if (EFLAG) {
    *evdwl += v::acc_reduce_add(vsevdwl);
  }
}


template<class flt_t, class acc_t, lmp_intel::CalculationMode mic, bool pack_i>
IntelKernelTersoff<flt_t,acc_t,mic,pack_i>::fvec IntelKernelTersoff<flt_t, acc_t, mic, pack_i>::zeta_vector(
    const c_inner_t * param,
    ivec xjw, bvec mask,
    fvec vrij, fvec rsq2,
    fvec vdijx, fvec vdijy, fvec vdijz,
    fvec dikx, fvec diky, fvec dikz
) {
  fvec v_1_0(1.0);
  fvec v_0_5(0.5);
  fvec vph = v::zero(), vpc2 = v::zero(), vpd2 = v::zero(), vpgamma = v::zero(), vplam3 = v::zero(), vppowermint = v::zero(), vpbigr = v::zero(), vpbigd = v::zero();
  // TDO: Specialize on number of species
  v::gather_8(xjw, mask, &param[0].lam3, &vplam3, &vppowermint, &vpbigr, &vpbigd, &vpc2, &vpd2, &vph, &vpgamma);
  fvec vrik = sqrt(rsq2);
  fvec vcostheta = (vdijx * dikx + vdijy * diky + vdijz * dikz) * v::recip(vrij * vrik);
  fvec vhcth = vph - vcostheta;
  fvec vgijk_a = vhcth * vhcth;
  fvec vgijk = vpgamma * (v_1_0 + vpc2 * vgijk_a * v::recip(vpd2 * (vpd2 + vgijk_a)));
  fvec varg1 = vplam3 * (vrij - vrik);
  fvec varg3 = varg1 * varg1 * varg1;
  bvec mask_ex = v::cmpeq(vppowermint, fvec(3.));
  fvec varg  = v::blend(mask_ex, varg1, varg3);
  fvec vex_delr = v::min(fvec(1.e30), exp(varg));
  bvec vmask_need_sine = v::cmpnle(vrik, vpbigr - vpbigd) & mask;
  fvec vfc = v_1_0;
  // Its kind of important to check the mask.
  // Some simulations never/rarely invoke this branch.
  if (! v::mask_testz(vmask_need_sine)) {
    vfc = v::blend(vmask_need_sine, vfc,
        v_0_5 * (v_1_0 - sin(fvec(MY_PI2) * (vrik - vpbigr) * v::recip(vpbigd))));
  }
  return vgijk * vex_delr * vfc;
}

template<class flt_t, class acc_t, lmp_intel::CalculationMode mic, bool pack_i>
void IntelKernelTersoff<flt_t, acc_t, mic, pack_i>::force_zeta_vector(
    const c_outer_t * param,
    ivec xjw,
    bvec mask,
    fvec vrij, fvec vzeta_ij,
    fvec *vfpair, fvec *vprefactor, int EVDWL, fvec *vevdwl,
    bvec vmask_repulsive
) {
  fvec v_0_0(0.0);
  fvec v_0_5(0.5);
  fvec v_m0_5(-0.5);
  fvec v_1_0(1.0);
  fvec v_m1_0(-1.0);
  fvec v_2_0(2.0);
  fvec vpbigr = v::zero(), vpbigd = v::zero(), vplam1 = v::zero(), vpbiga = v::zero(), vplam2 = v::zero(), vpbeta = v::zero(), vpbigb = v::zero(), vppowern = v::zero();
  v::gather_8(xjw, mask, &param[0].bigr, &vpbigr, &vpbigd, &vplam1, &vpbiga, &vplam2, &vpbeta, &vpbigb, &vppowern);
  fvec vfccos;

  // This is pretty much a literal translation.
  bvec vmask_need_sine = v::cmpnle(vrij, vpbigr - vpbigd) & mask;
  fvec vfc = v_1_0;
  fvec vfc_d = v_0_0;
  if (! v::mask_testz(vmask_need_sine)) {
    fvec vtmp = fvec(MY_PI2) * v::recip(vpbigd);
    vfc = v::blend(vmask_need_sine, vfc,
        v_0_5 * (v_1_0 - v::sincos(&vfccos, vtmp * (vrij - vpbigr))));
    vfc_d = v::blend(vmask_need_sine, vfc_d, v_m0_5 * vtmp * vfccos);
  }
  fvec vpminus_lam2 =  - vplam2;

  fvec vpminus_bigb = -vpbigb;
  fvec vexp = exp(vpminus_lam2 * vrij);
  fvec vfa = vpminus_bigb * vexp * vfc;
  fvec vfa_d = vpminus_lam2 * vfa + vpminus_bigb * vexp * vfc_d;

  fvec vpc1 = v::zero(), vpc2 = v::zero(), vpc3 = v::zero(), vpc4 = v::zero();
  v::gather_4(xjw, mask, &param[0].c1, &vpc1, &vpc2, &vpc3, &vpc4);
  fvec vpminus_powern = - vppowern;
  fvec vbij(0.), vbij_d(0.);
  fvec vtmp = vpbeta * vzeta_ij;
  bvec vmc1 = v::cmple(vpc1, vtmp) & mask;
  if (! v::mask_testz(vmc1)) {
    vbij = v::invsqrt(vtmp);
    vbij_d = vpbeta * v_m0_5 * vbij * v::recip(vtmp);
  }
  bvec vmc2 = v::cmple(vpc2, vtmp) & ~ vmc1 & mask;
  if (! v::mask_testz(vmc2)) {
    fvec vpowminus_powern = pow(vtmp, vpminus_powern);
    fvec vinvsqrt = v::invsqrt(vtmp);
    fvec vrcp2powern = v::recip(v_2_0 * vppowern);
    fvec va = (v_1_0 - vpowminus_powern * vrcp2powern) * vinvsqrt;
    fvec va_d = vpbeta * v_m0_5 * vinvsqrt * v::recip(vtmp) *
            (v_1_0 + v_m0_5 * vpowminus_powern * (v_1_0 + vrcp2powern));
    vbij = v::blend(vmc2, vbij, va);
    vbij_d = v::blend(vmc2, vbij_d, va_d);
  }
  bvec vmc3 = v::cmplt(vtmp, vpc4) & ~vmc2 & ~vmc1 & mask;
  if (! v::mask_testz(vmc3)) {
    vbij = v::blend(vmc3, vbij, v_1_0);
    vbij_d = v::blend(vmc3, vbij_d, v_0_0);
  }
  bvec vmc4 = v::cmple(vtmp, vpc3) & ~vmc3 & ~vmc2 & ~ vmc1 & mask;
  if (! v::mask_testz(vmc4)) {
    fvec vpowm1 = pow(vtmp, vppowern - v_1_0);
    fvec vrcp2powern = v::recip(v_2_0 * vppowern);
    fvec va = v_1_0 - vtmp * vrcp2powern * vpowm1;
    fvec va_d = v_m0_5 * vpbeta * vpowm1;
    vbij = v::blend(vmc4, vbij, va);
    vbij_d = v::blend(vmc4, vbij_d, va_d);
  }
  bvec vmc5 = mask & ~vmc1 & ~vmc2 & ~vmc3 & ~vmc4;
  if (! v::mask_testz(vmc5)) {
    fvec vtmp_n = pow(vtmp, vppowern);
    fvec vpow2 = pow(v_1_0 + vtmp_n, v_m1_0 - v::recip(v_2_0 * vppowern));
    fvec va = (v_1_0 + vtmp_n) * vpow2;
    fvec va_d = v_m0_5 * vpow2 * vtmp_n * v::recip(vzeta_ij);
    vbij = v::blend(vmc5, vbij, va);
    vbij_d = v::blend(vmc5, vbij_d, va_d);
  }
  fvec vtmp_exp = exp(-vplam1 * vrij);
  fvec vrep_fforce = vpbiga * vtmp_exp * (vfc_d - vfc * vplam1);
  fvec vfz_fforce = v_0_5 * vbij * vfa_d;

  *vfpair = v::mask_add(vfz_fforce, vmask_repulsive, vfz_fforce, vrep_fforce) * v::recip(vrij);
  *vprefactor = v_m0_5 * vfa * vbij_d;
  if (EVDWL) {
    fvec vrep_eng = vfc * vpbiga * vtmp_exp;
    fvec vfz_eng = v_0_5 * vfa * vbij;
    *vevdwl = v::mask_add(vfz_eng, vmask_repulsive, vfz_eng, vrep_eng);
  }
}

template<class flt_t, class acc_t, lmp_intel::CalculationMode mic, bool pack_i>
template<bool ZETA>
void IntelKernelTersoff<flt_t,acc_t,mic, pack_i>::attractive_vector(
    const c_inner_t * param,
    ivec xjw,
    bvec mask,
    fvec vprefactor,
    fvec vrij, fvec rsq2,
    fvec vdijx, fvec vdijy, fvec vdijz,
    fvec dikx, fvec diky, fvec dikz,
    fvec *fix, fvec *fiy, fvec *fiz,
    fvec *fjx, fvec *fjy, fvec *fjz,
    fvec *fkx, fvec *fky, fvec *fkz,
    fvec *zeta
) {
  fvec v_1_0 = fvec(1.0);

  fvec vph = v::zero(), vpc2 = v::zero(), vpd2 = fvec(1.0), vpgamma = v::zero(), vplam3 = v::zero(), vppowermint = v::zero(), vpbigr = v::zero(), vpbigd = fvec(1.0);
  v::gather_8(xjw, mask, &param[0].lam3, &vplam3, &vppowermint, &vpbigr, &vpbigd, &vpc2, &vpd2, &vph, &vpgamma);
  fvec vrijinv = v::recip(vrij);
  fvec vrij_hatx = vrijinv * vdijx;
  fvec vrij_haty = vrijinv * vdijy;
  fvec vrij_hatz = vrijinv * vdijz;
  fvec rikinv = v::invsqrt(rsq2);
  fvec rik_hatx = rikinv * dikx;
  fvec rik_haty = rikinv * diky;
  fvec rik_hatz = rikinv * dikz;

  fvec vrik = sqrt(rsq2);
  fvec vcostheta = (vdijx * dikx + vdijy * diky + vdijz * dikz) * v::recip(vrij * vrik);
  fvec vhcth = vph - vcostheta;
  fvec vdenominator = v::recip(vpd2 + vhcth * vhcth);
  fvec vgijk = vpgamma * (v_1_0 + vpc2 * v::recip(vpd2) - vpc2 * vdenominator);
  fvec vnumerator = fvec(-2.) * vpc2 * vhcth;
  fvec vgijk_d = vpgamma * vnumerator * vdenominator * vdenominator;
  fvec varg1 = vplam3 * (vrij - vrik);
  fvec varg3 = varg1 * varg1 * varg1;
  bvec mask_ex = v::cmpeq(vppowermint, fvec(3.));
  fvec varg  = v::blend(mask_ex, varg1, varg3);
  fvec vex_delr = min(fvec(1.e30), exp(varg));
  fvec vex_delr_d_factor = v::blend(mask_ex, v_1_0, fvec(3.0) * varg1 * varg1);
  fvec vex_delr_d = vplam3 * vex_delr_d_factor * vex_delr;
  bvec vmask_need_sine = v::cmpnle(vrik, vpbigr - vpbigd) & mask;
  fvec vfccos;
  fvec vfc = v_1_0;
  fvec vfc_d = v::zero();
  if (! v::mask_testz(vmask_need_sine)) {
    fvec vtmp = fvec(MY_PI2) * v::recip(vpbigd);
    vfc = v::blend(vmask_need_sine, vfc,
        fvec(0.5) * (v_1_0 - v::sincos(&vfccos, vtmp * (vrik - vpbigr))));
    vfc_d = v::blend(vmask_need_sine, vfc_d, fvec(-0.5) * vtmp * vfccos);
  }

  fvec vzeta_d_fc = vfc_d * vgijk * vex_delr;
  fvec vzeta_d_gijk = vfc * vgijk_d * vex_delr;
  fvec vzeta_d_ex_delr = vfc * vgijk * vex_delr_d;
  if (ZETA) *zeta = vfc * vgijk * vex_delr;

  fvec vminus_costheta = - vcostheta;
  fvec vdcosdrjx = vrijinv * fmadd(vminus_costheta, vrij_hatx, rik_hatx);
  fvec vdcosdrjy = vrijinv * fmadd(vminus_costheta, vrij_haty, rik_haty);
  fvec vdcosdrjz = vrijinv * fmadd(vminus_costheta, vrij_hatz, rik_hatz);
  fvec vdcosdrkx = rikinv * fmadd(vminus_costheta, rik_hatx, vrij_hatx);
  fvec vdcosdrky = rikinv * fmadd(vminus_costheta, rik_haty, vrij_haty);
  fvec vdcosdrkz = rikinv * fmadd(vminus_costheta, rik_hatz, vrij_hatz);
  fvec vdcosdrix = -(vdcosdrjx + vdcosdrkx);
  fvec vdcosdriy = -(vdcosdrjy + vdcosdrky);
  fvec vdcosdriz = -(vdcosdrjz + vdcosdrkz);

  *fix = vprefactor * (vzeta_d_gijk * vdcosdrix + vzeta_d_ex_delr * (rik_hatx - vrij_hatx) - vzeta_d_fc * rik_hatx);
  *fiy = vprefactor * (vzeta_d_gijk * vdcosdriy + vzeta_d_ex_delr * (rik_haty - vrij_haty) - vzeta_d_fc * rik_haty);
  *fiz = vprefactor * (vzeta_d_gijk * vdcosdriz + vzeta_d_ex_delr * (rik_hatz - vrij_hatz) - vzeta_d_fc * rik_hatz);
  *fjx = vprefactor * (vzeta_d_gijk * vdcosdrjx + vzeta_d_ex_delr * vrij_hatx);
  *fjy = vprefactor * (vzeta_d_gijk * vdcosdrjy + vzeta_d_ex_delr * vrij_haty);
  *fjz = vprefactor * (vzeta_d_gijk * vdcosdrjz + vzeta_d_ex_delr * vrij_hatz);
  *fkx = vprefactor * ((vzeta_d_fc - vzeta_d_ex_delr) * rik_hatx + vzeta_d_gijk * vdcosdrkx);
  *fky = vprefactor * ((vzeta_d_fc - vzeta_d_ex_delr) * rik_haty + vzeta_d_gijk * vdcosdrky);
  *fkz = vprefactor * ((vzeta_d_fc - vzeta_d_ex_delr) * rik_hatz + vzeta_d_gijk * vdcosdrkz);
}


#ifdef _LMP_INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif

#endif
