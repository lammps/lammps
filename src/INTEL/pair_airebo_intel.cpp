// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Markus Hohnerbach (RWTH)
------------------------------------------------------------------------- */

#include "pair_airebo_intel.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "suffix.h"

#include <cassert>
#include <cmath>
#include <cstring>

#include "intel_preprocess.h"
#include "intel_intrinsics_airebo.h"

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace LAMMPS_NS;
using namespace MathConst;

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

template<typename flt_t, typename acc_t>
struct LAMMPS_NS::PairAIREBOIntelParam {
  flt_t cutlj, cutljrebosq, cut3rebo;
  flt_t sigmin, sigcut;
  flt_t cutljsq[2][2];
  flt_t lj1[2][2], lj2[2][2], lj3[2][2], lj4[2][2];

  flt_t smin, Nmin, Nmax, NCmin, NCmax, thmin, thmax;
  flt_t rcmin[2][2], rcmax[2][2], rcmaxsq[2][2], rcmaxp[2][2];
  flt_t Q[2][2], alpha[2][2], A[2][2], rho[2][2], BIJc[2][2][3],
      Beta[2][2][3];
  flt_t rcLJmin[2][2], rcLJmax[2][2], rcLJmaxsq[2][2], bLJmin[2][2],
      bLJmax[2][2];
  flt_t epsilon[2][2], sigma[2][2], epsilonT[2][2];

  // spline coefficients

  flt_t gCdom[5], gC1[4][6], gC2[4][6], gHdom[4], gH[3][6];
  flt_t gDom[5+4];
  flt_t gVal[(4+4+3)*6];
  flt_t pCCdom[2][2], pCHdom[2][2], pCC[4][4][16], pCH[4][4][16];
  flt_t piCCdom[3][2], piCHdom[3][2], piHHdom[3][2];
  acc_t piCC[4][4][9][64], piCH[4][4][9][64], piHH[4][4][9][64];
  flt_t Tijdom[3][2];
  acc_t Tijc[4][4][9][64];

  // spline knot values

  flt_t PCCf[5][5], PCCdfdx[5][5], PCCdfdy[5][5], PCHf[5][5];
  flt_t PCHdfdx[5][5], PCHdfdy[5][5];
  flt_t piCCf[5][5][11], piCCdfdx[5][5][11];
  flt_t piCCdfdy[5][5][11], piCCdfdz[5][5][11];
  flt_t piCHf[5][5][11], piCHdfdx[5][5][11];
  flt_t piCHdfdy[5][5][11], piCHdfdz[5][5][11];
  flt_t piHHf[5][5][11], piHHdfdx[5][5][11];
  flt_t piHHdfdy[5][5][11], piHHdfdz[5][5][11];
  flt_t Tf[5][5][10], Tdfdx[5][5][10], Tdfdy[5][5][10], Tdfdz[5][5][10];
};

namespace {

struct NeighListLMPAIREBO {
  int * num; /* num_all */
  int * num_half; /* num_all */
  int * offset; /* num_all */
  int ** entries; /* num_all * num_neighs_per_atom */
};

struct NeighListAIREBO {
  int * num; /* num_all */
  int * num_half; /* num_all */
  int * offset; /* num_all */
  int * entries; /* num_all * num_neighs_per_atom */
};

template<typename flt_t>
struct AtomAIREBOT {
  flt_t x, y, z;
  int w;
};

template<typename acc_t>
struct ResultForceT {
  acc_t x, y, z, w;
};

template<typename flt_t, typename acc_t>
struct KernelArgsAIREBOT {
  int num_local;
  int num_all;
  int num_neighs_per_atom;
  int num_types;
  int frebo_from_atom, frebo_to_atom;
  int neigh_from_atom, neigh_to_atom;
  int rebuild_flag;
  flt_t skin;
  struct NeighListLMPAIREBO neigh_lmp;
  struct NeighListAIREBO neigh_rebo;
  PairAIREBOIntelParam<flt_t,acc_t> params;
  struct AtomAIREBOT<flt_t> * x; /* num_all */
  tagint * tag; /* num_all */
  flt_t * nC, * nH; /* num_all */
  int * map; /* num_types+1 */
  struct ResultForceT<acc_t> * result_f; /* num_all */
  acc_t result_eng;
};

template<typename flt_t, typename acc_t>
void aut_lennard_jones(KernelArgsAIREBOT<flt_t,acc_t> * ka, int morseflag);
template<typename flt_t, typename acc_t>
void aut_rebo_neigh(KernelArgsAIREBOT<flt_t,acc_t> * ka);
template<typename flt_t, typename acc_t>
void aut_frebo(KernelArgsAIREBOT<flt_t,acc_t> * ka, int torsion_flag);

}

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif

/* ---------------------------------------------------------------------- */

PairAIREBOIntel::PairAIREBOIntel(LAMMPS *lmp) : PairAIREBO(lmp)
{
  suffix_flag |= Suffix::INTEL;
  REBO_cnumneigh = nullptr;
  REBO_num_skin = nullptr;
  REBO_list_data = nullptr;
  fix = nullptr;
}

/* ---------------------------------------------------------------------- */

PairAIREBOIntel::~PairAIREBOIntel()
{
  memory->destroy(REBO_cnumneigh);
  memory->destroy(REBO_num_skin);
  memory->destroy(REBO_list_data);
}

/* ---------------------------------------------------------------------- */

void PairAIREBOIntel::init_style()
{
  PairAIREBO::init_style();
  neighbor->find_request(this)->intel = 1;

  if (utils::strmatch(force->pair_style,"^hybrid"))
    error->all(FLERR, "Cannot yet use airebo/intel with hybrid.");

  int ifix = modify->find_fix("package_intel");
  if (ifix < 0)
    error->all(FLERR,
               "The 'package intel' command is required for /intel styles");
  fix = static_cast<FixIntel *>(modify->fix[ifix]);

  fix->pair_init_check();
  #ifdef _LMP_INTEL_OFFLOAD
  _cop = fix->coprocessor_number();
  #endif

  if (fix->precision() == FixIntel::PREC_MODE_MIXED) {
    pack_force_const(fix->get_mixed_buffers());
    fix->get_mixed_buffers()->need_tag(1);
  } else if (fix->precision() == FixIntel::PREC_MODE_DOUBLE) {
    pack_force_const(fix->get_double_buffers());
    fix->get_double_buffers()->need_tag(1);
  } else {
    pack_force_const(fix->get_single_buffers());
    fix->get_single_buffers()->need_tag(1);
  }

  #ifdef _LMP_INTEL_OFFLOAD
  if (fix->offload_noghost())
    error->all(FLERR,"The 'ghost no' option cannot be used with airebo/intel.");
  #endif
}

/* ---------------------------------------------------------------------- */

template<typename T>
T * calloc_it(size_t size) {
  return static_cast<T*>(calloc(size, sizeof(T)));
}

void PairAIREBOIntel::compute(int eflag, int vflag)
{
  if (fix->precision()==FixIntel::PREC_MODE_MIXED)
    compute<float,double>(eflag, vflag, fix->get_mixed_buffers());
  else if (fix->precision()==FixIntel::PREC_MODE_DOUBLE)
    compute<double,double>(eflag, vflag, fix->get_double_buffers());
  else
    compute<float,float>(eflag, vflag, fix->get_single_buffers());

  fix->balance_stamp();
  vflag_fdotr = 0;
}

/* ---------------------------------------------------------------------- */

template<class flt_t, class acc_t>
PairAIREBOIntelParam<flt_t,acc_t> PairAIREBOIntel::get_param()
{
  PairAIREBOIntelParam<flt_t,acc_t> fc;

#define A(a)                                                           \
  for (int i = 0; i < sizeof(this->a)/sizeof(double); i++) {           \
    reinterpret_cast<flt_t*>(&fc.a)[i] =                               \
      reinterpret_cast<double*>(&this->a)[i];                          \
  }
#define A0(a)                                                           \
  for (int i = 0; i < sizeof(fc.a)/sizeof(flt_t); i++) {                \
    reinterpret_cast<flt_t*>(&fc.a)[i] =                                \
      reinterpret_cast<double*>(this->a[0])[i];                         \
  }
#define B(a)                                                            \
  for (int i = 0; i < sizeof(this->a)/sizeof(double); i++) {            \
    reinterpret_cast<acc_t*>(&fc.a)[i] =                                \
      reinterpret_cast<double*>(&this->a)[i];                           \
  }

  A(cutlj) A(cutljrebosq) A(cut3rebo) A(sigmin);
  A(sigcut) A0(cutljsq) A0(lj1) A0(lj2) A0(lj3);
  A0(lj4) A(smin) A(Nmin) A(Nmax) A(NCmin) A(NCmax) A(thmin) A(thmax);
  A(rcmin) A(rcmax) A(rcmaxsq) A(rcmaxp) A(Q) A(alpha) A(A) A(rho) A(BIJc);
  A(Beta) A(rcLJmin) A(rcLJmax) A(rcLJmaxsq) A(bLJmin) A(bLJmax) A(epsilon);
  A(sigma) A(epsilonT) A(gCdom) A(gC1) A(gC2) A(gHdom) A(gH) A(pCCdom);
  A(pCHdom) A(pCC) A(pCH) A(piCCdom) A(piCHdom) A(piHHdom) B(piCC);
  B(piCH) B(piHH) A(Tijdom) B(Tijc) A(PCCf) A(PCCdfdx) A(PCCdfdy) A(PCHf);
  A(PCHdfdx) A(PCHdfdy) A(piCCf) A(piCCdfdx) A(piCCdfdy) A(piCCdfdz);
  A(piCHf) A(piCHdfdx) A(piCHdfdy) A(piCHdfdz) A(piHHf) A(piHHdfdx);
  A(piHHdfdy) A(piHHdfdz) A(Tf) A(Tdfdx) A(Tdfdy) A(Tdfdz);

#undef A
#undef A0
#undef B
  for (int i = 0; i < 5; i++) fc.gDom[i] = fc.gCdom[i];
  for (int i = 0; i < 4; i++) fc.gDom[5+i] = fc.gHdom[i];
  for (int i = 0; i < 4; i++) for (int j = 0; j < 6; j++)
                                fc.gVal[6*i+j] = fc.gC1[i][j];
  for (int i = 0; i < 4; i++) for (int j = 0; j < 6; j++)
                                fc.gVal[4*6+6*i+j] = fc.gC2[i][j];
  for (int i = 0; i < 3; i++) for (int j = 0; j < 6; j++)
                                fc.gVal[8*6+6*i+j] = fc.gH[i][j];

  return fc;
}

/* ---------------------------------------------------------------------- */

template<class flt_t, class acc_t>
void PairAIREBOIntel::compute(
    int eflag, int vflag, IntelBuffers<flt_t,acc_t> * buffers
) {
  ev_init(eflag,vflag);
  if (vflag_atom)
    error->all(FLERR,"INTEL package does not support per-atom stress");
  if (vflag && !vflag_fdotr && force->newton_pair)
    error->all(FLERR,"INTEL package does not support pair_modify nofdotr "
               "with newton on");

  pvector[0] = pvector[1] = pvector[2] = 0.0;

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
    #pragma omp parallel if (packthreads > 1)
    #endif
    {
      int ifrom, ito, tid;
      IP_PRE_omp_range_id_align(ifrom, ito, tid, atom->nlocal + atom->nghost,
                                packthreads, sizeof(ATOM_T));
        buffers->thr_pack(ifrom,ito,ago);
    }
    fix->stop_watch(TIME_PACK);
  }

  if (atom->nmax > maxlocal) {
    #ifdef LMP_INTEL_OFFLOAD
    if (maxlocal > 0 && _cop >= 0) {
      int * const REBO_numneigh = this->REBO_numneigh;
      int * const REBO_num_skin = this->REBO_num_skin;
      int * const REBO_cnumneigh = this->REBO_cnumneigh;
      int * const REBO_list_data = this->REBO_list_data;
      double * const nC = this->nC;
      double * const nH = this->nH;
      #pragma offload_transfer target(mic:_cop) \
        nocopy(REBO_numneigh: alloc_if(0) free_if(1)) \
        nocopy(REBO_cnumneigh: alloc_if(0) free_if(1)) \
        nocopy(REBO_num_skin: alloc_if(0) free_if(1)) \
        nocopy(REBO_list_data: alloc_if(0) free_if(1)) \
        nocopy(nH: alloc_if(0) free_if(1)) \
        nocopy(nC: alloc_if(0) free_if(1))
    }
    #endif
    maxlocal = atom->nmax;
    memory->destroy(REBO_numneigh);
    memory->destroy(REBO_cnumneigh);
    memory->destroy(REBO_list_data);
    memory->sfree(REBO_firstneigh);
    memory->destroy(nC);
    memory->destroy(nH);
    memory->create(REBO_numneigh,maxlocal,"AIREBO:numneigh");
    memory->create(REBO_cnumneigh,maxlocal,"AIREBO:cnumneigh");
    memory->create(REBO_num_skin,maxlocal,"AIREBO:cnumneigh");
    int max_nbors = buffers->get_max_nbors();
    memory->create(REBO_list_data,maxlocal * max_nbors,"AIREBO:list_data");
    REBO_firstneigh = (int **) memory->smalloc(maxlocal*sizeof(int *),
                                               "AIREBO:firstneigh");
    memory->create(nC,maxlocal,"AIREBO:nC");
    memory->create(nH,maxlocal,"AIREBO:nH");
    #ifdef _LMP_INTEL_OFFLOAD
    if (_cop >= 0) {
      int * const REBO_numneigh = this->REBO_numneigh;
      int * const REBO_num_skin = this->REBO_num_skin;
      int * const REBO_cnumneigh = this->REBO_cnumneigh;
      int * const REBO_list_data = this->REBO_list_data;
      double * const nC = this->nC;
      double * const nH = this->nH;
      const int mnml = max_nbors * maxlocal;
      #pragma offload_transfer target(mic:_cop) \
        nocopy(REBO_numneigh: length(maxlocal) alloc_if(1) free_if(0)) \
        nocopy(REBO_cnumneigh:length(maxlocal) alloc_if(1) free_if(0)) \
        nocopy(REBO_num_skin: length(maxlocal) alloc_if(1) free_if(0)) \
        nocopy(REBO_list_data:length(mnml) alloc_if(1) free_if(0)) \
        nocopy(nH: length(maxlocal) alloc_if(1) free_if(0)) \
        nocopy(nC: length(maxlocal) alloc_if(1) free_if(0))
    }
    #endif
  }

  if (evflag || vflag_fdotr) {
    int ovflag = 0;
    if (vflag_fdotr) ovflag = 2;
    else if (vflag) ovflag = 1;
    if (eflag) {
      eval<1,1>(1, ovflag, buffers, 0, offload_end);
      eval<1,1>(0, ovflag, buffers, host_start, inum);
    } else {
      eval<1,0>(1, ovflag, buffers, 0, offload_end);
      eval<1,0>(0, ovflag, buffers, host_start, inum);
    }
  } else {
    eval<0,0>(1, 0, buffers, 0, offload_end);
    eval<0,0>(0, 0, buffers, host_start, inum);
  }
}

/* ---------------------------------------------------------------------- */

template<int EVFLAG, int EFLAG, class flt_t, class acc_t>
void PairAIREBOIntel::eval(
    const int offload, const int vflag,
    IntelBuffers<flt_t,acc_t> * buffers,
    const int astart, const int aend
) {
  const int inum = aend - astart;
  if (inum == 0) {
    return;
  }
  int nlocal, nall, minlocal;
  fix->get_buffern(offload, nlocal, nall, minlocal);

  const int ago = neighbor->ago;
  IP_PRE_pack_separate_buffers(fix, buffers, ago, offload, nlocal, nall);

  ATOM_T * _noalias const x = buffers->get_x(offload);
  const int * _noalias const numneighhalf = buffers->get_atombin();
  const int * _noalias const numneigh = list->numneigh;
  const int ** _noalias const firstneigh = (const int **)list->firstneigh;
  tagint * const tag = atom->tag;

  const int ntypes = atom->ntypes + 1;
  const int eatom = this->eflag_atom;

  int x_size, q_size, f_stride, ev_size, separate_flag;
  IP_PRE_get_transfern(ago, 1 /*NEWTON_PAIR*/, EFLAG, vflag,
                       buffers, offload, fix, separate_flag,
                       x_size, q_size, ev_size, f_stride);

  int tc;
  FORCE_T * _noalias f_start;
  acc_t * _noalias ev_global;
  IP_PRE_get_buffers(offload, buffers, fix, tc, f_start, ev_global);

  const int nthreads = tc;
  const double skin = neighbor->skin;
  const int max_nbor = buffers->get_max_nbors();
  const PairAIREBOIntelParam<flt_t,acc_t> param = get_param<flt_t,acc_t>();

  // offload here
  #ifdef _LMP_INTEL_OFFLOAD
  int *overflow = fix->get_off_overflow_flag();
  double *timer_compute = fix->off_watch_pair();

  int * const REBO_numneigh = this->REBO_numneigh;
  int * const REBO_num_skin = this->REBO_num_skin;
  int * const REBO_cnumneigh = this->REBO_cnumneigh;
  int * const REBO_list_data = this->REBO_list_data;
  double * const nC = this->nC;
  double * const nH = this->nH;
  const int torflag = this->torflag;
  const int ljflag = this->ljflag;
  const int morseflag = this->morseflag;
  int * const map = this->map;

  if (offload) fix->start_watch(TIME_OFFLOAD_LATENCY);

  #pragma offload target(mic:_cop) if (offload) \
    in(firstneigh:length(0) alloc_if(0) free_if(0)) \
    in(numneigh:length(0) alloc_if(0) free_if(0)) \
    in(numneighhalf:length(0) alloc_if(0) free_if(0)) \
    in(x:length(x_size) alloc_if(0) free_if(0)) \
    in(overflow:length(0) alloc_if(0) free_if(0)) \
    in(astart,nthreads,inum,nall,ntypes,vflag,eatom) \
    in(f_stride,nlocal,minlocal,separate_flag,offload) \
    out(f_start:length(f_stride) alloc_if(0) free_if(0)) \
    out(ev_global:length(ev_size) alloc_if(0) free_if(0)) \
    out(timer_compute:length(1) alloc_if(0) free_if(0)) \
    in(param,skin,max_nbor) \
    in(tag: length(0) alloc_if(0) free_if(0)) \
    in(torflag, ljflag, morseflag, ago) \
    in(nC: length(0) alloc_if(0) free_if(0)) \
    in(nH: length(0) alloc_if(0) free_if(0)) \
    in(REBO_numneigh: length(0) alloc_if(0) free_if(0)) \
    in(REBO_cnumneigh: length(0) alloc_if(0) free_if(0)) \
    in(REBO_num_skin: length(0) alloc_if(0) free_if(0)) \
    in(REBO_list_data: length(0) alloc_if(0) free_if(0)) \
    in(map: length(0) alloc_if(0) free_if(0)) \
    signal(f_start)
  #endif
  {
    #if defined(__MIC__) && defined(_LMP_INTEL_OFFLOAD)
    *timer_compute = MIC_Wtime();
    #endif

    IP_PRE_repack_for_offload(1 /*NEWTON_PAIR*/, separate_flag, nlocal, nall,
                              f_stride, x, 0/*q*/);

    acc_t oevdwl, oecoul, ov0, ov1, ov2, ov3, ov4, ov5;
    if (EVFLAG)
      oevdwl = oecoul = ov0 = ov1 = ov2 = ov3 = ov4 = ov5 = (acc_t)0;

    // loop over neighbors of my atoms
    #if defined(_OPENMP)
    #pragma omp parallel \
      shared(f_start,f_stride,nlocal,nall,minlocal)     \
      reduction(+:oevdwl,oecoul,ov0,ov1,ov2,ov3,ov4,ov5)
    #endif
    {
      int iifrom, iito, tid;
      IP_PRE_omp_range_id(iifrom, iito, tid, inum, nthreads);
      iifrom += astart;
      iito += astart;

      int neigh_iifrom, neigh_iito;
      IP_PRE_omp_range(neigh_iifrom, neigh_iito, tid, nall, nthreads);

      FORCE_T * _noalias const f = f_start - minlocal + (tid * f_stride);
      memset(f + minlocal, 0, f_stride * sizeof(FORCE_T));

      KernelArgsAIREBOT<flt_t,acc_t> args;
      args.num_local = nlocal;
      args.num_all = nall;
      args.num_neighs_per_atom = max_nbor;
      args.num_types = ntypes;
      args.frebo_from_atom = 0;
      args.frebo_to_atom = args.num_local;
      args.neigh_from_atom = 0;
      args.neigh_to_atom = args.num_all;
      args.rebuild_flag = ago == 0;
      args.skin = skin;
      args.neigh_lmp.num = const_cast<int*>(numneigh);
      args.neigh_lmp.num_half = const_cast<int*>(numneighhalf);
      args.neigh_lmp.entries = const_cast<int**>(firstneigh);
      args.neigh_rebo.num = REBO_numneigh;
      args.neigh_rebo.num_half = REBO_num_skin;
      args.neigh_rebo.offset = REBO_cnumneigh;
      args.neigh_rebo.entries = REBO_list_data;
      args.params = param;
      args.tag = tag;
      args.nC = reinterpret_cast<flt_t*>(nC);
      args.nH = reinterpret_cast<flt_t*>(nH);
      args.map = map;
      args.result_eng = 0;
      args.x = (AtomAIREBOT<flt_t>*) x;

      args.result_f = (ResultForceT<acc_t> *) f;
      args.neigh_from_atom = neigh_iifrom;
      args.neigh_to_atom = neigh_iito;
      args.frebo_from_atom = iifrom;
      args.frebo_to_atom = iito;

      aut_rebo_neigh(&args);
      #if defined(_OPENMP)
      #pragma omp barrier
      #endif
      aut_frebo(&args, torflag);
      if (ljflag) aut_lennard_jones(&args, morseflag);

      oevdwl += args.result_eng;

      IP_PRE_fdotr_reduce_omp(1, nall, minlocal, nthreads, f_start, f_stride, x,
                              offload, vflag, ov0, ov1, ov2, ov3, ov4, ov5);
    } // end of omp parallel region
    IP_PRE_fdotr_reduce(1, nall, nthreads, f_stride, vflag,
                        ov0, ov1, ov2, ov3, ov4, ov5);
    if (EVFLAG) {
      ev_global[0] = oevdwl;
      ev_global[1] = oecoul;
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
  } // end of offload region

  if (offload)
    fix->stop_watch(TIME_OFFLOAD_LATENCY);
  else
    fix->stop_watch(TIME_HOST_PAIR);

  if (EVFLAG)
    fix->add_result_array(f_start, ev_global, offload, eatom, 0, vflag);
  else
    fix->add_result_array(f_start, 0, offload);
}

/* ---------------------------------------------------------------------- */

template<class flt_t, class acc_t>
void PairAIREBOIntel::pack_force_const(IntelBuffers<flt_t,acc_t> * buffers) {
  int tp1 = atom->ntypes + 1;

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

  #ifdef _LMP_INTEL_OFFLOAD
  if (_cop < 0) return;
  size_t VL = 512 / 8 / sizeof(flt_t);
  int ntypes = tp1;
  int tp1sq = tp1 * tp1;
  // TODO the lifecycle of "map" is currently not 100% correct
  // it might not be freed if this method is called more than once
  int * map = this->map;
  #pragma offload_transfer target(mic:_cop) \
    in(map: length(tp1) alloc_if(1) free_if(0))
  #endif

}

/* ----------------------------------------------------------------------
    Implementation
   ---------------------------------------------------------------------- */

namespace {

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(push, target(mic))
#endif

namespace overloaded {
  double sqrt(double a) { return ::sqrt(a); }
  float sqrt(float a) { return ::sqrtf(a); }
  double sin(double a) { return ::sin(a); }
  float sin(float a) { return ::sinf(a); }
  double cos(double a) { return ::cos(a); }
  float cos(float a) { return ::cosf(a); }
  double exp(double a) { return ::exp(a); }
  float exp(float a) { return ::expf(a); }
  double pow(double a, double b) { return ::pow(a, b); }
  float pow(float a, float b) { return ::powf(a, b); }
}

/* ----------------------------------------------------------------------
    Scalar AIREBO implementation, standalone, with massive code reuse
    compared to original code.
   ---------------------------------------------------------------------- */

#define CARBON 0
#define HYDROGEN 1
#define TOL 1.0e-9

template<typename T>
inline T fmin_nonan(T a, T b) {
  return a < b ? a : b;
}
template<typename T>
inline T fmax_nonan(T a, T b) {
  return a > b ? a : b;
}

template<typename flt_t>
inline flt_t Sp(flt_t r, flt_t lo, flt_t hi, flt_t * del) {
  flt_t t = (r - lo) / (hi - lo);
  if (t <= 0) {
    if (del) *del = 0;
    return 1;
  } else if (t >= 1) {
    if (del) *del = 0;
    return 0;
  } else {
    t *= static_cast<flt_t>(MY_PI);
    if (del) *del = static_cast<flt_t>(-0.5 * MY_PI)
                  * overloaded::sin(t) / (hi - lo);
    return static_cast<flt_t>(0.5) * (1 + overloaded::cos(t));
  }
}

template<typename flt_t>
inline flt_t Sp2(flt_t r, flt_t lo, flt_t hi, flt_t * del) {
  flt_t t = (r - lo) / (hi - lo);
  if (t <= 0) {
    if (del) *del = 0;
    return 1;
  } else if (t >= 1) {
    if (del) *del = 0;
    return 0;
  } else {
    if (del) *del = 6 * (t * t - t) / (hi - lo);
    return 1 - t * t * (3 - 2 * t);
  }
}

template<typename flt_t>
inline flt_t eval_poly_lin(int n, flt_t * coeffs, flt_t x, flt_t * deriv) {
  flt_t result = coeffs[n - 1];
  *deriv = coeffs[n - 1] * (n - 1);
  for (int i = n - 2; i > 0; i--) {
    result = coeffs[i] + x * result;
    *deriv = coeffs[i] * i + x * (*deriv);
  }
  result = coeffs[0] + x * result;
  return result;
}

template<typename flt_t, typename acc_t>
inline flt_t gSpline(KernelArgsAIREBOT<flt_t,acc_t> * ka, int itype, flt_t cos, flt_t N, flt_t * dgdc, flt_t * dgdN) {
  flt_t NCmin = ka->params.NCmin;
  flt_t NCmax = ka->params.NCmax;
  int index = 0;
  flt_t * gDom = nullptr;
  int nDom = 0;
  int offs = 0;
  if (itype == 0) {
    nDom = 4;
    gDom = &ka->params.gCdom[0];
    if (N > NCmin) offs = 4 * 6;
  } else {
    nDom = 3;
    gDom = &ka->params.gHdom[0];
    offs = 8 * 6;
  }
  cos = fmax_nonan(gDom[0], fmin_nonan(gDom[nDom], cos));
  int i;
  for (i = 0; i < nDom; i++) {
    if (cos >= gDom[i] && cos <= gDom[i + 1]) {
      index = i;
    }
  }
  flt_t g = eval_poly_lin(6, &ka->params.gVal[offs+index*6], cos, dgdc);
  *dgdN = 0;
  if (itype == 0 && N > NCmin && N < NCmax) {
    flt_t dg1;
    flt_t g1 = eval_poly_lin(6, &ka->params.gVal[index*6], cos, &dg1);
    flt_t dS;
    flt_t cut = Sp(N, NCmin, NCmax, &dS);
    *dgdN = dS * (g1 - g);
    g = g + cut * (g1 - g);
    *dgdc = *dgdc + cut * (dg1 - *dgdc);
  }
  return g;
}

template<typename flt_t>
inline flt_t eval_poly_bi(int n, flt_t * coeffs, flt_t x, flt_t y,
                          flt_t * deriv) {
  flt_t dy;
  flt_t vy = eval_poly_lin(n, &coeffs[n * (n - 1)], y, &dy);
  flt_t result = vy;
  deriv[0] = vy * (n - 1);
  deriv[1] = dy;
  for (int i = n - 2; i > 0; i--) {
    vy = eval_poly_lin(n, &coeffs[n * i], y, &dy);
    result = vy + x * result;
    deriv[0] = vy * i + x * deriv[0];
    deriv[1] = dy + x * deriv[1];
  }
  result = eval_poly_lin(n, &coeffs[0], y, &dy) + x * result;
  deriv[1] = dy + x * deriv[1];
  return result;
}

template<typename flt_t>
inline flt_t eval_poly_tri(int n, flt_t * coeffs, flt_t x, flt_t y, flt_t z,
                           flt_t * deriv) {
  flt_t dyz[2];
  flt_t vyz = eval_poly_bi(n, &coeffs[n * n * (n - 1)], y, z, &dyz[0]);
  flt_t result = vyz;
  deriv[0] = vyz * (n - 1);
  deriv[1] = dyz[0];
  deriv[2] = dyz[1];
  for (int i = n - 2; i > 0; i--) {
    vyz = eval_poly_bi(n, &coeffs[n * n * i], y, z, &dyz[0]);
    result = vyz + x * result;
    deriv[0] = vyz * i + x * deriv[0];
    deriv[1] = dyz[0] + x * deriv[1];
    deriv[2] = dyz[1] + x * deriv[2];
  }
  result = eval_poly_bi(n, &coeffs[0], y, z, &dyz[0]) + x * result;
  deriv[1] = dyz[0] + x * deriv[1];
  deriv[2] = dyz[1] + x * deriv[2];
  return result;
}

template<typename flt_t, typename acc_t>
inline flt_t PijSpline(KernelArgsAIREBOT<flt_t,acc_t> * ka, int itype,
                       int jtype, flt_t NC, flt_t NH, flt_t * dN) {
  dN[0] = 0.0;
  dN[1] = 0.0;
  if (itype == HYDROGEN) return 0;
  flt_t *pCJdom = jtype == CARBON ? &ka->params.pCCdom[0][0] :
    &ka->params.pCHdom[0][0];
  NC = fmax_nonan(pCJdom[0], fmin_nonan(pCJdom[1], NC));
  NH = fmax_nonan(pCJdom[2], fmin_nonan(pCJdom[3], NH));
  int nC = floor(NC);
  int nH = floor(NH);
  #define PijSelect(a, b) (jtype == CARBON ? ka->params.a : ka->params.b)
  if (fabs(NC - nC) < TOL && fabs(NH - nH) < TOL) {
    dN[0] = PijSelect(PCCdfdx, PCHdfdx)[nC][nH];
    dN[1] = PijSelect(PCCdfdy, PCHdfdy)[nC][nH];
    return PijSelect(PCCf, PCHf)[nC][nH];
  }
  if (NC == pCJdom[1]) nC -= 1;
  if (NH == pCJdom[3]) nH -= 1;
  return eval_poly_bi(4, &PijSelect(pCC, pCH)[nC][nH][0], NC, NH, dN);
  #undef PijSelect
}

template<typename flt_t, typename acc_t>
inline flt_t TijSpline(KernelArgsAIREBOT<flt_t,acc_t> * ka, flt_t Nij,
    flt_t Nji, flt_t Nijconj, acc_t * dN3) {
  flt_t * Tijdom = &ka->params.Tijdom[0][0];
  Nij = fmax_nonan(Tijdom[0], fmin_nonan(Tijdom[1], Nij));
  Nji = fmax_nonan(Tijdom[2], fmin_nonan(Tijdom[3], Nji));
  Nijconj = fmax_nonan(Tijdom[4], fmin_nonan(Tijdom[5], Nijconj));
  int nij = floor(Nij);
  int nji = floor(Nji);
  int nijconj = floor(Nijconj);
  if (fabs(Nij - nij) < TOL && fabs(Nji - nji) <
                          TOL && fabs(Nijconj - nijconj) < TOL) {
    dN3[0] = ka->params.Tdfdx[nij][nji][nijconj];
    dN3[1] = ka->params.Tdfdy[nij][nji][nijconj];
    dN3[2] = ka->params.Tdfdz[nij][nji][nijconj];
    return ka->params.Tf[nij][nji][nijconj];
  }
  if (Nij == Tijdom[1]) nij -= 1;
  if (Nji == Tijdom[3]) nji -= 1;
  if (Nijconj == Tijdom[5]) nijconj -= 1;
  return eval_poly_tri<acc_t>(4, &ka->params.Tijc[nij][nji][nijconj][0], Nij,
    Nji, Nijconj, dN3);
}

template<typename flt_t, typename acc_t>
inline flt_t piRCSpline(KernelArgsAIREBOT<flt_t,acc_t> * ka, int itype,
    int jtype, flt_t Nij, flt_t Nji, flt_t Nijconj, acc_t * dN3) {
  const int HH = 2;
  const int CH = 1;
  /* const int CC = 0; */
  int select = itype + jtype;
  #define piRCSelect(a, b, c) (select == HH ? ka->params.a : select == CH ? \
                               ka->params.b : ka->params.c)
  flt_t * piIJdom = &piRCSelect(piHHdom, piCHdom, piCCdom)[0][0];
  if (select == HH) {
    if (Nij < piIJdom[0] || Nij > piIJdom[1] || Nji < piIJdom[2] ||
        Nji > piIJdom[3] || Nijconj < piIJdom[4] || Nijconj > piIJdom[5]) {
      Nij = 0;
      Nji = 0;
      Nijconj = 0;
    }
  }
  Nij = fmax_nonan(piIJdom[0], fmin_nonan(piIJdom[1], Nij));
  Nji = fmax_nonan(piIJdom[2], fmin_nonan(piIJdom[3], Nji));
  Nijconj = fmax_nonan(piIJdom[4], fmin_nonan(piIJdom[5], Nijconj));
  int nij = floor(Nij);
  int nji = floor(Nji);
  int nijconj = floor(Nijconj);
  if (fabs(Nij - nij) < TOL && fabs(Nji - nji) <
                          TOL && fabs(Nijconj - nijconj) < TOL) {
    dN3[0] = piRCSelect(piHHdfdx, piCHdfdx, piCCdfdx)[nij][nji][nijconj];
    dN3[1] = piRCSelect(piHHdfdy, piCHdfdy, piCCdfdy)[nij][nji][nijconj];
    dN3[2] = piRCSelect(piHHdfdz, piCHdfdz, piCCdfdz)[nij][nji][nijconj];
    return piRCSelect(piHHf, piCHf, piCCf)[nij][nji][nijconj];
  }
  if (Nij == piIJdom[1]) nij -= 1;
  if (Nji == piIJdom[3]) nji -= 1;
  if (Nijconj == piIJdom[5]) nijconj -= 1;
  return eval_poly_tri<acc_t>(4,
    &piRCSelect(piHH, piCH, piCC)[nij][nji][nijconj][0], Nij, Nji, Nijconj,
    dN3);
  #undef piRCSelect
}

/*
 * Implements the p_ij term in airebo, which occurs on 4 different occasions
 * in the original lammps code.
 */
template<typename flt_t, typename acc_t>
inline flt_t frebo_pij(KernelArgsAIREBOT<flt_t,acc_t> * ka, int i, int j,
    flt_t rijx, flt_t rijy, flt_t rijz, flt_t rijmag, flt_t wij, flt_t VA,
    flt_t * sum_N, acc_t fij[3]) {
  ResultForceT<acc_t> * result_f = ka->result_f;
  AtomAIREBOT<flt_t> * x = ka->x;
  int * map = ka->map;
  flt_t * nC = ka->nC;
  flt_t * nH = ka->nH;
  flt_t x_i = x[i].x;
  flt_t y_i = x[i].y;
  flt_t z_i = x[i].z;
  int itype = map[x[i].w];
  int jtype = map[x[j].w];
  flt_t invrijm = 1 / rijmag;
  flt_t invrijm2 = invrijm * invrijm;
  flt_t rcminij = ka->params.rcmin[itype][jtype];
  flt_t rcmaxij = ka->params.rcmax[itype][jtype];
  flt_t Nmin = ka->params.Nmin;
  flt_t Nmax = ka->params.Nmax;
  flt_t Nij = nC[i] + nH[i] - wij;
  flt_t NijC = nC[i] - wij * (1 - jtype);
  flt_t NijH = nH[i] - wij * jtype;
  flt_t sum_pij = 0;
  flt_t sum_dpij_dN = 0;
  flt_t dN2[2] = {0};
  flt_t pij = 0;
  *sum_N = 0;
  int * neighs = ka->neigh_rebo.entries + ka->neigh_rebo.offset[i];
  int pass;
  for (pass = 0; pass < 2; pass++) {
    int kk;
    int knum = ka->neigh_rebo.num[i];
    for (kk = 0; kk < knum; kk++) {
      int k = neighs[kk];
      if (k == j) continue;
      flt_t rikx = x_i - x[k].x;
      flt_t riky = y_i - x[k].y;
      flt_t rikz = z_i - x[k].z;
      int ktype = map[x[k].w];
      flt_t rikmag = overloaded::sqrt(rikx * rikx + riky * riky + rikz * rikz);
      flt_t rho_k = ka->params.rho[ktype][1];
      flt_t rho_j = ka->params.rho[jtype][1];
      flt_t lamdajik = 4 * itype * ((rho_k - rikmag) - (rho_j - rijmag));
      flt_t ex_lam = exp(lamdajik);
      flt_t rcminik = ka->params.rcmin[itype][ktype];
      flt_t rcmaxik = ka->params.rcmax[itype][ktype];
      flt_t dwik;
      flt_t wik = Sp(rikmag, rcminik, rcmaxik, &dwik);
      flt_t Nki = nC[k] + nH[k] - wik;
      flt_t cosjik = (rijx * rikx + rijy * riky + rijz * rikz) /
        (rijmag * rikmag);
      cosjik = fmin_nonan<flt_t>(1, fmax_nonan<flt_t>(-1, cosjik));
      flt_t dgdc, dgdN;
      flt_t g = gSpline(ka, itype, cosjik, Nij, &dgdc, &dgdN);
      if (pass == 0) {
        sum_pij += wik * g * ex_lam;
        sum_dpij_dN += wik * dgdN * ex_lam;
        flt_t cutN = Sp<flt_t>(Nki, Nmin, Nmax, nullptr);
        *sum_N += (1 - ktype) * wik * cutN;
      } else {
        flt_t tmp = -0.5 * pij * pij * pij;
        flt_t invrikm = 1 / rikmag;
        flt_t rjkx = rikx - rijx;
        flt_t rjky = riky - rijy;
        flt_t rjkz = rikz - rijz;
        flt_t rjkmag = sqrt(rjkx * rjkx + rjky * rjky + rjkz * rjkz);
        flt_t rijrik = 2 * rijmag * rikmag;
        flt_t rr = rijmag * rijmag - rikmag * rikmag;
        flt_t dctdjk = -2 / rijrik;
        flt_t dctdik = (-rr + rjkmag * rjkmag) / (rijrik * rikmag * rikmag);
        flt_t dctdij = (rr + rjkmag * rjkmag) / (rijrik * rijmag * rijmag);

        acc_t fi[3], fj[3], fk[3];
        flt_t pref = 0.5 * VA * tmp;
        flt_t tmp20 = pref * wik * dgdc * ex_lam;
        fj[0] = fj[1] = fj[2] = 0;
        fi[0] = -tmp20 * dctdik * rikx;
        fi[1] = -tmp20 * dctdik * riky;
        fi[2] = -tmp20 * dctdik * rikz;
        fk[0] =  tmp20 * dctdik * rikx;
        fk[1] =  tmp20 * dctdik * riky;
        fk[2] =  tmp20 * dctdik * rikz;

        fij[0] += -tmp20 * dctdij * rijx;
        fij[1] += -tmp20 * dctdij * rijy;
        fij[2] += -tmp20 * dctdij * rijz;

        fi[0] += -tmp20 * dctdjk * rjkx;
        fi[1] += -tmp20 * dctdjk * rjky;
        fi[2] += -tmp20 * dctdjk * rjkz;
        fk[0] +=  tmp20 * dctdjk * rjkx;
        fk[1] +=  tmp20 * dctdjk * rjky;
        fk[2] +=  tmp20 * dctdjk * rjkz;
        fij[0] -= -tmp20 * dctdjk * rjkx;
        fij[1] -= -tmp20 * dctdjk * rjky;
        fij[2] -= -tmp20 * dctdjk * rjkz;

        flt_t tmp21 = pref * (wik * g * ex_lam * 4 * itype);
        fij[0] -= 1 * tmp21 * rijx * invrijm;
        fij[1] -= 1 * tmp21 * rijy * invrijm;
        fij[2] -= 1 * tmp21 * rijz * invrijm;
        fi[0] -= tmp21 * (-rikx * invrikm);
        fi[1] -= tmp21 * (-riky * invrikm);
        fi[2] -= tmp21 * (-rikz * invrikm);
        fk[0] -= tmp21 * (rikx * invrikm);
        fk[1] -= tmp21 * (riky * invrikm);
        fk[2] -= tmp21 * (rikz * invrikm);

        // coordination forces

        // dwik forces
        flt_t tmp22 = pref * dwik * g * ex_lam * invrikm;
        fi[0] -= tmp22 * rikx;
        fi[1] -= tmp22 * riky;
        fi[2] -= tmp22 * rikz;
        fk[0] += tmp22 * rikx;
        fk[1] += tmp22 * riky;
        fk[2] += tmp22 * rikz;

        // PIJ forces
        flt_t tmp23 = pref * dN2[ktype] * dwik * invrikm;
        fi[0] -= tmp23 * rikx;
        fi[1] -= tmp23 * riky;
        fi[2] -= tmp23 * rikz;
        fk[0] += tmp23 * rikx;
        fk[1] += tmp23 * riky;
        fk[2] += tmp23 * rikz;

        // dgdN forces
        flt_t tmp24 = pref * sum_dpij_dN * dwik * invrikm;
        fi[0] -= tmp24 * rikx;
        fi[1] -= tmp24 * riky;
        fi[2] -= tmp24 * rikz;
        fk[0] += tmp24 * rikx;
        fk[1] += tmp24 * riky;
        fk[2] += tmp24 * rikz;

        result_f[i].x += fi[0];
        result_f[i].y += fi[1];
        result_f[i].z += fi[2];
        result_f[j].x += fj[0];
        result_f[j].y += fj[1];
        result_f[j].z += fj[2];
        result_f[k].x += fk[0];
        result_f[k].y += fk[1];
        result_f[k].z += fk[2];
      }
    }
    if (pass == 0) {
      flt_t PijS = PijSpline(ka, itype, jtype, NijC, NijH, dN2);
      pij = 1 / overloaded::sqrt(1 + sum_pij + PijS);
    }
  }
  return pij;
}

template<typename flt_t, typename acc_t>
inline flt_t frebo_pi_rc(KernelArgsAIREBOT<flt_t,acc_t> * ka, int itype,
    int jtype, flt_t Nij, flt_t Nji, flt_t Nijconj, flt_t * dN3) {
  acc_t dN3tmp[3] = {0};
  flt_t ret = piRCSpline(ka, itype, jtype, Nij, Nji, Nijconj, dN3tmp);
  dN3[0] = dN3tmp[0];
  dN3[1] = dN3tmp[1];
  dN3[2] = dN3tmp[2];
  return ret;
}

template<typename flt_t, typename acc_t>
inline flt_t frebo_Tij(KernelArgsAIREBOT<flt_t,acc_t> * ka, int itype,
    int jtype, flt_t Nij, flt_t Nji, flt_t Nijconj, flt_t * dN3) {
  dN3[0] = 0;
  dN3[1] = 0;
  dN3[2] = 0;
  if (itype == HYDROGEN || jtype == HYDROGEN) return 0;
  acc_t dN3tmp[3] = {0};
  flt_t ret = TijSpline(ka, Nij, Nji, Nijconj, dN3tmp);
  dN3[0] = dN3tmp[0];
  dN3[1] = dN3tmp[1];
  dN3[2] = dN3tmp[2];
  return ret;
}

/*
 * Implements a scalar version of the sum cos^1(omega) term used in pi^dh_ij.
 * Occurs in both bondorder and bondorderLJ.
 */
template<typename flt_t, typename acc_t>
inline flt_t frebo_sum_omega(KernelArgsAIREBOT<flt_t,acc_t> * ka, int i, int j,
 flt_t r23x, flt_t r23y, flt_t r23z, flt_t r23mag, flt_t VA, acc_t fij[3]) {
  ResultForceT<acc_t> * result_f = ka->result_f;
  acc_t sum_omega = 0;
  int a2 = i;
  int a3 = j;
  flt_t r32x = - r23x;
  flt_t r32y = - r23y;
  flt_t r32z = - r23z;
  int * map = ka->map;
  AtomAIREBOT<flt_t> * x = ka->x;
  flt_t thmin = ka->params.thmin;
  flt_t thmax = ka->params.thmax;
  int itype = map[x[i].w];
  int jtype = map[x[j].w];
  int * neighs_i = ka->neigh_rebo.entries + ka->neigh_rebo.offset[i];
  int * neighs_j = ka->neigh_rebo.entries + ka->neigh_rebo.offset[j];
  int num_i = ka->neigh_rebo.num[i];
  int num_j = ka->neigh_rebo.num[j];
  int kk;
  for (kk = 0; kk < num_i; kk++) {
    int k = neighs_i[kk];
    if (k == j) continue;
    int a1 = k;
    int ktype = map[x[k].w];
    flt_t r21x = x[a2].x - x[a1].x;
    flt_t r21y = x[a2].y - x[a1].y;
    flt_t r21z = x[a2].z - x[a1].z;
    flt_t r21mag = overloaded::sqrt(r21x * r21x + r21y * r21y + r21z * r21z);
    flt_t cos321 = (r23x * r21x + r23y * r21y + r23z * r21z) /
      (r23mag * r21mag);
    cos321 = fmin_nonan<flt_t>(1, fmax_nonan<flt_t>(-1, cos321));
    flt_t sin321 = overloaded::sqrt(1 - cos321 * cos321);
    if (sin321 == 0) continue;
    flt_t sink2i = 1 / (sin321 * sin321);
    flt_t rik2i = 1 / (r21mag * r21mag);
    flt_t rr = r23mag * r23mag - r21mag * r21mag;
    flt_t r31x = r21x - r23x;
    flt_t r31y = r21y - r23y;
    flt_t r31z = r21z - r23z;
    flt_t r31mag2 = r31x * r31x + r31y * r31y + r31z * r31z;
    flt_t rijrik = 2 * r23mag * r21mag;
    flt_t r21mag2 = r21mag * r21mag;
    flt_t dctik = (-rr + r31mag2) / (rijrik * r21mag2);
    flt_t dctij = (rr + r31mag2) / (rijrik * r23mag * r23mag);
    flt_t dctjk = -2 / rijrik;
    flt_t rcmin21  = ka->params.rcmin [itype][ktype];
    flt_t rcmaxp21 = ka->params.rcmaxp[itype][ktype];
    flt_t dw21;
    flt_t w21 = Sp(r21mag, rcmin21, rcmaxp21, &dw21);
    // why does this additional cutoff in the cosine exist?
    // the original code by stuart answers this:
    // it avoid issues when bonds in the dihedral are linear
    // by switching the dihedral off beforehand.
    // This is the reason for both the sin == 0 checks and the
    // tspjik = Sp2(..) calls.
    // Unfortunately, this is not exactly stated in the original paper.
    // It might be similar in purpose to the H(sin - s^min) term that
    // appears in that paper, but can not be found in original REBO papers.
    flt_t dtsjik;
    flt_t tspjik = Sp2(cos321, thmin, thmax, &dtsjik);
    dtsjik = - dtsjik;
    int ll;
    for (ll = 0; ll < num_j; ll++) {
      int l = neighs_j[ll];
      if (l == i || l == k) continue;
      int ltype = map[x[l].w];
      int a4 = l;
      flt_t r34x = x[a3].x - x[a4].x;
      flt_t r34y = x[a3].y - x[a4].y;
      flt_t r34z = x[a3].z - x[a4].z;
      flt_t r34mag = overloaded::sqrt(r34x * r34x + r34y * r34y + r34z * r34z);
      flt_t cos234 = (r32x * r34x + r32y * r34y + r32z * r34z) /
        (r23mag * r34mag);
      cos234 = fmin_nonan<flt_t>(1, fmax_nonan<flt_t>(-1, cos234));
      flt_t sin234 = overloaded::sqrt(1 - cos234 * cos234);
      if (sin234 == 0) continue;
      flt_t sinl2i = 1 / (sin234 * sin234);
      flt_t rjl2i = 1 / (r34mag * r34mag);

      flt_t rcminjl = ka->params.rcmin[jtype][ltype];
      flt_t rcmaxpjl = ka->params.rcmaxp[jtype][ltype];
      flt_t dw34;
      flt_t w34 = Sp(r34mag, rcminjl, rcmaxpjl, &dw34);
      flt_t rr = (r23mag * r23mag) - (r34mag * r34mag);
      flt_t r24x = r23x + r34x;
      flt_t r24y = r23y + r34y;
      flt_t r24z = r23z + r34z;
      flt_t r242 =
          (r24x * r24x) + (r24y * r24y) + (r24z * r24z);
      flt_t rijrjl = 2 * r23mag * r34mag;
      flt_t rjl2 = r34mag * r34mag;
      flt_t dctjl = (-rr + r242) / (rijrjl * rjl2);
      flt_t dctji = (rr + r242) / (rijrjl * r23mag * r23mag);
      flt_t dctil = -2 / rijrjl;
      flt_t dtsijl;
      flt_t tspijl = Sp2(cos234, thmin, thmax, &dtsijl);
      dtsijl = -dtsijl; // need minus sign
      flt_t prefactor = VA;

      flt_t cross321x = (r32y * r21z) - (r32z * r21y);
      flt_t cross321y = (r32z * r21x) - (r32x * r21z);
      flt_t cross321z = (r32x * r21y) - (r32y * r21x);
      flt_t cross234x = (r23y * r34z) - (r23z * r34y);
      flt_t cross234y = (r23z * r34x) - (r23x * r34z);
      flt_t cross234z = (r23x * r34y) - (r23y * r34x);

      flt_t cwnum = (cross321x * cross234x) +
              (cross321y * cross234y) +
              (cross321z * cross234z);
      flt_t cwnom = r21mag * r34mag * r23mag * r23mag * sin321 * sin234;
      flt_t om1234 = cwnum / cwnom;
      flt_t cw = om1234;
      sum_omega += ((1 - (om1234 * om1234)) * w21 * w34) *
              (1 - tspjik) * (1 - tspijl);
      if (VA == static_cast<flt_t>(0.0)) continue;

      flt_t dt1dik = (rik2i) - (dctik * sink2i * cos321);
      flt_t dt1djk = (-dctjk * sink2i * cos321);
      flt_t dt1djl = (rjl2i) - (dctjl * sinl2i * cos234);
      flt_t dt1dil = (-dctil * sinl2i * cos234);
      flt_t dt1dij = (2 / (r23mag * r23mag)) -
               (dctij * sink2i * cos321) -
               (dctji * sinl2i * cos234);

      flt_t dt2dikx = (-r23z * cross234y) + (r23y * cross234z);
      flt_t dt2diky = (-r23x * cross234z) + (r23z * cross234x);
      flt_t dt2dikz = (-r23y * cross234x) + (r23x * cross234y);

      flt_t dt2djlx = (-r23y * cross321z) + (r23z * cross321y);
      flt_t dt2djly = (-r23z * cross321x) + (r23x * cross321z);
      flt_t dt2djlz = (-r23x * cross321y) + (r23y * cross321x);

      flt_t dt2dijx = (r21z * cross234y) - (r34z * cross321y) -
      flt_t      (r21y * cross234z) + (r34y * cross321z);
      flt_t dt2dijy = (r21x * cross234z) - (r34x * cross321z) -
      flt_t      (r21z * cross234x) + (r34z * cross321x);
      flt_t dt2dijz = (r21y * cross234x) - (r34y * cross321x) -
      flt_t      (r21x * cross234y) + (r34x * cross321y);

      flt_t aa = (prefactor * 2 * cw / cwnom) * w21 * w34 *
           (1 - tspjik) * (1 - tspijl);
      flt_t aaa1 = -prefactor * (1 - (om1234 * om1234)) *
             (1 - tspjik) * (1 - tspijl);
      flt_t aaa2 = -prefactor * (1 - (om1234 * om1234)) * w21 * w34;
      flt_t at2 = aa * cwnum;

      flt_t fcijpc = (-dt1dij * at2) +
              (aaa2 * dtsjik * dctij * (1 - tspijl)) +
              (aaa2 * dtsijl * dctji * (1 - tspjik));
      flt_t fcikpc = (-dt1dik * at2) +
              (aaa2 * dtsjik * dctik * (1 - tspijl));
      flt_t fcjlpc = (-dt1djl * at2) +
              (aaa2 * dtsijl * dctjl * (1 - tspjik));
      flt_t fcjkpc = (-dt1djk * at2) +
              (aaa2 * dtsjik * dctjk * (1 - tspijl));
      flt_t fcilpc = (-dt1dil * at2) +
              (aaa2 * dtsijl * dctil * (1 - tspjik));

      flt_t F23x = (fcijpc * r23x) + (aa * dt2dijx);
      flt_t F23y = (fcijpc * r23y) + (aa * dt2dijy);
      flt_t F23z = (fcijpc * r23z) + (aa * dt2dijz);

      flt_t F12x = (fcikpc * r21x) + (aa * dt2dikx);
      flt_t F12y = (fcikpc * r21y) + (aa * dt2diky);
      flt_t F12z = (fcikpc * r21z) + (aa * dt2dikz);

      flt_t F34x = (fcjlpc * r34x) + (aa * dt2djlx);
      flt_t F34y = (fcjlpc * r34y) + (aa * dt2djly);
      flt_t F34z = (fcjlpc * r34z) + (aa * dt2djlz);

      flt_t F31x = (fcjkpc * r31x);
      flt_t F31y = (fcjkpc * r31y);
      flt_t F31z = (fcjkpc * r31z);

      flt_t F24x = (fcilpc * r24x);
      flt_t F24y = (fcilpc * r24y);
      flt_t F24z = (fcilpc * r24z);

      flt_t f1x = -F12x - F31x;
      flt_t f1y = -F12y - F31y;
      flt_t f1z = -F12z - F31z;
      flt_t f2x = F12x + F31x;
      flt_t f2y = F12y + F31y;
      flt_t f2z = F12z + F31z;
      flt_t f3x = F34x + F24x;
      flt_t f3y = F34y + F24y;
      flt_t f3z = F34z + F24z;
      flt_t f4x = -F34x - F24x;
      flt_t f4y = -F34y - F24y;
      flt_t f4z = -F34z - F24z;

      fij[0] += F23x + F24x - F31x;
      fij[1] += F23y + F24y - F31y;
      fij[2] += F23z + F24z - F31z;

      // coordination forces

      flt_t tmp20 = VA * ((1 - (om1234 * om1234))) *
             (1 - tspjik) * (1 - tspijl) * dw21 * w34 / r21mag;
      f2x -= tmp20 * r21x;
      f2y -= tmp20 * r21y;
      f2z -= tmp20 * r21z;
      f1x += tmp20 * r21x;
      f1y += tmp20 * r21y;
      f1z += tmp20 * r21z;

      flt_t tmp21 = VA * ((1 - (om1234 * om1234))) *
             (1 - tspjik) * (1 - tspijl) * w21 * dw34 / r34mag;
      f3x -= tmp21 * r34x;
      f3y -= tmp21 * r34y;
      f3z -= tmp21 * r34z;
      f4x += tmp21 * r34x;
      f4y += tmp21 * r34y;
      f4z += tmp21 * r34z;

      result_f[a1].x += f1x;
      result_f[a1].y += f1y;
      result_f[a1].z += f1z;
      result_f[a2].x += f2x;
      result_f[a2].y += f2y;
      result_f[a2].z += f2z;
      result_f[a3].x += f3x;
      result_f[a3].y += f3y;
      result_f[a3].z += f3z;
      result_f[a4].x += f4x;
      result_f[a4].y += f4y;
      result_f[a4].z += f4z;
    }
  }
  return sum_omega;
}

/*
 * Implements a scalar implementation the force update due to splines.
 * It is used for both pi^rc_ij and T_ij.
 * Occurs four times in each bondorder and bondorderLJ.
 */
template<typename flt_t, typename acc_t>
inline void frebo_N_spline_force(KernelArgsAIREBOT<flt_t,acc_t> * ka, int i,
    int j, flt_t VA, flt_t dN, flt_t dNconj, flt_t Nconj) {
  int * map = ka->map;
  AtomAIREBOT<flt_t> * x = ka->x;
  ResultForceT<acc_t> * result_f = ka->result_f;
  flt_t * nC = ka->nC;
  flt_t * nH = ka->nH;
  flt_t Nmin = ka->params.Nmin;
  flt_t Nmax = ka->params.Nmax;
  int itype = map[x[i].w];
  int * neighs = ka->neigh_rebo.entries + ka->neigh_rebo.offset[i];
  int knum = ka->neigh_rebo.num[i];
  int kk;
  for (kk = 0; kk < knum; kk++) {
    int k = neighs[kk];
    if (k == j) continue;
    flt_t rikx = x[i].x - x[k].x;
    flt_t riky = x[i].y - x[k].y;
    flt_t rikz = x[i].z - x[k].z;
    flt_t rikmag = overloaded::sqrt(rikx * rikx + riky * riky + rikz * rikz);
    int ktype = map[x[k].w];
    flt_t rcminik = ka->params.rcmin[itype][ktype];
    flt_t rcmaxik = ka->params.rcmax[itype][ktype];
    flt_t dwik;
    flt_t wik = Sp(rikmag, rcminik, rcmaxik, &dwik);
    flt_t Nki = nC[k] + nH[k] - wik;
    flt_t dNki;
    flt_t SpN = Sp(Nki, Nmin, Nmax, &dNki);
    flt_t fdN = VA * dN * dwik / rikmag;
    flt_t fdNconj = VA * dNconj * 2 * Nconj * dwik * SpN / rikmag;
    flt_t ffactor = fdN;
    if (ktype == 0) ffactor += fdNconj;
    flt_t fkx = ffactor * rikx;
    flt_t fky = ffactor * riky;
    flt_t fkz = ffactor * rikz;
    result_f[i].x -= fkx;
    result_f[i].y -= fky;
    result_f[i].z -= fkz;
    result_f[k].x += fkx;
    result_f[k].y += fky;
    result_f[k].z += fkz;
    if (ktype != 0 || fabs(dNki) <= TOL) continue;
    int * neighs_k = ka->neigh_rebo.entries + ka->neigh_rebo.offset[k];
    int nnum = ka->neigh_rebo.num[k];
    int nn;
    for (nn = 0; nn < nnum; nn++) {
      int n = neighs_k[nn];
      if (n == i) continue;
      flt_t rknx = x[k].x - x[n].x;
      flt_t rkny = x[k].y - x[n].y;
      flt_t rknz = x[k].z - x[n].z;
      flt_t rknmag = overloaded::sqrt(rknx * rknx + rkny * rkny + rknz * rknz);
      int ntype = map[x[n].w];
      flt_t rcminkn = ka->params.rcmin[ktype][ntype];
      flt_t rcmaxkn = ka->params.rcmax[ktype][ntype];
      flt_t dwkn;
      Sp(rknmag, rcminkn, rcmaxkn, &dwkn);
      flt_t ffactor = VA * dNconj * 2 * Nconj * wik * dNki * dwkn / rknmag;
      result_f[k].x -= ffactor * rknx;
      result_f[k].y -= ffactor * rkny;
      result_f[k].z -= ffactor * rknz;
      result_f[n].x += ffactor * rknx;
      result_f[n].y += ffactor * rkny;
      result_f[n].z += ffactor * rknz;
    }
  }
}

/*
 * This data-structure contains the result of a search through neighbor-lists.
 * It is used to calculate C_ij and the corresponding force updates.
 */
template<typename flt_t>
struct LennardJonesPathAIREBOT {
  AtomAIREBOT<flt_t> del[3];
  int num;
  flt_t w[3];
  flt_t dw[3];
  flt_t r[3];
  int idx[4];
};

/*
 * Checks a candidate path stored in idxs whether it is better than *path
 * and updates *path accordingly.
 */
template<typename flt_t, typename acc_t>
inline flt_t ref_lennard_jones_test_path_single(
 KernelArgsAIREBOT<flt_t,acc_t> * ka, flt_t best, int num, int * idxs,
 LennardJonesPathAIREBOT<flt_t> * path) {
  LennardJonesPathAIREBOT<flt_t> result;
  AtomAIREBOT<flt_t> * x = ka->x;
  int * map = ka->map;
  result.num = num;
  flt_t combined = 1;
  for (int i = num - 2; i >= 0; i--) {
    int a0 = idxs[i+0];
    int a1 = idxs[i+1];
    flt_t delx = x[a1].x - x[a0].x;
    flt_t dely = x[a1].y - x[a0].y;
    flt_t delz = x[a1].z - x[a0].z;
    flt_t rsq = delx * delx + dely * dely + delz * delz;
    int type0 = map[x[a0].w];
    int type1 = map[x[a1].w];
    if (rsq >= ka->params.rcmaxsq[type0][type1]) return best;
    flt_t r = overloaded::sqrt(rsq);
    flt_t dw, w = Sp<flt_t>(r, ka->params.rcmin[type0][type1],
                            ka->params.rcmax[type0][type1], &dw);
    if (w == 0) return best;
    combined *= w;
    if (combined <= best) return best;
    result.idx[i] = a0;
    result.del[i].x = delx;
    result.del[i].y = dely;
    result.del[i].z = delz;
    result.r[i] = r;
    result.w[i] = w;
    result.dw[i] = dw;
  }
  result.idx[num - 1] = idxs[num - 1];
  *path = result;
  return combined;
}

/*
 * Test through all paths surrounding i and j to find the corresponding
 * best path. Uses the same iteration ordering as FLJ() does.
 * Note that an optimization would use the j neighlist instead in the inner
 * loop.
 */
template<typename flt_t, typename acc_t>
inline flt_t ref_lennard_jones_test_path(KernelArgsAIREBOT<flt_t,acc_t> * ka,
    int i, int j, flt_t rij, flt_t rcmax,
    LennardJonesPathAIREBOT<flt_t> * path) {
  int idxs[4];
  idxs[0] = i;
  idxs[1] = j;
  flt_t best = 0;
  if (rij <= rcmax) {
    best = ref_lennard_jones_test_path_single(ka, best, 2, idxs, path);
    if (best == static_cast<flt_t>(1.0)) return 0;
  }
  for (int kk = 0; kk < ka->neigh_rebo.num[i]; kk++) {
    int k = ka->neigh_rebo.entries[ka->neigh_rebo.offset[i] + kk];
    if (k == j) continue;
    idxs[1] = k;
    idxs[2] = j;
    best = ref_lennard_jones_test_path_single(ka, best, 3, idxs, path);
    if (best == static_cast<flt_t>(1.0)) return 0;
    for (int mm = 0; mm < ka->neigh_rebo.num[k]; mm++) {
      int m = ka->neigh_rebo.entries[ka->neigh_rebo.offset[k] + mm];
      if (m == i || m == j) continue;
      idxs[2] = m;
      idxs[3] = j;
      best = ref_lennard_jones_test_path_single(ka, best, 4, idxs, path);
      if (best == static_cast<flt_t>(1.0)) return 0;
    }
  }
  return 1 - best;
}

/*
 * Conducts the force update due to C_ij, given the active path.
 */
template<typename flt_t, typename acc_t>
inline void ref_lennard_jones_force_path(KernelArgsAIREBOT<flt_t,acc_t> * ka,
    flt_t dC, LennardJonesPathAIREBOT<flt_t> * path) {
  AtomAIREBOT<flt_t> * x = ka->x;
  ResultForceT<acc_t> * result_f = ka->result_f;
  for (int i = 0; i < path->num - 1; i++) {
    flt_t fpair = dC * path->dw[i] / path->r[i];
    for (int j = 0; j < path->num - 1; j++) {
      if (i != j) fpair *= path->w[j];
    }
    result_f[path->idx[i+0]].x -= fpair * path->del[i].x;
    result_f[path->idx[i+0]].y -= fpair * path->del[i].y;
    result_f[path->idx[i+0]].z -= fpair * path->del[i].z;
    result_f[path->idx[i+1]].x += fpair * path->del[i].x;
    result_f[path->idx[i+1]].y += fpair * path->del[i].y;
    result_f[path->idx[i+1]].z += fpair * path->del[i].z;
  }
}

/*
 * Calculate the bondorderLJ term.
 */
template<typename flt_t, typename acc_t>
inline flt_t ref_lennard_jones_bondorder(KernelArgsAIREBOT<flt_t,acc_t> * ka,
    int i, int j, flt_t VA, acc_t fij[3]) {
  AtomAIREBOT<flt_t> * x = ka->x;
  int * map = ka->map;
  ResultForceT<acc_t> * result_f = ka->result_f;

  int itype = map[x[i].w];
  int jtype = map[x[j].w];

  flt_t delx = x[i].x - x[j].x;
  flt_t dely = x[i].y - x[j].y;
  flt_t delz = x[i].z - x[j].z;
  flt_t rsq = delx * delx + dely * dely + delz * delz;
  flt_t rij = overloaded::sqrt(rsq);

  flt_t rcminij = ka->params.rcmin[itype][jtype];
  flt_t rcmaxij = ka->params.rcmax[itype][jtype];
  flt_t dwij;
  flt_t wij = Sp(rij, rcminij, rcmaxij, &dwij);

  flt_t the_r = ka->params.rcmin[itype][jtype];
  flt_t scale = the_r / rij;
  flt_t Nij = ka->nH[i] + ka->nC[i] - wij;
  flt_t Nji = ka->nH[j] + ka->nC[j] - wij;
  flt_t NconjtmpI;
  acc_t fijc[3] = {0}, fjic[3] = {0};
  flt_t pij = frebo_pij<flt_t,acc_t>(ka, i, j, delx * scale, dely * scale,
    delz * scale, the_r, wij, 0.0, &NconjtmpI, fijc);
  flt_t NconjtmpJ;
  flt_t pji = frebo_pij<flt_t,acc_t>(ka, j, i, -delx * scale, -dely * scale,
    -delz * scale, the_r, wij, 0.0, &NconjtmpJ, fjic);
  flt_t Nijconj = 1.0 + (NconjtmpI * NconjtmpI) + (NconjtmpJ * NconjtmpJ);
  flt_t dN3_pi_rc[3];
  flt_t pi_rc = frebo_pi_rc<flt_t,acc_t>(ka, itype, jtype, Nij, Nji, Nijconj,
    dN3_pi_rc);
  flt_t dN3_Tij[3];
  flt_t Tij = frebo_Tij<flt_t,acc_t>(ka, itype, jtype, Nij, Nji, Nijconj,
    dN3_Tij);
  flt_t sum_omega = 0;
  if (fabs(Tij) > TOL) {
    sum_omega = frebo_sum_omega<flt_t,acc_t>(ka, i, j, delx * scale, dely *
                                             scale, delz * scale, the_r, 0.0,
                                             fijc);
  }
  flt_t pi_dh = Tij * sum_omega;
  flt_t bij = 0.5 * (pij + pji) + pi_rc + pi_dh;
  flt_t dStb;
  flt_t Stb = Sp2<flt_t>(bij, ka->params.bLJmin[itype][jtype],
    ka->params.bLJmax[itype][jtype], &dStb);
  if (dStb != 0) {
    flt_t pij_reverse = frebo_pij<flt_t,acc_t>(ka, i, j, delx * scale,
      dely * scale, delz * scale, the_r, wij, VA * dStb, &NconjtmpI, fijc);
    flt_t pji_reverse = frebo_pij<flt_t,acc_t>(ka, j, i, -delx * scale,
      -dely * scale, -delz * scale, the_r, wij, VA * dStb, &NconjtmpJ, fjic);
    fijc[0] -= fjic[0];
    fijc[1] -= fjic[1];
    fijc[2] -= fjic[2];
    frebo_N_spline_force<flt_t,acc_t>(ka, i, j, VA * dStb, dN3_pi_rc[0],
      dN3_pi_rc[2], NconjtmpI);
    frebo_N_spline_force<flt_t,acc_t>(ka, j, i, VA * dStb, dN3_pi_rc[1],
      dN3_pi_rc[2], NconjtmpJ);
    if (fabs(Tij) > TOL) {
      flt_t sum_omega_reverse = frebo_sum_omega<flt_t,acc_t>(ka, i, j,
        delx * scale, dely * scale, delz * scale, the_r, VA * dStb * Tij, fijc);
      frebo_N_spline_force(ka, i, j, VA * dStb * sum_omega, dN3_Tij[0],
        dN3_Tij[2], NconjtmpI);
      frebo_N_spline_force(ka, j, i, VA * dStb * sum_omega, dN3_Tij[1],
        dN3_Tij[2], NconjtmpJ);
    }
    assert(fij[0] == 0);
    assert(fij[1] == 0);
    assert(fij[2] == 0);
    fij[0] = scale * (fijc[0] - (delx * delx * fijc[0] + dely * delx *
                                 fijc[1] + delz * delx * fijc[2]) / rsq);
    fij[1] = scale * (fijc[1] - (delx * dely * fijc[0] + dely * dely *
                                 fijc[1] + delz * dely * fijc[2]) / rsq);
    fij[2] = scale * (fijc[2] - (delx * delz * fijc[0] + dely * delz *
                                 fijc[1] + delz * delz * fijc[2]) / rsq);
  }
  return Stb;
}

/*
 * Scalar reference implementation of neighbor routine.
 */
template<typename flt_t, typename acc_t>
void ref_rebo_neigh(KernelArgsAIREBOT<flt_t,acc_t> * ka) {
  int offset = ka->neigh_from_atom * ka->num_neighs_per_atom;
  for (int i = ka->neigh_from_atom; i < ka->neigh_to_atom; i++) {
    ka->neigh_rebo.offset[i] = offset;
    int itype = ka->map[ka->x[i].w];
    int n = 0;
    ka->nC[i] = 0;
    ka->nH[i] = 0;
    for (int j = 0; j < ka->neigh_lmp.num[i]; j++) {
      int ji = ka->neigh_lmp.entries[i][j];
      flt_t delx = ka->x[i].x - ka->x[ji].x;
      flt_t dely = ka->x[i].y - ka->x[ji].y;
      flt_t delz = ka->x[i].z - ka->x[ji].z;
      flt_t rsq = delx * delx + dely * dely + delz * delz;
      int jtype = ka->map[ka->x[ji].w];
      if (rsq < ka->params.rcmaxsq[itype][jtype]) {
        ka->neigh_rebo.entries[offset + n++] = ji;
        flt_t rcmin = ka->params.rcmin[itype][jtype];
        flt_t rcmax = ka->params.rcmax[itype][jtype];
        if (jtype == CARBON)
          ka->nC[i] += Sp<flt_t>(overloaded::sqrt(rsq), rcmin, rcmax, nullptr);
        else
          ka->nH[i] += Sp<flt_t>(overloaded::sqrt(rsq), rcmin, rcmax, nullptr);
      }
    }
    ka->neigh_rebo.num[i] = n;
    offset += n;
  }
}

template<typename flt_t, typename acc_t>
void ref_torsion_single_interaction(KernelArgsAIREBOT<flt_t,acc_t> * ka, int i,
                                    int j) {
  AtomAIREBOT<flt_t> * x = ka->x;
  int * map = ka->map;
  ResultForceT<acc_t> * f = ka->result_f;
  flt_t (*rcmin)[2] = ka->params.rcmin;
  flt_t (*rcmax)[2] = ka->params.rcmax;
  flt_t (*epsilonT)[2] = ka->params.epsilonT;
  flt_t thmin = ka->params.thmin;
  flt_t thmax = ka->params.thmax;
  int itype = map[x[i].w];
  flt_t xtmp = x[i].x;
  flt_t ytmp = x[i].y;
  flt_t ztmp = x[i].z;
  int * REBO_neighs_i = &ka->neigh_rebo.entries[ka->neigh_rebo.offset[i]];
  int jnum = ka->neigh_rebo.num[i];
  int jtype = map[x[j].w];

  flt_t del32x = x[j].x-x[i].x;
  flt_t del32y = x[j].y-x[i].y;
  flt_t del32z = x[j].z-x[i].z;
  flt_t rsq = del32x*del32x + del32y*del32y + del32z*del32z;
  flt_t r32 = overloaded::sqrt(rsq);
  flt_t del23x = -del32x;
  flt_t del23y = -del32y;
  flt_t del23z = -del32z;
  flt_t r23 = r32;
  flt_t dw23, w23 = Sp<flt_t>(r23,rcmin[itype][jtype],rcmax[itype][jtype],
    &dw23);

  assert(itype == 0);
  assert(jtype == 0);

  for (int kk = 0; kk < jnum; kk++) {
    int k = REBO_neighs_i[kk];
    int ktype = map[x[k].w];
    if (k == j) continue;
    flt_t del21x = x[i].x-x[k].x;
    flt_t del21y = x[i].y-x[k].y;
    flt_t del21z = x[i].z-x[k].z;
    flt_t rsq = del21x*del21x + del21y*del21y + del21z*del21z;
    flt_t r21 = overloaded::sqrt(rsq);
    flt_t cos321 = - ((del21x*del32x) + (del21y*del32y) +
                (del21z*del32z)) / (r21*r32);
    cos321 = fmin(cos321,1);
    cos321 = fmax(cos321,-1);
    flt_t sin321 = overloaded::sqrt(1 - cos321*cos321);
    if (sin321 < TOL) continue;

    flt_t deljkx = del21x-del23x;
    flt_t deljky = del21y-del23y;
    flt_t deljkz = del21z-del23z;
    flt_t rjk2 = deljkx*deljkx + deljky*deljky + deljkz*deljkz;
    flt_t rjk = overloaded::sqrt(rjk2);
    flt_t rik2 = r21*r21;
    flt_t dw21, w21 = Sp<flt_t>(r21,rcmin[itype][ktype],rcmax[itype][ktype],
      &dw21);

    flt_t rij = r32;
    flt_t rik = r21;
    flt_t rij2 = r32*r32;
    flt_t costmp = static_cast<flt_t>(0.5)*(rij2+rik2-rjk2)/rij/rik;
    flt_t dtsjik, tspjik = Sp2<flt_t>(costmp,thmin,thmax,&dtsjik);
    dtsjik = -dtsjik;

    int * REBO_neighs_j = &ka->neigh_rebo.entries[ka->neigh_rebo.offset[j]];
    int lnum = ka->neigh_rebo.num[j];
    for (int ll = 0; ll < lnum; ll++) {
      int l = REBO_neighs_j[ll];
      int ltype = map[x[l].w];
      if (l == i || l == k) continue;
      flt_t del34x = x[j].x-x[l].x;
      flt_t del34y = x[j].y-x[l].y;
      flt_t del34z = x[j].z-x[l].z;
      flt_t rsq = del34x*del34x + del34y*del34y + del34z*del34z;
      flt_t r34 = overloaded::sqrt(rsq);
      flt_t cos234 = (del32x*del34x + del32y*del34y +
                del32z*del34z) / (r32*r34);
      cos234 = fmin(cos234,1);
      cos234 = fmax(cos234,-1);
      flt_t sin234 = overloaded::sqrt(1 - cos234*cos234);
      if (sin234 < TOL) continue;
      flt_t dw34, w34 = Sp<flt_t>(r34,rcmin[jtype][ltype],rcmax[jtype][ltype],
        &dw34);
      flt_t delilx = del23x + del34x;
      flt_t delily = del23y + del34y;
      flt_t delilz = del23z + del34z;
      flt_t ril2 = delilx*delilx + delily*delily + delilz*delilz;
      flt_t ril = overloaded::sqrt(ril2);
      flt_t rjl2 = r34*r34;

      flt_t rjl = r34;
      flt_t costmp = static_cast<flt_t>(0.5)*(rij2+rjl2-ril2)/rij/rjl;
      flt_t dtsijl, tspijl = Sp2<flt_t>(costmp,thmin,thmax,&dtsijl);
      dtsijl = -dtsijl; //need minus sign
      flt_t cross321x = (del32y*del21z)-(del32z*del21y);
      flt_t cross321y = (del32z*del21x)-(del32x*del21z);
      flt_t cross321z = (del32x*del21y)-(del32y*del21x);
      flt_t cross321mag = overloaded::sqrt(cross321x*cross321x+
                         cross321y*cross321y + cross321z*cross321z);
      flt_t cross234x = (del23y*del34z)-(del23z*del34y);
      flt_t cross234y = (del23z*del34x)-(del23x*del34z);
      flt_t cross234z = (del23x*del34y)-(del23y*del34x);
      flt_t cross234mag = overloaded::sqrt(cross234x*cross234x+
                         cross234y*cross234y + cross234z*cross234z);
      flt_t cwnum = (cross321x*cross234x) +
        (cross321y*cross234y)+(cross321z*cross234z);
      flt_t cwnom = r21*r34*r32*r32*sin321*sin234;
      flt_t cw = cwnum/cwnom;

      flt_t cw2 = (static_cast<flt_t>(.5)*(1-cw));
      flt_t ekijl = epsilonT[ktype][ltype];
      flt_t Ec = 256*ekijl/405;
      flt_t Vtors = (Ec*(overloaded::pow(cw2,5)))-(ekijl/10);

      ka->result_eng += Vtors*w21*w23*w34*(1-tspjik)*(1-tspijl);

      flt_t dndijx = (cross234y*del21z)-(cross234z*del21y);
      flt_t dndijy = (cross234z*del21x)-(cross234x*del21z);
      flt_t dndijz = (cross234x*del21y)-(cross234y*del21x);

      flt_t tmpvecx = (del34y*cross321z)-(del34z*cross321y);
      flt_t tmpvecy = (del34z*cross321x)-(del34x*cross321z);
      flt_t tmpvecz = (del34x*cross321y)-(del34y*cross321x);

      dndijx = dndijx+tmpvecx;
      dndijy = dndijy+tmpvecy;
      dndijz = dndijz+tmpvecz;

      flt_t dndikx = (del23y*cross234z)-(del23z*cross234y);
      flt_t dndiky = (del23z*cross234x)-(del23x*cross234z);
      flt_t dndikz = (del23x*cross234y)-(del23y*cross234x);

      flt_t dndjlx = (cross321y*del23z)-(cross321z*del23y);
      flt_t dndjly = (cross321z*del23x)-(cross321x*del23z);
      flt_t dndjlz = (cross321x*del23y)-(cross321y*del23x);

      flt_t dcidij = ((r23*r23)-(r21*r21)+(rjk*rjk))/(2*r23*r23*r21);
      flt_t dcidik = ((r21*r21)-(r23*r23)+(rjk*rjk))/(2*r23*r21*r21);
      flt_t dcidjk = (-rjk)/(r23*r21);
      flt_t dcjdji = ((r23*r23)-(r34*r34)+(ril*ril))/(2*r23*r23*r34);
      flt_t dcjdjl = ((r34*r34)-(r23*r23)+(ril*ril))/(2*r23*r34*r34);
      flt_t dcjdil = (-ril)/(r23*r34);

      flt_t dsidij = (-cos321/sin321)*dcidij;
      flt_t dsidik = (-cos321/sin321)*dcidik;
      flt_t dsidjk = (-cos321/sin321)*dcidjk;

      flt_t dsjdji = (-cos234/sin234)*dcjdji;
      flt_t dsjdjl = (-cos234/sin234)*dcjdjl;
      flt_t dsjdil = (-cos234/sin234)*dcjdil;

      flt_t dxidij = (r21*sin321)+(r23*r21*dsidij);
      flt_t dxidik = (r23*sin321)+(r23*r21*dsidik);
      flt_t dxidjk = (r23*r21*dsidjk);

      flt_t dxjdji = (r34*sin234)+(r23*r34*dsjdji);
      flt_t dxjdjl = (r23*sin234)+(r23*r34*dsjdjl);
      flt_t dxjdil = (r23*r34*dsjdil);

      flt_t ddndij = (dxidij*cross234mag)+(cross321mag*dxjdji);
      flt_t ddndik = dxidik*cross234mag;
      flt_t ddndjk = dxidjk*cross234mag;
      flt_t ddndjl = cross321mag*dxjdjl;
      flt_t ddndil = cross321mag*dxjdil;
      flt_t dcwddn = -cwnum/(cwnom*cwnom);
      flt_t dcwdn = 1/cwnom;
      flt_t dvpdcw = (-1)*Ec*static_cast<flt_t>(-0.5)*5*overloaded::pow(cw2,4)*
                      w23*w21*w34*(1-tspjik)*(1-tspijl);

      flt_t Ftmpx = dvpdcw*((dcwdn*dndijx)+(dcwddn*ddndij*del23x/r23));
      flt_t Ftmpy = dvpdcw*((dcwdn*dndijy)+(dcwddn*ddndij*del23y/r23));
      flt_t Ftmpz = dvpdcw*((dcwdn*dndijz)+(dcwddn*ddndij*del23z/r23));
      flt_t fix = Ftmpx;
      flt_t fiy = Ftmpy;
      flt_t fiz = Ftmpz;
      flt_t fjx = -Ftmpx;
      flt_t fjy = -Ftmpy;
      flt_t fjz = -Ftmpz;

      Ftmpx = dvpdcw*((dcwdn*dndikx)+(dcwddn*ddndik*del21x/r21));
      Ftmpy = dvpdcw*((dcwdn*dndiky)+(dcwddn*ddndik*del21y/r21));
      Ftmpz = dvpdcw*((dcwdn*dndikz)+(dcwddn*ddndik*del21z/r21));
      fix += Ftmpx;
      fiy += Ftmpy;
      fiz += Ftmpz;
      flt_t fkx = -Ftmpx;
      flt_t fky = -Ftmpy;
      flt_t fkz = -Ftmpz;

      Ftmpx = (dvpdcw*dcwddn*ddndjk*deljkx)/rjk;
      Ftmpy = (dvpdcw*dcwddn*ddndjk*deljky)/rjk;
      Ftmpz = (dvpdcw*dcwddn*ddndjk*deljkz)/rjk;
      fjx += Ftmpx;
      fjy += Ftmpy;
      fjz += Ftmpz;
      fkx -= Ftmpx;
      fky -= Ftmpy;
      fkz -= Ftmpz;

      Ftmpx = dvpdcw*((dcwdn*dndjlx)+(dcwddn*ddndjl*del34x/r34));
      Ftmpy = dvpdcw*((dcwdn*dndjly)+(dcwddn*ddndjl*del34y/r34));
      Ftmpz = dvpdcw*((dcwdn*dndjlz)+(dcwddn*ddndjl*del34z/r34));
      fjx += Ftmpx;
      fjy += Ftmpy;
      fjz += Ftmpz;
      flt_t flx = -Ftmpx;
      flt_t fly = -Ftmpy;
      flt_t flz = -Ftmpz;

      Ftmpx = (dvpdcw*dcwddn*ddndil*delilx)/ril;
      Ftmpy = (dvpdcw*dcwddn*ddndil*delily)/ril;
      Ftmpz = (dvpdcw*dcwddn*ddndil*delilz)/ril;
      fix += Ftmpx;
      fiy += Ftmpy;
      fiz += Ftmpz;
      flx -= Ftmpx;
      fly -= Ftmpy;
      flz -= Ftmpz;

      // coordination forces

      flt_t fpair = Vtors*dw21*w23*w34*(1-tspjik)*(1-tspijl) / r21;
      fix -= del21x*fpair;
      fiy -= del21y*fpair;
      fiz -= del21z*fpair;
      fkx += del21x*fpair;
      fky += del21y*fpair;
      fkz += del21z*fpair;

      fpair = Vtors*w21*dw23*w34*(1-tspjik)*(1-tspijl) / r23;
      fix -= del23x*fpair;
      fiy -= del23y*fpair;
      fiz -= del23z*fpair;
      fjx += del23x*fpair;
      fjy += del23y*fpair;
      fjz += del23z*fpair;

      fpair = Vtors*w21*w23*dw34*(1-tspjik)*(1-tspijl) / r34;
      fjx -= del34x*fpair;
      fjy -= del34y*fpair;
      fjz -= del34z*fpair;
      flx += del34x*fpair;
      fly += del34y*fpair;
      flz += del34z*fpair;

      // additional cut off function forces

      flt_t fcpc = -Vtors*w21*w23*w34*dtsjik*(1-tspijl);
      fpair = fcpc*dcidij/rij;
      fix += fpair*del23x;
      fiy += fpair*del23y;
      fiz += fpair*del23z;
      fjx -= fpair*del23x;
      fjy -= fpair*del23y;
      fjz -= fpair*del23z;

      fpair = fcpc*dcidik/rik;
      fix += fpair*del21x;
      fiy += fpair*del21y;
      fiz += fpair*del21z;
      fkx -= fpair*del21x;
      fky -= fpair*del21y;
      fkz -= fpair*del21z;

      fpair = fcpc*dcidjk/rjk;
      fjx += fpair*deljkx;
      fjy += fpair*deljky;
      fjz += fpair*deljkz;
      fkx -= fpair*deljkx;
      fky -= fpair*deljky;
      fkz -= fpair*deljkz;

      fcpc = -Vtors*w21*w23*w34*(1-tspjik)*dtsijl;
      fpair = fcpc*dcjdji/rij;
      fix += fpair*del23x;
      fiy += fpair*del23y;
      fiz += fpair*del23z;
      fjx -= fpair*del23x;
      fjy -= fpair*del23y;
      fjz -= fpair*del23z;

      fpair = fcpc*dcjdjl/rjl;
      fjx += fpair*del34x;
      fjy += fpair*del34y;
      fjz += fpair*del34z;
      flx -= fpair*del34x;
      fly -= fpair*del34y;
      flz -= fpair*del34z;

      fpair = fcpc*dcjdil/ril;
      fix += fpair*delilx;
      fiy += fpair*delily;
      fiz += fpair*delilz;
      flx -= fpair*delilx;
      fly -= fpair*delily;
      flz -= fpair*delilz;

      // sum per-atom forces into atom force array

      f[i].x += fix; f[i].y += fiy; f[i].z += fiz;
      f[j].x += fjx; f[j].y += fjy; f[j].z += fjz;
      f[k].x += fkx; f[k].y += fky; f[k].z += fkz;
      f[l].x += flx; f[l].y += fly; f[l].z += flz;
    }
  }
}

template<typename flt_t, typename acc_t>
void ref_torsion(KernelArgsAIREBOT<flt_t,acc_t> * ka) {
  AtomAIREBOT<flt_t> * x = ka->x;
  int * map = ka->map;
  tagint * tag = ka->tag;
  for (int ii = ka->frebo_from_atom; ii < ka->frebo_to_atom; ii++) {
    int i = ii;
    tagint itag = tag[i];
    int itype = map[x[i].w];
    if (itype != 0) continue;
    flt_t xtmp = x[i].x;
    flt_t ytmp = x[i].y;
    flt_t ztmp = x[i].z;
    int * REBO_neighs_i = &ka->neigh_rebo.entries[ka->neigh_rebo.offset[i]];
    int jnum = ka->neigh_rebo.num[i];
    for (int jj = 0; jj < jnum; jj++) {
      int j = REBO_neighs_i[jj];
      tagint jtag = tag[j];

      if (itag > jtag) {
        if (((itag+jtag) & 1) == 0) continue;
      } else if (itag < jtag) {
        if (((itag+jtag) & 1) == 1) continue;
      } else {
        if (x[j].z < ztmp) continue;
        if (x[j].z == ztmp && x[j].y < ytmp) continue;
        if (x[j].z == ztmp && x[j].y == ytmp && x[j].x < xtmp) continue;
      }

      int jtype = map[x[j].w];
      if (jtype != 0) continue;
      ref_torsion_single_interaction(ka, i, j);
    }
  }
}

/*
 * Calculate single REBO interaction.
 * Corresponds to FREBO method. Note that the bondorder() function is
 * inlined.
 */
template<typename flt_t, typename acc_t>
void ref_frebo_single_interaction(KernelArgsAIREBOT<flt_t,acc_t> * ka, int i,
    int j) {
  AtomAIREBOT<flt_t> * x = ka->x;
  int * map = ka->map;
  ResultForceT<acc_t> * result_f = ka->result_f;
  int jj;
  int itype = map[x[i].w];
  flt_t x_i = x[i].x;
  flt_t y_i = x[i].y;
  flt_t z_i = x[i].z;
  int jtype = map[x[j].w];
  flt_t delx = x[i].x - x[j].x;
  flt_t dely = x[i].y - x[j].y;
  flt_t delz = x[i].z - x[j].z;
  flt_t rsq = delx * delx + dely * dely + delz * delz;
  flt_t rij = overloaded::sqrt(rsq);
  flt_t rcminij = ka->params.rcmin[itype][jtype];
  flt_t rcmaxij = ka->params.rcmax[itype][jtype];
  flt_t dwij;
  flt_t wij = Sp(rij, rcminij, rcmaxij, &dwij);
  if (wij <= TOL) return;

  flt_t Qij = ka->params.Q[itype][jtype];
  flt_t Aij = ka->params.A[itype][jtype];
  flt_t alphaij = ka->params.alpha[itype][jtype];

  flt_t exp_alphar = exp(-alphaij * rij);
  flt_t VR_by_wij = (1.0 + (Qij / rij)) * Aij * exp_alphar;
  flt_t VR = wij * VR_by_wij;
  flt_t pre = wij * Aij * exp_alphar;
  flt_t dVRdi = pre * ((-alphaij) - (Qij / rsq) - (Qij * alphaij / rij));
  dVRdi += VR_by_wij * dwij;

  flt_t VA_by_wij = 0, dVA = 0;
  for (int k = 0; k < 3; k++) {
    flt_t BIJc = ka->params.BIJc[itype][jtype][k];
    flt_t Betaij = ka->params.Beta[itype][jtype][k];
    flt_t term = -BIJc * overloaded::exp(-Betaij * rij);
    VA_by_wij += term;
    dVA += -Betaij * wij * term;
  }
  dVA += VA_by_wij * dwij;
  flt_t VA = VA_by_wij * wij;

  acc_t fij[3] = {0};
  flt_t Nij = ka->nH[i] + ka->nC[i] - wij;
  flt_t Nji = ka->nH[j] + ka->nC[j] - wij;
  flt_t NconjtmpI;
  flt_t pij = frebo_pij(ka, i, j, delx, dely, delz, rij, wij, VA, &NconjtmpI,
    fij);
  flt_t NconjtmpJ;
  acc_t fji[3] = {0};
  flt_t pji = frebo_pij(ka, j, i, -delx, -dely, -delz, rij, wij, VA,
    &NconjtmpJ, fji);
  fij[0] -= fji[0]; fij[1] -= fji[1]; fij[2] -= fji[2];
  flt_t Nijconj = 1.0 + (NconjtmpI * NconjtmpI) + (NconjtmpJ * NconjtmpJ);
  flt_t dN3[3];
  flt_t pi_rc = frebo_pi_rc(ka, itype, jtype, Nij, Nji, Nijconj, dN3);
  frebo_N_spline_force(ka, i, j, VA, dN3[0], dN3[2], NconjtmpI);
  frebo_N_spline_force(ka, j, i, VA, dN3[1], dN3[2], NconjtmpJ);
  flt_t Tij = frebo_Tij(ka, itype, jtype, Nij, Nji, Nijconj, dN3);
  flt_t sum_omega = 0.0;
  if (fabs(Tij) > TOL) {
    sum_omega = frebo_sum_omega(ka, i, j, delx, dely, delz, rij, VA * Tij, fij);
    frebo_N_spline_force(ka, i, j, VA * sum_omega, dN3[0], dN3[2], NconjtmpI);
    frebo_N_spline_force(ka, j, i, VA * sum_omega, dN3[1], dN3[2], NconjtmpJ);
  }
  flt_t pi_dh = Tij * sum_omega;
  flt_t bij = static_cast<flt_t>(0.5) * (pij + pji) + pi_rc + pi_dh;
  flt_t dVAdi = bij * dVA;
  flt_t fpair = -(dVRdi + dVAdi) / rij;

  result_f[i].x += fpair * delx + fij[0];
  result_f[i].y += fpair * dely + fij[1];
  result_f[i].z += fpair * delz + fij[2];
  result_f[j].x -= fpair * delx + fij[0];
  result_f[j].y -= fpair * dely + fij[1];
  result_f[j].z -= fpair * delz + fij[2];

  flt_t evdwl = VR + bij * VA;
  ka->result_eng += evdwl;
  result_f[i].w += 0.5 * evdwl;
  result_f[j].w += 0.5 * evdwl;
}


template<typename flt_t, typename acc_t>
inline void ref_frebo_single_atom(KernelArgsAIREBOT<flt_t,acc_t> * ka, int i) {
  AtomAIREBOT<flt_t> * x = ka->x;
  tagint * tag = ka->tag;
  int jj;
  tagint itag = tag[i];
  flt_t x_i = x[i].x;
  flt_t y_i = x[i].y;
  flt_t z_i = x[i].z;
  int * neighs = ka->neigh_rebo.entries + ka->neigh_rebo.offset[i];
  int jnum = ka->neigh_rebo.num[i];
  for (jj = 0; jj < jnum; jj++) {
    int j = neighs[jj];
    tagint jtag = tag[j];
    if (itag > jtag) {
      if (((itag + jtag) & 1) == 0)
        continue;
    } else if (itag < jtag) {
      if (((itag + jtag) & 1) == 1)
        continue;
    } else {
      if (x[j].z < z_i)
        continue;
      if (x[j].z == z_i && x[j].y < y_i)
        continue;
      if (x[j].z == z_i && x[j].y == y_i && x[j].x < x_i)
        continue;
    }
    ref_frebo_single_interaction(ka, i, j);
  }
}


template<typename flt_t, typename acc_t>
void ref_frebo(KernelArgsAIREBOT<flt_t,acc_t> * ka, int torflag) {
  for (int i = ka->frebo_from_atom; i < ka->frebo_to_atom; i++) {
    ref_frebo_single_atom(ka, i);
  }
  if (torflag) ref_torsion(ka);
}

template<typename flt_t, typename acc_t>
void ref_lennard_jones_single_interaction(KernelArgsAIREBOT<flt_t,acc_t> * ka,
    int i, int j, int morseflag) {
  AtomAIREBOT<flt_t> * x = ka->x;
  int * map = ka->map;
  ResultForceT<acc_t> * result_f = ka->result_f;

  int itype = map[x[i].w];
  int jtype = map[x[j].w];

  flt_t delx = x[i].x - x[j].x;
  flt_t dely = x[i].y - x[j].y;
  flt_t delz = x[i].z - x[j].z;
  flt_t rsq = delx * delx + dely * dely + delz * delz;

  if (rsq >= ka->params.cutljsq[itype][jtype]) { return; }
  flt_t rij = overloaded::sqrt(rsq);

  LennardJonesPathAIREBOT<flt_t> testpath;
  flt_t cij = 1.0;
  if (rij < ka->params.cut3rebo) {
    #pragma noinline
    cij = ref_lennard_jones_test_path<flt_t,acc_t>(ka, i, j, rij,
      ka->params.rcmax[itype][jtype], &testpath);
  }
  if (cij == 0) {
    return;
  }

  flt_t sigcut = ka->params.sigcut;
  flt_t sigmin = ka->params.sigmin;
  flt_t sigma = ka->params.sigma[itype][jtype];
  flt_t rljmax = sigcut * sigma;
  flt_t rljmin = sigmin * sigma;

  flt_t dslw, slw = Sp2(rij, rljmin, rljmax, &dslw);

  flt_t vdw, dvdw;
  if (morseflag) {
    const flt_t exr = exp(-rij * ka->params.lj4[itype][jtype]);
    vdw = ka->params.lj1[itype][jtype] * exr *
      (ka->params.lj2[itype][jtype]*exr - 2);
    dvdw = ka->params.lj3[itype][jtype] * exr *
      (1 - ka->params.lj2[itype][jtype]*exr);
  } else {
    flt_t r2inv = 1 / rsq;
    flt_t r6inv = r2inv * r2inv * r2inv;

    vdw = r6inv * (ka->params.lj3[itype][jtype]*r6inv -
                   ka->params.lj4[itype][jtype]);
    dvdw = -r6inv * (ka->params.lj1[itype][jtype]*r6inv -
                     ka->params.lj2[itype][jtype]) / rij;
  }

  flt_t VLJ = vdw * slw;
  flt_t dVLJ = dvdw * slw + vdw * dslw;

  flt_t dStr, Str = Sp2<flt_t>(rij, ka->params.rcLJmin[itype][jtype],
    ka->params.rcLJmax[itype][jtype], &dStr);
  flt_t VA = Str * cij * VLJ;
  flt_t Stb = 0;
  acc_t fij[3] = {0};
  if (Str > 0) {
    #pragma noinline
    Stb = ref_lennard_jones_bondorder(ka, i, j, VA, fij);
  }
  flt_t fpair = -(dStr * (Stb * cij * VLJ - cij * VLJ) +
                   dVLJ * (Str * Stb * cij + cij - Str * cij)) / rij;
  flt_t evdwl = VA * Stb + (1 - Str) * cij * VLJ;
  result_f[i].x += fpair * delx + fij[0];
  result_f[i].y += fpair * dely + fij[1];
  result_f[i].z += fpair * delz + fij[2];
  result_f[j].x -= fpair * delx + fij[0];
  result_f[j].y -= fpair * dely + fij[1];
  result_f[j].z -= fpair * delz + fij[2];
  ka->result_eng += evdwl;

  if (cij < 1) {
    #pragma noinline
    ref_lennard_jones_force_path(ka, Str * Stb * VLJ + (1 - Str) * VLJ,
      &testpath);
  }
}

template<typename flt_t, typename acc_t>
void ref_lennard_jones_single_atom(KernelArgsAIREBOT<flt_t,acc_t> * ka, int i,
                                   int morseflag) {
  AtomAIREBOT<flt_t> * x = ka->x;
  int jj;
  int * neighs = ka->neigh_lmp.entries[i];
  int jnum = ka->neigh_lmp.num_half[i];
  for (jj = 0; jj < jnum; jj++) {
    int j = neighs[jj];
    ref_lennard_jones_single_interaction(ka, i, j, morseflag);
  }
}

template<typename flt_t, typename acc_t>
void ref_lennard_jones(KernelArgsAIREBOT<flt_t,acc_t> * ka, int morseflag) {
  for (int i = ka->frebo_from_atom; i < ka->frebo_to_atom; i++) {
    #pragma noinline
    ref_lennard_jones_single_atom(ka, i, morseflag);
  }
}

/* ----------------------------------------------------------------------
    Vectorized AIREBO implementation, standalone, using caching to reduce
    memory access.
   ---------------------------------------------------------------------- */

template<typename flt_t, typename acc_t>
struct aut_wrap {

typedef typename intr_types<flt_t, acc_t>::fvec fvec;
typedef typename intr_types<flt_t, acc_t>::avec avec;
typedef typename intr_types<flt_t, acc_t>::ivec ivec;
typedef typename intr_types<flt_t, acc_t>::bvec bvec;

VEC_INLINE inline
static void aut_loadatoms_vec(
    AtomAIREBOT<flt_t> * atoms, ivec j_vec,
    fvec *x, fvec * y, fvec * z, bvec * type_mask, int * /*map*/, ivec map_i,
    ivec c_1
) {
  const ivec c_4 = ivec::set1(4);
  ivec j_vec_4 = ivec::mullo(c_4, j_vec);
  fvec w;
  fvec::gather_4_adjacent(j_vec_4, &atoms[0].x, sizeof(flt_t), x, y, z, &w);
  ivec jtype = fvec::unpackloepi32(w);
  jtype = ivec::srlv(map_i, jtype); //_mm512_castpd_si512(w));
  jtype = ivec::the_and(c_1, jtype);
  bvec jtype_mask = ivec::cmpneq(jtype, ivec::setzero());
  *type_mask = jtype_mask;
}

VEC_INLINE inline
static void aut_loadatoms_vec_notype(
    AtomAIREBOT<flt_t> * atoms, ivec j_vec,
    fvec *x, fvec * y, fvec * z
) {
  const ivec c_4 = ivec::set1(4);
  ivec j_vec_4 = ivec::mullo(c_4, j_vec);
  fvec::gather_3_adjacent(j_vec_4, &atoms[0].x, sizeof(flt_t), x, y, z);
}

static fvec aut_Sp2_deriv(fvec r, fvec lo, fvec hi, fvec * d) {
  fvec c_1 = fvec::set1(1);
  fvec c_2 = fvec::set1(2);
  fvec c_3 = fvec::set1(3);
  fvec c_6 = fvec::set1(6);
  bvec m_lo = fvec::cmple(r, lo);
  bvec m_hi = fvec::cmpnlt(r, hi); // nlt == ge
  bvec m_tr = bvec::kandn(m_lo, ~ m_hi);
  fvec ret = c_1;
  ret = fvec::mask_blend(m_hi, ret, fvec::setzero());
  fvec der = fvec::setzero();
  if (bvec::test_any_set(m_tr)) {
    fvec diff = hi -  lo;
    fvec rcp = fvec::recip(diff);
    fvec t = (r -  lo) *  rcp;
    fvec v = c_1 -  t *  t * ( c_3 -  c_2 *  t);
    ret = fvec::mask_blend(m_tr, ret, v);
    fvec dv = c_6 *  rcp * ( t *  t -  t);
    der = fvec::mask_blend(m_tr, der, dv);
  }
  *d = der;
  return ret;
}

static fvec aut_Sp_deriv(fvec r, fvec lo, fvec hi, fvec * d) {
  fvec c_1 = fvec::set1(1);
  fvec c_0_5 = fvec::set1(0.5);
  fvec c_m0_5 = fvec::set1(-0.5);
  fvec c_PI = fvec::set1(MY_PI);
  bvec m_lo = fvec::cmple(r, lo);
  bvec m_hi = fvec::cmpnlt(r, hi); // nlt == ge
  bvec m_tr = bvec::kandn(m_lo, ~ m_hi);
  fvec ret = c_1;
  ret = fvec::mask_blend(m_hi, ret, fvec::setzero());
  fvec der = fvec::setzero();
  if (bvec::test_any_set(m_tr)) {
    fvec diff = hi -  lo;
    fvec rcp = fvec::mask_recip(c_1, m_tr, diff);
    fvec t = (r -  lo) /  diff;
    fvec sinval, cosval;
    sinval = fvec::mask_sincos(&cosval, fvec::setzero(), c_1, m_tr, c_PI *  t);
    fvec v = c_0_5 * ( c_1 +  cosval);
    ret = fvec::mask_blend(m_tr, ret, v);
    fvec dv = c_PI *  c_m0_5 *  rcp *  sinval;
    der = fvec::mask_blend(m_tr, der, dv);
  }
  *d = der;
  return ret;
}

static fvec aut_mask_Sp(bvec mask, fvec r, fvec lo, fvec hi) {
  fvec c_1 = fvec::set1(1);
  fvec c_0_5 = fvec::set1(0.5);
  fvec c_PI = fvec::set1(MY_PI);
  bvec m_lo = fvec::mask_cmple(mask, r, lo);
  bvec m_hi = fvec::mask_cmpnlt(mask, r, hi); // nlt == ge
  bvec m_tr = bvec::kandn(m_lo, bvec::kandn(m_hi, mask));
  fvec ret = c_1;
  ret = fvec::mask_blend(m_hi, ret, fvec::setzero());
  if (bvec::test_any_set(m_tr)) {
    fvec rcp = fvec::mask_recip(c_1, m_tr, hi -  lo);
    fvec t = (r -  lo) *  rcp;
    fvec v = c_0_5 * ( c_1 +  fvec::mask_cos(c_1, m_tr, c_PI *  t));
    ret = fvec::mask_blend(m_tr, ret, v);
  }
  return ret;
}

static void aut_rebo_neigh(KernelArgsAIREBOT<flt_t,acc_t> * ka) {
  int offset = ka->neigh_from_atom * ka->num_neighs_per_atom;
  ivec c_CARBON = ivec::setzero();
  int map_i = 0;
  int i;
  for (i = 1; i < ka->num_types; i++) {
    if (ka->map[i])
      map_i |= (1 << i);
  }
  ivec c_i1 = ivec::set1(1);
  ivec c_im = ivec::set1(map_i);
  AtomAIREBOT<flt_t> * _noalias x = ka->x;

  for (i = ka->neigh_from_atom; i < ka->neigh_to_atom; i++) {

    fvec x_i = fvec::set1(x[i].x);
    fvec y_i = fvec::set1(x[i].y);
    fvec z_i = fvec::set1(x[i].z);
    int itype = ka->map[ka->x[i].w];

    fvec rcmaxsq0 = fvec::set1(ka->params.rcmaxsq[itype][0]);
    fvec rcmaxsq1 = fvec::set1(ka->params.rcmaxsq[itype][1]);
    fvec rcmax0 = fvec::set1(ka->params.rcmax[itype][0]);
    fvec rcmax1 = fvec::set1(ka->params.rcmax[itype][1]);
    fvec rcmin0 = fvec::set1(ka->params.rcmin[itype][0]);
    fvec rcmin1 = fvec::set1(ka->params.rcmin[itype][1]);
    fvec rcmaxskinsq0 = fvec::set1(
        (ka->params.rcmax[itype][0] + ka->skin) * (ka->params.rcmax[itype][0] +
                                                   ka->skin));
    fvec rcmaxskinsq1 = fvec::set1(
        (ka->params.rcmax[itype][1] + ka->skin) * (ka->params.rcmax[itype][1] +
                                                   ka->skin));
    fvec nC = fvec::setzero();
    fvec nH = fvec::setzero();

    ka->neigh_rebo.offset[i] = offset;

    int jnum = ka->rebuild_flag ? ka->neigh_lmp.num[i] :
      ka->neigh_rebo.num_half[i];
    int * neighs = ka->rebuild_flag ?
      ka->neigh_lmp.entries[i] :
      &ka->neigh_rebo.entries[ka->neigh_rebo.offset[i]+jnum];
    int * skin_target = &ka->neigh_rebo.entries[offset+ka->num_neighs_per_atom];
    int n = 0;
    int n_skin = 0;

    int lowest_idx;
    //#pragma unroll(4)
    for (lowest_idx = 0; lowest_idx < jnum; lowest_idx += fvec::VL) {
      bvec j_mask = bvec::full();
      if (lowest_idx + fvec::VL > jnum) j_mask = bvec::only(jnum - lowest_idx);

      int * _noalias neighs_l = neighs + lowest_idx;
      fvec x_j, y_j, z_j;
      bvec jtype_mask;
      ivec ji = ivec::maskz_loadu(j_mask, neighs_l);
      aut_loadatoms_vec(x, ji,
          &x_j, &y_j, &z_j, &jtype_mask, ka->map, c_im, c_i1);
      fvec delx = x_i -  x_j;
      fvec dely = y_i -  y_j;
      fvec delz = z_i -  z_j;
      fvec rsq = delx *  delx +  dely *  dely +  delz *  delz;
      if (ka->rebuild_flag) {
        fvec rcmaxskinsq = fvec::mask_blend(jtype_mask, rcmaxskinsq0,
                                            rcmaxskinsq1);
        bvec c_mask = fvec::mask_cmplt(j_mask, rsq, rcmaxskinsq);
        ivec::mask_compressstore(c_mask, &skin_target[n_skin], ji);
        n_skin += bvec::popcnt(c_mask);
      }
      fvec rcmaxsq = fvec::mask_blend(jtype_mask, rcmaxsq0, rcmaxsq1);
      bvec c_mask = fvec::mask_cmplt(j_mask, rsq, rcmaxsq);
      if (bvec::test_all_unset(c_mask)) continue;
      ivec::mask_compressstore(c_mask, &ka->neigh_rebo.entries[offset + n], ji);
      n += bvec::popcnt(c_mask);
      fvec rcmax = fvec::mask_blend(jtype_mask, rcmax0, rcmax1);
      fvec rcmin = fvec::mask_blend(jtype_mask, rcmin0, rcmin1);
      fvec sp = aut_mask_Sp(c_mask, fvec::sqrt(rsq), rcmin, rcmax);
      nC = fvec::mask_add(nC, bvec::kandn(jtype_mask, c_mask), nC, sp);
      nH = fvec::mask_add(nH, bvec::kand (jtype_mask, c_mask), nH, sp);
    }
    ka->neigh_rebo.num[i] = n;
    if (ka->rebuild_flag) {
      for (int i = 0; i < n_skin; i++) {
        ka->neigh_rebo.entries[offset+n_skin+i] = skin_target[i];
      }
    }
    if (ka->rebuild_flag) {
      assert(n <= n_skin);
      offset += 2 * n_skin;
      ka->neigh_rebo.num_half[i] = n_skin;
    } else {
      assert(n <= jnum);
      offset += 2 * jnum;
    }
    ka->nC[i] = fvec::reduce_add(nC);
    ka->nH[i] = fvec::reduce_add(nH);
  }
}


static fvec aut_eval_poly_lin_pd_2(int n, flt_t * vals, ivec idx, fvec x,
                                   fvec * deriv) {
  fvec c_1 = fvec::set1(1);
  fvec x_i = c_1;
  fvec x_im1 = fvec::setzero();
  fvec result = fvec::setzero();
  fvec i_v = fvec::setzero();
  *deriv = fvec::setzero();
  int i;
  for (i = 0; i < n; i++) {
    fvec coeff = fvec::gather(idx, vals + i, sizeof(flt_t));
    result = result +  coeff *  x_i;
    *deriv = *deriv +  coeff *  x_im1 *  i_v;
    x_im1 = x_i;
    x_i = x_i *  x;
    i_v = i_v +  c_1;
  }
  return result;
}

static fvec aut_mask_gSpline_pd_2(KernelArgsAIREBOT<flt_t,acc_t> * ka,
                                  bvec /*active_mask*/, int itype, fvec cosjik,
                                  fvec Nij, fvec *dgdc, fvec *dgdN) {
  int i;
  flt_t * gDom = nullptr;
  int nDom = 0;
  ivec offs = ivec::setzero();
  fvec NCmin = fvec::set1(ka->params.NCmin);
  bvec Ngt = fvec::cmpnle(Nij, NCmin); //gt
  if (itype == 0) {
    nDom = 4;
    gDom = &ka->params.gCdom[0];
    offs = ivec::mask_blend(Ngt, offs, ivec::set1(4*6));
  } else {
    nDom = 3;
    gDom = &ka->params.gHdom[0];
    offs = ivec::set1(8 * 6);
  }
  cosjik = fvec::max(fvec::set1(gDom[0]), fvec::min(fvec::set1(gDom[nDom]),
                                                    cosjik));
  ivec index6 = ivec::setzero();
  for (i = 0; i < nDom; i++) {
    bvec cosge = fvec::cmpnlt(cosjik, fvec::set1(gDom[i])); //ge
    bvec cosle = fvec::cmple(cosjik, fvec::set1(gDom[i+1]));
    index6 = ivec::mask_blend(cosge & cosle, index6, ivec::set1(6*i));
  }
  fvec g = aut_eval_poly_lin_pd_2(6, &ka->params.gVal[0], offs +  index6,
                                  cosjik, dgdc);
  *dgdN = fvec::setzero();
  if (itype == 0) {
    fvec NCmax = fvec::set1(ka->params.NCmax);
    bvec Nlt = fvec::cmplt(Nij, NCmax); //gt
    bvec Nmask = Ngt & Nlt;
    if (bvec::test_any_set(Nmask)) {
      fvec dg1;
      fvec g1 = aut_eval_poly_lin_pd_2(6, &ka->params.gVal[0], index6, cosjik,
                                       &dg1);
      fvec dS;
      fvec cut = aut_Sp_deriv(Nij, NCmin, NCmax, &dS);
      *dgdN = fvec::mask_mul(*dgdN, Nmask, dS, g1 -  g);
      g = fvec::mask_add(g, Nmask, g, cut * ( g1 -  g));
      *dgdc = fvec::mask_add(*dgdc, Nmask, *dgdc, cut * ( dg1 -  *dgdc));
    }
  }
  return g;
}

static fvec aut_PijSpline(KernelArgsAIREBOT<flt_t,acc_t> * ka, int itype,
                          int jtype, fvec NijC, fvec NijH, fvec *dN2) {
  flt_t ret[fvec::VL] __attribute__((aligned(64)));
  flt_t dN20[fvec::VL] __attribute__((aligned(64)));
  flt_t dN21[fvec::VL] __attribute__((aligned(64)));
  flt_t NijC_[fvec::VL] __attribute__((aligned(64)));
  flt_t NijH_[fvec::VL] __attribute__((aligned(64)));
  flt_t tmp_dN2[2];
  fvec::store(NijC_, NijC);
  fvec::store(NijH_, NijH);
  int i;
  for (i = 0; i < fvec::VL; i++) {
    ret[i] = PijSpline(ka, itype, jtype, NijC_[i], NijH_[i], tmp_dN2);
    dN20[i] = tmp_dN2[0];
    dN21[i] = tmp_dN2[1];
  }
  dN2[0] = fvec::load(dN20);
  dN2[1] = fvec::load(dN21);
  return fvec::load(ret);
}

/*
 * aut_frebo_data stores all the short-ranged coordinations
 * and intermediate values that get reused frequently during
 * bondorder calculations.
 * BUF_CAP should rarely exceed 4, so 8 is a very conservative
 * value.
 */
static const int BUF_CAP = 8;
struct aut_frebo_data {
  fvec rikx_buf[BUF_CAP];
  fvec riky_buf[BUF_CAP];
  fvec rikz_buf[BUF_CAP];
  fvec rikmag_buf[BUF_CAP];
  fvec cosjik_buf[BUF_CAP];
  ivec k_buf[BUF_CAP];
  fvec g_buf[BUF_CAP];
  fvec dgdc_buf[BUF_CAP];
  fvec ex_lam_buf[BUF_CAP];
  fvec wik_buf[BUF_CAP];
  fvec dwik_buf[BUF_CAP];
  fvec cutN_buf[BUF_CAP];
  fvec dcutN_buf[BUF_CAP];
  bvec ktype_buf[BUF_CAP];
  bvec mask_buf[BUF_CAP];
  fvec force_k_x_buf[BUF_CAP];
  fvec force_k_y_buf[BUF_CAP];
  fvec force_k_z_buf[BUF_CAP];
  int buf_len;
  fvec x_i;
  fvec y_i;
  fvec z_i;
  fvec x_j;
  fvec y_j;
  fvec z_j;
  fvec nCi;
  fvec nHi;
  fvec force_i_x;
  fvec force_i_y;
  fvec force_i_z;
  fvec force_j_x;
  fvec force_j_y;
  fvec force_j_z;
};

/*
 * Initialize values in aut_frebo_data and perform the calculations
 * for p_ij.
 */
static fvec aut_frebo_pij_pd_2(
    KernelArgsAIREBOT<flt_t,acc_t> * _noalias ka,
    struct aut_frebo_data * _noalias data,
    int itype, int jtype,
    ivec vi, ivec vj,
    fvec rijx, fvec rijy, fvec rijz, fvec rijmag,
    fvec wij, fvec VA, fvec * sum_N, fvec fij[3]
) {
  AtomAIREBOT<flt_t> * _noalias x = ka->x;
  int * _noalias map = ka->map;
  flt_t * _noalias nC = ka->nC;
  flt_t * _noalias nH = ka->nH;
  fvec x_i, y_i, z_i;
  fvec x_j, y_j, z_j;
  x_i = data->x_i;
  y_i = data->y_i;
  z_i = data->z_i;
  x_j = data->x_j;
  y_j = data->y_j;
  z_j = data->z_j;
  fvec invrijm = fvec::recip(rijmag);
  fvec invrijm2 = invrijm *  invrijm;
  fvec rcminij = fvec::set1(ka->params.rcmin[itype][jtype]);
  fvec rcmaxij = fvec::set1(ka->params.rcmax[itype][jtype]);
  fvec Nmin = fvec::set1(ka->params.Nmin);
  fvec Nmax = fvec::set1(ka->params.Nmax);
  int map_i_scalar = 0;
  {
    int i;
    for (i = 1; i < ka->num_types; i++) {
      if (ka->map[i])
        map_i_scalar |= (1 << i);
    }
  }
  ivec map_i = ivec::set1(map_i_scalar);
  fvec nCi = data->nCi;
  fvec nHi = data->nHi;
  fvec Nij = nHi +  nCi -  wij;
  fvec factor_jtype, factor_not_jtype;
  if (jtype) {
    factor_jtype = fvec::set1(1);
    factor_not_jtype = fvec::set1(0);
  } else {
    factor_jtype = fvec::set1(0);
    factor_not_jtype = fvec::set1(1);
  }
  fvec NijC = nCi -  wij *  factor_not_jtype;
  fvec NijH = nHi -  wij *  factor_jtype;
  fvec sum_pij = fvec::setzero();
  fvec sum_dpij_dN = fvec::setzero();
  fvec dN2[2];
  ivec offseti = ivec::mask_gather(ivec::setzero(), bvec::full(), vi,
                                   ka->neigh_rebo.offset, sizeof(int));
  int buf_len = 0;
  ivec knum = ivec::mask_gather(ivec::setzero(), bvec::full(), vi,
                                ka->neigh_rebo.num, sizeof(int));
  ivec kk = ivec::setzero();
  bvec active_mask = ivec::cmplt(kk, knum);
  ivec c_i1 = ivec::set1(1);
  fvec rho_j = fvec::set1(ka->params.rho[jtype][1]);
  fvec rho_k0 = fvec::set1(ka->params.rho[0][1]);
  fvec rho_k1 = fvec::set1(ka->params.rho[1][1]);
  fvec c_4 = fvec::set1(4);
  fvec c_2_0 = fvec::set1(2.0);
  fvec c_m2_0 = fvec::set1(-2.0);
  fvec c_4_0 = fvec::set1(4.0);
  fvec c_0_5 = fvec::set1(0.5);
  fvec c_m0_5 = fvec::set1(-0.5);
  fvec c_1 = fvec::set1(1);
  fvec c_m1 = fvec::set1(-1);
  fvec factor_itype = itype ? c_1 : fvec::setzero();
  fvec rcmax0 = fvec::set1(ka->params.rcmax[itype][0]);
  fvec rcmax1 = fvec::set1(ka->params.rcmax[itype][1]);
  fvec rcmin0 = fvec::set1(ka->params.rcmin[itype][0]);
  fvec rcmin1 = fvec::set1(ka->params.rcmin[itype][1]);
  fvec result_f_i_x = fvec::setzero();
  fvec result_f_i_y = fvec::setzero();
  fvec result_f_i_z = fvec::setzero();
  fvec result_f_j_x = fvec::setzero();
  fvec result_f_j_y = fvec::setzero();
  fvec result_f_j_z = fvec::setzero();
  *sum_N = fvec::setzero();
  {
    while (bvec::test_any_set(active_mask)) {
      ivec k = ivec::mask_gather(ivec::setzero(), active_mask, kk +  offseti,
                                 ka->neigh_rebo.entries, sizeof(int));
      bvec excluded_mask = ivec::cmpeq(k, vj) & active_mask;
      if (bvec::test_any_set(excluded_mask)) {
        kk = ivec::mask_add(kk, excluded_mask, kk, c_i1);
        active_mask = ivec::cmplt(kk, knum);
        continue;
      }
      fvec x_k, y_k, z_k;
      bvec ktype_mask;
      aut_loadatoms_vec(x, k, &x_k, &y_k, &z_k, &ktype_mask, ka->map, map_i,
                        c_i1);
      fvec rikx = x_i -  x_k;
      fvec riky = y_i -  y_k;
      fvec rikz = z_i -  z_k;
      fvec rikmag = fvec::sqrt(rikx *  rikx +  riky *  riky +  rikz *  rikz);
      fvec rho_k = fvec::mask_blend(ktype_mask, rho_k0, rho_k1);
      fvec lamdajik = c_4 *  factor_itype * ( rho_k -  rikmag - ( rho_j -
                                                                  rijmag));
      fvec ex_lam = fvec::exp(lamdajik);
      fvec rcmax = fvec::mask_blend(ktype_mask, rcmax0, rcmax1);
      fvec rcmin = fvec::mask_blend(ktype_mask, rcmin0, rcmin1);
      fvec dwik;
      fvec wik = aut_Sp_deriv(rikmag, rcmin, rcmax, &dwik);
      fvec Nki = fvec::gather(k, nC, sizeof(flt_t)) +
        fvec::gather(k, nH, sizeof(flt_t)) -  wik;
      fvec cosjik = (rijx *  rikx +  rijy *  riky +  rijz *  rikz) /
        ( rijmag *  rikmag);
      cosjik = fvec::min(c_1, fvec::max(c_m1, cosjik));
      fvec dgdc, dgdN;
      fvec g = aut_mask_gSpline_pd_2(ka, active_mask, itype, cosjik, Nij,
                                     &dgdc, &dgdN);
      sum_pij = fvec::mask_add(sum_pij, active_mask, sum_pij, wik * g * ex_lam);
      sum_dpij_dN = fvec::mask_add(sum_dpij_dN, active_mask, sum_dpij_dN,
                                   wik * ex_lam * dgdN);
      fvec dcutN;
      fvec cutN = aut_Sp_deriv(Nki, Nmin, Nmax, &dcutN);
      *sum_N = fvec::mask_add(*sum_N, active_mask, *sum_N,
                              fvec::mask_blend(ktype_mask, c_1,
                                               fvec::setzero()) * wik * cutN);
      if (buf_len == BUF_CAP) goto exceed_buffer;
      data->rikx_buf[buf_len] = rikx;
      data->riky_buf[buf_len] = riky;
      data->rikz_buf[buf_len] = rikz;
      data->rikmag_buf[buf_len] = rikmag;
      data->cosjik_buf[buf_len] = cosjik;
      data->ktype_buf[buf_len] = ktype_mask;
      data->k_buf[buf_len] = k;
      data->g_buf[buf_len] = g;
      data->dgdc_buf[buf_len] = dgdc;
      data->ex_lam_buf[buf_len] = ex_lam;
      data->wik_buf[buf_len] = wik;
      data->dwik_buf[buf_len] = dwik;
      data->mask_buf[buf_len] = active_mask;
      data->cutN_buf[buf_len] = cutN;
      data->dcutN_buf[buf_len] = dcutN;
      buf_len += 1;
      kk = ivec::mask_add(kk, active_mask, kk, c_i1);
      active_mask = ivec::cmplt(kk, knum);
    }
    data->buf_len = buf_len;
    fvec PijS = aut_PijSpline(ka, itype, jtype, NijC, NijH, &dN2[0]);
    fvec pij = fvec::invsqrt(c_1 + sum_pij + PijS);
    fvec tmp = c_m0_5 * pij * pij * pij;
    int buf_idx;
    for (buf_idx = 0; buf_idx < buf_len; buf_idx++) {
      fvec rikx = data->rikx_buf[buf_idx];
      fvec riky = data->riky_buf[buf_idx];
      fvec rikz = data->rikz_buf[buf_idx];
      fvec rikmag = data->rikmag_buf[buf_idx];
      fvec cosjik = data->cosjik_buf[buf_idx];
      bvec ktype_mask = data->ktype_buf[buf_idx];
      ivec k = data->k_buf[buf_idx];
      fvec g = data->g_buf[buf_idx];
      fvec dgdc = data->dgdc_buf[buf_idx];
      fvec ex_lam = data->ex_lam_buf[buf_idx];
      fvec wik = data->wik_buf[buf_idx];
      fvec dwik = data->dwik_buf[buf_idx];
      bvec mask = data->mask_buf[buf_idx];
      fvec invrikm = fvec::recip(rikmag);
      fvec rjkx = rikx -  rijx;
      fvec rjky = riky -  rijy;
      fvec rjkz = rikz -  rijz;
      fvec rjkmag = fvec::sqrt(
           rjkx *  rjkx +  rjky *  rjky +  rjkz *  rjkz);
      fvec rijrik = c_2_0 *  rijmag *  rikmag;
      fvec rr = rijmag *  rijmag -  rikmag *  rikmag;
      fvec dctdjk = c_m2_0 /  rijrik;
      fvec dctdik = (rjkmag *  rjkmag -  rr) / ( rijrik *  rikmag *  rikmag);
      fvec dctdij = (rjkmag *  rjkmag +  rr) / ( rijrik *  rijmag *  rijmag);
      fvec fi[3], fj[3], fk[3];
      fvec pref = c_0_5 *  VA *  tmp;
      fvec tmp20 = pref *  wik *  dgdc *  ex_lam;
      fj[0] = fj[1] = fj[2] = fvec::setzero();
      fvec tmpdik = tmp20 *  dctdik;
      fi[0] = fvec::setzero() -  tmpdik *  rikx;
      fi[1] = fvec::setzero() -  tmpdik *  riky;
      fi[2] = fvec::setzero() -  tmpdik *  rikz;
      fk[0] = tmpdik *  rikx;
      fk[1] = tmpdik *  riky;
      fk[2] = tmpdik *  rikz;

      fvec tmpdij = tmp20 *  dctdij;
      fij[0] = fvec::mask_sub(fij[0], mask, fij[0], tmpdij *  rijx);
      fij[1] = fvec::mask_sub(fij[1], mask, fij[1], tmpdij *  rijy);
      fij[2] = fvec::mask_sub(fij[2], mask, fij[2], tmpdij *  rijz);

      fvec tmpdjk = tmp20 *  dctdjk;
      fi[0] = fi[0] -  tmpdjk *  rjkx;
      fi[1] = fi[1] -  tmpdjk *  rjky;
      fi[2] = fi[2] -  tmpdjk *  rjkz;
      fk[0] = fk[0] +  tmpdjk *  rjkx;
      fk[1] = fk[1] +  tmpdjk *  rjky;
      fk[2] = fk[2] +  tmpdjk *  rjkz;
      fij[0] = fvec::mask_add(fij[0], mask, fij[0], tmpdjk *  rjkx);
      fij[1] = fvec::mask_add(fij[1], mask, fij[1], tmpdjk *  rjky);
      fij[2] = fvec::mask_add(fij[2], mask, fij[2], tmpdjk *  rjkz);

      if (itype) {
        fvec tmp21 = pref *  wik *  g *  ex_lam *  c_4_0;
        fvec tmp21ij = tmp21 *  invrijm;
        fij[0] = fvec::mask_sub(fij[0], mask, fij[0], tmp21ij * rijx);
        fij[1] = fvec::mask_sub(fij[1], mask, fij[1], tmp21ij * rijy);
        fij[2] = fvec::mask_sub(fij[2], mask, fij[2], tmp21ij * rijz);
        fvec tmp21ik = tmp21 * invrikm;
        fi[0] = fi[0] +  tmp21ik *  rikx;
        fi[1] = fi[1] +  tmp21ik *  riky;
        fi[2] = fi[2] +  tmp21ik *  rikz;
        fk[0] = fk[0] -  tmp21ik *  rikx;
        fk[1] = fk[1] -  tmp21ik *  riky;
        fk[2] = fk[2] -  tmp21ik *  rikz;
      }

      // coordination forces

      // dwik forces
      fvec tmp22 = pref *  dwik *  g *  ex_lam *  invrikm;
      fi[0] = fi[0] -  tmp22 *  rikx;
      fi[1] = fi[1] -  tmp22 *  riky;
      fi[2] = fi[2] -  tmp22 *  rikz;
      fk[0] = fk[0] +  tmp22 *  rikx;
      fk[1] = fk[1] +  tmp22 *  riky;
      fk[2] = fk[2] +  tmp22 *  rikz;

      // PIJ forces
      fvec dN2ktype = fvec::mask_blend(ktype_mask, dN2[0], dN2[1]);
      fvec tmp23 = pref *  dN2ktype *  dwik *  invrikm;
      fi[0] = fi[0] -  tmp23 *  rikx;
      fi[1] = fi[1] -  tmp23 *  riky;
      fi[2] = fi[2] -  tmp23 *  rikz;
      fk[0] = fk[0] +  tmp23 *  rikx;
      fk[1] = fk[1] +  tmp23 *  riky;
      fk[2] = fk[2] +  tmp23 *  rikz;

      // dgdN forces
      fvec tmp24 = pref *  sum_dpij_dN *  dwik *  invrikm;
      fi[0] = fi[0] -  tmp24 *  rikx;
      fi[1] = fi[1] -  tmp24 *  riky;
      fi[2] = fi[2] -  tmp24 *  rikz;
      fk[0] = fk[0] +  tmp24 *  rikx;
      fk[1] = fk[1] +  tmp24 *  riky;
      fk[2] = fk[2] +  tmp24 *  rikz;

      result_f_i_x = fvec::mask_add(result_f_i_x, mask, result_f_i_x, fi[0]);
      result_f_i_y = fvec::mask_add(result_f_i_y, mask, result_f_i_y, fi[1]);
      result_f_i_z = fvec::mask_add(result_f_i_z, mask, result_f_i_z, fi[2]);
      result_f_j_x = fvec::mask_add(result_f_j_x, mask, result_f_j_x, fj[0]);
      result_f_j_y = fvec::mask_add(result_f_j_y, mask, result_f_j_y, fj[1]);
      result_f_j_z = fvec::mask_add(result_f_j_z, mask, result_f_j_z, fj[2]);

      data->force_k_x_buf[buf_idx] = fk[0];
      data->force_k_y_buf[buf_idx] = fk[1];
      data->force_k_z_buf[buf_idx] = fk[2];
    }
    data->force_i_x = result_f_i_x;
    data->force_i_y = result_f_i_y;
    data->force_i_z = result_f_i_z;
    data->force_j_x = result_f_j_x;
    data->force_j_y = result_f_j_y;
    data->force_j_z = result_f_j_z;
    return pij;
  }
  exceed_buffer:
  data->buf_len = -1;
  return fvec::setzero();
}

/*
 * Apply the force values stored iin aut_frebo_data to
 * the respective neighbors.
 */
static void aut_frebo_data_writeback(
    KernelArgsAIREBOT<flt_t,acc_t> * _noalias ka,
    struct aut_frebo_data * _noalias data) {
  ResultForceT<acc_t> * _noalias result_f = ka->result_f;
  flt_t fk_x_buf[fvec::VL] __attribute__((aligned(64)));
  flt_t fk_y_buf[fvec::VL] __attribute__((aligned(64)));
  flt_t fk_z_buf[fvec::VL] __attribute__((aligned(64)));
  int fk_k_buf[ivec::VL] __attribute__((aligned(64)));
  int buf_idx;
  for (buf_idx = 0; buf_idx < data->buf_len; buf_idx++) {
    ivec k = data->k_buf[buf_idx];
    bvec active_mask = data->mask_buf[buf_idx];

    fvec::store(fk_x_buf, data->force_k_x_buf[buf_idx]);
    fvec::store(fk_y_buf, data->force_k_y_buf[buf_idx]);
    fvec::store(fk_z_buf, data->force_k_z_buf[buf_idx]);
    ivec::store(fk_k_buf, k);

    int lane;
    for (lane = 0; lane < fvec::VL; lane++) {
      if (bvec::test_at(active_mask, lane)) {} else continue;
      int kk = fk_k_buf[lane];
      result_f[kk].x += fk_x_buf[lane];
      result_f[kk].y += fk_y_buf[lane];
      result_f[kk].z += fk_z_buf[lane];
    }
  }
}

static void aut_frebo_N_spline_force(
     KernelArgsAIREBOT<flt_t,acc_t> * _noalias ka,
     struct aut_frebo_data * _noalias data, int itype, int jtype, ivec vi,
     ivec /*vj*/, fvec VA, fvec dN, fvec dNconj, fvec Nconj) {
  ivec c_i1 = ivec::set1(1);
  fvec c_2 = fvec::set1(2);
  fvec c_TOL = fvec::set1(TOL);
  ResultForceT<acc_t> * _noalias result_f = ka->result_f;
  AtomAIREBOT<flt_t> * _noalias x = ka->x;
  int * _noalias map = ka->map;
  flt_t * _noalias nC = ka->nC;
  flt_t * _noalias nH = ka->nH;
  fvec x_i, y_i, z_i;
  x_i = data->x_i;
  y_i = data->y_i;
  z_i = data->z_i;
  fvec Nmin = fvec::set1(ka->params.Nmin);
  fvec Nmax = fvec::set1(ka->params.Nmax);
  int map_i_scalar = 0;
  {
    int i;
    for (i = 1; i < ka->num_types; i++) {
      if (ka->map[i])
        map_i_scalar |= (1 << i);
    }
  }
  ivec map_i = ivec::set1(map_i_scalar);
  fvec dN2[2];
  ivec kk = ivec::setzero();
  fvec rcmax0 = fvec::set1(ka->params.rcmax[itype][0]);
  fvec rcmax1 = fvec::set1(ka->params.rcmax[itype][1]);
  fvec rcmin0 = fvec::set1(ka->params.rcmin[itype][0]);
  fvec rcmin1 = fvec::set1(ka->params.rcmin[itype][1]);
  fvec result_f_i_x = fvec::setzero();
  fvec result_f_i_y = fvec::setzero();
  fvec result_f_i_z = fvec::setzero();
  int buf_idx;
  for (buf_idx = 0; buf_idx < data->buf_len; buf_idx++) {
    ivec k = data->k_buf[buf_idx];
    bvec active_mask = data->mask_buf[buf_idx];
    fvec rikx = data->rikx_buf[buf_idx];
    fvec riky = data->riky_buf[buf_idx];
    fvec rikz = data->rikz_buf[buf_idx];
    fvec rikmag = data->rikmag_buf[buf_idx];
    bvec ktype_mask = data->ktype_buf[buf_idx];

    fvec dwik = data->dwik_buf[buf_idx];
    fvec wik = data->wik_buf[buf_idx];

    fvec dNki = data->dcutN_buf[buf_idx];
    fvec SpN = data->cutN_buf[buf_idx];

    fvec invrikmag = fvec::recip(rikmag);
    fvec pref = VA *  dwik *  invrikmag;
    fvec fdN = dN *  pref;
    fvec fdNconj = pref *  SpN *  c_2 *  dNconj *  Nconj;
    fvec ffactor = fdN;
    bvec ktype_is_C = ~ ktype_mask;
    ffactor = fvec::mask_add(ffactor, ktype_is_C, ffactor,  fdNconj);

    fvec fkx = ffactor *  rikx;
    fvec fky = ffactor *  riky;
    fvec fkz = ffactor *  rikz;

    data->force_k_x_buf[buf_idx] = data->force_k_x_buf[buf_idx] +  fkx;
    data->force_k_y_buf[buf_idx] = data->force_k_y_buf[buf_idx] +  fky;
    data->force_k_z_buf[buf_idx] = data->force_k_z_buf[buf_idx] +  fkz;

    result_f_i_x = fvec::mask_sub(result_f_i_x, active_mask, result_f_i_x, fkx);
    result_f_i_y = fvec::mask_sub(result_f_i_y, active_mask, result_f_i_y, fky);
    result_f_i_z = fvec::mask_sub(result_f_i_z, active_mask, result_f_i_z, fkz);

    bvec need_k_neighs = fvec::mask_cmpnle(active_mask, fvec::abs(dNki), c_TOL)
      & ktype_is_C;
    if (bvec::test_any_set(need_k_neighs)) {
      int lane;
      for (lane = 0; lane < fvec::VL; lane++) {
        if (! bvec::test_at(need_k_neighs, lane)) continue;
        int kk = ivec::at(k, lane);
        int k = kk;
        int ktype = map[x[k].w];
        int i = ivec::at(vi, lane);
        fvec oldVA = VA;
        double VA = fvec::at(oldVA, lane);
        fvec oldwik = wik;
        double wik = fvec::at(oldwik, lane);
        fvec olddNconj = dNconj;
        double dNconj = fvec::at(olddNconj, lane);
        fvec oldNconj = Nconj;
        double Nconj = fvec::at(oldNconj, lane);
        fvec olddNki = dNki;
        double dNki = fvec::at(olddNki, lane);
        int * neighs_k = ka->neigh_rebo.entries + ka->neigh_rebo.offset[k];
        int nnum = ka->neigh_rebo.num[k];
        int nn;
        for (nn = 0; nn < nnum; nn++) {
          int n = neighs_k[nn];
          if (n == i) continue;
          double rknx = x[k].x - x[n].x;
          double rkny = x[k].y - x[n].y;
          double rknz = x[k].z - x[n].z;
          double rknmag = sqrt(rknx * rknx + rkny * rkny + rknz * rknz);
          int ntype = map[x[n].w];
          double rcminkn = ka->params.rcmin[ktype][ntype];
          double rcmaxkn = ka->params.rcmax[ktype][ntype];
          double dwkn;
          Sp(rknmag, rcminkn, rcmaxkn, &dwkn);
          double ffactor = VA * dNconj * 2 * Nconj * wik * dNki * dwkn / rknmag;
          result_f[k].x -= ffactor * rknx;
          result_f[k].y -= ffactor * rkny;
          result_f[k].z -= ffactor * rknz;
          result_f[n].x += ffactor * rknx;
          result_f[n].y += ffactor * rkny;
          result_f[n].z += ffactor * rknz;
        }
      }
    }
  }
  data->force_i_x = data->force_i_x +  result_f_i_x;
  data->force_i_y = data->force_i_y +  result_f_i_y;
  data->force_i_z = data->force_i_z +  result_f_i_z;
}

static fvec aut_frebo_pi_rc_pd(KernelArgsAIREBOT<flt_t,acc_t> * ka, int itype,
                               int jtype, fvec Nij, fvec Nji, fvec Nijconj,
                               fvec * dN3) {
  flt_t ret[fvec::VL] __attribute__((aligned(64)));
  flt_t dN3ret[3][fvec::VL] __attribute__((aligned(64)));
  int i;
  for (i = 0; i < fvec::VL; i++) {
    flt_t dN3tmp[3];
    ret[i] = frebo_pi_rc(ka, itype, jtype, fvec::at(Nij, i), fvec::at(Nji, i),
                         fvec::at(Nijconj, i), &dN3tmp[0]);
    dN3ret[0][i] = dN3tmp[0];
    dN3ret[1][i] = dN3tmp[1];
    dN3ret[2][i] = dN3tmp[2];
  }
  dN3[0] = fvec::load(&dN3ret[0][0]);
  dN3[1] = fvec::load(&dN3ret[1][0]);
  dN3[2] = fvec::load(&dN3ret[2][0]);
  return fvec::load(&ret[0]);
}

static fvec aut_frebo_Tij(KernelArgsAIREBOT<flt_t,acc_t> * ka, int itype,
                          int jtype, fvec Nij, fvec Nji, fvec Nijconj,
                          fvec * dN3) {
  flt_t ret[fvec::VL] __attribute__((aligned(64)));
  flt_t dN3ret[3][fvec::VL] __attribute__((aligned(64)));
  int i;
  for (i = 0; i < fvec::VL; i++) {
    flt_t dN3tmp[3];
    ret[i] = frebo_Tij(ka, itype, jtype, fvec::at(Nij, i), fvec::at(Nji, i),
                       fvec::at(Nijconj, i), &dN3tmp[0]);
    dN3ret[0][i] = dN3tmp[0];
    dN3ret[1][i] = dN3tmp[1];
    dN3ret[2][i] = dN3tmp[2];
  }
  dN3[0] = fvec::load(&dN3ret[0][0]);
  dN3[1] = fvec::load(&dN3ret[1][0]);
  dN3[2] = fvec::load(&dN3ret[2][0]);
  return fvec::load(&ret[0]);
}

static fvec aut_frebo_sum_omega(
    KernelArgsAIREBOT<flt_t,acc_t> * _noalias ka,
    struct aut_frebo_data * _noalias i_data,
    struct aut_frebo_data * _noalias j_data,
    int /*itype*/, int /*jtype*/,
    ivec /*vi*/, ivec /*vj*/,
    fvec r23x, fvec r23y, fvec r23z, fvec r23mag,
    fvec VA, fvec fij[3]
) {
  fvec c_1 = fvec::set1(1);
  fvec c_m1 = fvec::set1(-1);
  fvec c_2 = fvec::set1(2);
  fvec c_m2 = fvec::set1(-2);
  fvec sum_omega = fvec::setzero();
  fvec thmin = fvec::set1(ka->params.thmin);
  fvec thmax = fvec::set1(ka->params.thmax);
  // 2 == i, 3 == j
  fvec r32x = fvec::setzero() -  r23x;
  fvec r32y = fvec::setzero() -  r23y;
  fvec r32z = fvec::setzero() -  r23z;
  int buf_idx_i, buf_idx_j;
  for (buf_idx_i = 0; buf_idx_i < i_data->buf_len; buf_idx_i++) {
    // a1 == k == buf_idx_i
    bvec mask_start = i_data->mask_buf[buf_idx_i];
    fvec r21x = i_data->rikx_buf[buf_idx_i]; // a2 - a1 -> i - k
    fvec r21y = i_data->riky_buf[buf_idx_i]; // a2 - a1 -> i - k
    fvec r21z = i_data->rikz_buf[buf_idx_i]; // a2 - a1 -> i - k
    fvec r21mag = i_data->rikmag_buf[buf_idx_i];
    // TODO use buffered cosjik
    fvec cos321 = (
        r23x *  r21x +  r23y *  r21y +  r23z *  r21z) / ( r23mag *  r21mag);
    cos321 = fvec::min(c_1, fvec::max(c_m1, cos321));
    fvec sin321 = fvec::sqrt(c_1 -  cos321 *  cos321);
    bvec mask_outer = fvec::cmpneq(fvec::setzero(), sin321) & mask_start;
    // add "continue"
    fvec sink2i = fvec::mask_recip(fvec::undefined(), mask_outer,
                                   sin321 * sin321);
    fvec rik2i = fvec::mask_recip(fvec::undefined(), mask_outer,
                                  r21mag * r21mag);
    fvec rr = r23mag *  r23mag -  r21mag *  r21mag;
    fvec r31x = r21x -  r23x;
    fvec r31y = r21y -  r23y;
    fvec r31z = r21z -  r23z;
    fvec r31mag2 = r31x *  r31x +  r31y *  r31y +  r31z *  r31z;
    fvec rijrik = c_2 *  r23mag *  r21mag;
    fvec r21mag2 = r21mag *  r21mag;
    fvec dctik = fvec::mask_div(fvec::undefined(), mask_outer, r31mag2 -  rr,
                                rijrik *  r21mag2);
    fvec dctij = fvec::mask_div(fvec::undefined(), mask_outer, r31mag2 +  rr,
                                rijrik *  r23mag *  r23mag);
    fvec dctjk = fvec::mask_div(fvec::undefined(), mask_outer, c_m2, rijrik);
    fvec dw21 = i_data->dwik_buf[buf_idx_i];
    fvec w21 = i_data->wik_buf[buf_idx_i];
    fvec dtsjik;
    fvec tspjik = aut_Sp2_deriv(cos321, thmin, thmax, &dtsjik);
    dtsjik = fvec::setzero() -  dtsjik; // todo replace by appropriate xor.
    ivec k = i_data->k_buf[buf_idx_i];
    for (buf_idx_j = 0; buf_idx_j < j_data->buf_len; buf_idx_j++) {
      // check l == k in second loop.
      // l == a4 == buf_idx_j
      ivec l = j_data->k_buf[buf_idx_j];
      bvec mask_inner_0 = ivec::mask_cmpneq(mask_outer, k, l) &
        j_data->mask_buf[buf_idx_j];
      // add "continue"
      fvec r34x = j_data->rikx_buf[buf_idx_j];
      fvec r34y = j_data->riky_buf[buf_idx_j];
      fvec r34z = j_data->rikz_buf[buf_idx_j];
      fvec r34mag = j_data->rikmag_buf[buf_idx_j];
      fvec cos234 = fvec::mask_div(fvec::undefined(), mask_inner_0,
                                   r32x * r34x + r32y * r34y + r32z * r34z,
                                   r23mag * r34mag);
      cos234 = fvec::min(c_1, fvec::max(c_m1, cos234));
      fvec sin234 = fvec::mask_sqrt(fvec::undefined(), mask_inner_0,
                                    c_1 - cos234 * cos234);
      bvec mask_inner_1 = fvec::mask_cmpneq(mask_inner_0, sin234,
                                            fvec::setzero());
      // add "continue"
      fvec sinl2i = fvec::mask_recip(fvec::undefined(), mask_inner_1,
                                     sin234 * sin234);
      fvec rjl2i = fvec::mask_recip(fvec::undefined(), mask_inner_1,
                                    r34mag * r34mag);
      fvec dw34 = j_data->dwik_buf[buf_idx_j];
      fvec w34 = j_data->wik_buf[buf_idx_j];
      fvec rr = r23mag *  r23mag - r34mag * r34mag;
      fvec r24x = r23x +  r34x;
      fvec r24y = r23y +  r34y;
      fvec r24z = r23z +  r34z;
      fvec r242 = r24x *  r24x +  r24y *  r24y +  r24z *  r24z;
      fvec rijrjl = c_2 *  r23mag *  r34mag;
      fvec rjl2 = r34mag *  r34mag;
      fvec dctjl = fvec::mask_div(fvec::undefined(), mask_inner_1, r242 -  rr,
                                  rijrjl *  rjl2);
      fvec dctji = fvec::mask_div(fvec::undefined(), mask_inner_1, r242 +  rr,
                                  rijrjl *  r23mag *  r23mag);
      fvec dctil = fvec::mask_div(fvec::undefined(), mask_inner_1, c_m2,
                                  rijrjl);
      fvec dtsijl;
      fvec tspijl = aut_Sp2_deriv(cos234, thmin, thmax, &dtsijl);
      dtsijl = fvec::setzero() -  dtsijl;
      fvec prefactor = VA;

      fvec cross321x = r32y *  r21z -  r32z *  r21y;
      fvec cross321y = r32z *  r21x -  r32x *  r21z;
      fvec cross321z = r32x *  r21y -  r32y *  r21x;
      fvec cross234x = r23y *  r34z -  r23z *  r34y;
      fvec cross234y = r23z *  r34x -  r23x *  r34z;
      fvec cross234z = r23x *  r34y -  r23y *  r34x;

      fvec cwnum = cross321x * cross234x + cross321y * cross234y + cross321z *
        cross234z;
      fvec cwnom = r21mag * r34mag * r23mag * r23mag * sin321 * sin234;
      fvec om1234 = fvec::mask_div(fvec::undefined(), mask_inner_1, cwnum,
                                   cwnom);
      fvec cw = om1234;
      fvec sum_omega_contrib = (c_1 -  om1234 *  om1234) *  w21 *  w34 *
        (c_1 -  tspjik) * ( c_1 -  tspijl);
      sum_omega = fvec::mask_add(sum_omega, mask_inner_1, sum_omega,
                                 sum_omega_contrib);
      fvec dt1dik = rik2i -  dctik *  sink2i *  cos321;
      fvec dt1djk = fvec::setzero() -  dctjk *  sink2i *  cos321;
      fvec dt1djl = rjl2i -  dctjl *  sinl2i *  cos234;
      fvec dt1dil = fvec::setzero() -  dctil *  sinl2i *  cos234;
      fvec dt1dij =   fvec::mask_div(fvec::undefined(), mask_inner_1, c_2,
                                     r23mag * r23mag) -
        dctij * sink2i * cos321 -  dctji *  sinl2i *  cos234;

      fvec dt2dikx = r23y *  cross234z -  r23z *  cross234y;
      fvec dt2diky = r23z *  cross234x -  r23x *  cross234z;
      fvec dt2dikz = r23x *  cross234y -  r23y *  cross234x;

      fvec dt2djlx = r23z *  cross321y -  r23y *  cross321z;
      fvec dt2djly = r23x *  cross321z -  r23z *  cross321x;
      fvec dt2djlz = r23y *  cross321x -  r23x *  cross321y;

      fvec dt2dijx = r21z *  cross234y +  r34y *  cross321z -
        ( r34z *  cross321y +  r21y *  cross234z);
      fvec dt2dijy = r21x *  cross234z +  r34z *  cross321x -
        ( r34x *  cross321z +  r21z *  cross234x);
      fvec dt2dijz = r21y *  cross234x +  r34x *  cross321y -
        ( r34y *  cross321x +  r21x *  cross234y);

      fvec aa = prefactor *  c_2 *  fvec::mask_div(fvec::undefined(),
                                                   mask_inner_1, cw, cwnom) *
        w21 *  w34 *  (c_1 -  tspjik) * ( c_1 -  tspijl);
      fvec aaa1 = (fvec::setzero() - prefactor) * (c_1 - om1234 * om1234) *
        (c_1 - tspjik) * (c_1 - tspijl);
      fvec aaa2 = (fvec::setzero() -  prefactor) * (c_1 -  om1234 *  om1234) *
        w21 * w34;
      fvec at2 = aa * cwnum;

      fvec fcijpc = aaa2 * dtsjik * dctij * (c_1 - tspijl) +  aaa2 * dtsijl *
        dctji * (c_1 - tspjik) - dt1dij * at2;
      fvec fcikpc =  aaa2 * dtsjik * dctik * (c_1 - tspijl) - dt1dik * at2;
      fvec fcjlpc =  aaa2 * dtsijl * dctjl * (c_1 - tspjik) - dt1djl * at2;
      fvec fcjkpc =  aaa2 * dtsjik * dctjk * (c_1 - tspijl) - dt1djk * at2;
      fvec fcilpc =  aaa2 * dtsijl * dctil * (c_1 - tspjik) - dt1dil * at2;

      fvec F23x = fcijpc *  r23x +  aa *  dt2dijx;
      fvec F23y = fcijpc *  r23y +  aa *  dt2dijy;
      fvec F23z = fcijpc *  r23z +  aa *  dt2dijz;

      fvec F12x = fcikpc *  r21x +  aa *  dt2dikx;
      fvec F12y = fcikpc *  r21y +  aa *  dt2diky;
      fvec F12z = fcikpc *  r21z +  aa *  dt2dikz;

      fvec F34x = fcjlpc *  r34x +  aa *  dt2djlx;
      fvec F34y = fcjlpc *  r34y +  aa *  dt2djly;
      fvec F34z = fcjlpc *  r34z +  aa *  dt2djlz;

      fvec F31x = fcjkpc *  r31x;
      fvec F31y = fcjkpc *  r31y;
      fvec F31z = fcjkpc *  r31z;

      fvec F24x = fcilpc *  r24x;
      fvec F24y = fcilpc *  r24y;
      fvec F24z = fcilpc *  r24z;

      fvec f1x = fvec::setzero() - ( F12x +  F31x);
      fvec f1y = fvec::setzero() - ( F12y +  F31y);
      fvec f1z = fvec::setzero() - ( F12z +  F31z);
      fvec f2x = F12x +  F31x;
      fvec f2y = F12y +  F31y;
      fvec f2z = F12z +  F31z;
      fvec f3x = F34x +  F24x;
      fvec f3y = F34y +  F24y;
      fvec f3z = F34z +  F24z;
      fvec f4x = fvec::setzero() - ( F34x +  F24x);
      fvec f4y = fvec::setzero() - ( F34y +  F24y);
      fvec f4z = fvec::setzero() - ( F34z +  F24z);

      fij[0] = fvec::mask_add(fij[0], mask_inner_1, fij[0],
          F23x +  F24x -  F31x);
      fij[1] = fvec::mask_add(fij[1], mask_inner_1, fij[1],
          F23y +  F24y -  F31y);
      fij[2] = fvec::mask_add(fij[2], mask_inner_1, fij[2],
          F23z +  F24z -  F31z);

      fvec tmp20 = VA * (c_1 - om1234 * om1234) * (c_1 - tspjik) *
        (c_1 - tspijl) * dw21 * w34 * fvec::mask_recip(fvec::undefined(),
                                                       mask_inner_1, r21mag);
      f2x = f2x -  tmp20 *  r21x;
      f2y = f2y -  tmp20 *  r21y;
      f2z = f2z -  tmp20 *  r21z;
      f1x = f1x +  tmp20 *  r21x;
      f1y = f1y +  tmp20 *  r21y;
      f1z = f1z +  tmp20 *  r21z;

      fvec tmp21 = VA * (c_1 - om1234 * om1234) * (c_1 - tspjik) *
        (c_1 - tspijl) * w21 * dw34 * fvec::mask_recip(fvec::undefined(),
                                                       mask_inner_1, r34mag);
      f3x = f3x -  tmp21 *  r34x;
      f3y = f3y -  tmp21 *  r34y;
      f3z = f3z -  tmp21 *  r34z;
      f4x = f4x +  tmp21 *  r34x;
      f4y = f4y +  tmp21 *  r34y;
      f4z = f4z +  tmp21 *  r34z;

      // 1 == buf_idx_i, 2 == i, 3 == j, 4 == buf_idx_j
      i_data->force_k_x_buf[buf_idx_i] =
        fvec::mask_add(i_data->force_k_x_buf[buf_idx_i],
                       mask_inner_1, i_data->force_k_x_buf[buf_idx_i], f1x);
      i_data->force_k_y_buf[buf_idx_i] =
        fvec::mask_add(i_data->force_k_y_buf[buf_idx_i], mask_inner_1,
                       i_data->force_k_y_buf[buf_idx_i], f1y);
      i_data->force_k_z_buf[buf_idx_i] =
        fvec::mask_add(i_data->force_k_z_buf[buf_idx_i], mask_inner_1,
                       i_data->force_k_z_buf[buf_idx_i], f1z);
      i_data->force_i_x =
        fvec::mask_add(i_data->force_i_x, mask_inner_1, i_data->force_i_x, f2x);
      i_data->force_i_y =
        fvec::mask_add(i_data->force_i_y, mask_inner_1, i_data->force_i_y, f2y);
      i_data->force_i_z =
        fvec::mask_add(i_data->force_i_z, mask_inner_1, i_data->force_i_z, f2z);
      j_data->force_i_x =
        fvec::mask_add(j_data->force_i_x, mask_inner_1, j_data->force_i_x, f3x);
      j_data->force_i_y =
        fvec::mask_add(j_data->force_i_y, mask_inner_1, j_data->force_i_y, f3y);
      j_data->force_i_z =
        fvec::mask_add(j_data->force_i_z, mask_inner_1, j_data->force_i_z, f3z);
      j_data->force_k_x_buf[buf_idx_j] =
        fvec::mask_add(j_data->force_k_x_buf[buf_idx_j], mask_inner_1,
                       j_data->force_k_x_buf[buf_idx_j], f4x);
      j_data->force_k_y_buf[buf_idx_j] =
        fvec::mask_add(j_data->force_k_y_buf[buf_idx_j], mask_inner_1,
                       j_data->force_k_y_buf[buf_idx_j], f4y);
      j_data->force_k_z_buf[buf_idx_j] =
        fvec::mask_add(j_data->force_k_z_buf[buf_idx_j], mask_inner_1,
                       j_data->force_k_z_buf[buf_idx_j], f4z);
    }
  }
  return sum_omega;
}

static fvec aut_frebo_pi_dh(
    KernelArgsAIREBOT<flt_t,acc_t> * _noalias ka,
    struct aut_frebo_data * _noalias i_data,
    struct aut_frebo_data * _noalias j_data,
    int itype, int jtype, ivec vi, ivec vj,
    fvec r23x, fvec r23y, fvec r23z, fvec r23mag,
    fvec VA,
    fvec Nij, fvec Nji, fvec Nijconj, fvec NconjtmpI, fvec NconjtmpJ,
    fvec fij[3]
) {
  fvec c_TOL = fvec::set1(TOL);
  fvec dN3[3];
  fvec Tij = aut_frebo_Tij(ka, itype, jtype, Nij, Nji, Nijconj, &dN3[0]);
  bvec TijgtTOLmask = fvec::cmpnle(fvec::abs(Tij), c_TOL);
  fvec sum_omega = fvec::setzero();
  if (bvec::test_any_set(TijgtTOLmask)) {
    sum_omega = aut_frebo_sum_omega(
        ka, i_data, j_data, itype, jtype, vi, vj,
        r23x, r23y, r23z, r23mag, VA *  Tij, fij);
    sum_omega = fvec::mask_blend(TijgtTOLmask, fvec::setzero(), sum_omega);
    aut_frebo_N_spline_force(ka, i_data, itype, jtype, vi, vj, VA * sum_omega,
                             dN3[0], dN3[2], NconjtmpI);
    aut_frebo_N_spline_force(ka, j_data, jtype, itype, vj, vi, VA * sum_omega,
                             dN3[1], dN3[2], NconjtmpJ);
  }
  return Tij *  sum_omega;
}

/*
 We can reuse the aut_frebo_data buffers here to do this calculation very
 cheaply.
*/
static void aut_torsion_vec(
    KernelArgsAIREBOT<flt_t,acc_t> * ka,
    struct aut_frebo_data * i_data,
    struct aut_frebo_data * j_data,
    ivec /*i*/, ivec /*j*/, fvec wij, fvec dwij
) {
  AtomAIREBOT<flt_t> * x = ka->x;
  int * map = ka->map;
  flt_t (*epsilonT)[2] = ka->params.epsilonT;
  fvec epsilonT00 = fvec::set1(epsilonT[0][0]);
  fvec epsilonT01 = fvec::set1(epsilonT[0][1]);
  fvec epsilonT10 = fvec::set1(epsilonT[1][0]);
  fvec epsilonT11 = fvec::set1(epsilonT[1][1]);
  fvec thmin = fvec::set1(ka->params.thmin);
  fvec thmax = fvec::set1(ka->params.thmax);

  const fvec c_1_0 = fvec::set1(1.0);
  const fvec c_0_5 = fvec::set1(0.5);
  const fvec c_0_1 = fvec::set1(0.1);
  const fvec c_2_0 = fvec::set1(2.0);
  const fvec c_2_5 = fvec::set1(2.5);
  const fvec c_256_405 = fvec::set1(256.0/405.0);

  fvec del32x = j_data->x_i -  i_data->x_i;
  fvec del32y = j_data->y_i -  i_data->y_i;
  fvec del32z = j_data->z_i -  i_data->z_i;
  fvec rsq = del32x * del32x +  del32y * del32y +  del32z * del32z;
  fvec r32 = fvec::sqrt(rsq);
  fvec del23x = fvec::setzero() -  del32x;
  fvec del23y = fvec::setzero() -  del32y;
  fvec del23z = fvec::setzero() -  del32z;
  fvec r23 = r32;
  fvec w23 = wij;
  fvec dw23 = dwij;

  for (int buf_idx_i = 0; buf_idx_i < i_data->buf_len; buf_idx_i++) {
    bvec mask_start = i_data->mask_buf[buf_idx_i];
    fvec del21x = i_data->rikx_buf[buf_idx_i]; // a2 - a1 -> i - k
    fvec del21y = i_data->riky_buf[buf_idx_i]; // a2 - a1 -> i - k
    fvec del21z = i_data->rikz_buf[buf_idx_i]; // a2 - a1 -> i - k
    fvec r21 = i_data->rikmag_buf[buf_idx_i];
    fvec cos321 = i_data->cosjik_buf[buf_idx_i];
    fvec sin321 = fvec::sqrt(c_1_0 -  cos321 *  cos321);
    // strictly equivalent to sin321 < TOL
    mask_start = fvec::mask_cmpneq(mask_start, fvec::setzero(), sin321);
    if (! bvec::test_any_set(mask_start)) continue;

    fvec deljkx = del21x -  del23x;
    fvec deljky = del21y -  del23y;
    fvec deljkz = del21z -  del23z;
    fvec rjk2 = deljkx * deljkx +  deljky * deljky + deljkz * deljkz;
    fvec rjk = fvec::sqrt(rjk2);
    fvec rik2 = r21 *  r21;
    fvec w21 = i_data->wik_buf[buf_idx_i];
    fvec dw21 = i_data->dwik_buf[buf_idx_i];

    fvec rij = r32;
    fvec rik = r21;
    fvec rij2 = r32 *  r32;
    fvec dtsjik;
    fvec tspjik = aut_Sp2_deriv(cos321, thmin, thmax, &dtsjik);
    dtsjik = fvec::setzero() -  dtsjik;

    bvec ktype_mask = i_data->ktype_buf[buf_idx_i];
    fvec epsilonT0 = fvec::mask_blend(ktype_mask, epsilonT00, epsilonT10);
    fvec epsilonT1 = fvec::mask_blend(ktype_mask, epsilonT01, epsilonT11);

    ivec k = i_data->k_buf[buf_idx_i];
    for (int buf_idx_j = 0; buf_idx_j < j_data->buf_len; buf_idx_j++) {
      ivec l = j_data->k_buf[buf_idx_j];
      bvec mask_inner_0 = ivec::mask_cmpneq(mask_start, k, l) &
        j_data->mask_buf[buf_idx_j];
      if (! bvec::test_any_set(mask_inner_0)) continue;
      fvec del34x = j_data->rikx_buf[buf_idx_j];
      fvec del34y = j_data->riky_buf[buf_idx_j];
      fvec del34z = j_data->rikz_buf[buf_idx_j];
      fvec r34 = j_data->rikmag_buf[buf_idx_j];
      bvec ltype_mask = j_data->ktype_buf[buf_idx_j];
      fvec cos234 = j_data->cosjik_buf[buf_idx_j];
      fvec sin234 = fvec::sqrt(c_1_0 -  cos234 *  cos234);
      // strictly equivalent to sin234 < TOL
      mask_inner_0 = fvec::mask_cmpneq(mask_inner_0, sin234, fvec::setzero());
      if (! bvec::test_any_set(mask_inner_0)) continue;
      fvec dw34 = j_data->dwik_buf[buf_idx_j];
      fvec w34 = j_data->wik_buf[buf_idx_j];
      fvec delilx = del23x +  del34x;
      fvec delily = del23y +  del34y;
      fvec delilz = del23z +  del34z;
      fvec ril2 = delilx * delilx +  delily * delily + delilz * delilz;
      fvec ril = fvec::sqrt(ril2);
      fvec rjl2 = r34 *  r34;

      fvec rjl = r34;
      fvec dtsijl;
      fvec tspijl = aut_Sp2_deriv(cos234, thmin, thmax, &dtsijl);
      dtsijl = fvec::setzero() -  dtsijl;
      fvec cross321x = del32y * del21z - del32z * del21y;
      fvec cross321y = del32z * del21x - del32x * del21z;
      fvec cross321z = del32x * del21y - del32y * del21x;
      fvec cross321mag = fvec::sqrt(cross321x * cross321x +
                                    cross321y * cross321y +
                                    cross321z * cross321z);
      fvec cross234x = del23y * del34z - del23z * del34y;
      fvec cross234y = del23z * del34x - del23x * del34z;
      fvec cross234z = del23x * del34y - del23y * del34x;
      fvec cross234mag = fvec::sqrt(cross234x * cross234x +
                                    cross234y * cross234y +
                                    cross234z * cross234z);
      fvec cwnum = cross321x * cross234x + cross321y * cross234y +
        cross321z * cross234z;
      fvec cwnom = r21 * r34 * r32 * r32 * sin321 * sin234;
      fvec cw = cwnum /  cwnom;

      fvec cw2 = c_0_5 * ( c_1_0 - cw);
      fvec ekijl = fvec::mask_blend(ltype_mask, epsilonT0, epsilonT1);
      fvec Ec = c_256_405 * ekijl;
      fvec cw2_5 = cw2 *  cw2 *  cw2 *  cw2 *  cw2;
      fvec Vtors = Ec *  cw2_5 -  ekijl *  c_0_1;

      fvec evdwl = Vtors * w21 * w23 * w34 * (c_1_0-tspjik) * (c_1_0-tspijl);
      ka->result_eng += fvec::mask_reduce_add(mask_inner_0, evdwl);

      fvec dndijx  = cross234y * del21z - cross234z * del21y;
      fvec dndijy  = cross234z * del21x - cross234x * del21z;
      fvec dndijz  = cross234x * del21y - cross234y * del21x;

      fvec tmpvecx = del34y * cross321z - del34z * cross321y;
      fvec tmpvecy = del34z * cross321x - del34x * cross321z;
      fvec tmpvecz = del34x * cross321y - del34y * cross321x;

      dndijx = dndijx + tmpvecx;
      dndijy = dndijy + tmpvecy;
      dndijz = dndijz + tmpvecz;

      fvec dndikx = del23y * cross234z - del23z * cross234y;
      fvec dndiky = del23z * cross234x - del23x * cross234z;
      fvec dndikz = del23x * cross234y - del23y * cross234x;

      fvec dndjlx = cross321y * del23z - cross321z * del23y;
      fvec dndjly = cross321z * del23x - cross321x * del23z;
      fvec dndjlz = cross321x * del23y - cross321y * del23x;

      fvec r23sq = r23 *  r23;
      fvec r21sq = r21 *  r21;
      fvec r34sq = r34 *  r34;
      fvec rjksq = rjk *  rjk;
      fvec rilsq = ril *  ril;
      fvec dcidij = (r23sq -  r21sq +  rjksq) / ( c_2_0 *  r23sq *  r21);
      fvec dcidik = (r21sq -  r23sq +  rjksq) / ( c_2_0 *  r21sq *  r23);
      fvec dcidjk = fvec::setzero() -  rjk / ( r23 *  r21);
      fvec dcjdji = (r23sq -  r34sq +  rilsq) / ( c_2_0 *  r23sq *  r34);
      fvec dcjdjl = (r34sq -  r23sq +  rilsq) / ( c_2_0 *  r34sq *  r23);
      fvec dcjdil = fvec::setzero() -  ril / ( r23 *  r34);

      fvec dsidij = fvec::setzero() -  cos321 / sin321 * dcidij;
      fvec dsidik = fvec::setzero() -  cos321 / sin321 * dcidik;
      fvec dsidjk = fvec::setzero() -  cos321 / sin321 * dcidjk;

      fvec dsjdji = fvec::setzero() -  cos234 / sin234 * dcjdji;
      fvec dsjdjl = fvec::setzero() -  cos234 / sin234 * dcjdjl;
      fvec dsjdil = fvec::setzero() -  cos234 / sin234 * dcjdil;

      fvec dxidij = r21 * sin321 + r23 * r21 * dsidij;
      fvec dxidik = r23 * sin321 + r23 * r21 * dsidik;
      fvec dxidjk = r23 * r21 * dsidjk;

      fvec dxjdji = r34 * sin234 + r23 * r34 * dsjdji;
      fvec dxjdjl = r23 * sin234 + r23 * r34 * dsjdjl;
      fvec dxjdil = r23 * r34 * dsjdil;

      fvec ddndij = dxidij * cross234mag + cross321mag * dxjdji;
      fvec ddndik = dxidik * cross234mag;
      fvec ddndjk = dxidjk * cross234mag;
      fvec ddndjl = cross321mag * dxjdjl;
      fvec ddndil = cross321mag * dxjdil;
      fvec dcwddn = fvec::setzero() -  cwnum / ( cwnom * cwnom);
      fvec dcwdn = fvec::recip(cwnom);
      fvec cw2_4 = cw2 *  cw2 *  cw2 *  cw2;
      fvec dvpdcw = c_2_5 * Ec * cw2_4 * w23 * w21 * w34 * (c_1_0 - tspjik) *
        (c_1_0 - tspijl);

      fvec Ftmpx = dvpdcw * (dcwdn * dndijx + dcwddn * ddndij * del23x / r23);
      fvec Ftmpy = dvpdcw * (dcwdn * dndijy + dcwddn * ddndij * del23y / r23);
      fvec Ftmpz = dvpdcw * (dcwdn * dndijz + dcwddn * ddndij * del23z / r23);
      fvec fix = Ftmpx;
      fvec fiy = Ftmpy;
      fvec fiz = Ftmpz;
      fvec fjx = fvec::setzero() - Ftmpx;
      fvec fjy = fvec::setzero() - Ftmpy;
      fvec fjz = fvec::setzero() - Ftmpz;

      Ftmpx = dvpdcw * (dcwdn * dndikx + dcwddn * ddndik * del21x / r21);
      Ftmpy = dvpdcw * (dcwdn * dndiky + dcwddn * ddndik * del21y / r21);
      Ftmpz = dvpdcw * (dcwdn * dndikz + dcwddn * ddndik * del21z / r21);
      fix = fix +  Ftmpx;
      fiy = fiy +  Ftmpy;
      fiz = fiz +  Ftmpz;
      fvec fkx = fvec::setzero() -  Ftmpx;
      fvec fky = fvec::setzero() -  Ftmpy;
      fvec fkz = fvec::setzero() -  Ftmpz;

      Ftmpx = dvpdcw * dcwddn * ddndjk * deljkx / rjk;
      Ftmpy = dvpdcw * dcwddn * ddndjk * deljky / rjk;
      Ftmpz = dvpdcw * dcwddn * ddndjk * deljkz / rjk;
      fjx = fjx +  Ftmpx;
      fjy = fjy +  Ftmpy;
      fjz = fjz +  Ftmpz;
      fkx = fkx -  Ftmpx;
      fky = fky -  Ftmpy;
      fkz = fkz -  Ftmpz;

      Ftmpx = dvpdcw * (dcwdn * dndjlx + dcwddn * ddndjl * del34x / r34);
      Ftmpy = dvpdcw * (dcwdn * dndjly + dcwddn * ddndjl * del34y / r34);
      Ftmpz = dvpdcw * (dcwdn * dndjlz + dcwddn * ddndjl * del34z / r34);
      fjx = fjx +  Ftmpx;
      fjy = fjy +  Ftmpy;
      fjz = fjz +  Ftmpz;
      fvec flx = fvec::setzero() -  Ftmpx;
      fvec fly = fvec::setzero() -  Ftmpy;
      fvec flz = fvec::setzero() -  Ftmpz;

      Ftmpx = dvpdcw * dcwddn * ddndil * delilx / ril;
      Ftmpy = dvpdcw * dcwddn * ddndil * delily / ril;
      Ftmpz = dvpdcw * dcwddn * ddndil * delilz / ril;
      fix = fix +  Ftmpx;
      fiy = fiy +  Ftmpy;
      fiz = fiz +  Ftmpz;
      flx = flx -  Ftmpx;
      fly = fly -  Ftmpy;
      flz = flz -  Ftmpz;

      // coordination forces

      fvec fpair = Vtors * dw21 * w23 * w34 * (c_1_0 - tspjik) *
        (c_1_0 - tspijl) /  r21;
      fix = fix -  del21x * fpair;
      fiy = fiy -  del21y * fpair;
      fiz = fiz -  del21z * fpair;
      fkx = fkx +  del21x * fpair;
      fky = fky +  del21y * fpair;
      fkz = fkz +  del21z * fpair;

      fpair = Vtors * w21 * dw23 * w34 * (c_1_0 - tspjik) * (c_1_0 - tspijl) /
        r23;
      fix = fix -  del23x * fpair;
      fiy = fiy -  del23y * fpair;
      fiz = fiz -  del23z * fpair;
      fjx = fjx +  del23x * fpair;
      fjy = fjy +  del23y * fpair;
      fjz = fjz +  del23z * fpair;

      fpair = Vtors * w21 * w23 * dw34 * (c_1_0 - tspjik) * (c_1_0 - tspijl) /
        r34;
      fjx = fjx -  del34x * fpair;
      fjy = fjy -  del34y * fpair;
      fjz = fjz -  del34z * fpair;
      flx = flx +  del34x * fpair;
      fly = fly +  del34y * fpair;
      flz = flz +  del34z * fpair;

      // additional cut off function forces

      fvec fcpc = fvec::setzero() - Vtors * w21 * w23 * w34 * dtsjik * (c_1_0 -
                                                                        tspijl);
      fpair = fcpc * dcidij / rij;
      fix = fix +  fpair * del23x;
      fiy = fiy +  fpair * del23y;
      fiz = fiz +  fpair * del23z;
      fjx = fjx -  fpair * del23x;
      fjy = fjy -  fpair * del23y;
      fjz = fjz -  fpair * del23z;

      fpair = fcpc * dcidik / rik;
      fix = fix +  fpair * del21x;
      fiy = fiy +  fpair * del21y;
      fiz = fiz +  fpair * del21z;
      fkx = fkx -  fpair * del21x;
      fky = fky -  fpair * del21y;
      fkz = fkz -  fpair * del21z;

      fpair = fcpc * dcidjk / rjk;
      fjx = fjx +  fpair * deljkx;
      fjy = fjy +  fpair * deljky;
      fjz = fjz +  fpair * deljkz;
      fkx = fkx -  fpair * deljkx;
      fky = fky -  fpair * deljky;
      fkz = fkz -  fpair * deljkz;

      fcpc = fvec::setzero() - Vtors * w21 * w23 * w34 * (c_1_0 - tspjik) *
        dtsijl;
      fpair = fcpc * dcjdji / rij;
      fix = fix +  fpair * del23x;
      fiy = fiy +  fpair * del23y;
      fiz = fiz +  fpair * del23z;
      fjx = fjx -  fpair * del23x;
      fjy = fjy -  fpair * del23y;
      fjz = fjz -  fpair * del23z;

      fpair = fcpc * dcjdjl / rjl;
      fjx = fjx +  fpair * del34x;
      fjy = fjy +  fpair * del34y;
      fjz = fjz +  fpair * del34z;
      flx = flx -  fpair * del34x;
      fly = fly -  fpair * del34y;
      flz = flz -  fpair * del34z;

      fpair = fcpc * dcjdil / ril;
      fix = fix +  fpair * delilx;
      fiy = fiy +  fpair * delily;
      fiz = fiz +  fpair * delilz;
      flx = flx -  fpair * delilx;
      fly = fly -  fpair * delily;
      flz = flz -  fpair * delilz;

      // sum per-atom forces into atom force array

      i_data->force_i_x = fvec::mask_add(i_data->force_i_x, mask_inner_0,
                                         i_data->force_i_x, fix);
      i_data->force_i_y = fvec::mask_add(i_data->force_i_y, mask_inner_0,
                                         i_data->force_i_y, fiy);
      i_data->force_i_z = fvec::mask_add(i_data->force_i_z, mask_inner_0,
                                         i_data->force_i_z, fiz);
      i_data->force_j_x = fvec::mask_add(i_data->force_j_x, mask_inner_0,
                                         i_data->force_j_x, fjx);
      i_data->force_j_y = fvec::mask_add(i_data->force_j_y, mask_inner_0,
                                         i_data->force_j_y, fjy);
      i_data->force_j_z = fvec::mask_add(i_data->force_j_z, mask_inner_0,
                                         i_data->force_j_z, fjz);
      i_data->force_k_x_buf[buf_idx_i] =
        fvec::mask_add(i_data->force_k_x_buf[buf_idx_i], mask_inner_0,
                       i_data->force_k_x_buf[buf_idx_i], fkx);
      i_data->force_k_y_buf[buf_idx_i] =
        fvec::mask_add(i_data->force_k_y_buf[buf_idx_i], mask_inner_0,
                       i_data->force_k_y_buf[buf_idx_i], fky);
      i_data->force_k_z_buf[buf_idx_i] =
        fvec::mask_add(i_data->force_k_z_buf[buf_idx_i], mask_inner_0,
                       i_data->force_k_z_buf[buf_idx_i], fkz);
      j_data->force_k_x_buf[buf_idx_j] =
        fvec::mask_add(j_data->force_k_x_buf[buf_idx_j], mask_inner_0,
                       j_data->force_k_x_buf[buf_idx_j], flx);
      j_data->force_k_y_buf[buf_idx_j] =
        fvec::mask_add(j_data->force_k_y_buf[buf_idx_j], mask_inner_0,
                       j_data->force_k_y_buf[buf_idx_j], fly);
      j_data->force_k_z_buf[buf_idx_j] =
        fvec::mask_add(j_data->force_k_z_buf[buf_idx_j], mask_inner_0,
                       j_data->force_k_z_buf[buf_idx_j], flz);
    }
  }
}

/*
 * Processes VL elements of the same type itype/jtype for REBO and TORSION
 * interactions. This allows us to reuse the aut_frebo_data buffes in the
 * torsion calculaltion.
 */
static void aut_frebo_batch_of_kind(KernelArgsAIREBOT<flt_t,acc_t> * ka,
                                    int torflag, int itype, int jtype,
                                    int * i_buf, int * j_buf) {
 { // jump-scope for exceed_limits
  AtomAIREBOT<flt_t> * x = ka->x;
  int * map = ka->map;
  ResultForceT<acc_t> * result_f = ka->result_f;
  flt_t rcminij = ka->params.rcmin[itype][jtype];
  flt_t rcmaxij = ka->params.rcmax[itype][jtype];
  flt_t Qij = ka->params.Q[itype][jtype];
  flt_t Aij = ka->params.A[itype][jtype];
  flt_t alphaij = ka->params.alpha[itype][jtype];
  fvec vrcminij = fvec::set1(ka->params.rcmin[itype][jtype]);
  fvec vrcmaxij = fvec::set1(ka->params.rcmax[itype][jtype]);
  fvec vQij = fvec::set1(ka->params.Q[itype][jtype]);
  fvec vAij = fvec::set1(ka->params.A[itype][jtype]);
  fvec malphaij = fvec::set1(-ka->params.alpha[itype][jtype]);
  fvec c_1_0 = fvec::set1(1);
  fvec c_0_5 = fvec::set1(0.5);
  fvec c_TOL = fvec::set1(1e-9);
  struct aut_frebo_data i_data, j_data;

  fvec evdwl_vacc = fvec::setzero();
  ivec vi = ivec::maskz_loadu(bvec::full(), i_buf);
  int tmp;
  ivec vj = ivec::maskz_loadu(bvec::full(), j_buf);
  fvec x_i, y_i, z_i;
  fvec x_j, y_j, z_j;
  aut_loadatoms_vec_notype(x, vi, &x_i, &y_i, &z_i);
  aut_loadatoms_vec_notype(x, vj, &x_j, &y_j, &z_j);
  i_data.x_i = x_i;
  i_data.y_i = y_i;
  i_data.z_i = z_i;
  i_data.x_j = x_j;
  i_data.y_j = y_j;
  i_data.z_j = z_j;
  j_data.x_i = x_j;
  j_data.y_i = y_j;
  j_data.z_i = z_j;
  j_data.x_j = x_i;
  j_data.y_j = y_i;
  j_data.z_j = z_i;
  fvec delx = x_i -  x_j;
  fvec dely = y_i -  y_j;
  fvec delz = z_i -  z_j;
  fvec rsq = delx *  delx +  dely *  dely +  delz *  delz;
  fvec rij = fvec::sqrt(rsq);
  fvec dwij;
  fvec wij = aut_Sp_deriv(rij, vrcminij, vrcmaxij, &dwij);

  fvec exp_alphar = fvec::exp(malphaij *  rij);
  fvec Qij_over_rij = vQij /  rij;
  fvec Qij_over_rsq = vQij /  rsq;
  fvec VR_by_wij = ( c_1_0 +  Qij_over_rij) *  vAij *  exp_alphar;
  fvec VR = wij * VR_by_wij;
  fvec pre = wij *  vAij *  exp_alphar;
  fvec dVRdi = pre * ( malphaij +  malphaij *  Qij_over_rij -  Qij_over_rsq);
  dVRdi = dVRdi + VR_by_wij *  dwij;

  fvec VA_by_wij = fvec::setzero();
  fvec dVA = fvec::setzero();

  int k;
  for (k = 0; k < 3; k++) {
    fvec mBIJc = fvec::set1(-ka->params.BIJc[itype][jtype][k]);
    fvec mBetaij = fvec::set1(-ka->params.Beta[itype][jtype][k]);
    fvec term = mBIJc *  fvec::exp(mBetaij *  rij);
    VA_by_wij = VA_by_wij +  term;
    dVA = dVA +  mBetaij * wij * term;
  }

  dVA = dVA +  dwij *  VA_by_wij;
  fvec VA = wij * VA_by_wij;

  bvec tol_check = fvec::cmplt(wij, c_TOL);
  VA = fvec::mask_blend(tol_check, VA, fvec::setzero());
  dVA = fvec::mask_blend(tol_check, dVA, fvec::setzero());
  VR = fvec::mask_blend(tol_check, VR, fvec::setzero());
  dVRdi = fvec::mask_blend(tol_check, dVRdi, fvec::setzero());

  fvec nHi = fvec::gather(vi, ka->nH, sizeof(flt_t));
  fvec nCi = fvec::gather(vi, ka->nC, sizeof(flt_t));
  fvec nHj = fvec::gather(vj, ka->nH, sizeof(flt_t));
  fvec nCj = fvec::gather(vj, ka->nC, sizeof(flt_t));
  fvec Nij = (nHi +  nCi) -  wij;
  fvec Nji = (nHj +  nCj) -  wij;
  i_data.nHi = nHi;
  i_data.nCi = nCi;
  j_data.nHi = nHj;
  j_data.nCi = nCj;
  fvec fij[3], fji[3];
  fij[0] = fvec::setzero(); fij[1] = fvec::setzero();
  fij[2] = fvec::setzero();
  fji[0] = fvec::setzero(); fji[1] = fvec::setzero();
  fji[2] = fvec::setzero();

  fvec NconjtmpI;
  fvec pij = aut_frebo_pij_pd_2(
      ka, &i_data, itype, jtype, vi, vj,
      delx, dely, delz, rij, wij, VA, &NconjtmpI, fij);

  if (i_data.buf_len < 0) goto exceed_limits;

  fvec NconjtmpJ;
  fvec rjix = fvec::setzero() -  delx;
  fvec rjiy = fvec::setzero() -  dely;
  fvec rjiz = fvec::setzero() -  delz;
  fvec pji = aut_frebo_pij_pd_2(
      ka, &j_data, jtype, itype, vj, vi,
      rjix, rjiy, rjiz, rij, wij, VA, &NconjtmpJ, fji);
  fij[0] = fij[0] -  fji[0];
  fij[1] = fij[1] -  fji[1];
  fij[2] = fij[2] -  fji[2];

  if (j_data.buf_len < 0) goto exceed_limits;

  if (torflag && itype == 0 && jtype == 0)
    aut_torsion_vec(ka, &i_data, &j_data, vi, vj, wij, dwij);

  fvec Nijconj = c_1_0 +  NconjtmpI *  NconjtmpI +  NconjtmpJ *  NconjtmpJ;
  fvec dN3[3];
  fvec pi_rc = aut_frebo_pi_rc_pd(ka, itype, jtype, Nij, Nji, Nijconj, dN3);
  aut_frebo_N_spline_force(ka, &i_data, itype, jtype, vi, vj, VA, dN3[0],
                           dN3[2], NconjtmpI);
  aut_frebo_N_spline_force(ka, &j_data, jtype, itype, vj, vi, VA, dN3[1],
                           dN3[2], NconjtmpJ);
  fvec pi_dh = aut_frebo_pi_dh(ka, &i_data, &j_data, itype, jtype, vi, vj,
                               delx, dely, delz, rij, VA, Nij, Nji, Nijconj,
                               NconjtmpI, NconjtmpJ, fij);

  fvec bij = c_0_5 * ( pij +  pji) +  pi_rc +  pi_dh;
  fvec dVAdi = bij *  dVA;
  fvec fpair = (dVAdi +  dVRdi) *  fvec::recip(rij);
  fvec result_f_j_x = fpair *  delx -  fij[0];
  fvec result_f_j_y = fpair *  dely -  fij[1];
  fvec result_f_j_z = fpair *  delz -  fij[2];
  fvec result_f_i_x = fvec::setzero() -  result_f_j_x;
  fvec result_f_i_y = fvec::setzero() -  result_f_j_y;
  fvec result_f_i_z = fvec::setzero() -  result_f_j_z;
  fvec evdwl = VR +  bij *  VA;
  evdwl_vacc = evdwl_vacc +  evdwl;

  aut_frebo_data_writeback(ka, &i_data);
  aut_frebo_data_writeback(ka, &j_data);

  flt_t fi_x_buf[fvec::VL] __attribute__((aligned(64)));
  flt_t fi_y_buf[fvec::VL] __attribute__((aligned(64)));
  flt_t fi_z_buf[fvec::VL] __attribute__((aligned(64)));
  int fi_i_buf[ivec::VL] __attribute__((aligned(64)));
  flt_t fj_x_buf[fvec::VL] __attribute__((aligned(64)));
  flt_t fj_y_buf[fvec::VL] __attribute__((aligned(64)));
  flt_t fj_z_buf[fvec::VL] __attribute__((aligned(64)));
  int fj_j_buf[ivec::VL] __attribute__((aligned(64)));
  flt_t evdwl_buf[fvec::VL] __attribute__((aligned(64)));

  result_f_i_x = i_data.force_i_x +  result_f_i_x;
  result_f_i_y = i_data.force_i_y +  result_f_i_y;
  result_f_i_z = i_data.force_i_z +  result_f_i_z;
  result_f_j_x = i_data.force_j_x +  result_f_j_x;
  result_f_j_y = i_data.force_j_y +  result_f_j_y;
  result_f_j_z = i_data.force_j_z +  result_f_j_z;

  result_f_i_x = j_data.force_j_x +  result_f_i_x;
  result_f_i_y = j_data.force_j_y +  result_f_i_y;
  result_f_i_z = j_data.force_j_z +  result_f_i_z;
  result_f_j_x = j_data.force_i_x +  result_f_j_x;
  result_f_j_y = j_data.force_i_y +  result_f_j_y;
  result_f_j_z = j_data.force_i_z +  result_f_j_z;

  fvec::store(fi_x_buf, result_f_i_x);
  fvec::store(fi_y_buf, result_f_i_y);
  fvec::store(fi_z_buf, result_f_i_z);
  ivec::store(fi_i_buf, vi);
  fvec::store(fj_x_buf, result_f_j_x);
  fvec::store(fj_y_buf, result_f_j_y);
  fvec::store(fj_z_buf, result_f_j_z);
  ivec::store(fj_j_buf, vj);
  fvec::store(evdwl_buf, evdwl);

  int lane;
  for (lane = 0; lane < fvec::VL; lane++) {
    int ii = fi_i_buf[lane];
    result_f[ii].x += fi_x_buf[lane];
    result_f[ii].y += fi_y_buf[lane];
    result_f[ii].z += fi_z_buf[lane];
    result_f[ii].w += 0.5 * evdwl_buf[lane];
    int jj = fj_j_buf[lane];
    result_f[jj].x += fj_x_buf[lane];
    result_f[jj].y += fj_y_buf[lane];
    result_f[jj].z += fj_z_buf[lane];
    result_f[jj].w += 0.5 * evdwl_buf[lane];
  }
  ka->result_eng += fvec::reduce_add(evdwl_vacc);
  return;
 }
exceed_limits:
  for (int l = 0; l < fvec::VL; l++) {
    int i = i_buf[l];
    int j = j_buf[l];
    ref_frebo_single_interaction(ka, i, j);
    if (torflag && itype == 0 && jtype == 0)
      ref_torsion_single_interaction(ka, i, j);
  }
}

/*
 Orders the interactions by itype and jtype and passes chunks to the above
 method.
*/
static void aut_frebo(KernelArgsAIREBOT<flt_t,acc_t> * ka, int torflag) {
  AtomAIREBOT<flt_t> * _noalias x = ka->x;
  tagint * _noalias tag = ka->tag;
  int * _noalias map = ka->map;
  int i_buf[2][2][fvec::VL];
  int j_buf[2][2][fvec::VL];
  int n_buf[2][2] = {0};
  for (int i = ka->frebo_from_atom; i < ka->frebo_to_atom; i++) {
    tagint itag = tag[i];
    int itype = map[x[i].w];
    flt_t x_i = x[i].x;
    flt_t y_i = x[i].y;
    flt_t z_i = x[i].z;
    int * neighs = ka->neigh_rebo.entries + ka->neigh_rebo.offset[i];
    int jnum = ka->neigh_rebo.num[i];
    for (int jj = 0; jj < jnum; jj++) {
      int j = neighs[jj];
      tagint jtag = tag[j];
      if (itag > jtag) {
        if (((itag + jtag) & 1) == 0)
          continue;
      } else if (itag < jtag) {
        if (((itag + jtag) & 1) == 1)
          continue;
      } else {
        if (x[j].z < z_i)
          continue;
        if (x[j].z == z_i && x[j].y < y_i)
          continue;
        if (x[j].z == z_i && x[j].y == y_i && x[j].x < x_i)
          continue;
      }
      int jtype = map[x[j].w];
      int ins = n_buf[itype][jtype];
      i_buf[itype][jtype][ins] = i;
      j_buf[itype][jtype][ins] = j;
      n_buf[itype][jtype] += 1;
      if (n_buf[itype][jtype] == fvec::VL) {
        aut_frebo_batch_of_kind(ka, torflag, itype, jtype,
            i_buf[itype][jtype], j_buf[itype][jtype]);
        n_buf[itype][jtype] = 0;
      }
    }
  }
  for (int itype = 0; itype < 2; itype++) {
    for (int jtype = 0; jtype < 2; jtype++) {
      for (int l = 0; l < n_buf[itype][jtype]; l++) {
        int i = i_buf[itype][jtype][l];
        int j = j_buf[itype][jtype][l];
        ref_frebo_single_interaction(ka, i, j);
        if (torflag && itype == 0 && jtype == 0)
          ref_torsion_single_interaction(ka, i, j);
      }
    }
  }
}

/*
 * Apply paths in scalar fashion, not crucial for performance.
 */
static void aut_airebo_lj_force_path(KernelArgsAIREBOT<flt_t,acc_t> * ka,
   bvec mask, fvec dC, LennardJonesPathAIREBOT<flt_t> path[fvec::VL]) {
  for (int i = 0; i < fvec::VL; i++) {
    if (bvec::test_at(mask, i)) {
      ref_lennard_jones_force_path(ka, fvec::at(dC, i), &path[i]);
    }
  }
}

/*
 * Hash-Map for efficient calculation of C_ij.
 * Can have up to ITEMS entries with associated paths, as well as
 * 1024 entries. Open addressing, invalidation by using a different i.
 * Only needs to be reset once per timestep.
 */
static const int OPT_TEST_PATH_SIZE = 1024;
static const int OPT_TEST_PATH_ITEMS = 128;
struct aut_airebo_lj_test_path_result_data {
  LennardJonesPathAIREBOT<flt_t> testpath[OPT_TEST_PATH_ITEMS];
  int i[OPT_TEST_PATH_SIZE];
  int j[OPT_TEST_PATH_SIZE];
  flt_t cij[OPT_TEST_PATH_SIZE];
  int testpath_idx[OPT_TEST_PATH_SIZE];
};
static const unsigned int OPT_TEST_PATH_HASH = 2654435761;

static int aut_lj_tap_hash_fn(int j, int attempt) {
  uint32_t result = j;
  result *= (uint32_t) OPT_TEST_PATH_HASH;
  result += (uint32_t) attempt;
  result %= (uint32_t) OPT_TEST_PATH_SIZE;
  return result;
}

static ivec aut_airebo_lj_tap_hash_fn_vec(ivec val, ivec attempt) {
  const ivec golden = ivec::set1(OPT_TEST_PATH_HASH);
  const ivec mask = ivec::set1(OPT_TEST_PATH_SIZE - 1);
  ivec a = ivec::mullo(golden, val);
  ivec b = a +  attempt;
  ivec c = ivec::the_and(b, mask);
  return c;
}

/*
 * Enter all those (potential) neighbors of i (including 2nd and 3rd degree)
 * into the hash-map. There is no good way to vectorize this, and it does not
 * seem time-critical.
 */
static bool aut_airebo_lj_test_all_paths(KernelArgsAIREBOT<flt_t,acc_t> * ka,
    int i, struct aut_airebo_lj_test_path_result_data * result) {
  AtomAIREBOT<flt_t> * x = ka->x;
  int * map = ka->map;
  flt_t (*rcmin)[2] = &ka->params.rcmin[0];
  flt_t (*rcmax)[2] = &ka->params.rcmax[0];
  flt_t rcminsq[2][2];
  rcminsq[0][0] = rcmin[0][0] * rcmin[0][0];
  rcminsq[0][1] = rcmin[0][1] * rcmin[0][1];
  rcminsq[1][0] = rcmin[1][0] * rcmin[1][0];
  rcminsq[1][1] = rcmin[1][1] * rcmin[1][1];
  int * neighs_i = &ka->neigh_rebo.entries[ka->neigh_rebo.offset[i]];
  int itype = map[x[i].w];
  int path_insert_pos = 0;
  for (int jj = 0; jj < ka->neigh_rebo.num[i]; jj++) {
    int j = neighs_i[jj];
    int jtype = map[x[j].w];
    flt_t dijx = x[j].x - x[i].x;
    flt_t dijy = x[j].y - x[i].y;
    flt_t dijz = x[j].z - x[i].z;
    flt_t rijsq = dijx * dijx + dijy * dijy + dijz * dijz;
    flt_t wj = 1, dwj = 0;
    flt_t rij = 0;
    if (rijsq >= rcminsq[itype][jtype]) {
      rij = overloaded::sqrt(rijsq);
      wj = Sp(rij, rcmin[itype][jtype], rcmax[itype][jtype], &dwj);
    }
    int attempt = 0;
    int start_hash_slot = aut_lj_tap_hash_fn(j, attempt);
    int hash_slot = start_hash_slot;
    while (result->i[hash_slot] == i && result->j[hash_slot] != j &&
           attempt < OPT_TEST_PATH_SIZE) {
      hash_slot = aut_lj_tap_hash_fn(j, ++attempt);
    }
    if (attempt >= OPT_TEST_PATH_SIZE) goto exceed_limits;
    bool init_slot = result->i[hash_slot] != i;
    if (init_slot || (1 - wj < result->cij[hash_slot])) {
      result->i[hash_slot] = i;
      result->j[hash_slot] = j;
      result->cij[hash_slot] = 1 - wj;
      if (wj != 1.0) {
        if (path_insert_pos >= OPT_TEST_PATH_ITEMS) goto exceed_limits;
        result->testpath_idx[hash_slot] = path_insert_pos;
        LennardJonesPathAIREBOT<flt_t> *path =
          &result->testpath[path_insert_pos++];
        path->num = 2;
        path->del[0].x = dijx;
        path->del[0].y = dijy;
        path->del[0].z = dijz;
        if (rij == 0) rij = sqrt(rijsq);
        path->r[0] = rij;
        path->w[0] = wj;
        path->dw[0] = dwj;
        path->idx[0] = i;
        path->idx[1] = j;
      }
    }
    int * neighs_j = &ka->neigh_rebo.entries[ka->neigh_rebo.offset[j]];
    for (int kk = 0; kk < ka->neigh_rebo.num[j]; kk++) {
      int k = neighs_j[kk];
      if (k == i) continue;
      int ktype = map[x[k].w];
      flt_t djkx = x[k].x - x[j].x;
      flt_t djky = x[k].y - x[j].y;
      flt_t djkz = x[k].z - x[j].z;
      flt_t rjksq = djkx * djkx + djky * djky + djkz * djkz;
      flt_t wk = 1, dwk = 0;
      flt_t rjk = 0;
      if (rjksq >= rcminsq[jtype][ktype]) {
        rjk = overloaded::sqrt(rjksq);
        wk = Sp(rjk, rcmin[jtype][ktype], rcmax[jtype][ktype], &dwk);
      }
      int attempt = 0;
      int start_hash_slot = aut_lj_tap_hash_fn(k, attempt);
      int hash_slot = start_hash_slot;
      while (result->i[hash_slot] == i && result->j[hash_slot] != k &&
             attempt < OPT_TEST_PATH_SIZE) {
        hash_slot = aut_lj_tap_hash_fn(k, ++attempt);
      }
      if (attempt >= OPT_TEST_PATH_SIZE) goto exceed_limits;
      bool init_slot = result->i[hash_slot] != i;
      if (init_slot || (1 - wj * wk < result->cij[hash_slot])) {
        result->i[hash_slot] = i;
        result->j[hash_slot] = k;
        result->cij[hash_slot] = 1 - wj * wk;
        if (wj * wk != 1.0) {
          if (path_insert_pos >= OPT_TEST_PATH_ITEMS) goto exceed_limits;
          result->testpath_idx[hash_slot] = path_insert_pos;
          LennardJonesPathAIREBOT<flt_t> *path =
            &result->testpath[path_insert_pos++];
          path->num = 3;
          path->del[0].x = dijx;
          path->del[0].y = dijy;
          path->del[0].z = dijz;
          if (rij == 0) rij = sqrt(rijsq);
          path->r[0] = rij;
          path->del[1].x = djkx;
          path->del[1].y = djky;
          path->del[1].z = djkz;
          if (rjk == 0) rjk = sqrt(rjksq);
          path->r[1] = rjk;
          path->w[0] = wj;
          path->dw[0] = dwj;
          path->w[1] = wk;
          path->dw[1] = dwk;
          path->idx[0] = i;
          path->idx[1] = j;
          path->idx[2] = k;
        }
      }
      int * neighs_k = &ka->neigh_rebo.entries[ka->neigh_rebo.offset[k]];
      for (int ll = 0; ll < ka->neigh_rebo.num[k]; ll++) {
        int l = neighs_k[ll];
        if ((l == i) || (l == j)) continue;
        int ltype = map[x[l].w];
        flt_t dklx = x[l].x - x[k].x;
        flt_t dkly = x[l].y - x[k].y;
        flt_t dklz = x[l].z - x[k].z;
        flt_t rklsq = dklx * dklx + dkly * dkly + dklz * dklz;
        flt_t wl = 1, dwl = 0;
        flt_t rkl = 0;
        if (rklsq >= rcminsq[ktype][ltype]) {
          rkl = overloaded::sqrt(rklsq);
          wl = Sp(rkl, rcmin[ktype][ltype], rcmax[ktype][ltype], &dwl);
        }
        int attempt = 0;
        int start_hash_slot = aut_lj_tap_hash_fn(l, attempt);
        int hash_slot = start_hash_slot;
        while (result->i[hash_slot] == i && result->j[hash_slot] != l &&
               attempt < OPT_TEST_PATH_SIZE) {
          hash_slot = aut_lj_tap_hash_fn(l, ++attempt);
        }
        if (attempt >= OPT_TEST_PATH_SIZE) goto exceed_limits;
        bool init_slot = result->i[hash_slot] != i;
        if (init_slot || (1 - wj * wk * wl < result->cij[hash_slot])) {
          result->i[hash_slot] = i;
          result->j[hash_slot] = l;
          result->cij[hash_slot] = 1 - wj * wk * wl;
          if (wj * wk * wl != 1.0) {
            if (path_insert_pos >= OPT_TEST_PATH_ITEMS) goto exceed_limits;
            result->testpath_idx[hash_slot] = path_insert_pos;
            LennardJonesPathAIREBOT<flt_t> *path =
              &result->testpath[path_insert_pos++];
            path->num = 4;
            path->del[0].x = dijx;
            path->del[0].y = dijy;
            path->del[0].z = dijz;
            if (rij == 0) rij = sqrt(rijsq);
            path->r[0] = rij;
            path->del[1].x = djkx;
            path->del[1].y = djky;
            path->del[1].z = djkz;
            if (rjk == 0) rjk = sqrt(rjksq);
            path->r[1] = rjk;
            path->del[2].x = dklx;
            path->del[2].y = dkly;
            path->del[2].z = dklz;
            if (rkl == 0) rkl = sqrt(rklsq);
            path->r[2] = rkl;
            path->w[0] = wj;
            path->dw[0] = dwj;
            path->w[1] = wk;
            path->dw[1] = dwk;
            path->w[2] = wl;
            path->dw[2] = dwl;
            path->idx[0] = i;
            path->idx[1] = j;
            path->idx[2] = k;
            path->idx[3] = l;
          }
        }
      }
    }
  }
  return true;
exceed_limits:
  return false;
}

/*
 * Attempt to look up an element in the hash-map.
 */
static fvec aut_airebo_lj_tap_test_path(KernelArgsAIREBOT<flt_t,acc_t> * /*ka*/,
  struct aut_airebo_lj_test_path_result_data * test_path_result,
  bvec need_search, ivec i_bc, ivec j,
  LennardJonesPathAIREBOT<flt_t> path[fvec::VL]
) {
  const ivec c_i1 = ivec::set1(1);
  fvec cij = fvec::set1(1.0);
  // first round: hash all j
  // lookup i/j in hash list.
  // if i matches and j matches: congrats
  // if i matches and j does not: look up attempts
  // if attempts > current_attempts:
  //   do another round of hashing
  // for all those found:

  //   fill in the path
  // -----------------------------------------------
  // find all the correct hash slots, and a mask of where found.
  ivec attempt = ivec::setzero();
  ivec hash_slot = aut_airebo_lj_tap_hash_fn_vec(j, attempt);
  ivec lookup_i = ivec::mask_gather(ivec::undefined(), need_search, hash_slot,
      &test_path_result->i[0], sizeof(int));
  bvec correct_i = ivec::mask_cmpeq(need_search, lookup_i, i_bc);
  ivec lookup_j = ivec::mask_gather(ivec::undefined(), correct_i, hash_slot,
      &test_path_result->j[0], sizeof(int));
  bvec found_items = ivec::mask_cmpeq(correct_i, lookup_j, j);
  bvec another_attempt = correct_i & ~ found_items;
  while (bvec::test_any_set(another_attempt)) {
    attempt = ivec::mask_add(attempt, another_attempt, attempt, c_i1);
    hash_slot = aut_airebo_lj_tap_hash_fn_vec(j, attempt);
    ivec lookup_i_2 = ivec::mask_gather(lookup_i, another_attempt, hash_slot,
        &test_path_result->i[0], sizeof(int));
    lookup_i = lookup_i_2;
    correct_i = ivec::mask_cmpeq(need_search, lookup_i, i_bc);
    lookup_j = ivec::mask_gather(lookup_j, another_attempt, hash_slot,
        &test_path_result->j[0], sizeof(int));
    found_items = ivec::mask_cmpeq(correct_i, lookup_j, j);
    another_attempt = correct_i & ~ found_items;
  }
  cij = fvec::mask_gather(cij, found_items, hash_slot,
                          &test_path_result->cij[0], sizeof(flt_t));
  bvec need_testpath = fvec::mask_cmplt(found_items, fvec::setzero(), cij);
  if (bvec::test_any_set(need_testpath)) {
    for (int i = 0; i < fvec::VL; i++) {
      if (bvec::test_at(need_testpath, i)) {
        int testpath_idx =
          test_path_result->testpath_idx[ivec::at(hash_slot, i)];
        path[i] = test_path_result->testpath[testpath_idx];
      }
    }
  }
  return cij;
}

/*
 * This function calculates the Lennard-Jones interaciton for those
 * elements that require a bond-order calculation.
 * It is similarly structured as the aut_frebo_batch_of_kind function.
 * The forces due to bondorders are calculated speculatively and later
 * updated with the correct outer derivative.
 */
template<int MORSEFLAG>
static void aut_lj_with_bo(
    KernelArgsAIREBOT<flt_t,acc_t> * ka,
    int itype, int jtype,
    ivec i, ivec j,
    fvec cij, LennardJonesPathAIREBOT<flt_t> testpath[fvec::VL]
) {
 { // jump-scope for exceed_limits
  AtomAIREBOT<flt_t> * _noalias x = ka->x;
  ResultForceT<acc_t> * result_f = ka->result_f;

  ivec c_i4 = ivec::set1(4);
  fvec c_1_0 = fvec::set1(1.0);
  fvec c_2_0 = fvec::set1(2.0);
  fvec c_0_5 = fvec::set1(0.5);

  fvec x_i, y_i, z_i;
  aut_loadatoms_vec_notype(x, i, &x_i, &y_i, &z_i);
  fvec x_j, y_j, z_j;
  aut_loadatoms_vec_notype(x, j, &x_j, &y_j, &z_j);
  fvec delx = x_i -  x_j;
  fvec dely = y_i -  y_j;
  fvec delz = z_i -  z_j;
  fvec rsq = delx *  delx +  dely *  dely +  delz *  delz;

  fvec rij = fvec::sqrt(rsq);
  bvec need_path_force = fvec::cmplt(cij, c_1_0);
  flt_t sigcut = ka->params.sigcut;
  flt_t sigmin = ka->params.sigmin;
  flt_t sigma = ka->params.sigma[itype][jtype];
  flt_t rljmax = sigcut * sigma;
  flt_t rljmin = sigmin * sigma;
  fvec p_rljmin = fvec::set1(rljmin);
  fvec p_rljmax = fvec::set1(rljmax);

  fvec dslw, slw = aut_Sp2_deriv(rij, p_rljmin, p_rljmax, &dslw);

  fvec p_lj1 = fvec::set1(ka->params.lj1[itype][jtype]);
  fvec p_lj2 = fvec::set1(ka->params.lj2[itype][jtype]);
  fvec p_lj3 = fvec::set1(ka->params.lj3[itype][jtype]);
  fvec p_lj4 = fvec::set1(ka->params.lj4[itype][jtype]);

  fvec r2inv = fvec::recip(rsq);

  fvec vdw, dvdw;
  if (MORSEFLAG) {
    fvec exr = fvec::exp(fvec::setzero() - rij * p_lj4);
    vdw = p_lj1 * exr * (p_lj2 * exr - c_2_0);
    dvdw = p_lj3 * exr * (c_1_0 - p_lj2 * exr);
  } else {
    fvec r6inv = r2inv *  r2inv *  r2inv;

    vdw = r6inv * ( p_lj3 *  r6inv -  p_lj4);
    fvec r7inv = r6inv *  rij *  r2inv;
    dvdw = r7inv * ( p_lj2 -  p_lj1 *  r6inv);
  }

  fvec VLJ = vdw *  slw;
  fvec dVLJ = dvdw *  slw +  vdw *  dslw;

  fvec p_rcLJmin = fvec::set1(ka->params.rcLJmin[itype][jtype]);
  fvec p_rcLJmax = fvec::set1(ka->params.rcLJmax[itype][jtype]);
  fvec dStr, Str = aut_Sp2_deriv(rij, p_rcLJmin, p_rcLJmax, &dStr);
  fvec VA = cij *  VLJ *  Str;

  fvec fij[3], fji[3];
  fij[0] = fvec::setzero(); fij[1] = fvec::setzero();
  fij[2] = fvec::setzero();
  fji[0] = fvec::setzero(); fji[1] = fvec::setzero();
  fji[2] = fvec::setzero();

  ivec vi = i;
  ivec vj = j;

  struct aut_frebo_data i_data, j_data;
  i_data.x_i = x_i;
  i_data.y_i = y_i;
  i_data.z_i = z_i;
  i_data.x_j = x_j;
  i_data.y_j = y_j;
  i_data.z_j = z_j;
  j_data.x_i = x_j;
  j_data.y_i = y_j;
  j_data.z_i = z_j;
  j_data.x_j = x_i;
  j_data.y_j = y_i;
  j_data.z_j = z_i;

  fvec p_rcmin = fvec::set1(ka->params.rcmin[itype][jtype]);
  fvec p_rcmax = fvec::set1(ka->params.rcmax[itype][jtype]);
  fvec dwij;
  fvec wij = aut_Sp_deriv(rij, p_rcmin, p_rcmax, &dwij);

  fvec nHi = fvec::gather(vi, ka->nH, sizeof(flt_t));
  fvec nCi = fvec::gather(vi, ka->nC, sizeof(flt_t));
  fvec nHj = fvec::gather(vj, ka->nH, sizeof(flt_t));
  fvec nCj = fvec::gather(vj, ka->nC, sizeof(flt_t));
  fvec Nij = nHi +  nCi -  wij;
  fvec Nji = nHj +  nCj -  wij;
  i_data.nHi = nHi;
  i_data.nCi = nCi;
  j_data.nHi = nHj;
  j_data.nCi = nCj;

  fvec the_r = fvec::set1(ka->params.rcmin[itype][jtype]);
  fvec scale = the_r / rij;

  fvec NconjtmpI;
  fvec pij = aut_frebo_pij_pd_2(ka, &i_data, itype, jtype, vi, vj,
                                delx * scale, dely * scale, delz * scale,
                                the_r, wij, VA, &NconjtmpI, fij);

  if (i_data.buf_len < 0) goto exceed_limits;

  fvec NconjtmpJ;
  fvec rjix = fvec::setzero() -  delx;
  fvec rjiy = fvec::setzero() -  dely;
  fvec rjiz = fvec::setzero() -  delz;
  fvec pji = aut_frebo_pij_pd_2(ka, &j_data, jtype, itype, vj, vi,
                                rjix * scale, rjiy * scale, rjiz * scale,
                                the_r, wij, VA, &NconjtmpJ, fji);
  fij[0] = fij[0] -  fji[0];
  fij[1] = fij[1] -  fji[1];
  fij[2] = fij[2] -  fji[2];

  if (j_data.buf_len < 0) goto exceed_limits;

  fvec Nijconj = c_1_0 +  NconjtmpI *  NconjtmpI +  NconjtmpJ *  NconjtmpJ;
  fvec dN3[3];
  fvec pi_rc = aut_frebo_pi_rc_pd(ka, itype, jtype, Nij, Nji, Nijconj, dN3);

  fvec c_TOL = fvec::set1(TOL);
  fvec dN3_dh[3];
  fvec Tij = aut_frebo_Tij(ka, itype, jtype, Nij, Nji, Nijconj, &dN3_dh[0]);
  bvec TijgtTOLmask = fvec::cmpnle(fvec::abs(Tij), c_TOL);
  fvec sum_omega = fvec::setzero();
  if (bvec::test_any_set(TijgtTOLmask)) {
    sum_omega = aut_frebo_sum_omega(
        ka, &i_data, &j_data, itype, jtype, vi, vj,
        delx * scale, dely * scale, delz * scale, the_r, VA *  Tij, fij);
    sum_omega = fvec::mask_blend(TijgtTOLmask, fvec::setzero(), sum_omega);
  }
  fvec pi_dh = Tij *  sum_omega;

  fvec bij = c_0_5 * ( pij +  pji) + pi_rc +  pi_dh;

  fvec p_bLJmin = fvec::set1(ka->params.bLJmin[itype][jtype]);
  fvec p_bLJmax = fvec::set1(ka->params.bLJmax[itype][jtype]);
  fvec dStb, Stb = aut_Sp2_deriv(bij, p_bLJmin, p_bLJmax, &dStb);

  bvec need_bo_deriv = fvec::cmpneq(dStb, fvec::setzero());
  // fix up j_data, i_data, fij:
  // multiply each by dStb
  if (bvec::test_any_set(need_bo_deriv)) {
    i_data.force_i_x = dStb * i_data.force_i_x;
    i_data.force_i_y = dStb * i_data.force_i_y;
    i_data.force_i_z = dStb * i_data.force_i_z;
    i_data.force_j_x = dStb * i_data.force_j_x;
    i_data.force_j_y = dStb * i_data.force_j_y;
    i_data.force_j_z = dStb * i_data.force_j_z;
    j_data.force_i_x = dStb * j_data.force_i_x;
    j_data.force_i_y = dStb * j_data.force_i_y;
    j_data.force_i_z = dStb * j_data.force_i_z;
    j_data.force_j_x = dStb * j_data.force_j_x;
    j_data.force_j_y = dStb * j_data.force_j_y;
    j_data.force_j_z = dStb * j_data.force_j_z;
    for (int k = 0; k < i_data.buf_len; k++) {
      i_data.force_k_x_buf[k] = dStb * i_data.force_k_x_buf[k];
      i_data.force_k_y_buf[k] = dStb * i_data.force_k_y_buf[k];
      i_data.force_k_z_buf[k] = dStb * i_data.force_k_z_buf[k];
    }
    for (int k = 0; k < j_data.buf_len; k++) {
      j_data.force_k_x_buf[k] = dStb * j_data.force_k_x_buf[k];
      j_data.force_k_y_buf[k] = dStb * j_data.force_k_y_buf[k];
      j_data.force_k_z_buf[k] = dStb * j_data.force_k_z_buf[k];
    }
    fvec fijc[3];
    fijc[0] = dStb * fij[0];
    fijc[1] = dStb * fij[1];
    fijc[2] = dStb * fij[2];
    fij[0] = scale * (fijc[0] - (delx * delx * fijc[0] + dely * delx *
                                 fijc[1] + delz * delx * fijc[2]) / rsq);
    fij[1] = scale * (fijc[1] - (delx * dely * fijc[0] + dely * dely *
                                 fijc[1] + delz * dely * fijc[2]) / rsq);
    fij[2] = scale * (fijc[2] - (delx * delz * fijc[0] + dely * delz *
                                 fijc[1] + delz * delz * fijc[2]) / rsq);

    aut_frebo_N_spline_force(ka, &i_data, itype, jtype, vi, vj, dStb * VA,
                             dN3[0], dN3[2], NconjtmpI);
    aut_frebo_N_spline_force(ka, &j_data, jtype, itype, vj, vi, dStb * VA,
                             dN3[1], dN3[2], NconjtmpJ);
    if (bvec::test_any_set(TijgtTOLmask)) {
      aut_frebo_N_spline_force(ka, &i_data, itype, jtype, vi, vj,
                               dStb * VA * sum_omega, dN3_dh[0], dN3_dh[2],
                               NconjtmpI);
      aut_frebo_N_spline_force(ka, &j_data, jtype, itype, vj, vi,
                               dStb * VA * sum_omega, dN3_dh[1], dN3_dh[2],
                               NconjtmpJ);
    }

    aut_frebo_data_writeback(ka, &i_data);
    aut_frebo_data_writeback(ka, &j_data);
  } else {
    fij[0] = fvec::setzero();
    fij[1] = fvec::setzero();
    fij[2] = fvec::setzero();
  }

  fvec fpdVLJ = cij *  dVLJ * ( c_1_0 +  Str * ( Stb -  c_1_0));
  fvec fpdStr = dStr *  cij * ( Stb *  VLJ -  VLJ);
  fvec fpair = r2inv *  rij * ( fvec::setzero() - ( fpdVLJ +  fpdStr));
  fvec evdwl = VA *  Stb +  cij *  VLJ * ( c_1_0 -  Str);

  fvec result_f_i_x = fpair *  delx +  fij[0];
  fvec result_f_i_y = fpair *  dely +  fij[1];
  fvec result_f_i_z = fpair *  delz +  fij[2];
  fvec result_f_j_x = fvec::setzero() -  result_f_i_x;
  fvec result_f_j_y = fvec::setzero() -  result_f_i_y;
  fvec result_f_j_z = fvec::setzero() -  result_f_i_z;

  flt_t fi_x_buf[fvec::VL] __attribute__((aligned(64)));
  flt_t fi_y_buf[fvec::VL] __attribute__((aligned(64)));
  flt_t fi_z_buf[fvec::VL] __attribute__((aligned(64)));
  int fi_i_buf[ivec::VL] __attribute__((aligned(64)));
  flt_t fj_x_buf[fvec::VL] __attribute__((aligned(64)));
  flt_t fj_y_buf[fvec::VL] __attribute__((aligned(64)));
  flt_t fj_z_buf[fvec::VL] __attribute__((aligned(64)));
  int fj_j_buf[ivec::VL] __attribute__((aligned(64)));
  flt_t evdwl_buf[fvec::VL] __attribute__((aligned(64)));

  if (bvec::test_any_set(need_bo_deriv)) {
    result_f_i_x = i_data.force_i_x +  result_f_i_x;
    result_f_i_y = i_data.force_i_y +  result_f_i_y;
    result_f_i_z = i_data.force_i_z +  result_f_i_z;
    result_f_j_x = i_data.force_j_x +  result_f_j_x;
    result_f_j_y = i_data.force_j_y +  result_f_j_y;
    result_f_j_z = i_data.force_j_z +  result_f_j_z;

    result_f_i_x = j_data.force_j_x +  result_f_i_x;
    result_f_i_y = j_data.force_j_y +  result_f_i_y;
    result_f_i_z = j_data.force_j_z +  result_f_i_z;
    result_f_j_x = j_data.force_i_x +  result_f_j_x;
    result_f_j_y = j_data.force_i_y +  result_f_j_y;
    result_f_j_z = j_data.force_i_z +  result_f_j_z;
  }

  fvec::store(fi_x_buf, result_f_i_x);
  fvec::store(fi_y_buf, result_f_i_y);
  fvec::store(fi_z_buf, result_f_i_z);
  ivec::store(fi_i_buf, vi);
  fvec::store(fj_x_buf, result_f_j_x);
  fvec::store(fj_y_buf, result_f_j_y);
  fvec::store(fj_z_buf, result_f_j_z);
  ivec::store(fj_j_buf, vj);
  fvec::store(evdwl_buf, evdwl);

  int lane;
  for (lane = 0; lane < fvec::VL; lane++) {
    int ii = fi_i_buf[lane];
    result_f[ii].x += fi_x_buf[lane];
    result_f[ii].y += fi_y_buf[lane];
    result_f[ii].z += fi_z_buf[lane];
    result_f[ii].w += 0.5 * evdwl_buf[lane];
    int jj = fj_j_buf[lane];
    result_f[jj].x += fj_x_buf[lane];
    result_f[jj].y += fj_y_buf[lane];
    result_f[jj].z += fj_z_buf[lane];
    result_f[jj].w += 0.5 * evdwl_buf[lane];
  }
  ka->result_eng += fvec::reduce_add(evdwl);

  if (bvec::test_any_set(need_path_force)) {
    fvec dC = VLJ * ( Str *  Stb +  c_1_0 -  Str);
    aut_airebo_lj_force_path(ka, need_path_force, dC, testpath);
  }
  return;
 }
exceed_limits:
  for (int l = 0; l < fvec::VL; l++) {
    ref_lennard_jones_single_interaction(ka, ivec::at(i, l), ivec::at(j, l),
                                         MORSEFLAG);
  }
  return;
}

/*
 * Calculate the lennard-jones interaction.
 * Uses the above hash-map, and outlines the calculation if the bondorder is
 *  needed.
 * Aggressively compresses to get the most values calculated.
 */
template<int MORSEFLAG>
static void aut_lennard_jones(KernelArgsAIREBOT<flt_t,acc_t> * ka) {
  AtomAIREBOT<flt_t> * x = ka->x;
  int * map = ka->map;
  ResultForceT<acc_t> * result_f = ka->result_f;
  ivec c_i1 = ivec::set1(1);
  ivec c_i4 = ivec::set1(4);
  fvec c_1_0 = fvec::set1(1.0);
  fvec c_2_0 = fvec::set1(2.0);
  fvec c_0_0 = fvec::set1(0.0);
  int map_i_scalar = 0;
  {
    int i;
    for (i = 1; i < ka->num_types; i++) {
      if (ka->map[i])
        map_i_scalar |= (1 << i);
    }
  }
  ivec map_i = ivec::set1(map_i_scalar);
  fvec result_eng = fvec::setzero();

  struct aut_airebo_lj_test_path_result_data test_path_result;
  for (int i = 0; i < OPT_TEST_PATH_SIZE; i++) {
    test_path_result.i[i] = -1;
  }

  ivec i_bo[2][2];
  ivec j_bo[2][2];
  fvec cij_bo[2][2];
  LennardJonesPathAIREBOT<flt_t> testpath_bo[2][2][fvec::VL];
  int num_bo[2][2] = {0};

  for (int i = ka->frebo_from_atom; i < ka->frebo_to_atom; i++) {
    int itype = map[x[i].w];
    fvec x_i = fvec::set1(x[i].x);
    fvec y_i = fvec::set1(x[i].y);
    fvec z_i = fvec::set1(x[i].z);
    ivec i_bc = ivec::set1(i);

    fvec cutljsq0 = fvec::set1(ka->params.cutljsq[itype][0]);
    fvec cutljsq1 = fvec::set1(ka->params.cutljsq[itype][1]);
    fvec p_rcmax0 = fvec::set1(ka->params.rcmax[itype][0]);
    fvec p_rcmax1 = fvec::set1(ka->params.rcmax[itype][1]);
    flt_t sigcut = ka->params.sigcut;
    flt_t sigmin = ka->params.sigmin;
    flt_t sigma0 = ka->params.sigma[itype][0];
    flt_t rljmax0 = sigcut * sigma0;
    flt_t rljmin0 = sigmin * sigma0;
    flt_t sigma1 = ka->params.sigma[itype][1];
    flt_t rljmax1 = sigcut * sigma1;
    flt_t rljmin1 = sigmin * sigma1;
    fvec p_rljmax0 = fvec::set1(rljmax0);
    fvec p_rljmax1 = fvec::set1(rljmax1);
    fvec p_rljmin0 = fvec::set1(rljmin0);
    fvec p_rljmin1 = fvec::set1(rljmin1);
    fvec p_rcLJmax0 = fvec::set1(ka->params.rcLJmax[itype][0]);
    fvec p_rcLJmax1 = fvec::set1(ka->params.rcLJmax[itype][1]);
    fvec p_rcLJmin0 = fvec::set1(ka->params.rcLJmin[itype][0]);
    fvec p_rcLJmin1 = fvec::set1(ka->params.rcLJmin[itype][1]);
    fvec p_lj10 = fvec::set1(ka->params.lj1[itype][0]);
    fvec p_lj11 = fvec::set1(ka->params.lj1[itype][1]);
    fvec p_lj20 = fvec::set1(ka->params.lj2[itype][0]);
    fvec p_lj21 = fvec::set1(ka->params.lj2[itype][1]);
    fvec p_lj30 = fvec::set1(ka->params.lj3[itype][0]);
    fvec p_lj31 = fvec::set1(ka->params.lj3[itype][1]);
    fvec p_lj40 = fvec::set1(ka->params.lj4[itype][0]);
    fvec p_lj41 = fvec::set1(ka->params.lj4[itype][1]);

    int * neighs = ka->neigh_lmp.entries[i];
    int jnum = ka->neigh_lmp.num_half[i];

    bool tap_success = aut_airebo_lj_test_all_paths(ka, i, &test_path_result);
    if (! tap_success) {
      for (int jj = 0; jj < jnum; jj++) {
        ref_lennard_jones_single_interaction(ka, i, neighs[jj], MORSEFLAG);
      }
      continue;
    }

    ivec j_2;
    fvec delx_2, dely_2, delz_2, rsq_2;
    bvec jtype_mask_2;
    int num_2 = 0;

    fvec result_f_i_x = fvec::setzero();
    fvec result_f_i_y = fvec::setzero();
    fvec result_f_i_z = fvec::setzero();

    int jj = 0;
    bool rest_j = jj < jnum;
    bool rest_2 = fvec::fast_compress();
    #pragma forceinline recursive
    while (rest_j || rest_2) {
      fvec delx, dely, delz, rsq;
      bvec jtype_mask, within_cutoff;
      ivec j;
      if (rest_j) {
        bvec mask_0 = bvec::full();
        //0xFF >> (8 - (jnum - jj));
        if (jj + (fvec::VL - 1) >= jnum) mask_0 = bvec::only(jnum - jj);
        j = ivec::maskz_loadu(mask_0, &neighs[jj]);
        fvec x_j, y_j, z_j;
        aut_loadatoms_vec(x, j, &x_j, &y_j, &z_j, &jtype_mask, map, map_i,
                          c_i1);
        fvec::gather_prefetch0(ivec::mullo(c_i4,
          ivec::maskz_loadu(bvec::full(), &neighs[jj + fvec::VL])), x);
        _mm_prefetch((const char*)&neighs[jj + 2 * fvec::VL], _MM_HINT_T0);
        delx = x_i -  x_j;
        dely = y_i -  y_j;
        delz = z_i -  z_j;
        rsq = delx *  delx +  dely *  dely +  delz *  delz;
        fvec cutoff_sq = fvec::mask_blend(jtype_mask, cutljsq0, cutljsq1);
        within_cutoff = fvec::mask_cmplt(mask_0, rsq, cutoff_sq);

        if (fvec::fast_compress()) {
          j = ivec::masku_compress(within_cutoff, j);
          delx = fvec::masku_compress(within_cutoff, delx);
          dely = fvec::masku_compress(within_cutoff, dely);
          delz = fvec::masku_compress(within_cutoff, delz);
          rsq = fvec::masku_compress(within_cutoff, rsq);
          jtype_mask = bvec::masku_compress(within_cutoff, jtype_mask);
          //within_cutoff = 0xFF >> (8 - _cc_popcnt(within_cutoff));

          bvec mask_2 = bvec::after(num_2);//0xFF << num_2;
          j_2 = ivec::mask_expand(j_2, mask_2, j);
          delx_2 = fvec::mask_expand(delx_2, mask_2, delx);
          dely_2 = fvec::mask_expand(dely_2, mask_2, dely);
          delz_2 = fvec::mask_expand(delz_2, mask_2, delz);
          rsq_2 = fvec::mask_expand(rsq_2, mask_2, rsq);
          jtype_mask_2 = bvec::mask_expand(jtype_mask_2, mask_2, jtype_mask);
          num_2 = num_2 + bvec::popcnt(within_cutoff);
          if (num_2 < fvec::VL) {
            jj += fvec::VL;
            rest_j = jj < jnum;
            continue;
          }

          num_2 -= fvec::VL;
          //(0xFF >> (8 - num_2)) << (_cc_popcnt(within_cutoff) - num_2);
          mask_2 = bvec::onlyafter(num_2, bvec::popcnt(within_cutoff) - num_2);
          {
            ivec tmp_j = j_2;
            j_2 = ivec::masku_compress(mask_2, j);
            j = tmp_j;
            fvec tmp_delx = delx_2;
            delx_2 = fvec::masku_compress(mask_2, delx);
            delx = tmp_delx;
            fvec tmp_dely = dely_2;
            dely_2 = fvec::masku_compress(mask_2, dely);
            dely = tmp_dely;
            fvec tmp_delz = delz_2;
            delz_2 = fvec::masku_compress(mask_2, delz);
            delz = tmp_delz;
            fvec tmp_rsq = rsq_2;
            rsq_2 = fvec::masku_compress(mask_2, rsq);
            rsq = tmp_rsq;
            bvec tmp_jtype_mask = jtype_mask_2;
            jtype_mask_2 = bvec::masku_compress(mask_2, jtype_mask);
            jtype_mask = tmp_jtype_mask;
            within_cutoff = bvec::full();
          }
        }
      } else if (rest_2) {
        rest_2 = false;
        j = j_2;
        delx = delx_2;
        dely = dely_2;
        delz = delz_2;
        rsq = rsq_2;
        jtype_mask = jtype_mask_2;
        within_cutoff = bvec::only(num_2);
        num_2 = 0;
      }

      bvec current_mask = within_cutoff;
      if (bvec::test_all_unset(current_mask)) {
        jj += fvec::VL;
        rest_j = jj < jnum;
        continue;
      }

      fvec rij = fvec::sqrt(rsq);
      LennardJonesPathAIREBOT<flt_t> testpath[fvec::VL];
      fvec cij = c_1_0;
      fvec p_cut3rebo = fvec::set1(ka->params.cut3rebo);
      bvec need_search = fvec::mask_cmplt(current_mask, rij, p_cut3rebo);
      if (bvec::test_any_set(need_search)) {
        fvec p_rcmax = fvec::mask_blend(jtype_mask, p_rcmax0, p_rcmax1);
        #pragma noinline
        cij = aut_airebo_lj_tap_test_path(ka, &test_path_result, need_search,
                                          i_bc, j, testpath);
      }
      current_mask = fvec::mask_cmplt(current_mask, c_0_0, cij);
      if (bvec::test_all_unset(current_mask)) {
        jj += fvec::VL;
        rest_j = jj < jnum;
        continue;
      }
      bvec need_path_force = fvec::mask_cmplt(current_mask, cij, c_1_0);

      fvec p_rljmax = fvec::mask_blend(jtype_mask, p_rljmax0, p_rljmax1);
      fvec p_rljmin = fvec::mask_blend(jtype_mask, p_rljmin0, p_rljmin1);

      fvec dslw, slw = aut_Sp2_deriv(rij, p_rljmin, p_rljmax, &dslw);

      fvec p_lj1 = fvec::mask_blend(jtype_mask, p_lj10, p_lj11);
      fvec p_lj2 = fvec::mask_blend(jtype_mask, p_lj20, p_lj21);
      fvec p_lj3 = fvec::mask_blend(jtype_mask, p_lj30, p_lj31);
      fvec p_lj4 = fvec::mask_blend(jtype_mask, p_lj40, p_lj41);

      fvec vdw, dvdw;

      fvec r2inv = fvec::recip(rsq);

      if (MORSEFLAG) {
        fvec exr = fvec::exp(fvec::setzero() - rij * p_lj4);
        vdw = p_lj1 * exr * (p_lj2 * exr - c_2_0);
        dvdw = p_lj3 * exr * (c_1_0 - p_lj2 * exr);
      } else {
        fvec r6inv = r2inv *  r2inv *  r2inv;

        vdw = r6inv * ( p_lj3 *  r6inv -  p_lj4);
        fvec r7inv = r6inv *  rij *  r2inv;
        dvdw = r7inv * ( p_lj2 -  p_lj1 *  r6inv);
      }

      fvec VLJ = vdw *  slw;
      fvec dVLJ = dvdw *  slw +  vdw *  dslw;

      fvec p_rcLJmin = fvec::mask_blend(jtype_mask, p_rcLJmin0, p_rcLJmin1);
      fvec p_rcLJmax = fvec::mask_blend(jtype_mask, p_rcLJmax0, p_rcLJmax1);
      fvec dStr, Str = aut_Sp2_deriv(rij, p_rcLJmin, p_rcLJmax, &dStr);
      fvec VA = cij *  VLJ *  Str;
      bvec need_bondorder = fvec::mask_cmplt(current_mask, c_0_0, Str);
      fvec Stb = fvec::setzero();
      fvec fij[3];
      fij[0] = fvec::setzero();
      fij[1] = fvec::setzero();
      fij[2] = fvec::setzero();
      if (bvec::test_any_set(need_bondorder)) {
        for (int jtype = 0; jtype < 2; jtype++) {
          bvec need_bo_with_jtype = need_bondorder;
          if (jtype) need_bo_with_jtype = need_bo_with_jtype & jtype_mask;
          else need_bo_with_jtype = need_bo_with_jtype & ~ jtype_mask;
          ivec jtmp = ivec::masku_compress(need_bo_with_jtype, j);
          ivec itmp = ivec::masku_compress(need_bo_with_jtype, ivec::set1(i));
          fvec cijtmp = fvec::masku_compress(need_bo_with_jtype, cij);
          bvec insert_mask = bvec::after(num_bo[itype][jtype]);
          i_bo[itype][jtype] = ivec::mask_expand(i_bo[itype][jtype],
                                                 insert_mask, itmp);
          j_bo[itype][jtype] = ivec::mask_expand(j_bo[itype][jtype],
                                                 insert_mask, jtmp);
          cij_bo[itype][jtype] = fvec::mask_expand(cij_bo[itype][jtype],
                                                   insert_mask, cijtmp);
          bvec need_path_force_with_jtype = need_bo_with_jtype &
            need_path_force;
          int testpath_end = fvec::VL;
          if (bvec::test_any_set(need_path_force_with_jtype)) {
            int pos = num_bo[itype][jtype];
            for (int l = 0; l < fvec::VL; l++) {
              if (pos >= fvec::VL) {
                testpath_end = l;
                break;
              }
              if (bvec::test_at(need_path_force_with_jtype, l)) {
                testpath_bo[itype][jtype][pos] = testpath[l];
              }
              if (bvec::test_at(need_bo_with_jtype, l)) {
                pos += 1;
              }
            }
          }
          num_bo[itype][jtype] = num_bo[itype][jtype] +
            bvec::popcnt(need_bo_with_jtype);
          if (num_bo[itype][jtype] >= fvec::VL) {
            #pragma noinline
            aut_lj_with_bo<MORSEFLAG>(ka, itype, jtype, i_bo[itype][jtype],
                                      j_bo[itype][jtype], cij_bo[itype][jtype],
                                      testpath_bo[itype][jtype]);
            num_bo[itype][jtype] -= fvec::VL;
            insert_mask = bvec::onlyafter(num_bo[itype][jtype],
                                          bvec::popcnt(need_bo_with_jtype) -
                                          num_bo[itype][jtype]);
            i_bo[itype][jtype] = ivec::masku_compress(insert_mask, itmp);
            j_bo[itype][jtype] = ivec::masku_compress(insert_mask, jtmp);
            cij_bo[itype][jtype] = fvec::masku_compress(insert_mask, cijtmp);
            if (bvec::test_any_set(need_path_force_with_jtype)) {
              int pos = 0;
              for (int l = testpath_end; l < fvec::VL; l++) {
                if (bvec::test_at(need_path_force_with_jtype, l)) {
                  testpath_bo[itype][jtype][pos] = testpath[l];
                }
                if (bvec::test_at(need_bo_with_jtype, l)) {
                  pos += 1;
                }
              }
            }
          }
        }
        current_mask = current_mask & ~ need_bondorder;
        need_path_force = need_path_force & ~ need_bondorder;
      }

      fvec fpdVLJ = cij *  dVLJ * ( c_1_0 +  Str * ( Stb -  c_1_0));
      fvec fpdStr = dStr *  cij * ( Stb *  VLJ -  VLJ);
      fvec fpair = r2inv *  rij * ( fvec::setzero() - ( fpdVLJ +  fpdStr));
      fvec evdwl = VA *  Stb +  cij *  VLJ * ( c_1_0 -  Str);

      fvec fix = fpair *  delx +  fij[0];
      fvec fiy = fpair *  dely +  fij[1];
      fvec fiz = fpair *  delz +  fij[2];
      result_f_i_x = fvec::mask_add(result_f_i_x, current_mask, result_f_i_x,
                                    fix);
      result_f_i_y = fvec::mask_add(result_f_i_y, current_mask, result_f_i_y,
                                    fiy);
      result_f_i_z = fvec::mask_add(result_f_i_z, current_mask, result_f_i_z,
                                    fiz);
      result_eng = fvec::mask_add(result_eng, current_mask, result_eng, evdwl);

      ivec j_dbl_idx = ivec::mullo(j, c_i4);
      avec fjx = avec::mask_gather(avec::undefined(), current_mask, j_dbl_idx,
                                   &ka->result_f[0].x, sizeof(acc_t));
      avec fjy = avec::mask_gather(avec::undefined(), current_mask, j_dbl_idx,
                                   &ka->result_f[0].y, sizeof(acc_t));
      avec fjz = avec::mask_gather(avec::undefined(), current_mask, j_dbl_idx,
                                   &ka->result_f[0].z, sizeof(acc_t));

      fjx = fjx -  fix;
      fjy = fjy -  fiy;
      fjz = fjz -  fiz;
      avec::mask_i32loscatter(&ka->result_f[0].x, current_mask, j_dbl_idx, fjx,
                              sizeof(acc_t));
      avec::mask_i32loscatter(&ka->result_f[0].y, current_mask, j_dbl_idx, fjy,
                              sizeof(acc_t));
      avec::mask_i32loscatter(&ka->result_f[0].z, current_mask, j_dbl_idx, fjz,
                              sizeof(acc_t));

      if (bvec::test_any_set(need_path_force)) {
        fvec dC = VLJ * ( Str *  Stb +  c_1_0 -  Str);
        #pragma noinline
        aut_airebo_lj_force_path(ka, need_path_force, dC, testpath);
      }
      jj += fvec::VL;
      rest_j = jj < jnum;
    }
    ka->result_f[i].x += fvec::reduce_add(result_f_i_x);
    ka->result_f[i].y += fvec::reduce_add(result_f_i_y);
    ka->result_f[i].z += fvec::reduce_add(result_f_i_z);
  }
  for (int itype = 0; itype < 2; itype++) {
    for (int jtype = 0; jtype < 2; jtype++) {
      for (int l = 0; l < num_bo[itype][jtype]; l++) {
        ref_lennard_jones_single_interaction(ka,ivec::at(i_bo[itype][jtype],l),
                                             ivec::at(j_bo[itype][jtype], l),
                                             MORSEFLAG);
      }
    }
  }
  ka->result_eng += fvec::reduce_add(result_eng);
}

};

template<typename flt_t, typename acc_t>
void aut_lennard_jones(KernelArgsAIREBOT<flt_t,acc_t> * ka, int morseflag) {
#ifdef LMP_INTEL_AIREBO_REF
  ref_lennard_jones(ka, morseflag);
#else
  if (morseflag) {
    aut_wrap<flt_t,acc_t>::template aut_lennard_jones<1>(ka);
  } else {
    aut_wrap<flt_t,acc_t>::template aut_lennard_jones<0>(ka);
  }
#endif
}

template<typename flt_t, typename acc_t>
void aut_rebo_neigh(KernelArgsAIREBOT<flt_t,acc_t> * ka) {
#ifdef LMP_INTEL_AIREBO_REF
  ref_rebo_neigh(ka);
#else
  aut_wrap<flt_t,acc_t>::aut_rebo_neigh(ka);
#endif
}

template<typename flt_t, typename acc_t>
void aut_frebo(KernelArgsAIREBOT<flt_t,acc_t> * ka, int torsion_flag) {
#ifdef LMP_INTEL_AIREBO_REF
  ref_frebo(ka, torsion_flag);
#else
  aut_wrap<flt_t,acc_t>::aut_frebo(ka, torsion_flag);
#endif
}

#ifdef __INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif

}

