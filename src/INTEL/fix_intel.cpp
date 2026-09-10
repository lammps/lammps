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
   Contributing author: W. Michael Brown (Intel)
                        Anupama Kurpad (Intel) - Host Affinitization
------------------------------------------------------------------------- */

#include "fix_intel.h"

#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "pair.h"
#include "pair_hybrid.h"
#include "update.h"

#include <cstring>

#include "suffix.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixIntel::FixIntel(LAMMPS *lmp, int narg, char **arg) :  Fix(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal package intel command");

  // the first argument (formerly the number of coprocessors) is accepted for
  // backward compatibility but ignored: offload to Intel(R) Xeon Phi(TM)
  // coprocessors is no longer supported.
  utils::inumeric(FLERR,arg[3],false,lmp);

  _nbor_pack_width = 1;
  _three_body_neighbor = 0;
  _pair_intel_count = 0;
  _hybrid_nonpair = 0;
  _print_pkg_info = 1;
  _nthreads = comm->nthreads;
  _torque_flag = 0;

  _precision_mode = PREC_MODE_MIXED;
  _overflow_flag[LMP_OVERFLOW] = 0;

  _force_array_s = nullptr;
  _force_array_m = nullptr;
  _force_array_d = nullptr;
  _ev_array_s = nullptr;
  _ev_array_d = nullptr;

  // optional keywords

  int nomp = 0;
  int offload_deprecated = 0;
  _lrt = 0;
  _p3m_table = 1;

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"omp") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal package intel command");
      nomp = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"mode") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal package intel command");
      if (strcmp(arg[iarg+1],"single") == 0)
        _precision_mode = PREC_MODE_SINGLE;
      else if (strcmp(arg[iarg+1],"mixed") == 0)
        _precision_mode = PREC_MODE_MIXED;
      else if (strcmp(arg[iarg+1],"double") == 0)
        _precision_mode = PREC_MODE_DOUBLE;
      else error->all(FLERR,"Illegal package intel command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "lrt") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal package intel command");
      _lrt = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "pppm_table") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal package intel command");
      _p3m_table = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;

    // offload-related keywords are deprecated and ignored: offload to
    // Intel(R) Xeon Phi(TM) coprocessors is no longer supported

    } else if (strcmp(arg[iarg],"balance") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal package intel command");
      offload_deprecated = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg], "ghost") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal package intel command");
      offload_deprecated = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg], "tpc") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal package intel command");
      offload_deprecated = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"tptask") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal package intel command");
      offload_deprecated = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"no_affinity") == 0) {
      offload_deprecated = 1;
      iarg++;

    // undocumented offload keywords (also deprecated and ignored)

    } else if (strcmp(arg[iarg],"offload_affinity_balanced") == 0) {
      offload_deprecated = 1;
      iarg++;
    } else if (strcmp(arg[iarg],"buffers") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal package intel command");
      offload_deprecated = 1;
      iarg += 2;
    } else error->all(FLERR,"Illegal package intel command");
  }

  if (offload_deprecated && comm->me == 0)
    error->warning(FLERR,"Offload to Intel(R) Xeon Phi(TM) coprocessors is no "
                   "longer supported; ignoring offload-related 'package intel' "
                   "keyword(s)");

  // if using LRT mode, create the integrate style
  if (_lrt) {
    char *cmd[1];
    cmd[0] = (char *) "verlet/lrt/intel";
    update->create_integrate(1,cmd,0);
  }

  // error check

  if (nomp < 0) error->all(FLERR,"Illegal package intel command");

  // set OpenMP threads
  // nomp is user setting, default = 0

  #if defined(_OPENMP)
  #if defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
  kmp_set_blocktime(0);
  #endif
  if (nomp != 0) {
    omp_set_num_threads(nomp);
    _nthreads = comm->nthreads = nomp;
  }
  #endif

  // set precision

  if (_precision_mode == PREC_MODE_SINGLE)
    _single_buffers = new IntelBuffers<float,float>(lmp);
  else if (_precision_mode == PREC_MODE_MIXED)
    _mixed_buffers = new IntelBuffers<float,double>(lmp);
  else
    _double_buffers = new IntelBuffers<double,double>(lmp);
}

/* ---------------------------------------------------------------------- */

FixIntel::~FixIntel()
{
  if (_precision_mode == PREC_MODE_SINGLE)
    delete _single_buffers;
  else if (_precision_mode == PREC_MODE_MIXED)
    delete _mixed_buffers;
  else
    delete _double_buffers;
}

/* ---------------------------------------------------------------------- */

int FixIntel::setmask()
{
  int mask = 0;
  mask |= PRE_REVERSE;
  mask |= MIN_PRE_REVERSE;
  mask |= POST_FORCE;
  mask |= MIN_POST_FORCE;
  mask |= POST_RUN;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixIntel::init()
{
  _torque_flag = 0;
  if (force->pair_match("gayberne/intel$", 0)) _torque_flag = 1;

  const int nstyles = _pair_intel_count;
  if (force->pair_match("^hybrid", 0)) {
    _pair_hybrid_flag = 1;

    // Check if need to handle torque
    auto hybrid = dynamic_cast<PairHybrid *>(force->pair);
    if (hybrid) {
      for (int i = 0; i < hybrid->nstyles; i++)
        if (utils::strmatch(hybrid->keywords[i],"/intel$") &&
            utils::strmatch(hybrid->keywords[i],"gayberne"))
          _torque_flag = 1;
    }
    if (force->newton_pair != 0 && force->pair->no_virial_fdotr_compute)
      error->all(FLERR,"INTEL package requires fdotr virial with newton on.");
  } else
    _pair_hybrid_flag = 0;

  if (_torque_flag && nstyles > 1)
    error->all(FLERR,"gayberne/intel style cannot yet be used with other "
               "intel pair styles.");

  if (nstyles > 1 && _pair_hybrid_flag) _pair_hybrid_flag = 2;
  else if (force->newton_pair == 0) _pair_hybrid_flag = 0;

  _pair_hybrid_zero = 0;
  _zero_master = 0;

  if (_pair_hybrid_flag && _hybrid_nonpair)
    _pair_hybrid_zero = 1;
  _hybrid_nonpair = 0;

  _pair_intel_count = 0;

  const int off_mode = 0;
  if (_precision_mode == PREC_MODE_SINGLE) {
    _single_buffers->set_torque_flag(_torque_flag);
    _single_buffers->zero_ev();
    _single_buffers->grow_ncache(off_mode, comm->nthreads);
    _single_buffers->free_list_ptrs();
  } else if (_precision_mode == PREC_MODE_MIXED) {
    _mixed_buffers->set_torque_flag(_torque_flag);
    _mixed_buffers->zero_ev();
    _mixed_buffers->grow_ncache(off_mode, comm->nthreads);
    _mixed_buffers->free_list_ptrs();
  } else {
    _double_buffers->set_torque_flag(_torque_flag);
    _double_buffers->zero_ev();
    _double_buffers->grow_ncache(off_mode, comm->nthreads);
    _double_buffers->free_list_ptrs();
  }

  _need_reduce = 0;
}

/* ---------------------------------------------------------------------- */

void FixIntel::setup(int vflag)
{
  if (neighbor->style != Neighbor::BIN)
    error->all(FLERR,"Currently, neighbor style BIN must be used with INTEL package.");
  if (vflag > 3)
   error->all(FLERR,"Cannot currently get per-atom virials with INTEL package.");
}

/* ---------------------------------------------------------------------- */

void FixIntel::setup_pre_reverse(int eflag, int vflag)
{
  pre_reverse(eflag,vflag);
}

/* ---------------------------------------------------------------------- */

bool FixIntel::pair_hybrid_check()
{
  auto ph = dynamic_cast<PairHybrid *>(force->pair);
  bool has_intel = false;
  int nstyles = ph->nstyles;

  for (int i = 0; i < nstyles; ++i)
    if (ph->styles[i]->suffix_flag & Suffix::INTEL) has_intel = true;

  return has_intel;
}

/* ---------------------------------------------------------------------- */

void FixIntel::pair_init_check(const bool cdmessage)
{
  #ifdef INTEL_VMASK
  if (atom->sortfreq) atom->sortfreq = 1;
  #endif

  _nbor_pack_width = 1;

  _nthreads = comm->nthreads;

  #ifndef LMP_INTEL_NBOR_COMPAT
  if (force->pair->manybody_flag && (atom->molecular != Atom::ATOMIC)) {
    int flag = 0;
    if (atom->nbonds > 0 && force->special_lj[1] == 0.0 &&
        force->special_coul[1] == 0.0) flag = 1;
    if (atom->nangles > 0 && force->special_lj[2] == 0.0 &&
        force->special_coul[2] == 0.0) flag = 1;
    if (atom->ndihedrals > 0 && force->special_lj[3] == 0.0 &&
        force->special_coul[3] == 0.0) flag = 1;
    if (flag)
      error->all(FLERR,"Add -DLMP_INTEL_NBOR_COMPAT to build for special_bond"
                 " exclusions with Intel");
  }
  #endif

  int need_tag = 0;
  if (atom->molecular != Atom::ATOMIC || three_body_neighbor()) need_tag = 1;
  if (domain->triclinic && force->newton_pair) need_tag = 1;

  // Clear buffers used for pair style
  char kmode[80];
  if (_precision_mode == PREC_MODE_SINGLE) {
    strcpy(kmode, "single");
    get_single_buffers()->need_tag(need_tag);
  } else if (_precision_mode == PREC_MODE_MIXED) {
    strcpy(kmode, "mixed");
    get_mixed_buffers()->need_tag(need_tag);
  } else {
    strcpy(kmode, "double");
    get_double_buffers()->need_tag(need_tag);
  }

  _pair_intel_count++;

  if (_print_pkg_info && comm->me == 0) {
    utils::logmesg(lmp, "----------------------------------------------------------\n");
    utils::logmesg(lmp,"Using INTEL Package.\n");
    utils::logmesg(lmp,"Compiler: {}\n",platform::compiler_info());
    #ifdef LMP_SIMD_COMPILER
    utils::logmesg(lmp,"SIMD compiler directives: Enabled\n");
    #else
    utils::logmesg(lmp,"SIMD compiler directives: Disabled\n");
    #endif
    utils::logmesg(lmp,"Precision: {}\n",kmode);
    if (cdmessage) {
      #ifdef LMP_USE_AVXCD
      utils::logmesg(lmp,"AVX512 CD Optimizations: Enabled\n");
      #else
      utils::logmesg(lmp,"AVX512 CD Optimizations: Disabled\n");
      #endif
    }
    utils::logmesg(lmp, "----------------------------------------------------------\n");
  }
  _print_pkg_info = 0;
}

/* ---------------------------------------------------------------------- */

void FixIntel::bond_init_check()
{
  int intel_pair = 0;
  if (force->pair_match("/intel$", 0) != nullptr)
    intel_pair = 1;
  else if (force->pair_match("^hybrid", 0) != nullptr) {
    _hybrid_nonpair = 1;
    if (pair_hybrid_check()) intel_pair = 1;
  }

  if (intel_pair == 0)
    error->all(FLERR,"Intel styles for bond/angle/dihedral/improper require intel pair style.");
}

/* ---------------------------------------------------------------------- */

void FixIntel::kspace_init_check()
{
  int intel_pair = 0;
  if (force->pair_match("/intel$", 0) != nullptr)
    intel_pair = 1;
  else if (force->pair_match("^hybrid", 0) != nullptr) {
    _hybrid_nonpair = 1;
    if (pair_hybrid_check()) intel_pair = 1;
  }

  if (intel_pair == 0)
    error->all(FLERR,"Intel styles for kspace require intel pair style.");

  if (utils::strmatch(update->integrate_style, "^verlet/split"))
    error->all(FLERR,"Intel styles for kspace are not compatible with run_style verlet/split");
}

/* ---------------------------------------------------------------------- */

void FixIntel::_sync_main_arrays(const int prereverse)
{
  if (!prereverse) _zero_master = 1;
  int done_this_step = prereverse;
  if (_pair_hybrid_zero == 0) done_this_step = 1;
  if (_force_array_m != nullptr) {
    if (_need_reduce) {
      reduce_results(&_force_array_m[0].x);
      _need_reduce = 0;
    }
    add_results(_force_array_m, _ev_array_d, _results_eatom);
    if (done_this_step) _force_array_m = nullptr;
    else _ev_array_d = nullptr;
  } else if (_force_array_d != nullptr) {
    if (_need_reduce) {
      reduce_results(&_force_array_d[0].x);
      _need_reduce = 0;
    }
    add_results(_force_array_d, _ev_array_d, _results_eatom);
    if (done_this_step) _force_array_d = nullptr;
    else _ev_array_d = nullptr;
  } else if (_force_array_s != nullptr) {
    if (_need_reduce) {
      reduce_results(&_force_array_s[0].x);
      _need_reduce = 0;
    }
    add_results(_force_array_s, _ev_array_s, _results_eatom);
    if (done_this_step) _force_array_s = nullptr;
    else _ev_array_s = nullptr;
  }
}

/* ---------------------------------------------------------------------- */

void FixIntel::pre_reverse(int /*eflag*/, int /*vflag*/)
{
  _sync_main_arrays(1);
}

/* ---------------------------------------------------------------------- */

void FixIntel::post_force(int /*vflag*/)
{
  // Redundant call to sync Intel data structs with native for methods that
  // call force compute but do not call prereverse
  _sync_main_arrays(1);
}

/* ---------------------------------------------------------------------- */

template <class acc_t>
void FixIntel::reduce_results(acc_t * _noalias const f_scalar)
{
  int o_range, f_stride;
  if (force->newton_pair)
    o_range = atom->nlocal + atom->nghost;
  else
    o_range = atom->nlocal;
  IP_PRE_get_stride(f_stride, o_range, (sizeof(acc_t)*4), _torque_flag);

  o_range *= 4;
  const int f_stride4 = f_stride * 4;

  if (_nthreads <= INTEL_HTHREADS) {
    acc_t *f_scalar2 = f_scalar + f_stride4;
    if (_nthreads == 4) {
      acc_t *f_scalar3 = f_scalar2 + f_stride4;
      acc_t *f_scalar4 = f_scalar3 + f_stride4;
      #if defined(USE_OMP_SIMD)
      #pragma omp simd aligned(f_scalar,f_scalar2,f_scalar3,f_scalar4:64)
      #elif defined(LMP_SIMD_COMPILER)
      #pragma vector aligned
      #pragma simd
      #endif
      for (int n = 0; n < o_range; n++)
        f_scalar[n] += f_scalar2[n] + f_scalar3[n] + f_scalar4[n];
    } else if (_nthreads == 2) {
      #if defined(USE_OMP_SIMD)
      #pragma omp simd aligned(f_scalar,f_scalar2:64)
      #elif defined(LMP_SIMD_COMPILER)
      #pragma vector aligned
      #pragma simd
      #endif
      for (int n = 0; n < o_range; n++)
        f_scalar[n] += f_scalar2[n];
    } else {
      acc_t *f_scalar3 = f_scalar2 + f_stride4;
      #if defined(USE_OMP_SIMD)
      #pragma omp simd aligned(f_scalar,f_scalar2,f_scalar3:64)
      #elif defined(LMP_SIMD_COMPILER)
      #pragma vector aligned
      #pragma simd
      #endif
      for (int n = 0; n < o_range; n++)
        f_scalar[n] += f_scalar2[n] + f_scalar3[n];
    }
  } else {
    #if defined(_OPENMP)
    #pragma omp parallel
    #endif
    {
      int iifrom, iito, tid;
      IP_PRE_omp_range_id_align(iifrom, iito, tid, o_range, _nthreads,
                                sizeof(acc_t));

      acc_t *f_scalar2 = f_scalar + f_stride4;
      for (int t = 1; t < _nthreads; t++) {
        #if defined(USE_OMP_SIMD)
        #pragma omp simd aligned(f_scalar,f_scalar2:64)
        #elif defined(LMP_SIMD_COMPILER)
        #pragma vector aligned
        #pragma simd
        #endif
        for (int n = iifrom; n < iito; n++)
          f_scalar[n] += f_scalar2[n];
        f_scalar2 += f_stride4;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

template <class ft, class acc_t>
void FixIntel::add_results(const ft * _noalias const f_in,
                           const acc_t * _noalias const ev_global,
                           const int eatom) {
  start_watch(TIME_PACK);
  int f_length;
  if (force->newton_pair)
    f_length = atom->nlocal + atom->nghost;
  else
    f_length = atom->nlocal;
  add_oresults(f_in, ev_global, eatom, 0, f_length);
  stop_watch(TIME_PACK);
}

/* ---------------------------------------------------------------------- */

template <class ft, class acc_t>
void FixIntel::add_oresults(const ft * _noalias const f_in,
                            const acc_t * _noalias const ev_global,
                            const int eatom,
                            const int out_offset, const int nall) {
  lmp_ft * _noalias const f = (lmp_ft *) lmp->atom->f[0] + out_offset;
  if (_torque_flag) {
    if (f_in[1].w) {
      if (f_in[1].w == 1)
        error->all(FLERR,"Bad matrix inversion in mldivide3");
      else
        error->all(FLERR,"Sphere particles not yet supported for gayberne/intel");
    }
  }

  int packthreads;
  if (_nthreads > INTEL_HTHREADS) packthreads = _nthreads;
  else packthreads = 1;
  #if defined(_OPENMP)
  #pragma omp parallel if (packthreads > 1)
  #endif
  {
    #if defined(_OPENMP)
    const int tid = omp_get_thread_num();
    #else
    const int tid = 0;
    #endif
    int ifrom, ito;
    IP_PRE_omp_range_align(ifrom, ito, tid, nall, packthreads, sizeof(acc_t));
    if (_torque_flag) {
      int ii = ifrom * 2;
      lmp_ft * _noalias const tor = (lmp_ft *) lmp->atom->torque[0] +
        out_offset;
      if (eatom) {
        double * _noalias const lmp_eatom = force->pair->eatom + out_offset;
        #if defined(LMP_SIMD_COMPILER)
        #pragma vector aligned
        #pragma ivdep
        #endif
        for (int i = ifrom; i < ito; i++) {
          f[i].x += f_in[ii].x;
          f[i].y += f_in[ii].y;
          f[i].z += f_in[ii].z;
          lmp_eatom[i] += f_in[ii].w;
          tor[i].x += f_in[ii+1].x;
          tor[i].y += f_in[ii+1].y;
          tor[i].z += f_in[ii+1].z;
          ii += 2;
        }
      } else {
        #if defined(LMP_SIMD_COMPILER)
        #pragma vector aligned
        #pragma ivdep
        #endif
        for (int i = ifrom; i < ito; i++) {
          f[i].x += f_in[ii].x;
          f[i].y += f_in[ii].y;
          f[i].z += f_in[ii].z;
          tor[i].x += f_in[ii+1].x;
          tor[i].y += f_in[ii+1].y;
          tor[i].z += f_in[ii+1].z;
          ii += 2;
        }
      }
    } else {
      if (eatom) {
        double * _noalias const lmp_eatom = force->pair->eatom + out_offset;
        #if defined(LMP_SIMD_COMPILER)
        #pragma vector aligned
        #pragma ivdep
        #endif
        for (int i = ifrom; i < ito; i++) {
          f[i].x += f_in[i].x;
          f[i].y += f_in[i].y;
          f[i].z += f_in[i].z;
          lmp_eatom[i] += f_in[i].w;
        }
      } else {
        #if defined(LMP_SIMD_COMPILER)
        #pragma vector aligned
        #pragma ivdep
        #endif
        for (int i = ifrom; i < ito; i++) {
          f[i].x += f_in[i].x;
          f[i].y += f_in[i].y;
          f[i].z += f_in[i].z;
        }
      }
    }
  }

  if (ev_global != nullptr) {
    force->pair->eng_vdwl += ev_global[0];
    force->pair->eng_coul += ev_global[1];
    force->pair->virial[0] += ev_global[2];
    force->pair->virial[1] += ev_global[3];
    force->pair->virial[2] += ev_global[4];
    force->pair->virial[3] += ev_global[5];
    force->pair->virial[4] += ev_global[6];
    force->pair->virial[5] += ev_global[7];
  }
}

/* ---------------------------------------------------------------------- */

double FixIntel::memory_usage()
{
  double bytes;
  if (_precision_mode == PREC_MODE_SINGLE)
    bytes = _single_buffers->memory_usage(_nthreads);
  else if (_precision_mode == PREC_MODE_MIXED)
    bytes = _mixed_buffers->memory_usage(_nthreads);
  else
    bytes = _double_buffers->memory_usage(_nthreads);

  return bytes;
}
