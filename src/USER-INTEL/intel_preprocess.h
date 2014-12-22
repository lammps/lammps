/* -*- c++ -*- -------------------------------------------------------------
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
   Contributing author: W. Michael Brown (Intel)
------------------------------------------------------------------------- */

#ifdef __INTEL_OFFLOAD
#ifdef LMP_INTEL_OFFLOAD
#define _LMP_INTEL_OFFLOAD
#endif
#endif

#ifndef LMP_INTEL_PREPROCESS_H
#define LMP_INTEL_PREPROCESS_H

#ifndef LAMMPS_MEMALIGN
#error Please set -DLAMMPS_MEMALIGN=64 in CCFLAGS for your LAMMPS makefile.
#else
#if (LAMMPS_MEMALIGN != 64)
#error Please set -DLAMMPS_MEMALIGN=64 in CCFLAGS for your LAMMPS makefile.
#endif
#endif

#if defined(_OPENMP)
#define _use_omp_pragma(txt) _Pragma(txt)
#else
#define _use_omp_pragma(txt)
#endif

#if defined(__INTEL_COMPILER)
#define _use_simd_pragma(txt) _Pragma(txt)
#else
#define _use_simd_pragma(txt)
#endif

namespace LAMMPS_NS {

enum {LMP_OVERFLOW, LMP_LOCAL_MIN, LMP_LOCAL_MAX, LMP_GHOST_MIN,
      LMP_GHOST_MAX};
enum {TIME_PACK, TIME_HOST_NEIGHBOR, TIME_HOST_PAIR, TIME_OFFLOAD_NEIGHBOR,
      TIME_OFFLOAD_PAIR, TIME_OFFLOAD_WAIT, TIME_OFFLOAD_LATENCY,
      TIME_IMBALANCE};
#define NUM_ITIMERS ( TIME_IMBALANCE + 1 )

#define INTEL_DATA_ALIGN 64
#define INTEL_ONEATOM_FACTOR 2
#define INTEL_MIC_VECTOR_WIDTH 16
#define INTEL_MIC_NBOR_PAD INTEL_MIC_VECTOR_WIDTH
#define INTEL_VECTOR_WIDTH 8
#define INTEL_NBOR_PAD INTEL_VECTOR_WIDTH
#define INTEL_LB_MEAN_WEIGHT 0.1
#define INTEL_BIGP 1e15
#define INTEL_MAX_HOST_CORE_COUNT 512
#define INTEL_MAX_COI_CORES 2

#define IP_PRE_get_stride(stride, n, datasize, torque)	\
  {								\
    int blength = n;						\
    if (torque) blength *= 2;					\
    const int bytes = blength * datasize;			\
    stride = INTEL_DATA_ALIGN - (bytes % INTEL_DATA_ALIGN);     \
    stride = blength + stride / datasize;			\
  }

#if defined(_OPENMP)

#define IP_PRE_omp_range(ifrom, ito, tid, inum, nthreads) 	\
  {								\
    const int idelta = 1 + inum/nthreads;			\
    ifrom = tid * idelta;					\
    ito = ((ifrom + idelta) > inum) ? inum : ifrom + idelta;	\
  }

#define IP_PRE_omp_range_id(ifrom, ito, tid, inum, nthreads)	\
  {								\
    tid = omp_get_thread_num();         			\
    IP_PRE_omp_range(ifrom, ito, tid, inum, nthreads);		\
  }

#define IP_PRE_omp_range_align(ifrom, ito, tid, inum, nthreads, \
                             datasize)                          \
{                                                               \
  int chunk_size = INTEL_DATA_ALIGN / datasize;                 \
  int idelta = static_cast<int>(static_cast<float>(inum)	\
				/chunk_size/nthreads) + 1;	\
  idelta *= chunk_size;						\
  ifrom = tid*idelta;                                           \
  ito = ifrom + idelta;                                         \
  if (ito > inum) ito = inum;                                   \
}

#define IP_PRE_omp_range_id_align(ifrom, ito, tid, inum,        \
				nthreads, datasize)		\
  {								\
    tid = omp_get_thread_num();         			\
    IP_PRE_omp_range_align(ifrom, ito, tid, inum, nthreads,     \
			   datasize);				\
  }

#else

#define IP_PRE_omp_range(ifrom, ito, tid, inum, nthreads)	\
  {								\
    ifrom = 0;							\
    ito = inum;						        \
  }

#define IP_PRE_omp_range_id(ifrom, ito, tid, inum, nthreads)	\
  {								\
    tid = 0;							\
    ifrom = 0;							\
    ito = inum;							\
  }

#define IP_PRE_omp_range_align(ifrom, ito, tid, inum, nthreads, \
                             datasize)                          \
{                                                               \
    ifrom = 0;							\
    ito = inum;						        \
}

#define IP_PRE_omp_range_id_align(ifrom, ito, tid, inum,        \
				nthreads, datasize)		\
{								\
  tid = 0;							\
  ifrom = 0;							\
  ito = inum;							\
}

#endif

#ifdef _LMP_INTEL_OFFLOAD
#include <sys/time.h>

__declspec( target (mic))
inline double MIC_Wtime() {
  double time;
  struct timeval tv;

  gettimeofday(&tv, NULL);
  time = 1.0 * tv.tv_sec + 1.0e-6 * tv.tv_usec;
  return time;
}

#define IP_PRE_pack_separate_buffers(fix, buffers, ago, offload,	\
				     nlocal, nall)			\
{									\
    if (fix->separate_buffers() && ago != 0) {				\
    fix->start_watch(TIME_PACK);					\
    if (offload) {							\
      _use_omp_pragma("omp parallel default(none) shared(buffers,nlocal,nall)")	\
      {									\
        int ifrom, ito, tid;						\
	int nthreads = comm->nthreads;					\
	IP_PRE_omp_range_id_align(ifrom, ito, tid, nlocal,		\
				nthreads, sizeof(flt_t));		\
	buffers->thr_pack_cop(ifrom, ito, 0);				\
	int nghost = nall - nlocal;					\
	if (nghost) {							\
	  IP_PRE_omp_range_align(ifrom, ito, tid, nall - nlocal,	\
				 nthreads, sizeof(flt_t));		\
	  buffers->thr_pack_cop(ifrom + nlocal, ito + nlocal,		\
				fix->offload_min_ghost() - nlocal,	\
				ago == 1);				\
	}								\
      }									\
    } else {								\
      buffers->thr_pack_host(fix->host_min_local(), nlocal, 0);		\
      buffers->thr_pack_host(nlocal, nall,				\
			     fix->host_min_ghost()-nlocal);		\
    }									\
    fix->stop_watch(TIME_PACK);						\
  }									\
}

#define IP_PRE_get_transfern(ago, newton, evflag, eflag, vflag, 	\
			     buffers, offload, fix, separate_flag,	\
			     x_size, q_size, ev_size, f_stride)		\
{									\
  separate_flag = 0;							\
  if (ago == 0) {							\
    x_size = 0;								\
    q_size = nall;							\
    if (offload) {							\
      if (fix->separate_buffers()) {					\
	if (lmp->atom->torque)						\
	  separate_flag = 2;						\
	else								\
	  separate_flag = 1;						\
      } else								\
	separate_flag = 3;						\
    }									\
  } else {								\
    x_size = nall;							\
    q_size = 0;								\
  }									\
  ev_size = 0;								\
  if (evflag) {								\
    if (eflag) ev_size = 2;						\
    if (vflag) ev_size = 8;						\
  }									\
  int f_length;								\
  if (newton)								\
    f_length = nall;							\
  else									\
    f_length = nlocal;							\
  f_length -= minlocal;							\
  f_stride = buffers->get_stride(f_length);				\
}

#define IP_PRE_get_buffers(offload, buffers, fix, tc, f_start,    	\
			   ev_global)					\
{									\
  if (offload) {							\
    tc = buffers->get_off_threads();					\
    f_start = buffers->get_off_f();					\
    ev_global = buffers->get_ev_global();				\
  } else {								\
    tc = comm->nthreads;						\
    f_start = buffers->get_f();						\
    fix->start_watch(TIME_HOST_PAIR);					\
    ev_global = buffers->get_ev_global_host();				\
  }									\
}

#define IP_PRE_repack_for_offload(newton, separate_flag, nlocal, nall,	\
				  f_stride, x, q)			\
{									\
  if (separate_flag) {							\
    if (separate_flag < 3) {						\
      int all_local = nlocal;						\
      int ghost_min = overflow[LMP_GHOST_MIN];				\
      nlocal = overflow[LMP_LOCAL_MAX] + 1;				\
      int nghost = overflow[LMP_GHOST_MAX] + 1 - ghost_min;		\
      if (nghost < 0) nghost = 0;					\
      nall = nlocal + nghost;						\
      separate_flag--;							\
      int flength;							\
      if (NEWTON_PAIR) flength = nall;					\
      else flength = nlocal;						\
      IP_PRE_get_stride(f_stride, flength, sizeof(FORCE_T),		\
			   separate_flag);				\
      if (nghost) {							\
	if (nlocal < all_local || ghost_min > all_local) {		\
	  memmove(x + nlocal, x + ghost_min,				\
		  (nall - nlocal) * sizeof(ATOM_T));			\
	  if (q != 0)							\
	    memmove((void *)(q + nlocal), (void *)(q + ghost_min),	\
		    (nall - nlocal) * sizeof(flt_t));			\
	}								\
      }									\
    }									\
    x[nall].x = INTEL_BIGP;						\
    x[nall].y = INTEL_BIGP;						\
    x[nall].z = INTEL_BIGP;						\
  }									\
}


#else

#define MIC_Wtime MPI_Wtime
#define IP_PRE_pack_separate_buffers(fix, buffers, ago, offload,        \
                                     nlocal, nall)			

#define IP_PRE_get_transfern(ago, newton, evflag, eflag, vflag, 	\
			     buffers, offload, fix, separate_flag,	\
			     x_size, q_size, ev_size, f_stride)		\
{                                                                       \
  separate_flag = 0;							\
  int f_length;                                                         \
  if (newton)                                                           \
    f_length = nall;                                                    \
  else                                                                  \
    f_length = nlocal;                                                  \
  f_stride = buffers->get_stride(f_length);				\
}

#define IP_PRE_get_buffers(offload, buffers, fix, tc, f_start,    	\
			   ev_global)					\
{									\
  tc = comm->nthreads;							\
  f_start = buffers->get_f();						\
  fix->start_watch(TIME_HOST_PAIR);					\
  ev_global = buffers->get_ev_global_host();				\
}

#define IP_PRE_repack_for_offload(newton, separate_flag, nlocal, nall,	\
				  f_stride, x, q)			


#endif

#define IP_PRE_ev_tally_nbor(vflag, ev_pre, fpair, delx, dely, delz)	\
{									\
  if (vflag == 1) {							\
    sv0 += ev_pre * delx * delx * fpair;				\
    sv1 += ev_pre * dely * dely * fpair;				\
    sv2 += ev_pre * delz * delz * fpair;				\
    sv3 += ev_pre * delx * dely * fpair;				\
    sv4 += ev_pre * delx * delz * fpair;				\
    sv5 += ev_pre * dely * delz * fpair;				\
  }									\
}

#define IP_PRE_ev_tally_atom(evflag, eflag, vflag, f, fwtmp)    	\
{									\
  if (evflag) {								\
    if (eflag) {							\
      f[i].w += fwtmp;							\
      oevdwl += sevdwl;							\
    }									\
    if (vflag == 1) {							\
      ov0 += sv0;							\
      ov1 += sv1;							\
      ov2 += sv2;							\
      ov3 += sv3;							\
      ov4 += sv4;							\
      ov5 += sv5;							\
    }									\
  }									\
}

#define IP_PRE_ev_tally_atomq(evflag, eflag, vflag, f, fwtmp)    	\
{									\
  if (evflag) {								\
    if (eflag) {							\
      f[i].w += fwtmp;							\
      oevdwl += sevdwl;							\
      oecoul += secoul;							\
    }									\
    if (vflag == 1) {							\
      ov0 += sv0;							\
      ov1 += sv1;							\
      ov2 += sv2;							\
      ov3 += sv3;							\
      ov4 += sv4;							\
      ov5 += sv5;							\
    }									\
  }									\
}

#define IP_PRE_fdotr_acc_force(newton, evflag, eflag, vflag, eatom,	\
			       nall, nlocal, minlocal, nthreads,	\
			       f_start, f_stride, x)			\
{									\
  int o_range;								\
  if (newton)								\
    o_range = nall;							\
  else									\
    o_range = nlocal;							\
  if (offload == 0) o_range -= minlocal;				\
    IP_PRE_omp_range_align(iifrom, iito, tid, o_range, nthreads,	\
			 sizeof(acc_t));				\
									\
  int t_off = f_stride;						        \
  if (eflag && eatom) {							\
    for (int t = 1; t < nthreads; t++) {				\
      _use_simd_pragma("vector nontemporal")				\
      _use_simd_pragma("novector")					\
      for (int n = iifrom; n < iito; n++) {				\
        f_start[n].x += f_start[n + t_off].x;				\
        f_start[n].y += f_start[n + t_off].y;				\
	f_start[n].z += f_start[n + t_off].z;				\
	f_start[n].w += f_start[n + t_off].w;				\
      }									\
      t_off += f_stride;						\
    }									\
  } else {								\
    for (int t = 1; t < nthreads; t++) {				\
      _use_simd_pragma("vector nontemporal")  				\
      _use_simd_pragma("novector")					\
      for (int n = iifrom; n < iito; n++) {                             \
	f_start[n].x += f_start[n + t_off].x;                  	        \
        f_start[n].y += f_start[n + t_off].y;				\
        f_start[n].z += f_start[n + t_off].z;				\
      }									\
      t_off += f_stride;						\
    }									\
  }									\
									\
  if (evflag) {								\
    if (vflag == 2) {							\
      const ATOM_T * _noalias const xo = x + minlocal;			\
      _use_simd_pragma("vector nontemporal")   				\
      _use_simd_pragma("novector")					\
      for (int n = iifrom; n < iito; n++) {				\
	ov0 += f_start[n].x * xo[n].x;					\
	ov1 += f_start[n].y * xo[n].y;					\
	ov2 += f_start[n].z * xo[n].z;					\
	ov3 += f_start[n].y * xo[n].x;					\
	ov4 += f_start[n].z * xo[n].x;					\
	ov5 += f_start[n].z * xo[n].y;					\
      }									\
    }									\
  }									\
}

}

#endif
