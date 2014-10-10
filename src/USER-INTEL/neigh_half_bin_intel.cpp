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
   Contributing author: W. Michael Brown (Intel)
------------------------------------------------------------------------- */

#include "neighbor.h"
#include "neigh_list.h"
#include "atom.h"
#include "comm.h"
#include "group.h"
#include "fix_intel.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace LAMMPS_NS;

#ifdef _LMP_INTEL_OFFLOAD
#pragma offload_attribute(push,target(mic))
#endif

#define ofind_special(which, special, nspecial, i, tag, special_flag) \
{                                                                     \
  which = 0;                                                          \
  const int n1 = nspecial[i * 3];                                     \
  const int n2 = nspecial[i * 3 + 1];                                 \
  const int n3 = nspecial[i * 3 + 2];                                 \
  const tagint *sptr = special + i * maxspecial;                      \
  for (int s = 0; s < n3; s++) {                                      \
    if (sptr[s] == tag) {                                             \
      if (s < n1) {                                                   \
        if (special_flag[1] == 0) which = -1;                         \
        else if (special_flag[1] == 1) which = 0;                     \
        else which = 1;                                               \
      } else if (s < n2) {                                            \
        if (special_flag[2] == 0) which = -1;                         \
        else if (special_flag[2] == 1) which = 0;                     \
        else which = 2;                                               \
      } else {                                                        \
        if (special_flag[3] == 0) which = -1;                         \
        else if (special_flag[3] == 1) which = 0;                     \
        else which = 3;                                               \
      }                                                               \
    }                                                                 \
  }                                                                   \
}

#ifdef _LMP_INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif

template <class flt_t, class acc_t>
void Neighbor::bin_atoms(void * xin, int * _noalias const atombin) {
  const ATOM_T * _noalias const x = (const ATOM_T * _noalias const)xin;
  int nlocal = atom->nlocal;
  const int nall = nlocal + atom->nghost;

  const double sboxlo0 = bboxlo[0] + mbinxlo/bininvx;
  const double sboxlo1 = bboxlo[1] + mbinylo/bininvy;
  const double sboxlo2 = bboxlo[2] + mbinzlo/bininvz;

  int i, ibin;

  for (i = 0; i < mbins; i++) binhead[i] = -1;

  int *mask = atom->mask;

  if (includegroup) {
    int bitmask = group->bitmask[includegroup];
    for (i = nall-1; i >= nlocal; i--) {
      if (mask[i] & bitmask) {
        ibin = coord2bin(atom->x[i]);
        bins[i] = binhead[ibin];
        binhead[ibin] = i;
      }
    }
    for (i = atom->nfirst-1; i >= 0; i--) {
      ibin = coord2bin(atom->x[i]);
      atombin[i] = ibin;
      bins[i] = binhead[ibin];
      binhead[ibin] = i;
    }
  } else {
    for (i = nall-1; i >= nlocal; i--) {
      ibin = coord2bin(atom->x[i]);
      bins[i] = binhead[ibin];
      binhead[ibin] = i;
    }
    for (i = nlocal-1; i >= 0; i--) {
      ibin = coord2bin(atom->x[i]);
      atombin[i]=ibin;
      bins[i] = binhead[ibin];
      binhead[ibin] = i;
    }
  }
}

/* ----------------------------------------------------------------------
   binned neighbor list construction with partial Newton's 3rd law
   each owned atom i checks own bin and other bins in stencil
   pair stored once if i,j are both owned and i < j
   pair stored by me if j is ghost (also stored by proc owning j)
------------------------------------------------------------------------- */

void Neighbor::half_bin_no_newton_intel(NeighList *list)
{
  const int nlocal = (includegroup) ? atom->nfirst : atom->nlocal;
  list->inum = nlocal;

  // Get fix for intel stuff
  FixIntel *fix = static_cast<FixIntel *>(fix_intel);

  const int off_end = fix->offload_end_neighbor();
  int host_start = off_end;;
  #ifdef _LMP_INTEL_OFFLOAD
  if (fix->full_host_list()) host_start = 0;
  if (exclude) 
    error->all(FLERR, "Exclusion lists not yet supported for Intel offload");
  #endif

  if (fix->precision() == FixIntel::PREC_MODE_MIXED) {
    hbnni<float,double>(1, list, fix->get_mixed_buffers(),
                        0, off_end, fix);
    hbnni<float,double>(0, list, fix->get_mixed_buffers(),
                        host_start, nlocal,fix);
  } else if (fix->precision() == FixIntel::PREC_MODE_DOUBLE) {
    hbnni<double,double>(1, list, fix->get_double_buffers(),
                         0, off_end, fix);
    hbnni<double,double>(0, list, fix->get_double_buffers(),
                         host_start, nlocal, fix);
  } else {
    hbnni<float,float>(1, list, fix->get_single_buffers(),
                       0, off_end, fix);
    hbnni<float,float>(0, list, fix->get_single_buffers(),
                       host_start, nlocal, fix);
  }
}

template <class flt_t, class acc_t>
void Neighbor::hbnni(const int offload, NeighList *list, void *buffers_in,
                     const int astart, const int aend, void *fix_in) {
  IntelBuffers<flt_t,acc_t> *buffers = (IntelBuffers<flt_t,acc_t> *)buffers_in;
  FixIntel *fix = (FixIntel *)fix_in;
  const int nall = atom->nlocal + atom->nghost;
  int pad = 1;

  if (offload) {
    fix->start_watch(TIME_PACK);
    buffers->grow(nall, atom->nlocal, comm->nthreads, aend);
    buffers->grow_nbor(list, atom->nlocal, aend);

    ATOM_T biga;
    biga.x = INTEL_BIGP;
    biga.y = INTEL_BIGP;
    biga.z = INTEL_BIGP;
    biga.w = 1;
    buffers->get_x()[nall] = biga;

    const int nthreads = comm->nthreads;
    #if defined(_OPENMP)
    #pragma omp parallel default(none) shared(buffers)
    #endif
    {
      int ifrom, ito, tid;
      IP_PRE_omp_range_id_align(ifrom, ito, tid, nall, nthreads, 
				sizeof(ATOM_T));
      buffers->thr_pack(ifrom, ito, 0);
    }
    fix->stop_watch(TIME_PACK);

    fix->start_watch(TIME_HOST_NEIGHBOR);
    bin_atoms<flt_t,acc_t>(buffers->get_x(), buffers->get_atombin());
    if (INTEL_MIC_NBOR_PAD > 1)
      pad = INTEL_MIC_NBOR_PAD * sizeof(float) / sizeof(flt_t);
  } else {
    fix->start_watch(TIME_HOST_NEIGHBOR);
    if (INTEL_NBOR_PAD > 1)
      pad = INTEL_NBOR_PAD * sizeof(float) / sizeof(flt_t);
  }
  const int pad_width = pad;

  if (aend-astart == 0) {
    fix->stop_watch(TIME_HOST_NEIGHBOR);
    return;
  }

  const ATOM_T * _noalias const x = buffers->get_x();
  int * _noalias const firstneigh = buffers->firstneigh(list);

  const int molecular = atom->molecular;
  int *ns = NULL;
  tagint *s = NULL;
  int tag_size, special_size;
  if (molecular) {
    s = atom->special[0];
    ns = atom->nspecial[0];
    tag_size = nall;
    special_size = aend;
  } else {
    s = &buffers->_special_holder;
    ns = &buffers->_nspecial_holder;
    tag_size = 0;
    special_size = 0;
  }
  const tagint * _noalias const special = s;
  const int * _noalias const nspecial = ns;
  const int maxspecial = atom->maxspecial;
  const tagint * _noalias const tag = atom->tag;

  int * _noalias const ilist = list->ilist;
  int * _noalias numneigh = list->numneigh;
  int * _noalias const cnumneigh = buffers->cnumneigh(list);
  const int nstencil = list->nstencil;
  const int * _noalias const stencil = list->stencil;
  const flt_t * _noalias const cutneighsq = buffers->get_cutneighsq()[0];
  const int ntypes = atom->ntypes + 1;
  const int nlocal = atom->nlocal;

  #ifndef _LMP_INTEL_OFFLOAD
  int * const mask = atom->mask;
  tagint * const molecule = atom->molecule;
  #endif

  int tnum;
  int *overflow;
  double *timer_compute;
  if (offload) {
    timer_compute = fix->off_watch_neighbor();
    tnum = buffers->get_off_threads();
    overflow = fix->get_off_overflow_flag();
    fix->stop_watch(TIME_HOST_NEIGHBOR);
    fix->start_watch(TIME_OFFLOAD_LATENCY);
  } else {
    tnum = comm->nthreads;
    overflow = fix->get_overflow_flag();
  }
  const int nthreads = tnum;
  const int maxnbors = buffers->get_max_nbors();
  int * _noalias const atombin = buffers->get_atombin();

  // Make sure dummy coordinates to eliminate loop remainder not within cutoff
  {
    const flt_t dx = (INTEL_BIGP - bboxhi[0]);
    const flt_t dy = (INTEL_BIGP - bboxhi[1]);
    const flt_t dz = (INTEL_BIGP - bboxhi[2]);
    if (dx * dx + dy * dy + dz * dz < static_cast<flt_t>(cutneighmaxsq))
      error->one(FLERR,
	"Intel package expects no atoms within cutoff of {1e15,1e15,1e15}.");
  }

  #ifdef _LMP_INTEL_OFFLOAD
  const int * _noalias const binhead = this->binhead;
  const int * _noalias const special_flag = this->special_flag;
  const int * _noalias const bins = this->bins;
  const int cop = fix->coprocessor_number();
  const int separate_buffers = fix->separate_buffers();
  #pragma offload target(mic:cop) if(offload) \
    in(x:length(nall+1) alloc_if(0) free_if(0)) \
    in(tag:length(tag_size) alloc_if(0) free_if(0)) \
    in(special:length(special_size*maxspecial) alloc_if(0) free_if(0)) \
    in(nspecial:length(special_size*3) alloc_if(0) free_if(0)) \
    in(bins:length(nall) alloc_if(0) free_if(0)) \
    in(binhead:length(mbins) alloc_if(0) free_if(0)) \
    in(cutneighsq:length(0) alloc_if(0) free_if(0)) \
    in(firstneigh:length(0) alloc_if(0) free_if(0)) \
    in(cnumneigh:length(0) alloc_if(0) free_if(0)) \
    out(numneigh:length(0) alloc_if(0) free_if(0)) \
    in(ilist:length(0) alloc_if(0) free_if(0)) \
    in(atombin:length(aend) alloc_if(0) free_if(0)) \
    in(stencil:length(nstencil) alloc_if(0) free_if(0)) \
    in(special_flag:length(0) alloc_if(0) free_if(0)) \
    in(maxnbors,nthreads,maxspecial,nstencil,pad_width,offload) \
    in(separate_buffers, astart, aend, nlocal, molecular, ntypes) \
    out(overflow:length(5) alloc_if(0) free_if(0)) \
    out(timer_compute:length(1) alloc_if(0) free_if(0)) \
    signal(numneigh)
  #endif
  {
    #ifdef __MIC__
    *timer_compute = MIC_Wtime();
    #endif

    #ifdef _LMP_INTEL_OFFLOAD
    overflow[LMP_LOCAL_MIN] = astart;
    overflow[LMP_LOCAL_MAX] = aend - 1;
    overflow[LMP_GHOST_MIN] = nall;
    overflow[LMP_GHOST_MAX] = -1;
    #endif

    #if defined(_OPENMP)
    #pragma omp parallel default(none) shared(numneigh,overflow)
    #endif
    {
      #ifdef _LMP_INTEL_OFFLOAD
      int lmin = nall, lmax = -1, gmin = nall, gmax = -1;
      #endif

      const int num = aend - astart;
      int tid, ifrom, ito;
      IP_PRE_omp_range_id(ifrom, ito, tid, num, nthreads);
      ifrom += astart;
      ito += astart;

      int which;

      const int list_size = (ito + tid + 1) * maxnbors;
      int ct = (ifrom + tid) * maxnbors;
      int *neighptr = firstneigh + ct;
      for (int i = ifrom; i < ito; i++) {
        int j, k, n, n2, itype, jtype, ibin;
        double xtmp, ytmp, ztmp, delx, dely, delz, rsq;

        n = 0;
        n2 = maxnbors;

        xtmp = x[i].x;
        ytmp = x[i].y;
        ztmp = x[i].z;
        itype = x[i].w;
        const int ioffset = ntypes*itype;

        // loop over all atoms in other bins in stencil including self
        // only store pair if i < j
        // stores own/own pairs only once
        // stores own/ghost pairs on both procs

        ibin = atombin[i];

        for (k = 0; k < nstencil; k++) {
          for (j = binhead[ibin + stencil[k]]; j >= 0; j = bins[j]) {
            if (j <= i) continue;

            jtype = x[j].w;
            #ifndef _LMP_INTEL_OFFLOAD
            if (exclude && exclusion(i,j,itype,jtype,mask,molecule)) continue;
            #endif

            delx = xtmp - x[j].x;
            dely = ytmp - x[j].y;
            delz = ztmp - x[j].z;
            rsq = delx * delx + dely * dely + delz * delz;
            if (rsq <= cutneighsq[ioffset + jtype]) {
              if (j < nlocal) {
                neighptr[n++] = j;
                #ifdef _LMP_INTEL_OFFLOAD
		if (j < lmin) lmin = j;
		if (j > lmax) lmax = j;
                #endif
              } else {
                neighptr[n2++] = j;
	        #ifdef _LMP_INTEL_OFFLOAD
		if (j < gmin) gmin = j;
		if (j > gmax) gmax = j;
                #endif
	      }
	    }
          }
        }
        ilist[i] = i;

        cnumneigh[i] = ct;
        if (n > maxnbors) *overflow = 1;
        for (k = maxnbors; k < n2; k++) neighptr[n++] = neighptr[k];
        while( (n % pad_width) != 0 ) neighptr[n++] = nall;
        numneigh[i] = n;
        while((n % (INTEL_DATA_ALIGN / sizeof(int))) != 0) n++;
        ct += n;
        neighptr += n;
        if (ct + n + maxnbors > list_size) {
          *overflow = 1;
	  ct = (ifrom + tid) * maxnbors;
        }
      }

      if (*overflow == 1)
	for (int i = ifrom; i < ito; i++)
	  numneigh[i] = 0;

      #ifdef _LMP_INTEL_OFFLOAD
      if (separate_buffers) {
        #if defined(_OPENMP)
        #pragma omp critical
        #endif
        {
          if (lmin < overflow[LMP_LOCAL_MIN]) overflow[LMP_LOCAL_MIN] = lmin;
          if (lmax > overflow[LMP_LOCAL_MAX]) overflow[LMP_LOCAL_MAX] = lmax;
          if (gmin < overflow[LMP_GHOST_MIN]) overflow[LMP_GHOST_MIN] = gmin;
          if (gmax > overflow[LMP_GHOST_MAX]) overflow[LMP_GHOST_MAX] = gmax;
        }
        #pragma omp barrier
      }

      int ghost_offset = 0, nall_offset = nall;
      if (separate_buffers) {
        int nghost = overflow[LMP_GHOST_MAX] + 1 - overflow[LMP_GHOST_MIN];
        if (nghost < 0) nghost = 0;
        if (offload) {
          ghost_offset = overflow[LMP_GHOST_MIN] - overflow[LMP_LOCAL_MAX] - 1;
          nall_offset = overflow[LMP_LOCAL_MAX] + 1 + nghost;
	} else {
          ghost_offset = overflow[LMP_GHOST_MIN] - nlocal;
          nall_offset = nlocal + nghost;
        }
      }
      #endif

      if (molecular) {
        for (int i = ifrom; i < ito; ++i) {
          int * _noalias jlist = firstneigh + cnumneigh[i];
          const int jnum = numneigh[i];
          for (int jj = 0; jj < jnum; jj++) {
            const int j = jlist[jj];
            ofind_special(which, special, nspecial, i, tag[j], special_flag);
            #ifdef _LMP_INTEL_OFFLOAD
	    if (j >= nlocal) {
	      if (j == nall) 
		jlist[jj] = nall_offset;
	      else if (which > 0) 
		jlist[jj] = (j-ghost_offset) ^ (which << SBBITS);
	      else jlist[jj]-=ghost_offset;
            } else
            #endif
	      if (which > 0) jlist[jj] = j ^ (which << SBBITS);
          }
        }
      }
      #ifdef _LMP_INTEL_OFFLOAD
      else if (separate_buffers) {
	for (int i = ifrom; i < ito; ++i) {
          int * _noalias jlist = firstneigh + cnumneigh[i];
          const int jnum = numneigh[i];
	  int jj = 0;
	  for (jj = 0; jj < jnum; jj++)
	    if (jlist[jj] >= nlocal) break;
	  while (jj < jnum) {
	    if (jlist[jj] == nall) jlist[jj] = nall_offset;
	    else jlist[jj] -= ghost_offset;
	    jj++;
	  }
	}
      }
      #endif
    } // end omp
    #ifdef __MIC__
    *timer_compute = MIC_Wtime() - *timer_compute;
    #endif
  } // end offload

  if (offload) {
    fix->stop_watch(TIME_OFFLOAD_LATENCY);
    #ifdef _LMP_INTEL_OFFLOAD
    for (int n = 0; n < aend; n++) {
      ilist[n] = n;
      numneigh[n] = 0;
    }
    #endif
  } else {
    for (int i = astart; i < aend; i++)
      list->firstneigh[i] = firstneigh + cnumneigh[i];
    fix->stop_watch(TIME_HOST_NEIGHBOR);
    #ifdef _LMP_INTEL_OFFLOAD
    if (separate_buffers) {
      fix->start_watch(TIME_PACK);
      fix->set_neighbor_host_sizes();
      buffers->pack_sep_from_single(fix->host_min_local(),
				    fix->host_used_local(),
				    fix->host_min_ghost(),
				    fix->host_used_ghost());
      fix->stop_watch(TIME_PACK);
    }
    #endif
  }
}

/* ----------------------------------------------------------------------
   binned neighbor list construction with full Newton's 3rd law
   each owned atom i checks its own bin and other bins in Newton stencil
   every pair stored exactly once by some processor
------------------------------------------------------------------------- */

void Neighbor::half_bin_newton_intel(NeighList *list)
{
  const int nlocal = (includegroup) ? atom->nfirst : atom->nlocal;
  list->inum = nlocal;

  // Get fix for intel stuff
  FixIntel *fix = static_cast<FixIntel *>(fix_intel);

  const int off_end = fix->offload_end_neighbor();
  int host_start = fix->host_start_neighbor();;
  int offload_noghost = 0;
  #ifdef _LMP_INTEL_OFFLOAD
  if (fix->full_host_list()) host_start = 0;
  offload_noghost = fix->offload_noghost();
  if (exclude) 
    error->all(FLERR, "Exclusion lists not yet supported for Intel offload");
  #endif

  if (fix->precision() == FixIntel::PREC_MODE_MIXED) {
    if (offload_noghost) {
      hbni<float,double,1>(1, list, fix->get_mixed_buffers(),
                           0, off_end, fix);
      hbni<float,double,1>(0, list, fix->get_mixed_buffers(),
                           host_start, nlocal, fix, off_end);
    } else {
      hbni<float,double,0>(1, list, fix->get_mixed_buffers(),
                           0, off_end, fix);
      hbni<float,double,0>(0, list, fix->get_mixed_buffers(),
                           host_start, nlocal, fix);
    }
  } else if (fix->precision() == FixIntel::PREC_MODE_DOUBLE) {
    if (offload_noghost) {
      hbni<double,double,1>(1, list, fix->get_double_buffers(),
                            0, off_end, fix);
      hbni<double,double,1>(0, list, fix->get_double_buffers(),
                            host_start, nlocal, fix, off_end);
    } else {
      hbni<double,double,0>(1, list, fix->get_double_buffers(),
                            0, off_end, fix);
      hbni<double,double,0>(0, list, fix->get_double_buffers(),
                            host_start, nlocal, fix);
    }
  } else {
    if (offload_noghost) {
      hbni<float,float,1>(1, list, fix->get_single_buffers(), 0, off_end, fix);
      hbni<float,float,1>(0, list, fix->get_single_buffers(),
                          host_start, nlocal, fix, off_end);
    } else {
      hbni<float,float,0>(1, list, fix->get_single_buffers(), 0, off_end, fix);
      hbni<float,float,0>(0, list, fix->get_single_buffers(),
                          host_start, nlocal, fix);
    }
  }
}

template <class flt_t, class acc_t, int offload_noghost>
void Neighbor::hbni(const int offload, NeighList *list, void *buffers_in,
                    const int astart, const int aend, void *fix_in,
                    const int offload_end) {
  IntelBuffers<flt_t,acc_t> *buffers = (IntelBuffers<flt_t,acc_t> *)buffers_in;
  FixIntel *fix = (FixIntel *)fix_in;
  const int nall = atom->nlocal + atom->nghost;
  int pad = 1;

  if (offload) {
    fix->start_watch(TIME_PACK);
    buffers->grow(nall, atom->nlocal, comm->nthreads, aend);
    buffers->grow_nbor(list, atom->nlocal, aend);

    ATOM_T biga;
    biga.x = INTEL_BIGP;
    biga.y = INTEL_BIGP;
    biga.z = INTEL_BIGP;
    biga.w = 1;
    buffers->get_x()[nall]=biga;

    const int nthreads = comm->nthreads;
    #if defined(_OPENMP)
    #pragma omp parallel default(none) shared(buffers)
    #endif
    {
      int ifrom, ito, tid;
      IP_PRE_omp_range_id_align(ifrom, ito, tid, nall, nthreads, 
				sizeof(ATOM_T));
      buffers->thr_pack(ifrom, ito, 0);
    }
    fix->stop_watch(TIME_PACK);

    fix->start_watch(TIME_HOST_NEIGHBOR);
    bin_atoms<flt_t,acc_t>(buffers->get_x(), buffers->get_atombin());
    if (INTEL_MIC_NBOR_PAD > 1)
      pad = INTEL_MIC_NBOR_PAD * sizeof(float) / sizeof(flt_t);
  } else {
    fix->start_watch(TIME_HOST_NEIGHBOR);
    if (INTEL_NBOR_PAD > 1)
      pad = INTEL_NBOR_PAD * sizeof(float) / sizeof(flt_t);
  }    
  const int pad_width = pad;

  if (aend-astart == 0) {
    fix->stop_watch(TIME_HOST_NEIGHBOR);
    return;
  }

  const ATOM_T * _noalias const x = buffers->get_x();
  int * _noalias const firstneigh = buffers->firstneigh(list);
  int nall_t = nall;
  if (offload_noghost && offload) nall_t = atom->nlocal;
  const int e_nall = nall_t;

  const int molecular = atom->molecular;
  int *ns = NULL;
  tagint *s = NULL;
  int tag_size, special_size;
  if (molecular) {
    s = atom->special[0];
    ns = atom->nspecial[0];
    tag_size = e_nall;
    special_size = aend;
  } else {
    s = &buffers->_special_holder;
    ns = &buffers->_nspecial_holder;
    tag_size = 0;
    special_size = 0;
  }
  const tagint * _noalias const special = s;
  const int * _noalias const nspecial = ns;
  const int maxspecial = atom->maxspecial;
  const tagint * _noalias const tag = atom->tag;

  int * _noalias const ilist = list->ilist;
  int * _noalias numneigh = list->numneigh;
  int * _noalias const cnumneigh = buffers->cnumneigh(list);
  const int nstencil = list->nstencil;
  const int * _noalias const stencil = list->stencil;
  const flt_t * _noalias const cutneighsq = buffers->get_cutneighsq()[0];
  const int ntypes = atom->ntypes + 1;
  const int nlocal = atom->nlocal;

  #ifndef _LMP_INTEL_OFFLOAD
  int * const mask = atom->mask;
  tagint * const molecule = atom->molecule;
  #endif

  int tnum;
  int *overflow;
  double *timer_compute;
  if (offload) {
    timer_compute = fix->off_watch_neighbor();
    tnum = buffers->get_off_threads();
    overflow = fix->get_off_overflow_flag();
    fix->stop_watch(TIME_HOST_NEIGHBOR);
    fix->start_watch(TIME_OFFLOAD_LATENCY);
  } else {
    tnum = comm->nthreads;
    overflow = fix->get_overflow_flag();
  }
  const int nthreads = tnum;
  const int maxnbors = buffers->get_max_nbors();
  int * _noalias const atombin = buffers->get_atombin();

  // Make sure dummy coordinates to eliminate loop remainder not within cutoff
  {
    const flt_t dx = (INTEL_BIGP - bboxhi[0]);
    const flt_t dy = (INTEL_BIGP - bboxhi[1]);
    const flt_t dz = (INTEL_BIGP - bboxhi[2]);
    if (dx * dx + dy * dy + dz * dz < static_cast<flt_t>(cutneighmaxsq))
      error->one(FLERR,
	"Intel package expects no atoms within cutoff of {1e15,1e15,1e15}.");
  }

  #ifdef _LMP_INTEL_OFFLOAD
  const int * _noalias const binhead = this->binhead;
  const int * _noalias const special_flag = this->special_flag;
  const int * _noalias const bins = this->bins;
  const int cop = fix->coprocessor_number();
  const int separate_buffers = fix->separate_buffers();
  #pragma offload target(mic:cop) if(offload) \
    in(x:length(e_nall+1) alloc_if(0) free_if(0)) \
    in(tag:length(tag_size) alloc_if(0) free_if(0)) \
    in(special:length(special_size*maxspecial) alloc_if(0) free_if(0)) \
    in(nspecial:length(special_size*3) alloc_if(0) free_if(0)) \
    in(bins:length(nall) alloc_if(0) free_if(0)) \
    in(binhead:length(mbins) alloc_if(0) free_if(0)) \
    in(cutneighsq:length(0) alloc_if(0) free_if(0)) \
    in(firstneigh:length(0) alloc_if(0) free_if(0)) \
    in(cnumneigh:length(0) alloc_if(0) free_if(0)) \
    out(numneigh:length(0) alloc_if(0) free_if(0)) \
    in(ilist:length(0) alloc_if(0) free_if(0)) \
    in(atombin:length(aend) alloc_if(0) free_if(0)) \
    in(stencil:length(nstencil) alloc_if(0) free_if(0)) \
    in(special_flag:length(0) alloc_if(0) free_if(0)) \
    in(maxnbors,nthreads,maxspecial,nstencil,e_nall,offload,pad_width) \
    in(offload_end,separate_buffers,astart, aend, nlocal, molecular, ntypes) \
    out(overflow:length(5) alloc_if(0) free_if(0)) \
    out(timer_compute:length(1) alloc_if(0) free_if(0)) \
    signal(numneigh)
  #endif
  {
    #ifdef __MIC__
    *timer_compute = MIC_Wtime();
    #endif

    #ifdef _LMP_INTEL_OFFLOAD
    overflow[LMP_LOCAL_MIN] = astart;
    overflow[LMP_LOCAL_MAX] = aend - 1;
    overflow[LMP_GHOST_MIN] = e_nall;
    overflow[LMP_GHOST_MAX] = -1;
    #endif

    #if defined(_OPENMP)
    #pragma omp parallel default(none) shared(numneigh, overflow)
    #endif
    {
      #ifdef _LMP_INTEL_OFFLOAD
      int lmin = e_nall, lmax = -1, gmin = e_nall, gmax = -1;
      #endif
      
      const int num = aend - astart;
      int tid, ifrom, ito;
      IP_PRE_omp_range_id(ifrom, ito, tid, num, nthreads);
      ifrom += astart;
      ito += astart;

      int which;

      const int list_size = (ito + tid + 1) * maxnbors;
      int ct = (ifrom + tid) * maxnbors;
      int *neighptr = firstneigh + ct;
      for (int i = ifrom; i < ito; i++) {
        int j, k, n, n2, itype, jtype, ibin;
        double xtmp, ytmp, ztmp, delx, dely, delz, rsq;

        n = 0;
	n2 = maxnbors;

        xtmp = x[i].x;
        ytmp = x[i].y;
        ztmp = x[i].z;
        itype = x[i].w;
        const int ioffset = ntypes * itype;

        // loop over rest of atoms in i's bin, ghosts are at end of linked list
        // if j is owned atom, store it, since j is beyond i in linked list
        // if j is ghost, only store if j coords are "above/to the right" of i

        for (j = bins[i]; j >= 0; j = bins[j]) {
          if (j >= nlocal) {
            if (offload_noghost && offload) continue;
            if (x[j].z < ztmp) continue;
            if (x[j].z == ztmp) {
              if (x[j].y < ytmp) continue;
              if (x[j].y == ytmp && x[j].x < xtmp) continue;
            }
          } else if (offload_noghost && i < offload_end) continue;

          jtype = x[j].w;
          #ifndef _LMP_INTEL_OFFLOAD
          if (exclude && exclusion(i,j,itype,jtype,mask,molecule)) continue;
	  #endif

          delx = xtmp - x[j].x;
          dely = ytmp - x[j].y;
          delz = ztmp - x[j].z;
          rsq = delx * delx + dely * dely + delz * delz;

          if (rsq <= cutneighsq[ioffset + jtype]) {
	    if (j < nlocal) {
	      neighptr[n++] = j;
	      #ifdef _LMP_INTEL_OFFLOAD
	      if (j < lmin) lmin = j;
	      if (j > lmax) lmax = j;
              #endif
	    } else {
	      neighptr[n2++] = j;
	      #ifdef _LMP_INTEL_OFFLOAD
	      if (j < gmin) gmin = j;
	      if (j > gmax) gmax = j;
              #endif
            }
	  }
        }
        // loop over all atoms in other bins in stencil, store every pair

        ibin = atombin[i];

        for (k = 0; k < nstencil; k++) {
          for (j = binhead[ibin + stencil[k]]; j >= 0; j = bins[j]) {
            if (offload_noghost) {
              if (j < nlocal) {
                if (i < offload_end) continue;
              } else if (offload) continue;
            }

            jtype = x[j].w;
            #ifndef _LMP_INTEL_OFFLOAD
            if (exclude && exclusion(i,j,itype,jtype,mask,molecule)) continue;
	    #endif

            delx = xtmp - x[j].x;
            dely = ytmp - x[j].y;
            delz = ztmp - x[j].z;
            rsq = delx * delx + dely * dely + delz * delz;
            if (rsq <= cutneighsq[ioffset + jtype]) {
	      if (j < nlocal) {
		neighptr[n++] = j;
                #ifdef _LMP_INTEL_OFFLOAD
		if (j < lmin) lmin = j;
		if (j > lmax) lmax = j;
                #endif
	      } else {
		neighptr[n2++] = j;
	        #ifdef _LMP_INTEL_OFFLOAD
		if (j < gmin) gmin = j;
		if (j > gmax) gmax = j;
                #endif
	      }
	    }
          }
        }
        ilist[i] = i;

        cnumneigh[i] = ct;
        if (n > maxnbors) *overflow = 1;
	for (k = maxnbors; k < n2; k++) neighptr[n++] = neighptr[k];
        while( (n % pad_width) != 0 ) neighptr[n++] = e_nall;
        numneigh[i] = n;
	while((n % (INTEL_DATA_ALIGN / sizeof(int))) != 0) n++;
        ct += n;
        neighptr += n;
	if (ct + n + maxnbors > list_size) {
	  *overflow = 1;
	  ct = (ifrom + tid) * maxnbors;
	}
      }

      if (*overflow == 1)
        for (int i = ifrom; i < ito; i++)
          numneigh[i] = 0;

      #ifdef _LMP_INTEL_OFFLOAD
      if (separate_buffers) {
        #if defined(_OPENMP)
        #pragma omp critical
        #endif
        {
  	  if (lmin < overflow[LMP_LOCAL_MIN]) overflow[LMP_LOCAL_MIN] = lmin;
	  if (lmax > overflow[LMP_LOCAL_MAX]) overflow[LMP_LOCAL_MAX] = lmax;
	  if (gmin < overflow[LMP_GHOST_MIN]) overflow[LMP_GHOST_MIN] = gmin;
	  if (gmax > overflow[LMP_GHOST_MAX]) overflow[LMP_GHOST_MAX] = gmax;
        }
	#pragma omp barrier
      }

      int ghost_offset = 0, nall_offset = e_nall;
      if (separate_buffers) {
	int nghost = overflow[LMP_GHOST_MAX] + 1 - overflow[LMP_GHOST_MIN];
	if (nghost < 0) nghost = 0;
	if (offload) {
	  ghost_offset = overflow[LMP_GHOST_MIN] - overflow[LMP_LOCAL_MAX] - 1;
	  nall_offset = overflow[LMP_LOCAL_MAX] + 1 + nghost;
	} else {
	  ghost_offset = overflow[LMP_GHOST_MIN] - nlocal;
	  nall_offset = nlocal + nghost;
	}
      }
      #endif

      if (molecular) {
        for (int i = ifrom; i < ito; ++i) {
          int * _noalias jlist = firstneigh + cnumneigh[i];
          const int jnum = numneigh[i];
          for (int jj = 0; jj < jnum; jj++) {
            const int j = jlist[jj];
            ofind_special(which, special, nspecial, i, tag[j],
                          special_flag);
	    #ifdef _LMP_INTEL_OFFLOAD
	    if (j >= nlocal) {
	      if (j == e_nall)
		jlist[jj] = nall_offset;
	      else if (which > 0) 
		jlist[jj] = (j-ghost_offset) ^ (which << SBBITS);
	      else jlist[jj]-=ghost_offset;
            } else
	    #endif
            if (which > 0) jlist[jj] = j ^ (which << SBBITS);
          }
        }
      }
      #ifdef _LMP_INTEL_OFFLOAD
      else if (separate_buffers) {
	for (int i = ifrom; i < ito; ++i) {
          int * _noalias jlist = firstneigh + cnumneigh[i];
          const int jnum = numneigh[i];
	  int jj = 0;
	  for (jj = 0; jj < jnum; jj++)
	    if (jlist[jj] >= nlocal) break;
	  while (jj < jnum) {
	    if (jlist[jj] == e_nall) jlist[jj] = nall_offset;
	    else jlist[jj] -= ghost_offset;
	    jj++;
	  }
	}
      }
      #endif
    } // end omp
    #ifdef __MIC__
    *timer_compute = MIC_Wtime() - *timer_compute;
    #endif
  } // end offload

  if (offload) {
    fix->stop_watch(TIME_OFFLOAD_LATENCY);
    #ifdef _LMP_INTEL_OFFLOAD
    for (int n = 0; n < aend; n++) {
      ilist[n] = n;
      numneigh[n] = 0;
    }
    #endif
  } else {
    for (int i = astart; i < aend; i++)
      list->firstneigh[i] = firstneigh + cnumneigh[i];
    fix->stop_watch(TIME_HOST_NEIGHBOR);
    #ifdef _LMP_INTEL_OFFLOAD
    if (separate_buffers) {
      fix->start_watch(TIME_PACK);
      fix->set_neighbor_host_sizes();
      buffers->pack_sep_from_single(fix->host_min_local(),
				    fix->host_used_local(),
				    fix->host_min_ghost(),
				    fix->host_used_ghost());
      fix->stop_watch(TIME_PACK);
    }
    #endif
  }
}

/* ----------------------------------------------------------------------
   binned neighbor list construction with Newton's 3rd law for triclinic
   each owned atom i checks its own bin and other bins in triclinic stencil
   every pair stored exactly once by some processor
------------------------------------------------------------------------- */

void Neighbor::half_bin_newton_tri_intel(NeighList *list)
{
  const int nlocal = (includegroup) ? atom->nfirst : atom->nlocal;
  list->inum = nlocal;

  // Get fix for intel stuff
  FixIntel *fix = static_cast<FixIntel *>(fix_intel);

  const int off_end = fix->offload_end_neighbor();
  int host_start = fix->host_start_neighbor();
  int offload_noghost = 0;
  #ifdef _LMP_INTEL_OFFLOAD
  if (fix->full_host_list()) host_start = 0;
  offload_noghost = fix->offload_noghost();
  if (exclude) 
    error->all(FLERR, "Exclusion lists not yet supported for Intel offload");
  #endif

  if (fix->precision() == FixIntel::PREC_MODE_MIXED) {
    if (offload_noghost) {
      hbnti<float,double,1>(1, list, fix->get_mixed_buffers(),
			    0, off_end, fix);
      hbnti<float,double,1>(0, list, fix->get_mixed_buffers(),
			    host_start, nlocal, fix, off_end);
    } else {
      hbnti<float,double,0>(1, list, fix->get_mixed_buffers(),
			    0, off_end, fix);
      hbnti<float,double,0>(0, list, fix->get_mixed_buffers(),
			    host_start, nlocal, fix);
    }
  } else if (fix->precision() == FixIntel::PREC_MODE_DOUBLE) {
    if (offload_noghost) {
      hbnti<double,double,1>(1, list, fix->get_double_buffers(),
			     0, off_end, fix);
      hbnti<double,double,1>(0, list, fix->get_double_buffers(),
			     host_start, nlocal, fix, off_end);
    } else {
      hbnti<double,double,0>(1, list, fix->get_double_buffers(),
			     0, off_end, fix);
      hbnti<double,double,0>(0, list, fix->get_double_buffers(),
			     host_start, nlocal, fix);
    }
  } else {
    if (offload_noghost) {
      hbnti<float,float,1>(1, list, fix->get_single_buffers(),
			   0, off_end, fix);
      hbnti<float,float,1>(0, list, fix->get_single_buffers(),
			   host_start, nlocal, fix, off_end);
    } else {
      hbnti<float,float,0>(1, list, fix->get_single_buffers(),
			   0, off_end, fix);
      hbnti<float,float,0>(0, list, fix->get_single_buffers(),
			   host_start, nlocal, fix);
    }
  }
}

template <class flt_t, class acc_t, int offload_noghost>
void Neighbor::hbnti(const int offload, NeighList *list, void *buffers_in,
                     const int astart, const int aend, void *fix_in,
		     const int offload_end) {
  IntelBuffers<flt_t,acc_t> *buffers = (IntelBuffers<flt_t,acc_t> *)buffers_in;
  FixIntel *fix = (FixIntel *)fix_in;
  const int nall = atom->nlocal + atom->nghost;
  int pad = 1;

  if (offload) {
    fix->start_watch(TIME_PACK);
    buffers->grow(nall, atom->nlocal, comm->nthreads, aend);
    buffers->grow_nbor(list, atom->nlocal, aend);

    ATOM_T biga;
    biga.x = INTEL_BIGP;
    biga.y = INTEL_BIGP;
    biga.z = INTEL_BIGP;
    biga.w = 1;
    buffers->get_x()[nall]=biga;

    const int nthreads = comm->nthreads;
    #if defined(_OPENMP)
    #pragma omp parallel default(none) shared(buffers)
    #endif
    {
      int ifrom, ito, tid;
      IP_PRE_omp_range_id_align(ifrom, ito, tid, nall, nthreads, 
				sizeof(ATOM_T));
      buffers->thr_pack(ifrom, ito, 0);
    }
    fix->stop_watch(TIME_PACK);

    fix->start_watch(TIME_HOST_NEIGHBOR);
    bin_atoms<flt_t,acc_t>(buffers->get_x(), buffers->get_atombin());
    if (INTEL_MIC_NBOR_PAD > 1)
      pad = INTEL_MIC_NBOR_PAD * sizeof(float) / sizeof(flt_t);
  } else {
    fix->start_watch(TIME_HOST_NEIGHBOR);
    if (INTEL_NBOR_PAD > 1)
      pad = INTEL_NBOR_PAD * sizeof(float) / sizeof(flt_t);
  }
  const int pad_width = pad;

  if (aend-astart == 0) {
    fix->stop_watch(TIME_HOST_NEIGHBOR);
    return;
  }

  const ATOM_T * _noalias const x = buffers->get_x();
  int * _noalias const firstneigh = buffers->firstneigh(list);
  int nall_t = nall;
  if (offload_noghost && offload) nall_t = atom->nlocal;
  const int e_nall = nall_t;

  const int molecular = atom->molecular;
  int *ns = NULL;
  tagint *s = NULL;
  int tag_size, special_size;
  if (molecular) {
    s = atom->special[0];
    ns = atom->nspecial[0];
    tag_size = e_nall;
    special_size = aend;
  } else {
    s = &buffers->_special_holder;
    ns = &buffers->_nspecial_holder;
    tag_size = 0;
    special_size = 0;
  }
  const tagint * _noalias const special = s;
  const int * _noalias const nspecial = ns;
  const int maxspecial = atom->maxspecial;
  const tagint * _noalias const tag = atom->tag;

  int * _noalias const ilist = list->ilist;
  int * _noalias numneigh = list->numneigh;
  int * _noalias const cnumneigh = buffers->cnumneigh(list);
  const int nstencil = list->nstencil;
  const int * _noalias const stencil = list->stencil;
  const flt_t * _noalias const cutneighsq = buffers->get_cutneighsq()[0];
  const int ntypes = atom->ntypes + 1;
  const int nlocal = atom->nlocal;

  #ifndef _LMP_INTEL_OFFLOAD
  int * const mask = atom->mask;
  tagint * const molecule = atom->molecule;
  #endif

  int tnum;
  int *overflow;
  double *timer_compute;
  if (offload) {
    timer_compute = fix->off_watch_neighbor();
    tnum = buffers->get_off_threads();
    overflow = fix->get_off_overflow_flag();
    fix->stop_watch(TIME_HOST_NEIGHBOR);
    fix->start_watch(TIME_OFFLOAD_LATENCY);
  } else {
    tnum = comm->nthreads;
    overflow = fix->get_overflow_flag();
  }
  const int nthreads = tnum;
  const int maxnbors = buffers->get_max_nbors();
  int * _noalias const atombin = buffers->get_atombin();

  // Make sure dummy coordinates to eliminate loop remainder not within cutoff
  {
    const flt_t dx = (INTEL_BIGP - bboxhi[0]);
    const flt_t dy = (INTEL_BIGP - bboxhi[1]);
    const flt_t dz = (INTEL_BIGP - bboxhi[2]);
    if (dx * dx + dy * dy + dz * dz < static_cast<flt_t>(cutneighmaxsq))
      error->one(FLERR,
	"Intel package expects no atoms within cutoff of {1e15,1e15,1e15}.");
  }

  #ifdef _LMP_INTEL_OFFLOAD
  const int * _noalias const binhead = this->binhead;
  const int * _noalias const special_flag = this->special_flag;
  const int * _noalias const bins = this->bins;
  const int cop = fix->coprocessor_number();
  const int separate_buffers = fix->separate_buffers();
  #pragma offload target(mic:cop) if(offload) \
    in(x:length(e_nall+1) alloc_if(0) free_if(0)) \
    in(tag:length(tag_size) alloc_if(0) free_if(0)) \
    in(special:length(special_size*maxspecial) alloc_if(0) free_if(0)) \
    in(nspecial:length(special_size*3) alloc_if(0) free_if(0)) \
    in(bins:length(nall) alloc_if(0) free_if(0)) \
    in(binhead:length(mbins) alloc_if(0) free_if(0)) \
    in(cutneighsq:length(0) alloc_if(0) free_if(0)) \
    in(firstneigh:length(0) alloc_if(0) free_if(0)) \
    in(cnumneigh:length(0) alloc_if(0) free_if(0)) \
    out(numneigh:length(0) alloc_if(0) free_if(0)) \
    in(ilist:length(0) alloc_if(0) free_if(0)) \
    in(atombin:length(aend) alloc_if(0) free_if(0)) \
    in(stencil:length(nstencil) alloc_if(0) free_if(0)) \
    in(special_flag:length(0) alloc_if(0) free_if(0)) \
    in(maxnbors,nthreads,maxspecial,nstencil,offload_end,pad_width,e_nall) \
    in(offload,separate_buffers, astart, aend, nlocal, molecular, ntypes) \
    out(overflow:length(5) alloc_if(0) free_if(0)) \
    out(timer_compute:length(1) alloc_if(0) free_if(0)) \
    signal(numneigh)
  #endif
  {
    #ifdef __MIC__
    *timer_compute = MIC_Wtime();
    #endif

    #ifdef _LMP_INTEL_OFFLOAD
    overflow[LMP_LOCAL_MIN] = astart;
    overflow[LMP_LOCAL_MAX] = aend - 1;
    overflow[LMP_GHOST_MIN] = e_nall;
    overflow[LMP_GHOST_MAX] = -1;
    #endif

    #if defined(_OPENMP)
    #pragma omp parallel default(none) shared(numneigh, overflow)
    #endif
    {
      #ifdef _LMP_INTEL_OFFLOAD
      int lmin = e_nall, lmax = -1, gmin = e_nall, gmax = -1;
      #endif

      const int num = aend-astart;
      int tid, ifrom, ito;
      IP_PRE_omp_range_id(ifrom,ito,tid,num,nthreads);
      ifrom += astart;
      ito += astart;

      int which;

      const int list_size = (ito + tid + 1) * maxnbors;
      int ct = (ifrom + tid) * maxnbors;
      int *neighptr = firstneigh + ct;
      for (int i = ifrom; i < ito; i++) {
        int j, k, n, n2, itype, jtype, ibin;
        double xtmp, ytmp, ztmp, delx, dely, delz, rsq;

        n = 0;
        n2 = maxnbors;

        xtmp = x[i].x;
        ytmp = x[i].y;
        ztmp = x[i].z;
        itype = x[i].w;
        const int ioffset = ntypes * itype;

        // loop over all atoms in bins in stencil
        // pairs for atoms j "below" i are excluded
        // below = lower z or (equal z and lower y) or (equal zy and lower x)
        //         (equal zyx and j <= i)
        // latter excludes self-self interaction but allows superposed atoms

        ibin = atombin[i];

        for (k = 0; k < nstencil; k++) {
          for (j = binhead[ibin + stencil[k]]; j >= 0; j = bins[j]) {
	    if (offload_noghost) {
              if (j < nlocal) {
                if (i < offload_end) continue;
              } else if (offload) continue;
            }

            if (x[j].z < ztmp) continue;
            if (x[j].z == ztmp) {
              if (x[j].y < ytmp) continue;
              if (x[j].y == ytmp) {
                if (x[j].x < xtmp) continue;
                if (x[j].x == xtmp && j <= i) continue;
              }
            }

            jtype = x[j].w;
            #ifndef _LMP_INTEL_OFFLOAD
            if (exclude && exclusion(i,j,itype,jtype,mask,molecule)) continue;
	    #endif

            delx = xtmp - x[j].x;
            dely = ytmp - x[j].y;
            delz = ztmp - x[j].z;
            rsq = delx * delx + dely * dely + delz * delz;
            if (rsq <= cutneighsq[ioffset + jtype]) {
              if (j < nlocal) {
                neighptr[n++] = j;
                #ifdef _LMP_INTEL_OFFLOAD
		if (j < lmin) lmin = j;
		if (j > lmax) lmax = j;
                #endif
	      }  else {
                neighptr[n2++] = j;
  	        #ifdef _LMP_INTEL_OFFLOAD
		if (j < gmin) gmin = j;
		if (j > gmax) gmax = j;
                #endif
	      }
            }
          }
        }
        ilist[i] = i;

        cnumneigh[i] = ct;
        if (n > maxnbors) *overflow = 1;
        for (k = maxnbors; k < n2; k++) neighptr[n++] = neighptr[k];
        while( (n % pad_width) != 0 ) neighptr[n++] = e_nall;
        numneigh[i] = n;
        while((n % (INTEL_DATA_ALIGN / sizeof(int))) != 0) n++;
        ct += n;
        neighptr += n;
        if (ct + n + maxnbors > list_size) {
          *overflow = 1;
          ct = (ifrom + tid) * maxnbors;
        }
      }

      if (*overflow == 1)
        for (int i = ifrom; i < ito; i++)
          numneigh[i] = 0;

      #ifdef _LMP_INTEL_OFFLOAD
      if (separate_buffers) {
        #if defined(_OPENMP)
        #pragma omp critical
        #endif
        {
          if (lmin < overflow[LMP_LOCAL_MIN]) overflow[LMP_LOCAL_MIN] = lmin;
          if (lmax > overflow[LMP_LOCAL_MAX]) overflow[LMP_LOCAL_MAX] = lmax;
          if (gmin < overflow[LMP_GHOST_MIN]) overflow[LMP_GHOST_MIN] = gmin;
          if (gmax > overflow[LMP_GHOST_MAX]) overflow[LMP_GHOST_MAX] = gmax;
        }
        #pragma omp barrier
      }

      int ghost_offset = 0, nall_offset = e_nall;
      if (separate_buffers) {
        int nghost = overflow[LMP_GHOST_MAX] + 1 - overflow[LMP_GHOST_MIN];
        if (nghost < 0) nghost = 0;
        if (offload) {
          ghost_offset = overflow[LMP_GHOST_MIN] - overflow[LMP_LOCAL_MAX] - 1;
          nall_offset = overflow[LMP_LOCAL_MAX] + 1 + nghost;
	} else {
          ghost_offset = overflow[LMP_GHOST_MIN] - nlocal;
          nall_offset = nlocal + nghost;
        }
      }
      #endif

      if (molecular) {
        for (int i = ifrom; i < ito; ++i) {
          int * _noalias jlist = firstneigh + cnumneigh[i];
          const int jnum = numneigh[i];
          for (int jj = 0; jj < jnum; jj++) {
            const int j = jlist[jj];
            ofind_special(which, special, nspecial, i, tag[j], special_flag);
            #ifdef _LMP_INTEL_OFFLOAD
	    if (j >= nlocal) {
	      if (j == e_nall) 
		jlist[jj] = nall_offset;
	      else if (which > 0) 
		jlist[jj] = (j-ghost_offset) ^ (which << SBBITS);
	      else jlist[jj]-=ghost_offset;
            } else
            #endif
	      if (which > 0) jlist[jj] = j ^ (which << SBBITS);
          }
        }
      }
      #ifdef _LMP_INTEL_OFFLOAD
      else if (separate_buffers) {
	for (int i = ifrom; i < ito; ++i) {
          int * _noalias jlist = firstneigh + cnumneigh[i];
          const int jnum = numneigh[i];
	  int jj = 0;
	  for (jj = 0; jj < jnum; jj++)
	    if (jlist[jj] >= nlocal) break;
	  while (jj < jnum) {
	    if (jlist[jj] == e_nall) jlist[jj] = nall_offset;
	    else jlist[jj] -= ghost_offset;
	    jj++;
	  }
	}
      }
      #endif
    } // end omp
    #ifdef __MIC__
    *timer_compute = MIC_Wtime() - *timer_compute;
    #endif
  } // end offload

  if (offload) {
    fix->stop_watch(TIME_OFFLOAD_LATENCY);
    #ifdef _LMP_INTEL_OFFLOAD
    for (int n = 0; n < aend; n++) {
      ilist[n] = n;
      numneigh[n] = 0;
    }
    #endif
  } else {
    for (int i = astart; i < aend; i++)
      list->firstneigh[i] = firstneigh + cnumneigh[i];
    fix->stop_watch(TIME_HOST_NEIGHBOR);
    #ifdef _LMP_INTEL_OFFLOAD
    if (separate_buffers) {
      fix->start_watch(TIME_PACK);
      fix->set_neighbor_host_sizes();
      buffers->pack_sep_from_single(fix->host_min_local(),
				    fix->host_used_local(),
				    fix->host_min_ghost(),
				    fix->host_used_ghost());
      fix->stop_watch(TIME_PACK);
    }
    #endif
  }
}
