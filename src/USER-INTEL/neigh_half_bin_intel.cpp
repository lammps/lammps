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

//#define OUTER_CHUNK 1

#include "neighbor.h"
#include "neigh_list.h"
#include "atom.h"
#include "comm.h"
#include "group.h"
#include "fix_intel.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

#ifdef LMP_USE_AVXCD
#include "intel_simd.h"
#endif

#ifdef OUTER_CHUNK
#include "intel_simd.h"
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
        which = 1;						      \
      } else if (s < n2) {                                            \
        which = 2;						      \
      } else {                                                        \
        which = 3;						      \
      }                                                               \
    }                                                                 \
  }                                                                   \
}

#define ominimum_image_check(answer, dx, dy, dz)                      \
{								      \
  answer = 0;							      \
  if (xperiodic && fabs(dx) > xprd_half) answer = 1;		      \
  if (yperiodic && fabs(dy) > yprd_half) answer = 1;		      \
  if (zperiodic && fabs(dz) > zprd_half) answer = 1;		      \
}

#define dminimum_image_check(answer, dx, dy, dz)                      \
{								      \
  answer = 0;							      \
  if (domain->xperiodic && fabs(dx) > domain->xprd_half) answer = 1;  \
  if (domain->yperiodic && fabs(dy) > domain->yprd_half) answer = 1;  \
  if (domain->zperiodic && fabs(dz) > domain->zprd_half) answer = 1;  \
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

  int need_ic = 0;
  if (atom->molecular)
    dminimum_image_check(need_ic, cutneighmax, cutneighmax, cutneighmax);

  if (need_ic) {
    if (fix->precision() == FixIntel::PREC_MODE_MIXED) {
      hbnni<float,double,1>(1, list, fix->get_mixed_buffers(),
			    0, off_end, fix);
      hbnni<float,double,1>(0, list, fix->get_mixed_buffers(),
			    host_start, nlocal,fix);
    } else if (fix->precision() == FixIntel::PREC_MODE_DOUBLE) {
      hbnni<double,double,1>(1, list, fix->get_double_buffers(),
			     0, off_end, fix);
      hbnni<double,double,1>(0, list, fix->get_double_buffers(),
			     host_start, nlocal, fix);
    } else {
      hbnni<float,float,1>(1, list, fix->get_single_buffers(),
			   0, off_end, fix);
      hbnni<float,float,1>(0, list, fix->get_single_buffers(),
			   host_start, nlocal, fix);
    }
  } else {
    if (fix->precision() == FixIntel::PREC_MODE_MIXED) {
      hbnni<float,double,0>(1, list, fix->get_mixed_buffers(),
			    0, off_end, fix);
      hbnni<float,double,0>(0, list, fix->get_mixed_buffers(),
			    host_start, nlocal,fix);
    } else if (fix->precision() == FixIntel::PREC_MODE_DOUBLE) {
      hbnni<double,double,0>(1, list, fix->get_double_buffers(),
			     0, off_end, fix);
      hbnni<double,double,0>(0, list, fix->get_double_buffers(),
			     host_start, nlocal, fix);
    } else {
      hbnni<float,float,0>(1, list, fix->get_single_buffers(),
			   0, off_end, fix);
      hbnni<float,float,0>(0, list, fix->get_single_buffers(),
			   host_start, nlocal, fix);
    }
  }
}

template <class flt_t, class acc_t, int need_ic>
void Neighbor::hbnni(const int offload, NeighList *list, void *buffers_in,
                     const int astart, const int aend, void *fix_in) {
  IntelBuffers<flt_t,acc_t> *buffers = (IntelBuffers<flt_t,acc_t> *)buffers_in;
  FixIntel *fix = (FixIntel *)fix_in;
  const int nall = atom->nlocal + atom->nghost;
  int pad = 1;

  if (offload) {
    fix->start_watch(TIME_PACK);
    buffers->grow(nall, atom->nlocal, comm->nthreads, aend);
    buffers->grow_nbor(list, atom->nlocal, comm->nthreads, aend);

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
  int tag_size = 0, special_size;
  if (buffers->need_tag()) tag_size = nall;
  if (molecular) {
    s = atom->special[0];
    ns = atom->nspecial[0];
    special_size = aend;
  } else {
    s = &buffers->_special_holder;
    ns = &buffers->_nspecial_holder;
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

  const int xperiodic = domain->xperiodic;
  const int yperiodic = domain->yperiodic;
  const int zperiodic = domain->zperiodic;
  const flt_t xprd_half = domain->xprd_half;
  const flt_t yprd_half = domain->yprd_half;
  const flt_t zprd_half = domain->zprd_half;

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
    in(maxnbors,nthreads,maxspecial,nstencil,pad_width,offload)  \
    in(separate_buffers, astart, aend, nlocal, molecular, ntypes) \
    in(xperiodic, yperiodic, zperiodic, xprd_half, yprd_half, zprd_half) \
    out(overflow:length(5) alloc_if(0) free_if(0)) \
    out(timer_compute:length(1) alloc_if(0) free_if(0)) \
    signal(tag)
  #endif
  {
    #if defined(__MIC__) && defined(_LMP_INTEL_OFFLOAD)
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
                if (need_ic) {
                  int no_special;
                  ominimum_image_check(no_special, delx, dely, delz);
                  if (no_special)
                    neighptr[n++] = -j - 1;
		  else
                    neighptr[n++] = j;
                } else
                  neighptr[n++] = j;
                #ifdef _LMP_INTEL_OFFLOAD
		if (j < lmin) lmin = j;
		if (j > lmax) lmax = j;
                #endif
              } else {
                if (need_ic) {
                  int no_special;
                  ominimum_image_check(no_special, delx, dely, delz);
                  if (no_special)
                    neighptr[n2++] = -j - 1;
		  else
                    neighptr[n2++] = j;
                } else
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
	    if (need_ic && j < 0) {
	      which = 0;
	      jlist[jj] = -j - 1;
	    } else
              ofind_special(which, special, nspecial, i, tag[j], special_flag);
            #ifdef _LMP_INTEL_OFFLOAD
	    if (j >= nlocal) {
	      if (j == nall)
		jlist[jj] = nall_offset;
	      else if (which) 
		jlist[jj] = (j-ghost_offset) ^ (which << SBBITS);
	      else jlist[jj]-=ghost_offset;
            } else
            #endif
	      if (which) jlist[jj] = j ^ (which << SBBITS);
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
    #if defined(__MIC__) && defined(_LMP_INTEL_OFFLOAD)
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

  int need_ic = 0;
  if (atom->molecular)
    dminimum_image_check(need_ic, cutneighmax, cutneighmax, cutneighmax);

  if (need_ic) {
    if (fix->precision() == FixIntel::PREC_MODE_MIXED) {
      #ifdef _LMP_INTEL_OFFLOAD
      if (offload_noghost) {
	hbni<float,double,1,1>(1, list, fix->get_mixed_buffers(),
	                       0, off_end, fix);
	hbni<float,double,1,1>(0, list, fix->get_mixed_buffers(),
	                       host_start, nlocal, fix, off_end);
      } else
      #endif
      {
	hbni<float,double,0,1>(1, list, fix->get_mixed_buffers(),
	                       0, off_end, fix);
	hbni<float,double,0,1>(0, list, fix->get_mixed_buffers(),
	                       host_start, nlocal, fix);
      }
    } else if (fix->precision() == FixIntel::PREC_MODE_DOUBLE) {
      #ifdef _LMP_INTEL_OFFLOAD
      if (offload_noghost) {
        hbni<double,double,1,1>(1, list, fix->get_double_buffers(),
                                0, off_end, fix);
        hbni<double,double,1,1>(0, list, fix->get_double_buffers(),
                                host_start, nlocal, fix, off_end);
      } else
      #endif
      {
        hbni<double,double,0,1>(1, list, fix->get_double_buffers(),
                                0, off_end, fix);
        hbni<double,double,0,1>(0, list, fix->get_double_buffers(),
                                host_start, nlocal, fix);
      }
    } else {
      #ifdef _LMP_INTEL_OFFLOAD
      if (offload_noghost) {
        hbni<float,float,1,1>(1, list, fix->get_single_buffers(), 0, off_end,
	                      fix);
        hbni<float,float,1,1>(0, list, fix->get_single_buffers(),
                              host_start, nlocal, fix, off_end);
      } else
      #endif
      {
        hbni<float,float,0,1>(1, list, fix->get_single_buffers(), 0, off_end,
	                      fix);
        hbni<float,float,0,1>(0, list, fix->get_single_buffers(),
                              host_start, nlocal, fix);
      }
    }
  } else {
    if (fix->precision() == FixIntel::PREC_MODE_MIXED) {
      #ifdef _LMP_INTEL_OFFLOAD
      if (offload_noghost) {
	hbni<float,double,1,0>(1, list, fix->get_mixed_buffers(),
	                       0, off_end, fix);
	hbni<float,double,1,0>(0, list, fix->get_mixed_buffers(),
	                       host_start, nlocal, fix, off_end);
      } else
      #endif
      {
	hbni<float,double,0,0>(1, list, fix->get_mixed_buffers(),
	                       0, off_end, fix);
	hbni<float,double,0,0>(0, list, fix->get_mixed_buffers(),
	                       host_start, nlocal, fix);
      }
    } else if (fix->precision() == FixIntel::PREC_MODE_DOUBLE) {
      #ifdef _LMP_INTEL_OFFLOAD
      if (offload_noghost) {
        hbni<double,double,1,0>(1, list, fix->get_double_buffers(),
                                0, off_end, fix);
        hbni<double,double,1,0>(0, list, fix->get_double_buffers(),
                                host_start, nlocal, fix, off_end);
      } else
      #endif
      {
        hbni<double,double,0,0>(1, list, fix->get_double_buffers(),
                                0, off_end, fix);
        hbni<double,double,0,0>(0, list, fix->get_double_buffers(),
                                host_start, nlocal, fix);
      }
    } else {
      #ifdef _LMP_INTEL_OFFLOAD
      if (offload_noghost) {
        hbni<float,float,1,0>(1, list, fix->get_single_buffers(), 0, off_end,
	                      fix);
        hbni<float,float,1,0>(0, list, fix->get_single_buffers(),
                              host_start, nlocal, fix, off_end);
      } else
      #endif
      {
        hbni<float,float,0,0>(1, list, fix->get_single_buffers(), 0, off_end,
	                      fix);
        hbni<float,float,0,0>(0, list, fix->get_single_buffers(),
                              host_start, nlocal, fix);
      }
    }
  }
}

template <class flt_t, class acc_t, int offload_noghost, int need_ic>
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
    buffers->grow_nbor(list, atom->nlocal, comm->nthreads, aend);

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
  int tag_size = 0, special_size;
  if (buffers->need_tag()) tag_size = e_nall;
  if (molecular) {
    s = atom->special[0];
    ns = atom->nspecial[0];
    special_size = aend;
  } else {
    s = &buffers->_special_holder;
    ns = &buffers->_nspecial_holder;
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

  const int xperiodic = domain->xperiodic;
  const int yperiodic = domain->yperiodic;
  const int zperiodic = domain->zperiodic;
  const flt_t xprd_half = domain->xprd_half;
  const flt_t yprd_half = domain->yprd_half;
  const flt_t zprd_half = domain->zprd_half;

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
    in(xperiodic, yperiodic, zperiodic, xprd_half, yprd_half, zprd_half) \
    out(overflow:length(5) alloc_if(0) free_if(0)) \
    out(timer_compute:length(1) alloc_if(0) free_if(0)) \
    signal(tag)
  #endif
  {
    #if defined(__MIC__) && defined(_LMP_INTEL_OFFLOAD)
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

      #ifdef OUTER_CHUNK
      const int swidth = ip_simd::SIMD_type<flt_t>::width();
      IP_PRE_omp_range_id_vec(ifrom, ito, tid, num, nthreads, swidth);
      ifrom += astart;
      ito += astart;
      int e_ito = ito;
      if (ito == num) {
	int imod = ito % swidth;
	if (imod) e_ito += swidth - imod;
      }
      const int list_size = (e_ito + tid * 2 + 2) * maxnbors;
      #else
      const int swidth = 1;
      IP_PRE_omp_range_id(ifrom, ito, tid, num, nthreads);
      ifrom += astart;
      ito += astart;
      const int list_size = (ito + tid * 2 + 2) * maxnbors;
      #endif

      int which;

      int pack_offset = maxnbors * swidth;
      int ct = (ifrom + tid * 2) * maxnbors;
      int *neighptr = firstneigh + ct;
      const int obound = pack_offset + maxnbors * 2;

      int max_chunk = 0;
      int lane = 0;
      for (int i = ifrom; i < ito; i++) {
        const flt_t xtmp = x[i].x;
        const flt_t ytmp = x[i].y;
        const flt_t ztmp = x[i].z;
        const int itype = x[i].w;
        const int ioffset = ntypes * itype;

        // loop over rest of atoms in i's bin, ghosts are at end of linked list
        // if j is owned atom, store it, since j is beyond i in linked list
        // if j is ghost, only store if j coords are "above/to the right" of i

	int raw_count = pack_offset;
        for (int j = bins[i]; j >= 0; j = bins[j]) {
          if (j >= nlocal) {
            if (offload_noghost && offload) continue;
            if (x[j].z < ztmp) continue;
            if (x[j].z == ztmp) {
              if (x[j].y < ytmp) continue;
              if (x[j].y == ytmp && x[j].x < xtmp) continue;
            }
          } else if (offload_noghost && i < offload_end) continue;

          #ifndef _LMP_INTEL_OFFLOAD
          if (exclude) {
	    const int jtype = x[j].w;
	    if (exclusion(i,j,itype,jtype,mask,molecule)) continue;
	  }
	  #endif

	  neighptr[raw_count++] = j;
	}

        // loop over all atoms in other bins in stencil, store every pair

        const int ibin = atombin[i];
        for (int k = 0; k < nstencil; k++) {
          for (int j = binhead[ibin + stencil[k]]; j >= 0; j = bins[j]) {
            if (offload_noghost) {
              if (j < nlocal) {
                if (i < offload_end) continue;
              } else if (offload) continue;
            }

            #ifndef _LMP_INTEL_OFFLOAD
            if (exclude) {
	      const int jtype = x[j].w;
	      if (exclusion(i,j,itype,jtype,mask,molecule)) continue;
	    }
	    #endif

	    neighptr[raw_count++] = j;
	  }
	}

	if (raw_count > obound) *overflow = 1;

        #if defined(LMP_SIMD_COMPILER)
	#ifdef _LMP_INTEL_OFFLOAD
	int vlmin = lmin, vlmax = lmax, vgmin = gmin, vgmax = gmax;
	#pragma vector aligned
        #pragma simd reduction(max:vlmax,vgmax) reduction(min:vlmin, vgmin)
	#else
	#pragma vector aligned
        #pragma simd
	#endif
	#endif
	for (int u = pack_offset; u < raw_count; u++) {
          int j = neighptr[u];
          const flt_t delx = xtmp - x[j].x;
          const flt_t dely = ytmp - x[j].y;
          const flt_t delz = ztmp - x[j].z;
	  const int jtype = x[j].w;
          const flt_t rsq = delx * delx + dely * dely + delz * delz;
          if (rsq > cutneighsq[ioffset + jtype]) 
	    neighptr[u] = e_nall;
	  else {
	    if (need_ic) {
	      int no_special;
	      ominimum_image_check(no_special, delx, dely, delz);
	      if (no_special)
		neighptr[u] = -j - 1;
	    }
            #ifdef _LMP_INTEL_OFFLOAD
            if (j < nlocal) {
              if (j < vlmin) vlmin = j;
              if (j > vlmax) vlmax = j;
            } else {
              if (j < vgmin) vgmin = j;
              if (j > vgmax) vgmax = j;
            }
            #endif
	  }
	}
        #ifdef _LMP_INTEL_OFFLOAD
	lmin = MIN(lmin,vlmin);
	gmin = MIN(gmin,vgmin);
	lmax = MAX(lmax,vlmax);
	gmax = MAX(gmax,vgmax);
        #endif

        int n = lane, n2 = pack_offset;
	for (int u = pack_offset; u < raw_count; u++) {
	  const int j = neighptr[u];
	  int pj = j;
	  if (pj < e_nall) {
	    if (need_ic)
	      if (pj < 0) pj = -pj - 1;

	    if (pj < nlocal) {
	      neighptr[n] = j;
	      n += swidth;
	    } else
	      neighptr[n2++] = j;
	  }
	}
	int ns = (n - lane) / swidth;
	for (int u = pack_offset; u < n2; u++) {
	  neighptr[n] = neighptr[u];
	  n += swidth;
	}

        ilist[i] = i;
        cnumneigh[i] = ct + lane;
	ns += n2 - pack_offset;
	#ifndef OUTER_CHUNK
        while( (ns % pad_width) != 0 ) neighptr[ns++] = e_nall;
	#endif
        numneigh[i] = ns;

	#ifdef OUTER_CHUNK
	if (ns > max_chunk) max_chunk = ns;
	lane++;
	if (lane == swidth) {
	  ct += max_chunk * swidth;
	  const int alignb = (INTEL_DATA_ALIGN / sizeof(int));
	  const int edge = (ct % alignb);
	  if (edge) ct += alignb - edge;
	  neighptr = firstneigh + ct;
	  max_chunk = 0;
	  pack_offset = maxnbors * swidth;
	  lane = 0;
	  if (ct + obound > list_size) {
            if (i < ito - 1) {
	      *overflow = 1;
	      ct = (ifrom + tid * 2) * maxnbors;
            }
	  }
	}
	#else
	ct += ns;
	const int alignb = (INTEL_DATA_ALIGN / sizeof(int));
	const int edge = (ct % alignb);
	if (edge) ct += alignb - edge;
	neighptr = firstneigh + ct;
	if (ct + obound > list_size) {
	  if (i < ito - 1) {
	    *overflow = 1;
	    ct = (ifrom + tid * 2) * maxnbors;
	  }
	}
	#endif
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
	  #ifndef OUTER_CHUNK
          #if defined(LMP_SIMD_COMPILER)
	  #pragma vector aligned
          #pragma simd
	  #endif
          for (int jj = 0; jj < jnum; jj++) {
	  #else
	  const int trip = jnum * swidth;
          for (int jj = 0; jj < trip; jj+= swidth) {
	  #endif
            const int j = jlist[jj];
	    if (need_ic && j < 0) {
	      which = 0;
	      jlist[jj] = -j - 1;
            } else
              ofind_special(which, special, nspecial, i, tag[j],
                            special_flag);
	    #ifdef _LMP_INTEL_OFFLOAD
	    if (j >= nlocal) {
	      if (j == e_nall)
		jlist[jj] = nall_offset;
	      else if (which) 
		jlist[jj] = (j-ghost_offset) ^ (which << SBBITS);
	      else jlist[jj]-=ghost_offset;
            } else
	    #endif
            if (which) jlist[jj] = j ^ (which << SBBITS);
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
    #if defined(__MIC__) && defined(_LMP_INTEL_OFFLOAD)
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

  int need_ic = 0;
  if (atom->molecular)
    dminimum_image_check(need_ic, cutneighmax, cutneighmax, cutneighmax);

  if (need_ic) {
    if (fix->precision() == FixIntel::PREC_MODE_MIXED) {
      #ifdef _LMP_INTEL_OFFLOAD
      if (offload_noghost) {
        hbnti<float,double,1,1>(1, list, fix->get_mixed_buffers(),
  			        0, off_end, fix);
        hbnti<float,double,1,1>(0, list, fix->get_mixed_buffers(),
  			        host_start, nlocal, fix, off_end);
      } else
      #endif
      {
        hbnti<float,double,0,1>(1, list, fix->get_mixed_buffers(),
	  		        0, off_end, fix);
        hbnti<float,double,0,1>(0, list, fix->get_mixed_buffers(),
			        host_start, nlocal, fix);
      }
    } else if (fix->precision() == FixIntel::PREC_MODE_DOUBLE) {
      #ifdef _LMP_INTEL_OFFLOAD
      if (offload_noghost) {
        hbnti<double,double,1,1>(1, list, fix->get_double_buffers(),
	  		         0, off_end, fix);
        hbnti<double,double,1,1>(0, list, fix->get_double_buffers(),
			         host_start, nlocal, fix, off_end);
      } else
      #endif
      {
        hbnti<double,double,0,1>(1, list, fix->get_double_buffers(),
			         0, off_end, fix);
        hbnti<double,double,0,1>(0, list, fix->get_double_buffers(),
			         host_start, nlocal, fix);
      }
    } else {
      #ifdef _LMP_INTEL_OFFLOAD
      if (offload_noghost) {
        hbnti<float,float,1,1>(1, list, fix->get_single_buffers(),
	  		       0, off_end, fix);
        hbnti<float,float,1,1>(0, list, fix->get_single_buffers(),
			       host_start, nlocal, fix, off_end);
      } else
      #endif
      {
        hbnti<float,float,0,1>(1, list, fix->get_single_buffers(),
			       0, off_end, fix);
        hbnti<float,float,0,1>(0, list, fix->get_single_buffers(),
			       host_start, nlocal, fix);
      }
    }
  } else {
    if (fix->precision() == FixIntel::PREC_MODE_MIXED) {
      #ifdef _LMP_INTEL_OFFLOAD
      if (offload_noghost) {
        hbnti<float,double,1,0>(1, list, fix->get_mixed_buffers(),
  			        0, off_end, fix);
        hbnti<float,double,1,0>(0, list, fix->get_mixed_buffers(),
  			        host_start, nlocal, fix, off_end);
      } else
      #endif
      {
        hbnti<float,double,0,0>(1, list, fix->get_mixed_buffers(),
	  		        0, off_end, fix);
        hbnti<float,double,0,0>(0, list, fix->get_mixed_buffers(),
			        host_start, nlocal, fix);
      }
    } else if (fix->precision() == FixIntel::PREC_MODE_DOUBLE) {
      #ifdef _LMP_INTEL_OFFLOAD
      if (offload_noghost) {
        hbnti<double,double,1,0>(1, list, fix->get_double_buffers(),
	  		         0, off_end, fix);
        hbnti<double,double,1,0>(0, list, fix->get_double_buffers(),
			         host_start, nlocal, fix, off_end);
      } else
      #endif
      {
        hbnti<double,double,0,0>(1, list, fix->get_double_buffers(),
			         0, off_end, fix);
        hbnti<double,double,0,0>(0, list, fix->get_double_buffers(),
			         host_start, nlocal, fix);
      }
    } else {
      #ifdef _LMP_INTEL_OFFLOAD
      if (offload_noghost) {
        hbnti<float,float,1,0>(1, list, fix->get_single_buffers(),
	  		       0, off_end, fix);
        hbnti<float,float,1,0>(0, list, fix->get_single_buffers(),
			       host_start, nlocal, fix, off_end);
      } else
      #endif
      {
        hbnti<float,float,0,0>(1, list, fix->get_single_buffers(),
			       0, off_end, fix);
        hbnti<float,float,0,0>(0, list, fix->get_single_buffers(),
			       host_start, nlocal, fix);
      }
    }
  }
}

template <class flt_t, class acc_t, int offload_noghost, int need_ic>
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
    buffers->grow_nbor(list, atom->nlocal, comm->nthreads, aend);

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
  int tag_size = 0, special_size;
  if (buffers->need_tag()) tag_size = e_nall;
  if (molecular) {
    s = atom->special[0];
    ns = atom->nspecial[0];
    special_size = aend;
  } else {
    s = &buffers->_special_holder;
    ns = &buffers->_nspecial_holder;
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

  const int xperiodic = domain->xperiodic;
  const int yperiodic = domain->yperiodic;
  const int zperiodic = domain->zperiodic;
  const flt_t xprd_half = domain->xprd_half;
  const flt_t yprd_half = domain->yprd_half;
  const flt_t zprd_half = domain->zprd_half;

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
    in(xperiodic, yperiodic, zperiodic, xprd_half, yprd_half, zprd_half) \
    out(overflow:length(5) alloc_if(0) free_if(0)) \
    out(timer_compute:length(1) alloc_if(0) free_if(0)) \
    signal(tag)
  #endif
  {
    #if defined(__MIC__) && defined(_LMP_INTEL_OFFLOAD)
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

      const int list_size = (ito + tid * 2 + 2) * maxnbors;
      int ct = (ifrom + tid * 2) * maxnbors;
      int *neighptr = firstneigh + ct;
      const int obound = maxnbors * 3;

      for (int i = ifrom; i < ito; i++) {
        const flt_t xtmp = x[i].x;
        const flt_t ytmp = x[i].y;
        const flt_t ztmp = x[i].z;
        const int itype = x[i].w;
        const int ioffset = ntypes * itype;

        // loop over all atoms in bins in stencil
        // pairs for atoms j "below" i are excluded
        // below = lower z or (equal z and lower y) or (equal zy and lower x)
        //         (equal zyx and j <= i)
        // latter excludes self-self interaction but allows superposed atoms

        const int ibin = atombin[i];

	int raw_count = maxnbors;
        for (int k = 0; k < nstencil; k++) {
          for (int j = binhead[ibin + stencil[k]]; j >= 0; j = bins[j]) {
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

            #ifndef _LMP_INTEL_OFFLOAD
            if (exclude) {
	      const int jtype = x[j].w;
	      if (exclusion(i,j,itype,jtype,mask,molecule)) continue;
	    }
	    #endif

	    neighptr[raw_count++] = j;
	  }
	}
	if (raw_count > obound)
	  *overflow = 1;

        #if defined(LMP_SIMD_COMPILER)
	#ifdef _LMP_INTEL_OFFLOAD
	int vlmin = lmin, vlmax = lmax, vgmin = gmin, vgmax = gmax;
	#pragma vector aligned
        #pragma simd reduction(max:vlmax,vgmax) reduction(min:vlmin, vgmin)
	#else
	#pragma vector aligned
        #pragma simd
	#endif
	#endif
	for (int u = maxnbors; u < raw_count; u++) {
          int j = neighptr[u];
          const flt_t delx = xtmp - x[j].x;
          const flt_t dely = ytmp - x[j].y;
          const flt_t delz = ztmp - x[j].z;
	  const int jtype = x[j].w;
          const flt_t rsq = delx * delx + dely * dely + delz * delz;
          if (rsq > cutneighsq[ioffset + jtype]) 
	    neighptr[u] = e_nall;
	  else {
            if (need_ic) {
	      int no_special;
	      ominimum_image_check(no_special, delx, dely, delz);
	      if (no_special)
	        neighptr[u] = -j - 1;
	    }

            #ifdef _LMP_INTEL_OFFLOAD
            if (j < nlocal) {
              if (j < vlmin) vlmin = j;
              if (j > vlmax) vlmax = j;
            } else {
              if (j < vgmin) vgmin = j;
              if (j > vgmax) vgmax = j;
            }
            #endif
          }
        }

        int n = 0, n2 = maxnbors;
	for (int u = maxnbors; u < raw_count; u++) {
	  const int j = neighptr[u];
	  int pj = j;
	  if (pj < e_nall) {
	    if (need_ic)
	      if (pj < 0) pj = -pj - 1;

	    if (pj < nlocal)
	      neighptr[n++] = j;
	    else
	      neighptr[n2++] = j;
	  }
	}
	int ns = n;
	for (int u = maxnbors; u < n2; u++)
	  neighptr[n++] = neighptr[u];

        ilist[i] = i;
        cnumneigh[i] = ct;
	ns += n2 - maxnbors;
        while( (ns % pad_width) != 0 ) neighptr[ns++] = e_nall;
        numneigh[i] = ns;

	ct += ns;
	const int alignb = (INTEL_DATA_ALIGN / sizeof(int));
	const int edge = (ct % alignb);
	if (edge) ct += alignb - edge;
	neighptr = firstneigh + ct;
	if (ct + obound > list_size) {
	  if (i < ito - 1) {
	    *overflow = 1;
	    ct = (ifrom + tid * 2) * maxnbors;
	  }
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
          #if defined(LMP_SIMD_COMPILER)
	  #pragma vector aligned
          #pragma simd
	  #endif
          for (int jj = 0; jj < jnum; jj++) {
            const int j = jlist[jj];
	    if (need_ic && j < 0) {
	      which = 0;
	      jlist[jj] = -j - 1;
	    } else
              ofind_special(which, special, nspecial, i, tag[j], special_flag);
            #ifdef _LMP_INTEL_OFFLOAD
	    if (j >= nlocal) {
	      if (j == e_nall)
		jlist[jj] = nall_offset;
	      else if (which) 
		jlist[jj] = (j-ghost_offset) ^ (which << SBBITS);
	      else jlist[jj]-=ghost_offset;
            } else
            #endif
	      if (which) jlist[jj] = j ^ (which << SBBITS);
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
    #if defined(__MIC__) && defined(_LMP_INTEL_OFFLOAD)
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
   binned neighbor list construction for all neighbors
   every neighbor pair appears in list of both atoms i and j
------------------------------------------------------------------------- */

void Neighbor::full_bin_intel(NeighList *list)
{
  const int nlocal = (includegroup) ? atom->nfirst : atom->nlocal;
  list->inum = nlocal;
  list->gnum = 0;

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

  int need_ic = 0;
  if (atom->molecular)
    dminimum_image_check(need_ic, cutneighmax, cutneighmax, cutneighmax);

  if (need_ic) {
    if (fix->precision() == FixIntel::PREC_MODE_MIXED) {
      #ifdef _LMP_INTEL_OFFLOAD
      if (offload_noghost) {
        fbi<float,double,1,1>(1, list, fix->get_mixed_buffers(),
                              0, off_end, fix);
        fbi<float,double,1,1>(0, list, fix->get_mixed_buffers(),
                              host_start, nlocal, fix, off_end);
      } else
      #endif
      {
        fbi<float,double,0,1>(1, list, fix->get_mixed_buffers(),
                              0, off_end, fix);
        fbi<float,double,0,1>(0, list, fix->get_mixed_buffers(),
                              host_start, nlocal, fix);
      }
    } else if (fix->precision() == FixIntel::PREC_MODE_DOUBLE) {
      #ifdef _LMP_INTEL_OFFLOAD
      if (offload_noghost) {
        fbi<double,double,1,1>(1, list, fix->get_double_buffers(),
                               0, off_end, fix);
        fbi<double,double,1,1>(0, list, fix->get_double_buffers(),
                               host_start, nlocal, fix, off_end);
      } else
      #endif
      {
        fbi<double,double,0,1>(1, list, fix->get_double_buffers(),
                               0, off_end, fix);
        fbi<double,double,0,1>(0, list, fix->get_double_buffers(),
                               host_start, nlocal, fix);
      }
    } else {
      #ifdef _LMP_INTEL_OFFLOAD
      if (offload_noghost) {
        fbi<float,float,1,1>(1, list, fix->get_single_buffers(), 0, off_end,
			     fix);
        fbi<float,float,1,1>(0, list, fix->get_single_buffers(),
                             host_start, nlocal, fix, off_end);
      } else
      #endif
      {
        fbi<float,float,0,1>(1, list, fix->get_single_buffers(), 0, off_end,
			     fix);
        fbi<float,float,0,1>(0, list, fix->get_single_buffers(),
                             host_start, nlocal, fix);
      }
    }
  } else {
    if (fix->precision() == FixIntel::PREC_MODE_MIXED) {
      #ifdef _LMP_INTEL_OFFLOAD
      if (offload_noghost) {
        fbi<float,double,1,0>(1, list, fix->get_mixed_buffers(),
                              0, off_end, fix);
        fbi<float,double,1,0>(0, list, fix->get_mixed_buffers(),
                              host_start, nlocal, fix, off_end);
      } else
      #endif
      {
        fbi<float,double,0,0>(1, list, fix->get_mixed_buffers(),
                              0, off_end, fix);
        fbi<float,double,0,0>(0, list, fix->get_mixed_buffers(),
                              host_start, nlocal, fix);
      }
    } else if (fix->precision() == FixIntel::PREC_MODE_DOUBLE) {
      #ifdef _LMP_INTEL_OFFLOAD
      if (offload_noghost) {
        fbi<double,double,1,0>(1, list, fix->get_double_buffers(),
                               0, off_end, fix);
        fbi<double,double,1,0>(0, list, fix->get_double_buffers(),
                               host_start, nlocal, fix, off_end);
      } else
      #endif
      {
        fbi<double,double,0,0>(1, list, fix->get_double_buffers(),
                               0, off_end, fix);
        fbi<double,double,0,0>(0, list, fix->get_double_buffers(),
                               host_start, nlocal, fix);
      }
    } else {
      #ifdef _LMP_INTEL_OFFLOAD
      if (offload_noghost) {
        fbi<float,float,1,0>(1, list, fix->get_single_buffers(), 0, off_end,
			     fix);
        fbi<float,float,1,0>(0, list, fix->get_single_buffers(),
                             host_start, nlocal, fix, off_end);
      } else
      #endif
      {
        fbi<float,float,0,0>(1, list, fix->get_single_buffers(), 0, off_end,
			     fix);
        fbi<float,float,0,0>(0, list, fix->get_single_buffers(),
                             host_start, nlocal, fix);
      }
    }
  }
}

template <class flt_t, class acc_t, int offload_noghost, int need_ic>
void Neighbor::fbi(const int offload, NeighList *list, void *buffers_in,
                    const int astart, const int aend, void *fix_in,
                    const int offload_end) {
  IntelBuffers<flt_t,acc_t> *buffers = (IntelBuffers<flt_t,acc_t> *)buffers_in;
  FixIntel *fix = (FixIntel *)fix_in;
  const int nall = atom->nlocal + atom->nghost;
  int pad = 1;

  const int pack_width = fix->nbor_pack_width();

  if (offload) {
    fix->start_watch(TIME_PACK);
    buffers->grow(nall, atom->nlocal, comm->nthreads, aend);
    buffers->grow_nbor(list, atom->nlocal, comm->nthreads, aend, pack_width);

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
  } else {
    fix->start_watch(TIME_HOST_NEIGHBOR);
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
  int tag_size = 0, special_size;
  if (buffers->need_tag()) tag_size = e_nall;
  if (molecular) {
    s = atom->special[0];
    ns = atom->nspecial[0];
    special_size = aend;
  } else {
    s = &buffers->_special_holder;
    ns = &buffers->_nspecial_holder;
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

  const int xperiodic = domain->xperiodic;
  const int yperiodic = domain->yperiodic;
  const int zperiodic = domain->zperiodic;
  const flt_t xprd_half = domain->xprd_half;
  const flt_t yprd_half = domain->yprd_half;
  const flt_t zprd_half = domain->zprd_half;

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
    in(maxnbors,nthreads,maxspecial,nstencil,e_nall,offload,pack_width)	\
    in(offload_end,separate_buffers,astart, aend, nlocal, molecular, ntypes) \
    in(xperiodic, yperiodic, zperiodic, xprd_half, yprd_half, zprd_half) \
    out(overflow:length(5) alloc_if(0) free_if(0)) \
    out(timer_compute:length(1) alloc_if(0) free_if(0)) \
    signal(tag)
  #endif
  {
    #if defined(__MIC__) && defined(_LMP_INTEL_OFFLOAD)
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

      IP_PRE_omp_range_id_vec(ifrom, ito, tid, num, nthreads, pack_width);
      ifrom += astart;
      ito += astart;
      int e_ito = ito;
      if (ito == num) {
	int imod = ito % pack_width;
	if (imod) e_ito += pack_width - imod;
      }
      const int list_size = (e_ito + tid * 2 + 2) * maxnbors;
      int which;
      int pack_offset = maxnbors * pack_width;
      int ct = (ifrom + tid * 2) * maxnbors;
      int *neighptr = firstneigh + ct;
      const int obound = pack_offset + maxnbors * 2;

      int max_chunk = 0;
      int lane = 0;
      for (int i = ifrom; i < ito; i++) {
        const flt_t xtmp = x[i].x;
        const flt_t ytmp = x[i].y;
        const flt_t ztmp = x[i].z;
        const int itype = x[i].w;
        const tagint itag = tag[i];
        const int ioffset = ntypes * itype;

        const int ibin = atombin[i];
	int raw_count = pack_offset;

        // loop over all atoms in surrounding bins in stencil including self
        // skip i = j
        for (int k = 0; k < nstencil; k++) {
          for (int j = binhead[ibin + stencil[k]]; j >= 0; j = bins[j]) {
            if (i == j) continue;

            if (offload_noghost) {
              if (j < nlocal) {
                if (i < offload_end) continue;
              } else if (offload) continue;
            }

            #ifndef _LMP_INTEL_OFFLOAD
            if (exclude) {
              const int jtype = x[j].w;
              if (exclusion(i,j,itype,jtype,mask,molecule)) continue;
            }
            #endif
	    
	    neighptr[raw_count++] = j;
          }
        }

	if (raw_count > obound) *overflow = 1;

        #if defined(LMP_SIMD_COMPILER)
	#ifdef _LMP_INTEL_OFFLOAD
	int vlmin = lmin, vlmax = lmax, vgmin = gmin, vgmax = gmax;
	#pragma vector aligned
        #pragma simd reduction(max:vlmax,vgmax) reduction(min:vlmin, vgmin)
	#else
	#pragma vector aligned
        #pragma simd
	#endif
	#endif
	for (int u = pack_offset; u < raw_count; u++) {
          int j = neighptr[u];
          const flt_t delx = xtmp - x[j].x;
          const flt_t dely = ytmp - x[j].y;
          const flt_t delz = ztmp - x[j].z;
          const int jtype = x[j].w;
          const flt_t rsq = delx * delx + dely * dely + delz * delz;
          if (rsq > cutneighsq[ioffset + jtype]) 
	    neighptr[u] = e_nall;
	  else {
            if (need_ic) {
              int no_special;
              ominimum_image_check(no_special, delx, dely, delz);
              if (no_special)
                neighptr[u] = -j - 1;
            }
            #ifdef _LMP_INTEL_OFFLOAD
            if (j < nlocal) {
              if (j < vlmin) vlmin = j;
              if (j > vlmax) vlmax = j;
            } else {
              if (j < vgmin) vgmin = j;
              if (j > vgmax) vgmax = j;
            }
            #endif
	  }
	}
        #ifdef _LMP_INTEL_OFFLOAD
	lmin = MIN(lmin,vlmin);
	gmin = MIN(gmin,vgmin);
	lmax = MAX(lmax,vlmax);
	gmax = MAX(gmax,vgmax);
        #endif

        int n = lane, n2 = pack_offset;
	for (int u = pack_offset; u < raw_count; u++) {
	  const int j = neighptr[u];
	  int pj = j;
	  if (pj < e_nall) {
	    if (need_ic)
	      if (pj < 0) pj = -pj - 1;
	
	    const int jtag = tag[pj];
	    int flist = 0;
	    if (itag > jtag) {
	      if ((itag+jtag) % 2 == 0) flist = 1;
	    } else if (itag < jtag) {
	      if ((itag+jtag) % 2 == 1) flist = 1;
	    } else {
              if (x[pj].z < ztmp) flist = 1;
	      else if (x[pj].z == ztmp && x[pj].y < ytmp) flist = 1;
	      else if (x[pj].z == ztmp && x[pj].y == ytmp && x[pj].x < xtmp) 
	      flist = 1;
	    }
	    if (flist) {
	      neighptr[n2++] = j;
	    } else {
	      neighptr[n] = j;
	      n += pack_width;
	    }
          }
        }
	int ns = (n - lane) / pack_width;
	atombin[i] = ns;
	for (int u = pack_offset; u < n2; u++) {
	  neighptr[n] = neighptr[u];
	  n += pack_width;
	}

        ilist[i] = i;
        cnumneigh[i] = ct + lane;
	ns += n2 - pack_offset;
        numneigh[i] = ns;

	if (ns > max_chunk) max_chunk = ns;
	lane++;
	if (lane == pack_width) {
	  ct += max_chunk * pack_width;
	  const int alignb = (INTEL_DATA_ALIGN / sizeof(int));
	  const int edge = (ct % alignb);
	  if (edge) ct += alignb - edge;
	  neighptr = firstneigh + ct;
	  max_chunk = 0;
	  pack_offset = maxnbors * pack_width;
	  lane = 0;
	  if (ct + obound > list_size) {
  	    if (i < ito - 1) {
	      *overflow = 1;
	      ct = (ifrom + tid * 2) * maxnbors;
	    }
          }
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

	  const int trip = jnum * pack_width;
          for (int jj = 0; jj < trip; jj+=pack_width) {
            const int j = jlist[jj];
	    if (need_ic && j < 0) {
	      which = 0;
	      jlist[jj] = -j - 1;
	    } else
              ofind_special(which, special, nspecial, i, tag[j],
                            special_flag);
            #ifdef _LMP_INTEL_OFFLOAD
            if (j >= nlocal) {
              if (j == e_nall)
                jlist[jj] = nall_offset;
              else if (which)
                jlist[jj] = (j-ghost_offset) ^ (which << SBBITS);
              else jlist[jj]-=ghost_offset;
            } else
            #endif
            if (which) jlist[jj] = j ^ (which << SBBITS);
          }
        }
      }
      #ifdef _LMP_INTEL_OFFLOAD
      else if (separate_buffers) {
        for (int i = ifrom; i < ito; ++i) {
          int * _noalias jlist = firstneigh + cnumneigh[i];
          const int jnum = numneigh[i];
          int jj = 0;
          for (jj = 0; jj < jnum; jj++) {
            if (jlist[jj] >= nlocal) {
	      if (jlist[jj] == e_nall) jlist[jj] = nall_offset;
	      else jlist[jj] -= ghost_offset;
	    }
          }
        }
      }
      #endif
    } // end omp
    #if defined(__MIC__) && defined(_LMP_INTEL_OFFLOAD)
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

