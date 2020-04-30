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

#include "omp_compat.h"
#include "comm.h"
#include "domain.h"
#include "timer.h"
#include "modify.h"
#include "npair_intel.h"
#include "nstencil.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NPairIntel::NPairIntel(LAMMPS *lmp) : NPair(lmp) {
  int ifix = modify->find_fix("package_intel");
  if (ifix < 0)
    error->all(FLERR,
               "The 'package intel' command is required for /intel styles");
  _fix = static_cast<FixIntel *>(modify->fix[ifix]);
  #ifdef _LMP_INTEL_OFFLOAD
  _cop = _fix->coprocessor_number();
  _off_map_stencil = 0;
  #endif
}

/* ---------------------------------------------------------------------- */

NPairIntel::~NPairIntel() {
  #ifdef _LMP_INTEL_OFFLOAD
  if (_off_map_stencil) {
    const int * stencil = this->stencil;
    #pragma offload_transfer target(mic:_cop)   \
      nocopy(stencil:alloc_if(0) free_if(1))
  }
  #endif
}

/* ---------------------------------------------------------------------- */

void NPairIntel::copy_neighbor_info()
{
  NPair::copy_neighbor_info();
  if (_fix->precision() == FixIntel::PREC_MODE_MIXED)
    copy_cutsq_info(_fix->get_mixed_buffers());
  else if (_fix->precision() == FixIntel::PREC_MODE_DOUBLE)
    copy_cutsq_info(_fix->get_double_buffers());
  else
    copy_cutsq_info(_fix->get_single_buffers());
}

/* ---------------------------------------------------------------------- */

template <class flt_t, class acc_t>
void NPairIntel::copy_cutsq_info(IntelBuffers<flt_t,acc_t> *buffers) {
  int tp1 = atom->ntypes + 1;
  int use_ghost_cut = 0;
  if (cutneighghostsq)
    use_ghost_cut = 1;
  buffers->set_ntypes(tp1, use_ghost_cut);

  flt_t **cutneighsqb = buffers->get_cutneighsq();
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = 1; j <= atom->ntypes; j++)
      cutneighsqb[i][j] = cutneighsq[i][j];

  flt_t **cutneighghostsqb;
  if (use_ghost_cut) {
    cutneighghostsqb = buffers->get_cutneighghostsq();
    for (int i = 1; i <= atom->ntypes; i++)
      for (int j = 1; j <= atom->ntypes; j++)
        cutneighghostsqb[i][j] = cutneighghostsq[i][j];
  }

  #ifdef _LMP_INTEL_OFFLOAD
  if (_cop < 0) return;
  int tp1sq = tp1 * tp1;
  flt_t * ocutneighsq = cutneighsqb[0];
  #pragma offload_transfer target(mic:_cop) in(ocutneighsq: length(tp1sq))
  if (use_ghost_cut) {
    flt_t * ocutneighghostsq = cutneighghostsqb[0];
    #pragma offload_transfer target(mic:_cop) \
      in(ocutneighghostsq: length(tp1sq))
  }
  #endif
}

/* ---------------------------------------------------------------------- */

template <class flt_t, class acc_t, int offload_noghost, int need_ic,
          int FULL, int TRI, int THREE>
void NPairIntel::bin_newton(const int offload, NeighList *list,
                            IntelBuffers<flt_t,acc_t> *buffers,
                            const int astart, const int aend,
                            const int offload_end) {

  if (aend-astart == 0) return;

  const int nall = atom->nlocal + atom->nghost;
  int nall_t = nall;

  #ifdef _LMP_INTEL_OFFLOAD
  if (offload_noghost && offload) nall_t = atom->nlocal;
  #endif

  const int pack_width = _fix->nbor_pack_width();

  const ATOM_T * _noalias const x = buffers->get_x();
  int * _noalias const intel_list = buffers->intel_list(list);
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
  int ** _noalias const firstneigh = list->firstneigh;
  int * _noalias const numneigh = list->numneigh;
  int * _noalias const cnumneigh = buffers->cnumneigh();

  const int nstencil = this->nstencil;
  const int * _noalias const stencil = this->stencil;
  const flt_t * _noalias const cutneighsq = buffers->get_cutneighsq()[0];
  const int ntypes = atom->ntypes + 1;
  const int nlocal = atom->nlocal;

  #ifndef _LMP_INTEL_OFFLOAD
  int * _noalias const mask = atom->mask;
  tagint * _noalias const molecule = atom->molecule;
  #endif

  int tnum;
  int * _noalias overflow;
  #ifdef _LMP_INTEL_OFFLOAD
  double *timer_compute;
  if (offload) {
    timer_compute = _fix->off_watch_neighbor();
    tnum = buffers->get_off_threads();
    overflow = _fix->get_off_overflow_flag();
    _fix->stop_watch(TIME_HOST_NEIGHBOR);
    _fix->start_watch(TIME_OFFLOAD_LATENCY);
  } else
  #endif
  {
    tnum = comm->nthreads;
    overflow = _fix->get_overflow_flag();
  }
  const int nthreads = tnum;
  const int maxnbors = buffers->get_max_nbors();
  int * _noalias const atombin = buffers->get_atombin();
  const int * _noalias const binpacked = buffers->get_binpacked();

  const int xperiodic = domain->xperiodic;
  const int yperiodic = domain->yperiodic;
  const int zperiodic = domain->zperiodic;
  const flt_t xprd_half = domain->xprd_half;
  const flt_t yprd_half = domain->yprd_half;
  const flt_t zprd_half = domain->zprd_half;

  flt_t * _noalias const ncachex = buffers->get_ncachex();
  flt_t * _noalias const ncachey = buffers->get_ncachey();
  flt_t * _noalias const ncachez = buffers->get_ncachez();
  int * _noalias const ncachej = buffers->get_ncachej();
  int * _noalias const ncachejtype = buffers->get_ncachejtype();
  tagint * _noalias const ncachetag = buffers->get_ncachetag();
  const int ncache_stride = buffers->ncache_stride();

  int sb = 1;
  if (special_flag[1] == 0) {
    sb = 2;
    if (special_flag[2] == 0) {
      sb = 3;
      if (special_flag[3] == 0)
        sb = 4;
    }
  }
  const int special_bound = sb;

  #ifdef _LMP_INTEL_OFFLOAD
  const int * _noalias const binhead = this->binhead;
  const int * _noalias const bins = this->bins;
  const int cop = _fix->coprocessor_number();
  const int separate_buffers = _fix->separate_buffers();
  #pragma offload target(mic:cop) if(offload) \
    in(x:length(e_nall+1) alloc_if(0) free_if(0)) \
    in(tag:length(tag_size) alloc_if(0) free_if(0)) \
    in(special:length(special_size*maxspecial) alloc_if(0) free_if(0)) \
    in(nspecial:length(special_size*3) alloc_if(0) free_if(0)) \
    in(bins,binpacked:length(nall) alloc_if(0) free_if(0)) \
    in(binhead:length(mbins+1) alloc_if(0) free_if(0)) \
    in(cutneighsq:length(0) alloc_if(0) free_if(0)) \
    in(firstneigh:length(0) alloc_if(0) free_if(0)) \
    in(intel_list:length(0) alloc_if(0) free_if(0)) \
    in(cnumneigh:length(0) alloc_if(0) free_if(0)) \
    out(numneigh:length(0) alloc_if(0) free_if(0)) \
    in(ilist:length(0) alloc_if(0) free_if(0)) \
    in(atombin:length(aend) alloc_if(0) free_if(0)) \
    in(stencil:length(nstencil) alloc_if(0) free_if(0)) \
    in(ncachex,ncachey,ncachez,ncachej:length(0) alloc_if(0) free_if(0)) \
    in(ncachejtype,ncachetag:length(0) alloc_if(0) free_if(0)) \
    in(ncache_stride,maxnbors,nthreads,maxspecial,nstencil,e_nall,offload) \
    in(offload_end,separate_buffers,astart,aend,nlocal,molecular) \
    in(ntypes,xperiodic,yperiodic,zperiodic,xprd_half,yprd_half,zprd_half) \
    in(pack_width,special_bound)                                        \
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

    int nstencilp = 0;
    int binstart[INTEL_MAX_STENCIL], binend[INTEL_MAX_STENCIL];
    for (int k = 0; k < nstencil; k++) {
      binstart[nstencilp] = stencil[k];
      int end = stencil[k] + 1;
      for (int kk = k + 1; kk < nstencil; kk++) {
        if (stencil[kk-1]+1 == stencil[kk]) {
          end++;
          k++;
        } else break;
      }
      binend[nstencilp] = end;
      nstencilp++;
    }

    #if defined(_OPENMP)
    #pragma omp parallel LMP_DEFAULT_NONE \
      shared(overflow, nstencilp, binstart, binend)
    #endif
    {
      #ifdef _LMP_INTEL_OFFLOAD
      int lmin = e_nall, lmax = -1, gmin = e_nall, gmax = -1;
      #endif

      const int num = aend - astart;
      int tid, ifrom, ito;

      if (THREE) {
        IP_PRE_omp_range_id_vec(ifrom, ito, tid, num, nthreads, pack_width);
      } else {
        IP_PRE_omp_range_id(ifrom, ito, tid, num, nthreads);
      }
      ifrom += astart;
      ito += astart;
      int e_ito = ito;
      #ifdef LMP_INTEL_3BODY_FAST
      if (THREE && ito == num) {
        int imod = ito & (pack_width - 1);
        if (imod) e_ito += pack_width - imod;
      }
      #endif
      const int list_size = (e_ito + tid * 2 + 2) * maxnbors;

      #ifdef LMP_INTEL_3BODY_FAST
      const int pack_offset = maxnbors * pack_width;
      const int obound = pack_offset + maxnbors * 2;
      #else
      const int pack_offset = 0;
      const int obound = maxnbors * 3;
      #endif
      int ct = (ifrom + tid * 2) * maxnbors;
      int * _noalias neighptr = intel_list + ct;
      int * _noalias neighptr2;
      if (THREE) neighptr2 = neighptr;

      const int toffs = tid * ncache_stride;
      flt_t * _noalias const tx = ncachex + toffs;
      flt_t * _noalias const ty = ncachey + toffs;
      flt_t * _noalias const tz = ncachez + toffs;
      int * _noalias const tj = ncachej + toffs;
      int * _noalias const tjtype = ncachejtype + toffs;
      tagint * _noalias const ttag = ncachetag + toffs;

      flt_t * _noalias itx;
      flt_t * _noalias ity;
      flt_t * _noalias itz;
      int * _noalias itj;
      int * _noalias itjtype;

      // loop over all atoms in other bins in stencil, store every pair
      int istart, icount, ncount, oldbin = -9999999;
      #ifdef LMP_INTEL_3BODY_FAST
      int lane, max_chunk;
      if (THREE) {
        lane = 0;
        max_chunk = 0;
      }
      #endif
      for (int i = ifrom; i < ito; i++) {
        const flt_t xtmp = x[i].x;
        const flt_t ytmp = x[i].y;
        const flt_t ztmp = x[i].z;
        const int itype = x[i].w;
        tagint itag;
        if (THREE) itag = tag[i];
        const int ioffset = ntypes * itype;

        const int ibin = atombin[i];
        if (ibin != oldbin) {
          oldbin = ibin;
          ncount = 0;
          for (int k = 0; k < nstencilp; k++) {
            const int bstart = binhead[ibin + binstart[k]];
            const int bend = binhead[ibin + binend[k]];
            #if defined(LMP_SIMD_COMPILER)
            #pragma simd
            #endif
            for (int jj = bstart; jj < bend; jj++)
              tj[ncount++] = binpacked[jj];
          }
          #if defined(LMP_SIMD_COMPILER)
          #pragma vector aligned
          #pragma simd
          #endif
          for (int u = 0; u < ncount; u++) {
            const int j = tj[u];
            tx[u] = x[j].x;
            ty[u] = x[j].y;
            tz[u] = x[j].z;
            tjtype[u] = x[j].w;
            if (THREE) ttag[u] = tag[j];
          }

          if (FULL == 0 && TRI != 1) {
            icount = 0;
            istart = ncount;
            IP_PRE_edge_align(istart, sizeof(int));
            itx = tx + istart;
            ity = ty + istart;
            itz = tz + istart;
            itj = tj + istart;
            itjtype = tjtype + istart;

            const int bstart = binhead[ibin];
            const int bend = binhead[ibin + 1];
            #if defined(LMP_SIMD_COMPILER)
            #pragma simd
            #endif
            for (int jj = bstart; jj < bend; jj++) {
              const int j = binpacked[jj];
              itj[icount] = j;
              itx[icount] = x[j].x;
              ity[icount] = x[j].y;
              itz[icount] = x[j].z;
              itjtype[icount] = x[j].w;
              icount++;
            }
            if (icount + istart > obound) *overflow = 1;
          } else
            if (ncount > obound) *overflow = 1;
        }

        // ---------------------- Loop over i bin

        int n = 0;
        if (FULL == 0 && TRI != 1) {
          #if defined(LMP_SIMD_COMPILER)
          #pragma vector aligned
          #pragma ivdep
          #endif
          for (int u = 0; u < icount; u++) {
            int addme = 1;
            int j = itj[u];

            // Cutoff Check
            const flt_t delx = xtmp - itx[u];
            const flt_t dely = ytmp - ity[u];
            const flt_t delz = ztmp - itz[u];
            const int jtype = itjtype[u];
            const flt_t rsq = delx * delx + dely * dely + delz * delz;
            if (rsq > cutneighsq[ioffset + jtype]) addme = 0;

            // i bin (half) check and offload ghost check
            if (j < nlocal) {
              const int ijmod = (i + j) & 1;
              if (i > j) {
                if (ijmod == 0) addme = 0;
              } else if (i < j) {
                if (ijmod == 1) addme = 0;
              } else
                addme = 0;
              #ifdef _LMP_INTEL_OFFLOAD
              if (offload_noghost && i < offload_end) addme = 0;
              #endif
            } else {
              #ifdef _LMP_INTEL_OFFLOAD
              if (offload_noghost && offload) addme = 0;
              #endif
              if (itz[u] < ztmp) addme = 0;
              if (itz[u] == ztmp) {
                if (ity[u] < ytmp) addme = 0;
                if (ity[u] == ytmp && itx[u] < xtmp) addme = 0;
              }
            }

            if (need_ic) {
              int no_special;
              ominimum_image_check(no_special, delx, dely, delz);
              if (no_special)
                j = -j - 1;
            }

            if (addme)
              neighptr[n++] = j;
          }
        } // if FULL==0

        // ---------------------- Loop over other bins

        int n2;
        if (THREE) {
          #ifdef LMP_INTEL_3BODY_FAST
          n = pack_offset;
          #endif
          n2 = pack_offset + maxnbors;
        }
        #if defined(LMP_SIMD_COMPILER)
        #pragma vector aligned
        #pragma ivdep
        #endif
        for (int u = 0; u < ncount; u++) {
          int addme = 1;
          int j = tj[u];

          if (FULL)
            if (i == j) addme = 0;

          // Cutoff Check
          const flt_t delx = xtmp - tx[u];
          const flt_t dely = ytmp - ty[u];
          const flt_t delz = ztmp - tz[u];
          const int jtype = tjtype[u];
          const flt_t rsq = delx * delx + dely * dely + delz * delz;
          if (rsq > cutneighsq[ioffset + jtype]) addme = 0;

          // Triclinic
          if (TRI) {
            if (tz[u] < ztmp) addme = 0;
            if (tz[u] == ztmp) {
              if (ty[u] < ytmp) addme = 0;
              if (ty[u] == ytmp) {
                if (tx[u] < xtmp) addme = 0;
                if (tx[u] == xtmp && j <= i) addme = 0;
              }
            }
          }

          // offload ghost check
          #ifdef _LMP_INTEL_OFFLOAD
          if (offload_noghost) {
            if (j < nlocal) {
              if (i < offload_end) addme = 0;
            } else if (offload) addme = 0;
          }
          #endif

          if (need_ic) {
            int no_special;
            ominimum_image_check(no_special, delx, dely, delz);
            if (no_special)
              j = -j - 1;
          }

          if (THREE) {
            const tagint jtag = ttag[u];
            int flist = 0;
            if (itag > jtag) {
              if (((itag+jtag) & 1) == 0) flist = 1;
            } else if (itag < jtag) {
              if (((itag+jtag) & 1) == 1) flist = 1;
            } else {
              if (tz[u] < ztmp) flist = 1;
              else if (tz[u] == ztmp && ty[u] < ytmp) flist = 1;
              else if (tz[u] == ztmp && ty[u] == ytmp && tx[u] < xtmp)
                flist = 1;
            }
            if (addme) {
              if (flist)
                neighptr2[n2++] = j;
              else
                neighptr[n++] = j;
            }
          } else {
            if (addme)
              neighptr[n++] = j;
          }
        } // for u

        if (molecular) {
          if (!THREE) neighptr2 = neighptr;
          int alln = n;

          n = pack_offset;
          #if defined(LMP_SIMD_COMPILER)
          #pragma vector aligned
          #ifdef LMP_INTEL_NBOR_COMPAT
          #pragma ivdep
          #else
          #pragma simd
          #endif
          #endif
          for (int u = n; u < alln; u++) {
            int which;
            int addme = 1;
            int j = neighptr[u];
            if (need_ic && j < 0) {
              which = 0;
              j = -j - 1;
            } else
              ofind_special(which, special, nspecial, i, tag[j]);

            if (which) {
              j = j ^ (which << SBBITS);
              if (which < special_bound) addme = 0;
            }
            #ifdef LMP_INTEL_NBOR_COMPAT
            if (addme) neighptr2[n++] = j;
            #else
            neighptr2[n++]=j;
            #endif
          }

          if (THREE) {
            alln = n2;
            n2 = pack_offset + maxnbors;

            #if defined(LMP_SIMD_COMPILER)
            #pragma vector aligned
            #ifdef LMP_INTEL_NBOR_COMPAT
            #pragma ivdep
            #else
            #pragma simd
            #endif
            #endif
            for (int u = n2; u < alln; u++) {
              int which;
              int addme = 1;
              int j = neighptr[u];
              if (need_ic && j < 0) {
                which = 0;
                j = -j - 1;
              } else
                ofind_special(which, special, nspecial, i, tag[j]);
              if (which) {
                j = j ^ (which << SBBITS);
                if (which < special_bound) addme = 0;
              }
              #ifdef LMP_INTEL_NBOR_COMPAT
              if (addme) neighptr2[n2++] = j;
              #else
              neighptr2[n2++]=j;
              #endif
            }
          }
        }

        #ifndef _LMP_INTEL_OFFLOAD
        if (exclude) {
          neighptr2 = neighptr;
          int alln = n;
          n = pack_offset;

          #if defined(LMP_SIMD_COMPILER)
          #pragma vector aligned
          #pragma ivdep
          #endif
          for (int u = n; u < alln; u++) {
            int addme = 1;
            const int js = neighptr[u];
            const int j = js & NEIGHMASK;
            const int jtype = x[j].w;
            if (exclusion(i,j,itype,jtype,mask,molecule)) addme = 0;
            if (addme) neighptr2[n++] = js;
          }
          if (THREE) {
            alln = n2;
            n2 = pack_offset + maxnbors;
            #if defined(LMP_SIMD_COMPILER)
            #pragma vector aligned
            #pragma ivdep
            #endif
            for (int u = n2; u < alln; u++) {
              int addme = 1;
              const int js = neighptr[u];
              const int j = js & NEIGHMASK;
              const int jtype = x[j].w;
              if (exclusion(i,j,itype,jtype,mask,molecule)) addme = 0;
              if (addme) neighptr2[n2++] = js;
            }
          }
        }
        #endif

        int ns;
        if (THREE) {
          ns = n - pack_offset;
          atombin[i] = ns;
          ns += n2 - pack_offset - maxnbors;

          #ifdef LMP_INTEL_3BODY_FAST
          n = lane;
          for (int u = pack_offset; u < alln; u++) {
            neighptr[n] = neighptr2[u];
            n += pack_width;
          }
          #endif

          for (int u = pack_offset + maxnbors; u < n2; u++) {
            #ifdef LMP_INTEL_3BODY_FAST
            neighptr[n] = neighptr2[u];
            n += pack_width;
            #else
            neighptr[n++] = neighptr2[u];
            #endif
          }
          if (ns > maxnbors) *overflow = 1;
        } else
          if (n > maxnbors) *overflow = 1;

        ilist[i] = i;
        firstneigh[i] = intel_list + ct;
        if (THREE) {
          numneigh[i] = ns;
          cnumneigh[i] = ct;
          #ifdef LMP_INTEL_3BODY_FAST
          cnumneigh[i] += lane;
          #else
          // Pad anyways just in case we have hybrid with 2-body and newton off
          int pad_end = ns;
          IP_PRE_neighbor_pad(pad_end, offload);
          #if defined(LMP_SIMD_COMPILER)
          #pragma vector aligned
          #pragma loop_count min=1, max=INTEL_COMPILE_WIDTH-1, \
                  avg=INTEL_COMPILE_WIDTH/2
          #endif
          for ( ; ns < pad_end; ns++)
            neighptr[n++] = e_nall;
          #endif
        } else {
          numneigh[i] = n;
          int pad_end = n;
          IP_PRE_neighbor_pad(pad_end, offload);
          #if defined(LMP_SIMD_COMPILER)
          #pragma vector aligned
          #pragma loop_count min=1, max=INTEL_COMPILE_WIDTH-1, \
                  avg=INTEL_COMPILE_WIDTH/2
          #endif
          for ( ; n < pad_end; n++)
            neighptr[n] = e_nall;
        }

        #ifdef LMP_INTEL_3BODY_FAST
        if (THREE) {
          if (ns > max_chunk) max_chunk = ns;
          lane++;
          if (lane == pack_width) {
            ct += max_chunk * pack_width;
            IP_PRE_edge_align(ct, sizeof(int));
            neighptr = intel_list + ct;
            neighptr2 = neighptr;
            max_chunk = 0;
            lane = 0;
            if (ct + obound > list_size) {
              if (i < ito - 1) {
                *overflow = 1;
                ct = (ifrom + tid * 2) * maxnbors;
              }
            }
          }
        } else
        #endif
        {
          ct += n;
          //IP_PRE_edge_align(ct, sizeof(int));
          neighptr = intel_list + ct;
          if (THREE) neighptr2 = neighptr;
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
      int vlmin = lmin, vlmax = lmax, vgmin = gmin, vgmax = gmax;
      int ghost_offset = 0, nall_offset = e_nall;
      if (separate_buffers) {
        for (int i = ifrom; i < ito; ++i) {
          int * _noalias jlist = firstneigh[i];
          int jnum = numneigh[i];
          if (!THREE) IP_PRE_neighbor_pad(jnum, offload);
          #if __INTEL_COMPILER+0 > 1499
          #pragma vector aligned
          #pragma simd reduction(max:vlmax,vgmax) reduction(min:vlmin, vgmin)
          #endif
          for (int jj = 0; jj < jnum; jj++) {
            const int j = jlist[jj] & NEIGHMASK;
            if (j < nlocal) {
              if (j < vlmin) vlmin = j;
              if (j > vlmax) vlmax = j;
            } else {
              if (j < vgmin) vgmin = j;
              if (j > vgmax) vgmax = j;
            }
          }
        }
        lmin = MIN(lmin,vlmin);
        gmin = MIN(gmin,vgmin);
        lmax = MAX(lmax,vlmax);
        gmax = MAX(gmax,vgmax);

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

        int nghost = overflow[LMP_GHOST_MAX] + 1 - overflow[LMP_GHOST_MIN];
        if (nghost < 0) nghost = 0;
        if (offload) {
          ghost_offset = overflow[LMP_GHOST_MIN] - overflow[LMP_LOCAL_MAX] - 1;
          nall_offset = overflow[LMP_LOCAL_MAX] + 1 + nghost;
        } else {
          ghost_offset = overflow[LMP_GHOST_MIN] - nlocal;
          nall_offset = nlocal + nghost;
        }

        for (int i = ifrom; i < ito; ++i) {
          int * _noalias jlist = firstneigh[i];
          int jnum = numneigh[i];
          if (!THREE) IP_PRE_neighbor_pad(jnum, offload);
          int jj = 0;
          #pragma vector aligned
          #pragma simd
          for (jj = 0; jj < jnum; jj++) {
            const int which = jlist[jj] >> SBBITS & 3;
            const int j = jlist[jj] & NEIGHMASK;
            if (j >= nlocal) {
              if (j == e_nall) jlist[jj] = nall_offset;
              else jlist[jj] = (j - ghost_offset) ^ (which << SBBITS);
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

  #ifdef _LMP_INTEL_OFFLOAD
  if (offload) {
    _fix->stop_watch(TIME_OFFLOAD_LATENCY);
    _fix->start_watch(TIME_HOST_NEIGHBOR);
    firstneigh[0] = intel_list;
    for (int n = 0; n < aend; n++) {
      ilist[n] = n;
      numneigh[n] = 0;
    }
  } else {
    if (separate_buffers) {
      _fix->start_watch(TIME_PACK);
      _fix->set_neighbor_host_sizes();
      buffers->pack_sep_from_single(_fix->host_min_local(),
                                    _fix->host_used_local(),
                                    _fix->host_min_ghost(),
                                    _fix->host_used_ghost());
      _fix->stop_watch(TIME_PACK);
    }
  }
  #endif
}

/* ---------------------------------------------------------------------- */

#ifdef _LMP_INTEL_OFFLOAD
void NPairIntel::grow_stencil()
{
  if (_off_map_stencil != stencil) {
    if (_off_map_stencil) {
      const int * stencil = _off_map_stencil;
      #pragma offload_transfer target(mic:_cop) \
        nocopy(stencil:alloc_if(0) free_if(1))
    }
    _off_map_stencil = stencil;
    const int * stencil = _off_map_stencil;
    const int maxstencil = ns->get_maxstencil();
    #pragma offload_transfer target(mic:_cop)   \
      in(stencil:length(maxstencil) alloc_if(1) free_if(0))
  }
}
#endif

/* ---------------------------------------------------------------------- */

// ---- Half, no IC

template void NPairIntel::bin_newton<float, float, 0, 0, 0, 0, 0>
  (const int, NeighList *, IntelBuffers<float,float> *, const int, const int,
   const int);
template void NPairIntel::bin_newton<float, double, 0, 0, 0, 0, 0>
  (const int, NeighList *, IntelBuffers<float,double> *, const int, const int,
   const int);
template void NPairIntel::bin_newton<double, double, 0, 0, 0, 0, 0>
  (const int, NeighList *, IntelBuffers<double,double> *, const int, const int,
   const int);

// ---- Half, IC

template void NPairIntel::bin_newton<float, float, 0, 1, 0, 0, 0>
  (const int, NeighList *, IntelBuffers<float,float> *, const int, const int,
   const int);
template void NPairIntel::bin_newton<float, double, 0, 1, 0, 0, 0>
  (const int, NeighList *, IntelBuffers<float,double> *, const int, const int,
   const int);
template void NPairIntel::bin_newton<double, double, 0, 1, 0, 0, 0>
  (const int, NeighList *, IntelBuffers<double,double> *, const int, const int,
   const int);

// ---- Tri, no IC

template void NPairIntel::bin_newton<float, float, 0, 0, 0, 1, 0>
  (const int, NeighList *, IntelBuffers<float,float> *, const int, const int,
   const int);
template void NPairIntel::bin_newton<float, double, 0, 0, 0, 1, 0>
  (const int, NeighList *, IntelBuffers<float,double> *, const int, const int,
   const int);
template void NPairIntel::bin_newton<double, double, 0, 0, 0, 1, 0>
  (const int, NeighList *, IntelBuffers<double,double> *, const int, const int,
   const int);

// ---- Tri, IC

template void NPairIntel::bin_newton<float, float, 0, 1, 0, 1, 0>
  (const int, NeighList *, IntelBuffers<float,float> *, const int, const int,
   const int);
template void NPairIntel::bin_newton<float, double, 0, 1, 0, 1, 0>
  (const int, NeighList *, IntelBuffers<float,double> *, const int, const int,
   const int);
template void NPairIntel::bin_newton<double, double, 0, 1, 0, 1, 0>
  (const int, NeighList *, IntelBuffers<double,double> *, const int, const int,
   const int);

// ---- Full, no IC

template void NPairIntel::bin_newton<float, float, 0, 0, 1, 0, 0>
  (const int, NeighList *, IntelBuffers<float,float> *, const int, const int,
   const int);
template void NPairIntel::bin_newton<float, double, 0, 0, 1, 0, 0>
  (const int, NeighList *, IntelBuffers<float,double> *, const int, const int,
   const int);
template void NPairIntel::bin_newton<double, double, 0, 0, 1, 0, 0>
  (const int, NeighList *, IntelBuffers<double,double> *, const int, const int,
   const int);

// ---- Full, IC

template void NPairIntel::bin_newton<float, float, 0, 1, 1, 0, 0>
  (const int, NeighList *, IntelBuffers<float,float> *, const int, const int,
   const int);
template void NPairIntel::bin_newton<float, double, 0, 1, 1, 0, 0>
  (const int, NeighList *, IntelBuffers<float,double> *, const int, const int,
   const int);
template void NPairIntel::bin_newton<double, double, 0, 1, 1, 0, 0>
  (const int, NeighList *, IntelBuffers<double,double> *, const int, const int,
   const int);

// ---- 3-body, no IC

template void NPairIntel::bin_newton<float, float, 0, 0, 1, 0, 1>
  (const int, NeighList *, IntelBuffers<float,float> *, const int, const int,
   const int);
template void NPairIntel::bin_newton<float, double, 0, 0, 1, 0, 1>
  (const int, NeighList *, IntelBuffers<float,double> *, const int, const int,
   const int);
template void NPairIntel::bin_newton<double, double, 0, 0, 1, 0, 1>
  (const int, NeighList *, IntelBuffers<double,double> *, const int, const int,
   const int);

// ---- 3-body, IC

template void NPairIntel::bin_newton<float, float, 0, 1, 1, 0, 1>
  (const int, NeighList *, IntelBuffers<float,float> *, const int, const int,
   const int);
template void NPairIntel::bin_newton<float, double, 0, 1, 1, 0, 1>
  (const int, NeighList *, IntelBuffers<float,double> *, const int, const int,
   const int);
template void NPairIntel::bin_newton<double, double, 0, 1, 1, 0, 1>
  (const int, NeighList *, IntelBuffers<double,double> *, const int, const int,
   const int);

#ifdef _LMP_INTEL_OFFLOAD

// ---- Half, no IC, no ghost

template void NPairIntel::bin_newton<float, float, 1, 0, 0, 0, 0>
  (const int, NeighList *, IntelBuffers<float,float> *, const int, const int,
   const int);
template void NPairIntel::bin_newton<float, double, 1, 0, 0, 0, 0>
  (const int, NeighList *, IntelBuffers<float,double> *, const int, const int,
   const int);
template void NPairIntel::bin_newton<double, double, 1, 0, 0, 0, 0>
  (const int, NeighList *, IntelBuffers<double,double> *, const int, const int,
   const int);

// ---- Half, IC, no ghost

template void NPairIntel::bin_newton<float, float, 1, 1, 0, 0, 0>
  (const int, NeighList *, IntelBuffers<float,float> *, const int, const int,
   const int);
template void NPairIntel::bin_newton<float, double, 1, 1, 0, 0, 0>
  (const int, NeighList *, IntelBuffers<float,double> *, const int, const int,
   const int);
template void NPairIntel::bin_newton<double, double, 1, 1, 0, 0, 0>
  (const int, NeighList *, IntelBuffers<double,double> *, const int, const int,
   const int);

// ---- Tri, no IC, no ghost

template void NPairIntel::bin_newton<float, float, 1, 0, 0, 1, 0>
  (const int, NeighList *, IntelBuffers<float,float> *, const int, const int,
   const int);
template void NPairIntel::bin_newton<float, double, 1, 0, 0, 1, 0>
  (const int, NeighList *, IntelBuffers<float,double> *, const int, const int,
   const int);
template void NPairIntel::bin_newton<double, double, 1, 0, 0, 1, 0>
  (const int, NeighList *, IntelBuffers<double,double> *, const int, const int,
   const int);

// ---- Tri, IC, no ghost

template void NPairIntel::bin_newton<float, float, 1, 1, 0, 1, 0>
  (const int, NeighList *, IntelBuffers<float,float> *, const int, const int,
   const int);
template void NPairIntel::bin_newton<float, double, 1, 1, 0, 1, 0>
  (const int, NeighList *, IntelBuffers<float,double> *, const int, const int,
   const int);
template void NPairIntel::bin_newton<double, double, 1, 1, 0, 1, 0>
  (const int, NeighList *, IntelBuffers<double,double> *, const int, const int,
   const int);

// ---- Full, no IC, no ghost

template void NPairIntel::bin_newton<float, float, 1, 0, 1, 0, 0>
  (const int, NeighList *, IntelBuffers<float,float> *, const int, const int,
   const int);
template void NPairIntel::bin_newton<float, double, 1, 0, 1, 0, 0>
  (const int, NeighList *, IntelBuffers<float,double> *, const int, const int,
   const int);
template void NPairIntel::bin_newton<double, double, 1, 0, 1, 0, 0>
  (const int, NeighList *, IntelBuffers<double,double> *, const int, const int,
   const int);

// ---- Full, IC, no ghost

template void NPairIntel::bin_newton<float, float, 1, 1, 1, 0, 0>
  (const int, NeighList *, IntelBuffers<float,float> *, const int, const int,
   const int);
template void NPairIntel::bin_newton<float, double, 1, 1, 1, 0, 0>
  (const int, NeighList *, IntelBuffers<float,double> *, const int, const int,
   const int);
template void NPairIntel::bin_newton<double, double, 1, 1, 1, 0, 0>
  (const int, NeighList *, IntelBuffers<double,double> *, const int, const int,
   const int);

// ---- 3-body, no IC, no ghost

template void NPairIntel::bin_newton<float, float, 1, 0, 1, 0, 1>
  (const int, NeighList *, IntelBuffers<float,float> *, const int, const int,
   const int);
template void NPairIntel::bin_newton<float, double, 1, 0, 1, 0, 1>
  (const int, NeighList *, IntelBuffers<float,double> *, const int, const int,
   const int);
template void NPairIntel::bin_newton<double, double, 1, 0, 1, 0, 1>
  (const int, NeighList *, IntelBuffers<double,double> *, const int, const int,
   const int);

// ---- 3-body, IC, no ghost

template void NPairIntel::bin_newton<float, float, 1, 1, 1, 0, 1>
  (const int, NeighList *, IntelBuffers<float,float> *, const int, const int,
   const int);
template void NPairIntel::bin_newton<float, double, 1, 1, 1, 0, 1>
  (const int, NeighList *, IntelBuffers<float,double> *, const int, const int,
   const int);
template void NPairIntel::bin_newton<double, double, 1, 1, 1, 0, 1>
  (const int, NeighList *, IntelBuffers<double,double> *, const int, const int,
   const int);

#endif
