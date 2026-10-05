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
------------------------------------------------------------------------- */

#include "npair_intel.h"

#include "comm.h"
#include "domain.h"
#include "modify.h"

#include "omp_compat.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NPairIntel::NPairIntel(LAMMPS *lmp) : NPair(lmp) {
  _fix = static_cast<FixIntel *>(modify->get_fix_by_id("package_intel"));
  if (!_fix) error->all(FLERR, "The 'package intel' command is required for /intel styles");
}

/* ---------------------------------------------------------------------- */


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

}

/* ---------------------------------------------------------------------- */

// NOTE: the "offload_noghost" template parameter no longer has any effect
// (it is a leftover from the removed Xeon Phi offload support and is always
// instantiated as 0).  It is kept as a placeholder so the many explicit
// template instantiations below do not need to be renumbered; it could be
// renamed and repurposed for a future feature.

template <class flt_t, class acc_t, int offload_noghost, int need_ic,
          int FULL, int TRI, int THREE>
void NPairIntel::bin_newton(NeighList *list,
                            IntelBuffers<flt_t,acc_t> *buffers,
                            const int astart, const int aend) {

  if (aend-astart == 0) return;

  const int nall = atom->nlocal + atom->nghost;
  int nall_t = nall;


  const int pack_width = _fix->nbor_pack_width();

  const ATOM_T * _noalias const x = buffers->get_x();
  int * _noalias const intel_list = buffers->intel_list(list);
  const int e_nall = nall_t;

  const int molecular = atom->molecular;
  int *ns = nullptr;
  tagint *s = nullptr;
  int tag_size = 0, special_size;
  if (buffers->need_tag()) tag_size = e_nall;
  if (molecular != Atom::ATOMIC) {
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

  int * _noalias const mask = atom->mask;
  tagint * _noalias const molecule = atom->molecule;

  int tnum;
  int * _noalias overflow;
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

  const double delta = 0.01 * force->angstrom;

  {


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
      const bigint list_size = (bigint)(e_ito + tid * 2 + 2) *
        (bigint)maxnbors;

      #ifdef LMP_INTEL_3BODY_FAST
      const int pack_offset = maxnbors * pack_width;
      const bigint obound = pack_offset + maxnbors * 2;
      #else
      const int pack_offset = 0;
      const bigint obound = maxnbors * 3;
      #endif
      bigint ct = (bigint)(ifrom + tid * 2) * (bigint)maxnbors;
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
        if (THREE || (TRI && !FULL)) itag = tag[i];
        const int ioffset = ntypes * itype;

        const int ibin = atombin[i];
        if (ibin != oldbin) {
          oldbin = ibin;
          ncount = 0;
          for (int k = 0; k < nstencilp; k++) {
            const int bstart = binhead[ibin + binstart[k]];
            const int bend = binhead[ibin + binend[k]];
            #if defined(LMP_SIMD_COMPILER)
#if defined(USE_OMP_SIMD)
            #pragma omp simd
#else
            #pragma simd
#endif
            #endif
            for (int jj = bstart; jj < bend; jj++)
              tj[ncount++] = binpacked[jj];
          }
          #if defined(LMP_SIMD_COMPILER)
#if defined(USE_OMP_SIMD)
          #pragma omp simd
#else
          #pragma simd
#endif
          #pragma vector aligned
          #endif
          for (int u = 0; u < ncount; u++) {
            const int j = IP_PRE_dword_index(tj[u]);
            tx[u] = x[j].x;
            ty[u] = x[j].y;
            tz[u] = x[j].z;
            tjtype[u] = x[j].w;
            if (THREE || (TRI && !FULL)) ttag[u] = tag[j];
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
#if defined(USE_OMP_SIMD)
            #pragma omp simd
#else
            #pragma simd
#endif
            #endif
            for (int jj = bstart; jj < bend; jj++) {
              const int j = IP_PRE_dword_index(binpacked[jj]);
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

            // i bin (half) check and ghost check
            if (j < nlocal) {
              const int ijmod = (i + j) & 1;
              if (i > j) {
                if (ijmod == 0) addme = 0;
              } else if (i < j) {
                if (ijmod == 1) addme = 0;
              } else
                addme = 0;
            } else {
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
            if (FULL) {
              if (tz[u] < ztmp) addme = 0;
              if (tz[u] == ztmp) {
                if (ty[u] < ytmp) addme = 0;
                if (ty[u] == ytmp) {
                  if (tx[u] < xtmp) addme = 0;
                  if (tx[u] == xtmp && j <= i) addme = 0;
                }
              }
            } else {
              if (j <= i) addme = 0;
              if (j >= nlocal) {
                const tagint jtag = ttag[u];
                if (itag > jtag) {
                  if ((itag+jtag) % 2 == 0) addme = 0;
                } else if (itag < jtag) {
                  if ((itag+jtag) % 2 == 1) addme = 0;
                } else {
                  if (fabs(tz[u]-ztmp) > delta) {
                    if (tz[u] < ztmp) addme = 0;
                  } else if (fabs(ty[u]-ytmp) > delta) {
                    if (ty[u] < ytmp) addme = 0;
                  } else {
                    if (tx[u] < xtmp) addme = 0;
                  }
                }
              }
            }
          }

          // ghost check

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

        if (molecular != Atom::ATOMIC) {
          if (!THREE) neighptr2 = neighptr;
          int alln = n;

          n = pack_offset;
          #if defined(LMP_SIMD_COMPILER)
          #ifdef LMP_INTEL_NBOR_COMPAT
          #pragma ivdep
          #else
#if defined(USE_OMP_SIMD)
          #pragma omp simd
#else
          #pragma simd
#endif
          #endif
          #pragma vector aligned
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
            #ifdef LMP_INTEL_NBOR_COMPAT
            #pragma ivdep
            #else
#if defined(USE_OMP_SIMD)
            #pragma omp simd
#else
            #pragma simd
#endif
            #endif
            #pragma vector aligned
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

        int ns;
        if (THREE) {
          ns = n - pack_offset;
          atombin[i] = ns;
          ns += n2 - pack_offset - maxnbors;

          #ifdef LMP_INTEL_3BODY_FAST
          int alln = n;
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
          IP_PRE_neighbor_pad(pad_end);
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
          IP_PRE_neighbor_pad(pad_end);
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
                ct = (bigint)(ifrom + tid * 2) * (bigint)maxnbors;
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
              ct = (bigint)(ifrom + tid * 2) * (bigint)maxnbors;
            }
          }
        }
      }

      if (*overflow == 1)
        for (int i = ifrom; i < ito; i++)
          numneigh[i] = 0;

    } // end omp
  }

}

/* ---------------------------------------------------------------------- */


/* ---------------------------------------------------------------------- */

// ---- Half, no IC

template void NPairIntel::bin_newton<float, float, 0, 0, 0, 0, 0>
  (NeighList *, IntelBuffers<float,float> *, const int, const int);
template void NPairIntel::bin_newton<float, double, 0, 0, 0, 0, 0>
  (NeighList *, IntelBuffers<float,double> *, const int, const int);
template void NPairIntel::bin_newton<double, double, 0, 0, 0, 0, 0>
  (NeighList *, IntelBuffers<double,double> *, const int, const int);

// ---- Half, IC

template void NPairIntel::bin_newton<float, float, 0, 1, 0, 0, 0>
  (NeighList *, IntelBuffers<float,float> *, const int, const int);
template void NPairIntel::bin_newton<float, double, 0, 1, 0, 0, 0>
  (NeighList *, IntelBuffers<float,double> *, const int, const int);
template void NPairIntel::bin_newton<double, double, 0, 1, 0, 0, 0>
  (NeighList *, IntelBuffers<double,double> *, const int, const int);

// ---- Tri, no IC

template void NPairIntel::bin_newton<float, float, 0, 0, 0, 1, 0>
  (NeighList *, IntelBuffers<float,float> *, const int, const int);
template void NPairIntel::bin_newton<float, double, 0, 0, 0, 1, 0>
  (NeighList *, IntelBuffers<float,double> *, const int, const int);
template void NPairIntel::bin_newton<double, double, 0, 0, 0, 1, 0>
  (NeighList *, IntelBuffers<double,double> *, const int, const int);

// ---- Tri, IC

template void NPairIntel::bin_newton<float, float, 0, 1, 0, 1, 0>
  (NeighList *, IntelBuffers<float,float> *, const int, const int);
template void NPairIntel::bin_newton<float, double, 0, 1, 0, 1, 0>
  (NeighList *, IntelBuffers<float,double> *, const int, const int);
template void NPairIntel::bin_newton<double, double, 0, 1, 0, 1, 0>
  (NeighList *, IntelBuffers<double,double> *, const int, const int);

// ---- Full, no IC

template void NPairIntel::bin_newton<float, float, 0, 0, 1, 0, 0>
  (NeighList *, IntelBuffers<float,float> *, const int, const int);
template void NPairIntel::bin_newton<float, double, 0, 0, 1, 0, 0>
  (NeighList *, IntelBuffers<float,double> *, const int, const int);
template void NPairIntel::bin_newton<double, double, 0, 0, 1, 0, 0>
  (NeighList *, IntelBuffers<double,double> *, const int, const int);

// ---- Full, IC

template void NPairIntel::bin_newton<float, float, 0, 1, 1, 0, 0>
  (NeighList *, IntelBuffers<float,float> *, const int, const int);
template void NPairIntel::bin_newton<float, double, 0, 1, 1, 0, 0>
  (NeighList *, IntelBuffers<float,double> *, const int, const int);
template void NPairIntel::bin_newton<double, double, 0, 1, 1, 0, 0>
  (NeighList *, IntelBuffers<double,double> *, const int, const int);

// ---- 3-body, no IC

template void NPairIntel::bin_newton<float, float, 0, 0, 1, 0, 1>
  (NeighList *, IntelBuffers<float,float> *, const int, const int);
template void NPairIntel::bin_newton<float, double, 0, 0, 1, 0, 1>
  (NeighList *, IntelBuffers<float,double> *, const int, const int);
template void NPairIntel::bin_newton<double, double, 0, 0, 1, 0, 1>
  (NeighList *, IntelBuffers<double,double> *, const int, const int);

// ---- 3-body, IC

template void NPairIntel::bin_newton<float, float, 0, 1, 1, 0, 1>
  (NeighList *, IntelBuffers<float,float> *, const int, const int);
template void NPairIntel::bin_newton<float, double, 0, 1, 1, 0, 1>
  (NeighList *, IntelBuffers<float,double> *, const int, const int);
template void NPairIntel::bin_newton<double, double, 0, 1, 1, 0, 1>
  (NeighList *, IntelBuffers<double,double> *, const int, const int);

