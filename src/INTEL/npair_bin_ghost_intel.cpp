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
   Contributing authors: W. Michael Brown (Intel)
------------------------------------------------------------------------- */

#include "npair_bin_ghost_intel.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "neigh_list.h"
#include "neighbor.h"

#include <cmath>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NPairFullBinGhostIntel::NPairFullBinGhostIntel(LAMMPS *lmp) : NPairIntel(lmp) {}

/* ----------------------------------------------------------------------
   binned neighbor list construction for all neighbors
   include neighbors of ghost atoms, but no "special neighbors" for ghosts
   every neighbor pair appears in list of both atoms i and j
------------------------------------------------------------------------- */

void NPairFullBinGhostIntel::build(NeighList *list)
{

  if (nstencil > INTEL_MAX_STENCIL_CHECK)
    error->all(FLERR, "Too many neighbor bins for INTEL package" + utils::errorurl(9));


  if (_fix->precision() == FixIntel::PREC_MODE_MIXED)
    fbi(list, _fix->get_mixed_buffers());
  else if (_fix->precision() == FixIntel::PREC_MODE_DOUBLE)
    fbi(list, _fix->get_double_buffers());
  else
    fbi(list, _fix->get_single_buffers());

  _fix->stop_watch(TIME_HOST_NEIGHBOR);
}

/* ---------------------------------------------------------------------- */

template<class flt_t, class acc_t>
void NPairFullBinGhostIntel::fbi(NeighList * list,
                                 IntelBuffers<flt_t,acc_t> * buffers)
{
  const int nlocal = atom->nlocal;
  const int nall = atom->nlocal + atom->nghost;
  list->inum = atom->nlocal;
  list->gnum = atom->nghost;

  buffers->grow_list(list, nall, comm->nthreads, 0,
                     _fix->nbor_pack_width());

  int need_ic = 0;
  if (atom->molecular != Atom::ATOMIC)
    dminimum_image_check(need_ic, neighbor->cutneighmax, neighbor->cutneighmax,
                         neighbor->cutneighmax);

  if (need_ic) {
    fbi<flt_t,acc_t,1>(list, buffers, 0, nlocal);
  } else {
    fbi<flt_t,acc_t,0>(list, buffers, 0, nlocal);
  }
}

/* ---------------------------------------------------------------------- */

template<class flt_t, class acc_t, int need_ic>
void NPairFullBinGhostIntel::fbi(NeighList * list,
                                 IntelBuffers<flt_t,acc_t> * buffers,
                                 const int pstart, const int pend) {
  if (pend-pstart == 0) return;

  const int nall = atom->nlocal + atom->nghost;
  int nall_t = nall;
  const int aend = nall;

  const ATOM_T * _noalias const x = buffers->get_x();
  int * _noalias const intel_list = buffers->intel_list(list);
  int ** _noalias const firstneigh = list->firstneigh;  // NOLINT
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
  int * _noalias numneigh = list->numneigh;
  const int nstencil = this->nstencil;
  const int * _noalias const stencil = this->stencil;
  const flt_t * _noalias const cutneighsq = buffers->get_cutneighsq()[0];
  const flt_t * _noalias const cutneighghostsq =
    buffers->get_cutneighghostsq()[0];
  const int ntypes = atom->ntypes + 1;
  const int nlocal = atom->nlocal;

  int * _noalias const mask = atom->mask;
  tagint * _noalias const molecule = atom->molecule;

  int moltemplate;
  if (molecular == Atom::TEMPLATE) moltemplate = 1;
  else moltemplate = 0;
  if (moltemplate)
    error->all(FLERR,
               "Can't use moltemplate with npair style full/bin/ghost/intel.");

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

  const int mbinx = this->mbinx;
  const int mbiny = this->mbiny;
  const int mbinz = this->mbinz;
  const int * _noalias const stencilxyz = &this->stencilxyz[0][0];

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

    const int mbinyx = mbiny * mbinx;

    #if defined(_OPENMP)
    #pragma omp parallel
    #endif
    {
      const int num = aend;
      int tid, ifrom, ito;

      const double balance_factor = 2.0;
      const double ibalance_factor = 1.0 / balance_factor;
      const int gnum = num - nlocal;
      const int wlocal = static_cast<int>(ceil(balance_factor * nlocal));
      const int snum = wlocal + gnum;
      IP_PRE_omp_range_id(ifrom, ito, tid, snum, nthreads);
      if (ifrom < wlocal) ifrom = static_cast<int>(ibalance_factor * ifrom);
      else ifrom -= wlocal - nlocal;
      if (ito < wlocal) ito = static_cast<int>(ibalance_factor * ito);
      else ito -= wlocal - nlocal;

      int e_ito = ito;
      const bigint list_size = (bigint)(e_ito + tid * 2 + 2) *
        (bigint)maxnbors;

      int pack_offset = maxnbors;
      bigint ct = (bigint)(ifrom + tid * 2) * (bigint)maxnbors;
      int * _noalias neighptr = intel_list + ct;
      const bigint obound = pack_offset + maxnbors * 2;

      const int toffs = tid * ncache_stride;
      flt_t * _noalias const tx = ncachex + toffs;
      flt_t * _noalias const ty = ncachey + toffs;
      flt_t * _noalias const tz = ncachez + toffs;
      int * _noalias const tj = ncachej + toffs;
      int * _noalias const tjtype = ncachejtype + toffs;
      tagint * _noalias const ttag = ncachetag + toffs;

      // loop over all atoms in other bins in stencil, store every pair
      int ncount, oldbin = -9999999;
      for (int i = ifrom; i < ito; i++) {
        const flt_t xtmp = x[i].x;
        const flt_t ytmp = x[i].y;
        const flt_t ztmp = x[i].z;
        const int itype = x[i].w;
        const tagint itag = tag[i];
        const int ioffset = ntypes * itype;

        const int ibin = atombin[i];
        if (ibin != oldbin) {
          oldbin = ibin;
          ncount = 0;
          if (i < nlocal) {
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
          } else {
            const int zbin = ibin / mbinyx;
            const int zrem = ibin % mbinyx;
            const int ybin = zrem / mbinx;
            const int xbin = zrem % mbinx;
            for (int k = 0; k < nstencil; k++) {
              const int xbin2 = xbin + stencilxyz[3 * k + 0];
              const int ybin2 = ybin + stencilxyz[3 * k + 1];
              const int zbin2 = zbin + stencilxyz[3 * k + 2];
              if (xbin2 < 0 || xbin2 >= mbinx ||
                  ybin2 < 0 || ybin2 >= mbiny ||
                  zbin2 < 0 || zbin2 >= mbinz) continue;

              const int bstart = binhead[ibin + stencil[k]];
              const int bend = binhead[ibin + stencil[k] + 1];
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
          } // if i < nlocal
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
            ttag[u] = tag[j];
          }
        } // if ibin != oldbin

        // ---------------------- Loop over other bins

        int n = maxnbors;
        int n2 = n * 2;
        int * _noalias neighptr2 = neighptr;
        const flt_t * _noalias cutsq;
        if (i < nlocal) cutsq = cutneighsq;
        else cutsq = cutneighghostsq;

        const int icp = i;

        #if defined(LMP_SIMD_COMPILER)
        #pragma vector aligned
        #pragma ivdep
        #endif
        for (int u = 0; u < ncount; u++) {
          int addme = 1;
          int j = tj[u];

          if (i == j) addme = 0;

          // Cutoff Check
          const flt_t delx = xtmp - tx[u];
          const flt_t dely = ytmp - ty[u];
          const flt_t delz = ztmp - tz[u];
          const int jtype = tjtype[u];
          const tagint jtag = ttag[u];
          const flt_t rsq = delx * delx + dely * dely + delz * delz;
          if (rsq > cutsq[ioffset + jtype]) addme = 0;

          if (need_ic && icp < nlocal) {
            int no_special;
            ominimum_image_check(no_special, delx, dely, delz);
            if (no_special)
              j = -j - 1;
          }

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
        } // for u

        if ((molecular != Atom::ATOMIC) && (i < nlocal)) {
          int alln = n;
          n = 0;
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
          for (int u = 0; u < alln; u++) {
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
            neighptr2[n++] = j;
            #endif
          }
          alln = n2;
          n2 = maxnbors * 2;
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
            neighptr2[n2++] = j;
            #endif
          }
        }

        if (exclude) {
          int alln = n;
          n = maxnbors;
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
          alln = n2;
          n2 = maxnbors * 2;
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

        int ns = n - maxnbors;
        int alln = n;
        atombin[i] = ns;
        n = 0;
        for (int u = maxnbors; u < alln; u++)
          neighptr[n++] = neighptr2[u];
        ns += n2 - maxnbors * 2;
        for (int u = maxnbors * 2; u < n2; u++)
          neighptr[n++] = neighptr2[u];
        if (ns > maxnbors) *overflow = 1;

        ilist[i] = i;
        firstneigh[i] = intel_list + ct;
        numneigh[i] = ns;

        ct += ns;
        IP_PRE_edge_align(ct, sizeof(int));
        neighptr = intel_list + ct;
        if (ct + obound > list_size) {
          if (i < ito - 1) {
            *overflow = 1;
            ct = (bigint)(ifrom + tid * 2) * (bigint)maxnbors;
          }
        }
      }

      if (*overflow == 1)
        for (int i = ifrom; i < ito; i++)
          numneigh[i] = 0;

    } // end omp
  }

}
