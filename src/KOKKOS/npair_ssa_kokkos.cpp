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
   Contributing authors:
   James Larentzos and Timothy I. Mattox (Engility Corporation)
------------------------------------------------------------------------- */

#include "npair_ssa_kokkos.h"
#include "neigh_list.h"
#include "atom_kokkos.h"
#include "atom_masks.h"
#include "domain_kokkos.h"
#include "neighbor_kokkos.h"
#include "nbin_ssa_kokkos.h"
#include "nstencil_ssa.h"
#include "error.h"
#include "comm.h"

namespace LAMMPS_NS {

/* ---------------------------------------------------------------------- */

template<class DeviceType>
NPairSSAKokkos<DeviceType>::NPairSSAKokkos(LAMMPS *lmp) : NPair(lmp), ssa_phaseCt(27), ssa_gphaseCt(7)
{
  const int gphaseLenEstimate = 1; //FIXME make this 4 eventually
  k_ssa_gphaseLen = DAT::tdual_int_1d("NPairSSAKokkos:ssa_gphaseLen",ssa_gphaseCt);
  ssa_gphaseLen = k_ssa_gphaseLen.view<DeviceType>();

  k_ssa_gitemLoc = DAT::tdual_int_2d("NPairSSAKokkos::ssa_gitemLoc",ssa_gphaseCt,gphaseLenEstimate);
  ssa_gitemLoc = k_ssa_gitemLoc.view<DeviceType>();
  k_ssa_gitemLen = DAT::tdual_int_2d("NPairSSAKokkos::ssa_gitemLen",ssa_gphaseCt,gphaseLenEstimate);
  ssa_gitemLen = k_ssa_gitemLen.view<DeviceType>();
}

/* ----------------------------------------------------------------------
   copy needed info from Neighbor class to this build class
   ------------------------------------------------------------------------- */

template<class DeviceType>
void NPairSSAKokkos<DeviceType>::copy_neighbor_info()
{
  NPair::copy_neighbor_info();

  NeighborKokkos* neighborKK = (NeighborKokkos*) neighbor;

  // general params

  k_cutneighsq = neighborKK->k_cutneighsq;

  // exclusion info

  k_ex1_type = neighborKK->k_ex1_type;
  k_ex2_type = neighborKK->k_ex2_type;
  k_ex_type = neighborKK->k_ex_type;
  k_ex1_bit = neighborKK->k_ex1_bit;
  k_ex2_bit = neighborKK->k_ex2_bit;
  k_ex_mol_group = neighborKK->k_ex_mol_group;
  k_ex_mol_bit = neighborKK->k_ex_mol_bit;
  k_ex_mol_intra = neighborKK->k_ex_mol_intra;
}

/* ----------------------------------------------------------------------
 copy per-atom and per-bin vectors from NBinSSAKokkos class to this build class
 ------------------------------------------------------------------------- */

template<class DeviceType>
void NPairSSAKokkos<DeviceType>::copy_bin_info()
{
  NPair::copy_bin_info();

  NBinSSAKokkos<DeviceType>* nbKK = dynamic_cast<NBinSSAKokkos<DeviceType>*>(nb);
  if (!nbKK) error->one(FLERR, "NBin wasn't a NBinSSAKokkos object");

  atoms_per_bin = nbKK->atoms_per_bin;
  k_bincount = nbKK->k_bincount;
  k_bins = nbKK->k_bins;

  ghosts_per_gbin = nbKK->ghosts_per_gbin;
  k_gbincount = nbKK->k_gbincount;
  k_gbins = nbKK->k_gbins;

  lbinxlo = nbKK->h_lbinxlo();
  lbinxhi = nbKK->h_lbinxhi();
  lbinylo = nbKK->h_lbinylo();
  lbinyhi = nbKK->h_lbinyhi();
  lbinzlo = nbKK->h_lbinzlo();
  lbinzhi = nbKK->h_lbinzhi();
}

/* ----------------------------------------------------------------------
 copy needed info from NStencil class to this build class
 ------------------------------------------------------------------------- */

template<class DeviceType>
void NPairSSAKokkos<DeviceType>::copy_stencil_info()
{
  NPair::copy_stencil_info();

  nstencil = ns->nstencil;

  int maxstencil = ns->get_maxstencil();

  k_stencil = DAT::tdual_int_1d("NPairSSAKokkos:stencil",maxstencil);
  for (int k = 0; k < maxstencil; k++) {
    k_stencil.h_view(k) = ns->stencil[k];
  }
  k_stencil.modify<LMPHostType>();
  k_stencil.sync<DeviceType>();
  k_stencilxyz = DAT::tdual_int_1d_3("NPairSSAKokkos:stencilxyz",maxstencil);
  for (int k = 0; k < maxstencil; k++) {
    k_stencilxyz.h_view(k,0) = ns->stencilxyz[k][0];
    k_stencilxyz.h_view(k,1) = ns->stencilxyz[k][1];
    k_stencilxyz.h_view(k,2) = ns->stencilxyz[k][2];
  }
  k_stencilxyz.modify<LMPHostType>();
  k_stencilxyz.sync<DeviceType>();

  NStencilSSA *ns_ssa = dynamic_cast<NStencilSSA*>(ns);
  if (!ns_ssa) error->one(FLERR, "NStencil wasn't a NStencilSSA object");

  k_nstencil_ssa = DAT::tdual_int_1d("NPairSSAKokkos:nstencil_ssa",5);
  for (int k = 0; k < 5; ++k) {
    k_nstencil_ssa.h_view(k) = ns_ssa->nstencil_ssa[k];
  }
  k_nstencil_ssa.modify<LMPHostType>();
  k_nstencil_ssa.sync<DeviceType>();
  sx1 = ns_ssa->sx + 1;
  sy1 = ns_ssa->sy + 1;
  sz1 = ns_ssa->sz + 1;

  // Setup the phases of the workplan for locals
  ssa_phaseCt = sz1*sy1*sx1;
  if (ssa_phaseCt > (int) k_ssa_phaseLen.extent(0)) {
    k_ssa_phaseLen = DAT::tdual_int_1d("NPairSSAKokkos:ssa_phaseLen",ssa_phaseCt);
    ssa_phaseLen = k_ssa_phaseLen.view<DeviceType>();
    k_ssa_phaseOff = DAT::tdual_int_1d_3("NPairSSAKokkos:ssa_phaseOff",ssa_phaseCt);
    ssa_phaseOff = k_ssa_phaseOff.view<DeviceType>();
  }
  auto h_ssa_phaseOff = k_ssa_phaseOff.h_view;
  k_ssa_phaseOff.sync<LMPHostType>();
  int workPhase = 0;
  for (int zoff = sz1 - 1; zoff >= 0; --zoff) {
    for (int yoff = sy1 - 1; yoff >= 0; --yoff) {
      for (int xoff = sx1 - 1; xoff >= 0; --xoff) {
        h_ssa_phaseOff(workPhase, 0) = xoff;
        h_ssa_phaseOff(workPhase, 1) = yoff;
        h_ssa_phaseOff(workPhase, 2) = zoff;
        workPhase++;
      }
    }
  }
  k_ssa_phaseOff.modify<LMPHostType>();
  k_ssa_phaseOff.sync<DeviceType>();

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
int NPairSSAKokkosExecute<DeviceType>::find_special(const int &i, const int &j) const
{
  const int n1 = nspecial(i,0);
  const int n2 = nspecial(i,1);
  const int n3 = nspecial(i,2);

  for (int k = 0; k < n3; k++) {
    if (special(i,k) == tag(j)) {
      if (k < n1) {
        if (special_flag[1] == 0) return -1;
        else if (special_flag[1] == 1) return 0;
        else return 1;
      } else if (k < n2) {
        if (special_flag[2] == 0) return -1;
        else if (special_flag[2] == 1) return 0;
        else return 2;
      } else {
        if (special_flag[3] == 0) return -1;
        else if (special_flag[3] == 1) return 0;
        else return 3;
      }
    }
  }
  return 0;
};

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
int NPairSSAKokkosExecute<DeviceType>::exclusion(const int &i,const int &j,
                                             const int &itype,const int &jtype) const
{
  int m;

  if (nex_type && ex_type(itype,jtype)) return 1;

  if (nex_group) {
    for (m = 0; m < nex_group; m++) {
      if (mask(i) & ex1_bit(m) && mask(j) & ex2_bit(m)) return 1;
      if (mask(i) & ex2_bit(m) && mask(j) & ex1_bit(m)) return 1;
    }
  }

  if (nex_mol) {
    for (m = 0; m < nex_mol; m++)
      if (ex_mol_intra[m]) { // intra-chain: exclude i-j pair if on same molecule
        if (mask[i] & ex_mol_bit[m] && mask[j] & ex_mol_bit[m] &&
            molecule[i] == molecule[j]) return 1;
      } else                 // exclude i-j pair if on different molecules
        if (mask[i] & ex_mol_bit[m] && mask[j] & ex_mol_bit[m] &&
            molecule[i] != molecule[j]) return 1;
  }

  return 0;
}

/* ---------------------------------------------------------------------- */


/* ----------------------------------------------------------------------
   binned neighbor list construction with full Newton's 3rd law
   for use by Shardlow Splitting Algorithm
   each owned atom i checks its own bin and other bins in Newton stencil
   every pair stored exactly once by some processor
------------------------------------------------------------------------- */

template<class DeviceType>
void NPairSSAKokkos<DeviceType>::build(NeighList *list_)
{
  NeighListKokkos<DeviceType>* list = (NeighListKokkos<DeviceType>*) list_;
  const int nlocal = includegroup?atom->nfirst:atom->nlocal;
  int nl_size;

  int xbinCt = (lbinxhi - lbinxlo + sx1 - 1) / sx1 + 1;
  int ybinCt = (lbinyhi - lbinylo + sy1 - 1) / sy1 + 1;
  int zbinCt = (lbinzhi - lbinzlo + sz1 - 1) / sz1 + 1;
  int phaseLenEstimate = xbinCt*ybinCt*zbinCt;

  if ((ssa_phaseCt > (int) k_ssa_itemLoc.extent(0)) ||
      (phaseLenEstimate > (int) k_ssa_itemLoc.extent(1))) {
    k_ssa_itemLoc = DAT::tdual_int_2d("NPairSSAKokkos::ssa_itemLoc",ssa_phaseCt,phaseLenEstimate);
    ssa_itemLoc = k_ssa_itemLoc.view<DeviceType>();
    k_ssa_itemLen = DAT::tdual_int_2d("NPairSSAKokkos::ssa_itemLen",ssa_phaseCt,phaseLenEstimate);
    ssa_itemLen = k_ssa_itemLen.view<DeviceType>();
  }

  k_ssa_itemLoc.sync<LMPHostType>();
  k_ssa_itemLen.sync<LMPHostType>();
  k_ssa_gitemLoc.sync<LMPHostType>();
  k_ssa_gitemLen.sync<LMPHostType>();
  k_ssa_phaseOff.sync<LMPHostType>();
  k_ssa_phaseLen.sync<LMPHostType>();
  auto h_ssa_itemLoc = k_ssa_itemLoc.h_view;
  auto h_ssa_itemLen = k_ssa_itemLen.h_view;
  auto h_ssa_gitemLoc = k_ssa_gitemLoc.h_view;
  auto h_ssa_gitemLen = k_ssa_gitemLen.h_view;
  auto h_ssa_phaseOff = k_ssa_phaseOff.h_view;
  auto h_ssa_phaseLen = k_ssa_phaseLen.h_view;

{ // Preflight the neighbor list workplan
  k_bincount.sync<LMPHostType>();
  auto h_bincount = k_bincount.h_view;
  k_stencil.sync<LMPHostType>();
  auto h_stencil = k_stencil.h_view;
  k_nstencil_ssa.sync<LMPHostType>();
  auto h_nstencil_ssa = k_nstencil_ssa.h_view;
  int inum = 0;

  // loop over bins with local atoms, counting half of the neighbors
  for (int workPhase = 0; workPhase < ssa_phaseCt; ++workPhase) {
    int zoff = h_ssa_phaseOff(workPhase, 2);
    int yoff = h_ssa_phaseOff(workPhase, 1);
    int xoff = h_ssa_phaseOff(workPhase, 0);
    int workItem = 0;
  for (int zbin = lbinzlo + zoff; zbin < lbinzhi; zbin += sz1) {
  for (int ybin = lbinylo + yoff - sy1 + 1; ybin < lbinyhi; ybin += sy1) {
  for (int xbin = lbinxlo + xoff - sx1 + 1; xbin < lbinxhi; xbin += sx1) {
    int inum_start = inum;
//    if (workItem >= phaseLenEstimate) error->one(FLERR,"phaseLenEstimate was too small");

    for (int subphase = 0; subphase < 4; subphase++) {
      int s_ybin = ybin + ((subphase & 0x2) ? sy1 - 1 : 0);
      int s_xbin = xbin + ((subphase & 0x1) ? sx1 - 1 : 0);
      if ((s_ybin < lbinylo) || (s_ybin >= lbinyhi)) continue;
      if ((s_xbin < lbinxlo) || (s_xbin >= lbinxhi)) continue;

      const int ibin = zbin*mbiny*mbinx + s_ybin*mbinx + s_xbin;
      const int ibinCt = h_bincount(ibin);
      if (ibinCt > 0) {
        int base_n = 0;
        bool include_same = false;
        // count all local atoms in the current stencil "subphase" as potential neighbors
        for (int k = h_nstencil_ssa(subphase); k < h_nstencil_ssa(subphase+1); k++) {
          const int jbin = ibin+h_stencil(k);
          if (jbin != ibin) base_n += h_bincount(jbin);
          else include_same = true;
        }
        // Calculate how many ibin particles would have had some neighbors
        if (base_n > 0) inum += ibinCt;
        else if (include_same) inum += ibinCt - 1;
      }
    }
    h_ssa_itemLoc(workPhase,workItem) = inum_start; // record where workItem starts in ilist
    h_ssa_itemLen(workPhase,workItem) = inum - inum_start; // record workItem length
#ifdef DEBUG_SSA_BUILD_LOCALS
if (h_ssa_itemLen(workPhase,workItem) < 0) fprintf(stdout, "undr%03d phase (%3d,%3d) inum %d - inum_start %d UNDERFLOW\n"
  ,comm->me
  ,workPhase
  ,workItem
  ,inum
  ,inum_start
);
#endif
    workItem++;
  }
  }
  }

#ifdef DEBUG_SSA_BUILD_LOCALS
fprintf(stdout, "phas%03d phase %3d could use %6d inums, expected %6d inums. maxworkItems = %3d, inums/workItems = %g\n"
  ,comm->me
  ,workPhase
  ,inum - h_ssa_itemLoc(workPhase, 0)
  ,(nlocal*4 + ssa_phaseCt - 1) / ssa_phaseCt
  ,workItem
  ,(inum - h_ssa_itemLoc(workPhase, 0)) / (double) workItem
);
#endif
    // record where workPhase ends
    h_ssa_phaseLen(workPhase) = workItem;
  }
#ifdef DEBUG_SSA_BUILD_LOCALS
fprintf(stdout, "tota%03d total %3d could use %6d inums, expected %6d inums. inums/phase = %g\n"
  ,comm->me
  ,workPhase
  ,inum
  ,nlocal*4
  ,inum / (double) workPhase
);
#endif
  nl_size = inum; // record how much space is needed for the local work plan
}

  // count how many ghosts might have neighbors, and increase the work plan storage
  k_gbincount.sync<LMPHostType>();
  for (int workPhase = 0; workPhase < ssa_gphaseCt; workPhase++) {
    int len = k_gbincount.h_view(workPhase + 1);
    h_ssa_gitemLoc(workPhase,0) = nl_size; // record where workItem starts in ilist
    h_ssa_gitemLen(workPhase,0) = len;
    nl_size += len;
  }
  list->grow(nl_size); // Make special larger SSA neighbor list

  k_ssa_itemLoc.modify<LMPHostType>();
  k_ssa_itemLen.modify<LMPHostType>();
  k_ssa_gitemLoc.modify<LMPHostType>();
  k_ssa_gitemLen.modify<LMPHostType>();
  k_ssa_phaseLen.modify<LMPHostType>();
  k_ssa_itemLoc.sync<DeviceType>();
  k_ssa_itemLen.sync<DeviceType>();
  k_ssa_gitemLen.sync<DeviceType>();
  k_ssa_gitemLoc.sync<DeviceType>();
  k_ssa_phaseOff.sync<DeviceType>();
  k_ssa_phaseLen.sync<DeviceType>();
  k_ssa_gphaseLen.sync<DeviceType>();

  NPairSSAKokkosExecute<DeviceType>
    data(*list,
         k_cutneighsq.view<DeviceType>(),
         k_bincount.view<DeviceType>(),
         k_bins.view<DeviceType>(),
         k_gbincount.view<DeviceType>(),
         k_gbins.view<DeviceType>(),
         lbinxlo, lbinxhi, lbinylo, lbinyhi, lbinzlo, lbinzhi,
         nstencil, sx1, sy1, sz1,
         k_stencil.view<DeviceType>(),
         k_stencilxyz.view<DeviceType>(),
         k_nstencil_ssa.view<DeviceType>(),
         ssa_phaseCt,
         k_ssa_phaseLen.view<DeviceType>(),
         k_ssa_phaseOff.view<DeviceType>(),
         k_ssa_itemLoc.view<DeviceType>(),
         k_ssa_itemLen.view<DeviceType>(),
         ssa_gphaseCt,
         k_ssa_gphaseLen.view<DeviceType>(),
         k_ssa_gitemLoc.view<DeviceType>(),
         k_ssa_gitemLen.view<DeviceType>(),
         nlocal,
         atomKK->k_x.view<DeviceType>(),
         atomKK->k_type.view<DeviceType>(),
         atomKK->k_mask.view<DeviceType>(),
         atomKK->k_molecule.view<DeviceType>(),
         atomKK->k_tag.view<DeviceType>(),
         atomKK->k_special.view<DeviceType>(),
         atomKK->k_nspecial.view<DeviceType>(),
         atomKK->molecular,
         nbinx,nbiny,nbinz,mbinx,mbiny,mbinz,mbinxlo,mbinylo,mbinzlo,
         bininvx,bininvy,bininvz,
         exclude, nex_type,
         k_ex1_type.view<DeviceType>(),
         k_ex2_type.view<DeviceType>(),
         k_ex_type.view<DeviceType>(),
         nex_group,
         k_ex1_bit.view<DeviceType>(),
         k_ex2_bit.view<DeviceType>(),
         nex_mol,
         k_ex_mol_group.view<DeviceType>(),
         k_ex_mol_bit.view<DeviceType>(),
         k_ex_mol_intra.view<DeviceType>(),
         bboxhi,bboxlo,
         domain->xperiodic,domain->yperiodic,domain->zperiodic,
         domain->xprd_half,domain->yprd_half,domain->zprd_half);

  k_cutneighsq.sync<DeviceType>();
  k_ex1_type.sync<DeviceType>();
  k_ex2_type.sync<DeviceType>();
  k_ex_type.sync<DeviceType>();
  k_ex1_bit.sync<DeviceType>();
  k_ex2_bit.sync<DeviceType>();
  k_ex_mol_group.sync<DeviceType>();
  k_ex_mol_bit.sync<DeviceType>();
  k_ex_mol_intra.sync<DeviceType>();
  k_bincount.sync<DeviceType>();
  k_bins.sync<DeviceType>();
  k_gbincount.sync<DeviceType>();
  k_gbins.sync<DeviceType>();
  atomKK->sync(Device,X_MASK|TYPE_MASK|MASK_MASK|MOLECULE_MASK|TAG_MASK|SPECIAL_MASK);

  data.special_flag[0] = special_flag[0];
  data.special_flag[1] = special_flag[1];
  data.special_flag[2] = special_flag[2];
  data.special_flag[3] = special_flag[3];

  bool firstTry = true;
  data.h_resize()=1;
  while (data.h_resize()) {
    data.h_new_maxneighs() = list->maxneighs;
    data.h_resize() = 0;

    Kokkos::deep_copy(data.resize, data.h_resize);
    Kokkos::deep_copy(data.new_maxneighs, data.h_new_maxneighs);

    // loop over bins with local atoms, storing half of the neighbors
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0,ssa_phaseCt),
     LAMMPS_LAMBDA (const int workPhase) {
      data.build_locals_onePhase(firstTry, workPhase);
    });
    k_ssa_itemLoc.modify<DeviceType>();
    k_ssa_itemLen.modify<DeviceType>();
    k_ssa_phaseLen.modify<DeviceType>();
    k_ssa_itemLoc.sync<LMPHostType>();
    k_ssa_itemLen.sync<LMPHostType>();
    k_ssa_phaseLen.sync<LMPHostType>();
    data.neigh_list.inum = h_ssa_itemLoc(ssa_phaseCt-1,h_ssa_phaseLen(ssa_phaseCt-1)-1) +
      h_ssa_itemLen(ssa_phaseCt-1,h_ssa_phaseLen(ssa_phaseCt-1)-1);

    // loop over AIR ghost atoms, storing their local neighbors
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0,ssa_gphaseCt),
     LAMMPS_LAMBDA (const int workPhase) {
      data.build_ghosts_onePhase(workPhase);
    });
    k_ssa_gitemLoc.modify<DeviceType>();
    k_ssa_gitemLen.modify<DeviceType>();
    k_ssa_gphaseLen.modify<DeviceType>();
    k_ssa_gitemLoc.sync<LMPHostType>();
    k_ssa_gitemLen.sync<LMPHostType>();
    k_ssa_gphaseLen.sync<LMPHostType>();
    auto h_ssa_gphaseLen = k_ssa_gphaseLen.h_view;
    data.neigh_list.gnum = h_ssa_gitemLoc(ssa_gphaseCt-1,h_ssa_gphaseLen(ssa_gphaseCt-1)-1) +
      h_ssa_gitemLen(ssa_gphaseCt-1,h_ssa_gphaseLen(ssa_gphaseCt-1)-1) - data.neigh_list.inum;
    firstTry = false;

    Kokkos::deep_copy(data.h_resize, data.resize);

    if (data.h_resize()) {
      Kokkos::deep_copy(data.h_new_maxneighs, data.new_maxneighs);
      list->maxneighs = data.h_new_maxneighs() * 1.2;
      list->d_neighbors = typename ArrayTypes<DeviceType>::t_neighbors_2d("neighbors", list->d_neighbors.extent(0), list->maxneighs);
      data.neigh_list.d_neighbors = list->d_neighbors;
      data.neigh_list.maxneighs = list->maxneighs;
    }
  }

  //k_ssa_phaseLen.modify<DeviceType>();
  //k_ssa_itemLoc.modify<DeviceType>();
  //k_ssa_itemLen.modify<DeviceType>();
  //k_ssa_gphaseLen.modify<DeviceType>();
  //k_ssa_gitemLoc.modify<DeviceType>();
  //k_ssa_gitemLen.modify<DeviceType>();

  list->inum = data.neigh_list.inum; //FIXME once the above is in a parallel_for
  list->gnum = data.neigh_list.gnum; // it will need a deep_copy or something

#ifdef DEBUG_SSA_BUILD_LOCALS
fprintf(stdout, "Fina%03d %6d inum %6d gnum, total used %6d, allocated %6d\n"
  ,comm->me
  ,list->inum
  ,list->gnum
  ,list->inum + list->gnum
  ,nl_size
);
#endif

  list->k_ilist.template modify<DeviceType>();
}


template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void NPairSSAKokkosExecute<DeviceType>::build_locals_onePhase(const bool firstTry, int workPhase) const
{
  const typename ArrayTypes<DeviceType>::t_int_1d_const_um stencil = d_stencil;
  int which = 0;

  int zoff = d_ssa_phaseOff(workPhase, 2);
  int yoff = d_ssa_phaseOff(workPhase, 1);
  int xoff = d_ssa_phaseOff(workPhase, 0);
  int workItem = 0;
  int skippedItems = 0;
  for (int zbin = lbinzlo + zoff; zbin < lbinzhi; zbin += sz1) {
  for (int ybin = lbinylo + yoff - sy1 + 1; ybin < lbinyhi; ybin += sy1) {
  for (int xbin = lbinxlo + xoff - sx1 + 1; xbin < lbinxhi; xbin += sx1) {
    if (d_ssa_itemLen(workPhase, workItem + skippedItems) == 0) {
      if (firstTry) ++skippedItems;
      else ++workItem; // phase is done,should break out of three loops here if we could...
      continue;
    }
    int inum_start = d_ssa_itemLoc(workPhase, workItem + skippedItems);
    int inum = inum_start;

    for (int subphase = 0; subphase < 4; subphase++) {
      int s_ybin = ybin + ((subphase & 0x2) ? sy1 - 1 : 0);
      int s_xbin = xbin + ((subphase & 0x1) ? sx1 - 1 : 0);
      if ((s_ybin < lbinylo) || (s_ybin >= lbinyhi)) continue;
      if ((s_xbin < lbinxlo) || (s_xbin >= lbinxhi)) continue;

      int ibin = zbin*mbiny*mbinx + s_ybin*mbinx + s_xbin;
      for (int il = 0; il < c_bincount(ibin); ++il) {
        const int i = c_bins(ibin, il);
        int n = 0;

        const AtomNeighbors neighbors_i = neigh_list.get_neighbors(inum);
        const X_FLOAT xtmp = x(i, 0);
        const X_FLOAT ytmp = x(i, 1);
        const X_FLOAT ztmp = x(i, 2);
        const int itype = type(i);

        // loop over all local atoms in the current stencil "subphase"
        for (int k = d_nstencil_ssa(subphase); k < d_nstencil_ssa(subphase+1); k++) {
          const int jbin = ibin+stencil(k);
          int jl;
          if (jbin != ibin) jl = 0;
          else jl = il + 1; // same bin as i, so start just past i in the bin
          for (; jl < c_bincount(jbin); ++jl) {
            const int j = c_bins(jbin, jl);
            const int jtype = type(j);
            if (exclude && exclusion(i,j,itype,jtype)) continue;

            const X_FLOAT delx = xtmp - x(j, 0);
            const X_FLOAT dely = ytmp - x(j, 1);
            const X_FLOAT delz = ztmp - x(j, 2);
            const X_FLOAT rsq = delx*delx + dely*dely + delz*delz;
            if (rsq <= cutneighsq(itype,jtype)) {
              if (molecular != Atom::ATOMIC) {
                if (!moltemplate)
                  which = find_special(i,j);
                    /* else if (imol >= 0) */
                    /*   which = find_special(onemols[imol]->special[iatom], */
                    /*                        onemols[imol]->nspecial[iatom], */
                    /*                        tag[j]-tagprev); */
                    /* else which = 0; */
                if (which == 0) {
                  if (n<neigh_list.maxneighs) neighbors_i(n++) = j;
                  else n++;
                } else if (minimum_image_check(delx,dely,delz)) {
                  if (n<neigh_list.maxneighs) neighbors_i(n++) = j;
                  else n++;
                }
                else if (which > 0) {
                  if (n<neigh_list.maxneighs) neighbors_i(n++) = j ^ (which << SBBITS);
                  else n++;
                }
              } else {
                if (n<neigh_list.maxneighs) neighbors_i(n++) = j;
                else n++;
              }
            }
          }
        }

        if (n > 0) {
          neigh_list.d_numneigh(inum) = n;
          neigh_list.d_ilist(inum++) = i;
          if (n > neigh_list.maxneighs) {
            resize() = 1;
            if (n > new_maxneighs()) Kokkos::atomic_max(&new_maxneighs(),n);
          }
        }
      }
    }
#ifdef DEBUG_SSA_BUILD_LOCALS
    int len = inum - inum_start;
    if (len != d_ssa_itemLen(workPhase, workItem + skippedItems)) {
fprintf(stdout, "Leng%03d workphase (%2d,%3d,%3d): len  = %4d, but ssa_itemLen = %4d%s\n"
  ,me
  ,workPhase
  ,workItem
  ,workItem + skippedItems
  ,len
  ,d_ssa_itemLen(workPhase, workItem + skippedItems)
  ,(len > d_ssa_itemLen(workPhase, workItem + skippedItems)) ? " OVERFLOW" : ""
);
    }
#endif
    if (inum > inum_start) {
      d_ssa_itemLoc(workPhase,workItem) = inum_start; // record where workItem starts in ilist
      d_ssa_itemLen(workPhase,workItem) = inum - inum_start; // record actual workItem length
      workItem++;
    } else if (firstTry) ++skippedItems;
  }
  }
  }

#ifdef DEBUG_SSA_BUILD_LOCALS
fprintf(stdout, "Phas%03d phase %3d used %6d inums, workItems = %3d, skipped = %3d, inums/workItems = %g\n"
  ,me
  ,workPhase
  ,inum - d_ssa_itemLoc(workPhase, 0)
  ,workItem
  ,skippedItems
  ,(inum - d_ssa_itemLoc(workPhase, 0)) / (double) workItem
);
#endif
    // record where workPhase actually ends
    if (firstTry) {
      d_ssa_phaseLen(workPhase) = workItem;
      while (workItem < (int) d_ssa_itemLen.extent(1)) {
        d_ssa_itemLen(workPhase,workItem++) = 0;
      }
    }

}


template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void NPairSSAKokkosExecute<DeviceType>::build_ghosts_onePhase(int workPhase) const
{
  const typename ArrayTypes<DeviceType>::t_int_1d_const_um stencil = d_stencil;
  int which = 0;

  // since these are ghosts, must check if stencil bin is out of bounds
    int airnum = workPhase + 1;
    //FIXME for now, there is only 1 workItem for each ghost AIR
    int workItem;
    for (workItem = 0; workItem < 1; ++workItem) {
      int gNdx = d_ssa_gitemLoc(workPhase, workItem); // record where workItem starts in ilist
      for (int il = 0; il < c_gbincount(airnum); ++il) {
        const int i = c_gbins(airnum, il);
        int n = 0;

        const AtomNeighbors neighbors_i = neigh_list.get_neighbors(gNdx);
        const X_FLOAT xtmp = x(i, 0);
        const X_FLOAT ytmp = x(i, 1);
        const X_FLOAT ztmp = x(i, 2);
        const int itype = type(i);

        int loc[3];
        const int ibin = coord2bin(x(i, 0), x(i, 1), x(i, 2), &(loc[0]));

        // loop over AIR ghost atoms in all bins in "full" stencil
        // Note: the non-AIR ghost atoms have already been filtered out
        for (int k = 0; k < nstencil; k++) {
          int xbin2 = loc[0] + d_stencilxyz(k,0);
          int ybin2 = loc[1] + d_stencilxyz(k,1);
          int zbin2 = loc[2] + d_stencilxyz(k,2);
          // Skip it if this bin is outside the extent of local bins
          if (xbin2 < lbinxlo || xbin2 >= lbinxhi ||
              ybin2 < lbinylo || ybin2 >= lbinyhi ||
              zbin2 < lbinzlo || zbin2 >= lbinzhi) continue;
          const int jbin = ibin+stencil(k);
          for (int jl = 0; jl < c_bincount(jbin); ++jl) {
            const int j = c_bins(jbin, jl);
            const int jtype = type(j);
            if (exclude && exclusion(i,j,itype,jtype)) continue;

            const X_FLOAT delx = xtmp - x(j, 0);
            const X_FLOAT dely = ytmp - x(j, 1);
            const X_FLOAT delz = ztmp - x(j, 2);
            const X_FLOAT rsq = delx*delx + dely*dely + delz*delz;
            if (rsq <= cutneighsq(itype,jtype)) {
              if (molecular != Atom::ATOMIC) {
                if (!moltemplate)
                  which = find_special(j,i);
                    /* else if (jmol >= 0) */
                    /*   which = find_special(onemols[jmol]->special[jatom], */
                    /*                        onemols[jmol]->nspecial[jatom], */
                    /*                        tag[i]-jtagprev); */
                    /* else which = 0; */
                if (which == 0) {
                  if (n<neigh_list.maxneighs) neighbors_i(n++) = j;
                  else n++;
                } else if (minimum_image_check(delx,dely,delz)) {
                  if (n<neigh_list.maxneighs) neighbors_i(n++) = j;
                  else n++;
                }
                else if (which > 0) {
                  if (n<neigh_list.maxneighs) neighbors_i(n++) = j ^ (which << SBBITS);
                  else n++;
                }
              } else {
                if (n<neigh_list.maxneighs) neighbors_i(n++) = j;
                else n++;
              }
            }
          }
        }

        if (n > 0) {
          neigh_list.d_numneigh(gNdx) = n;
          neigh_list.d_ilist(gNdx++) = i;
          if (n > neigh_list.maxneighs) {
            resize() = 1;
            if (n > new_maxneighs()) Kokkos::atomic_max(&new_maxneighs(),n);
          }
        }
      }
      // record where workItem ends in ilist
      d_ssa_gitemLen(workPhase,workItem) = gNdx - d_ssa_gitemLoc(workPhase,workItem);
      // if (d_ssa_gitemLen(workPhase,workItem) > 0) workItem++;
    }
    d_ssa_gphaseLen(workPhase) = workItem;
}

}

namespace LAMMPS_NS {
template class NPairSSAKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class NPairSSAKokkos<LMPHostType>;
#endif
}
