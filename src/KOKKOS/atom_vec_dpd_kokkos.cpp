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

#include "atom_vec_dpd_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm_kokkos.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "modify.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecDPDKokkos::AtomVecDPDKokkos(LAMMPS *lmp) : AtomVec(lmp),
AtomVecKokkos(lmp), AtomVecDPD(lmp)
{
}

/* ----------------------------------------------------------------------
   process field strings to initialize data structs for all other methods
------------------------------------------------------------------------- */

void AtomVecDPDKokkos::init()
{
  AtomVecDPD::init();

  set_atom_masks();
}

/* ----------------------------------------------------------------------
   grow atom arrays
   n = 0 grows arrays by DELTA
   n > 0 allocates arrays to size n
------------------------------------------------------------------------- */

void AtomVecDPDKokkos::grow(int n)
{
  auto DELTA = LMP_KOKKOS_AV_DELTA;
  int step = MAX(DELTA,nmax*0.01);
  if (n == 0) nmax += step;
  else nmax = n;
  atomKK->nmax = nmax;
  if (nmax < 0 || nmax > MAXSMALLINT)
    error->one(FLERR,"Per-processor system is too big");

  atomKK->sync(Device,ALL_MASK);
  atomKK->modified(Device,ALL_MASK);

  memoryKK->grow_kokkos(atomKK->k_tag,atomKK->tag,nmax,"atom:tag");
  memoryKK->grow_kokkos(atomKK->k_type,atomKK->type,nmax,"atom:type");
  memoryKK->grow_kokkos(atomKK->k_mask,atomKK->mask,nmax,"atom:mask");
  memoryKK->grow_kokkos(atomKK->k_image,atomKK->image,nmax,"atom:image");

  memoryKK->grow_kokkos(atomKK->k_x,atomKK->x,nmax,"atom:x");
  memoryKK->grow_kokkos(atomKK->k_v,atomKK->v,nmax,"atom:v");
  memoryKK->grow_kokkos(atomKK->k_f,atomKK->f,nmax,"atom:f");


  memoryKK->grow_kokkos(atomKK->k_rho,atomKK->rho,nmax,"atom:rho");
  memoryKK->grow_kokkos(atomKK->k_dpdTheta,atomKK->dpdTheta,nmax,"atom:dpdTheta");
  memoryKK->grow_kokkos(atomKK->k_uCond,atomKK->uCond,nmax,"atom:uCond");
  memoryKK->grow_kokkos(atomKK->k_uMech,atomKK->uMech,nmax,"atom:uMech");
  memoryKK->grow_kokkos(atomKK->k_uChem,atomKK->uChem,nmax,"atom:uChem");
  memoryKK->grow_kokkos(atomKK->k_uCG,atomKK->uCG,nmax,"atom:uCG");
  memoryKK->grow_kokkos(atomKK->k_uCGnew,atomKK->uCGnew,nmax,"atom:uCGnew");
  memoryKK->grow_kokkos(atomKK->k_duChem,atomKK->duChem,nmax,"atom:duChem");

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->grow_arrays(nmax);

  grow_pointers();
  atomKK->sync(Host,ALL_MASK);
}

/* ----------------------------------------------------------------------
   reset local array ptrs
------------------------------------------------------------------------- */

void AtomVecDPDKokkos::grow_pointers()
{
  tag = atomKK->tag;
  d_tag = atomKK->k_tag.view_device();
  h_tag = atomKK->k_tag.view_host();

  type = atomKK->type;
  d_type = atomKK->k_type.view_device();
  h_type = atomKK->k_type.view_host();
  mask = atomKK->mask;
  d_mask = atomKK->k_mask.view_device();
  h_mask = atomKK->k_mask.view_host();
  image = atomKK->image;
  d_image = atomKK->k_image.view_device();
  h_image = atomKK->k_image.view_host();

  x = atomKK->x;
  d_x = atomKK->k_x.view_device();
  h_x = atomKK->k_x.view_hostkk();
  v = atomKK->v;
  d_v = atomKK->k_v.view_device();
  h_v = atomKK->k_v.view_hostkk();
  f = atomKK->f;
  d_f = atomKK->k_f.view_device();
  h_f = atomKK->k_f.view_hostkk();

  rho = atomKK->rho;
  d_rho = atomKK->k_rho.view_device();
  h_rho = atomKK->k_rho.view_hostkk();
  dpdTheta = atomKK->dpdTheta;
  d_dpdTheta = atomKK->k_dpdTheta.view_device();
  h_dpdTheta = atomKK->k_dpdTheta.view_hostkk();
  uCond = atomKK->uCond;
  d_uCond = atomKK->k_uCond.view_device();
  h_uCond = atomKK->k_uCond.view_hostkk();
  uMech = atomKK->uMech;
  d_uMech = atomKK->k_uMech.view_device();
  h_uMech = atomKK->k_uMech.view_hostkk();
  uChem = atomKK->uChem;
  d_uChem = atomKK->k_uChem.view_device();
  h_uChem = atomKK->k_uChem.view_hostkk();
  uCG = atomKK->uCG;
  d_uCG = atomKK->k_uCG.view_device();
  h_uCG = atomKK->k_uCG.view_hostkk();
  uCGnew = atomKK->uCGnew;
  d_uCGnew = atomKK->k_uCGnew.view_device();
  h_uCGnew = atomKK->k_uCGnew.view_hostkk();
  duChem = atomKK->duChem;
  d_duChem = atomKK->k_duChem.view_device();
  h_duChem = atomKK->k_duChem.view_hostkk();
}

/* ----------------------------------------------------------------------
   sort atom arrays on device
------------------------------------------------------------------------- */

void AtomVecDPDKokkos::sort_kokkos(Kokkos::BinSort<KeyViewType, BinOp> &Sorter)
{
  atomKK->sync(Device, ALL_MASK & ~F_MASK & ~DPDRHO_MASK & ~DUCHEM_MASK & ~DVECTOR_MASK);

  Sorter.sort(LMPDeviceType(), d_tag);
  Sorter.sort(LMPDeviceType(), d_type);
  Sorter.sort(LMPDeviceType(), d_mask);
  Sorter.sort(LMPDeviceType(), d_image);
  Sorter.sort(LMPDeviceType(), d_x);
  Sorter.sort(LMPDeviceType(), d_v);
  Sorter.sort(LMPDeviceType(), d_dpdTheta);
  Sorter.sort(LMPDeviceType(), d_uCond);
  Sorter.sort(LMPDeviceType(), d_uMech);
  Sorter.sort(LMPDeviceType(), d_uChem);
  Sorter.sort(LMPDeviceType(), d_uCG);
  Sorter.sort(LMPDeviceType(), d_uCGnew);

  atomKK->modified(Device, ALL_MASK & ~F_MASK & ~DPDRHO_MASK & ~DUCHEM_MASK & ~DVECTOR_MASK);
}

/* ---------------------------------------------------------------------- */

void AtomVecDPDKokkos::sync(ExecutionSpace space, uint64_t mask)
{
  if (space == Device) {
    if (mask & X_MASK) atomKK->k_x.sync_device();
    if (mask & V_MASK) atomKK->k_v.sync_device();
    if (mask & F_MASK) atomKK->k_f.sync_device();
    if (mask & TAG_MASK) atomKK->k_tag.sync_device();
    if (mask & TYPE_MASK) atomKK->k_type.sync_device();
    if (mask & MASK_MASK) atomKK->k_mask.sync_device();
    if (mask & IMAGE_MASK) atomKK->k_image.sync_device();
    if (mask & DPDRHO_MASK) atomKK->k_rho.sync_device();
    if (mask & DPDTHETA_MASK) atomKK->k_dpdTheta.sync_device();
    if (mask & UCOND_MASK) atomKK->k_uCond.sync_device();
    if (mask & UMECH_MASK) atomKK->k_uMech.sync_device();
    if (mask & UCHEM_MASK) atomKK->k_uChem.sync_device();
    if (mask & UCG_MASK) atomKK->k_uCG.sync_device();
    if (mask & UCGNEW_MASK) atomKK->k_uCGnew.sync_device();
    if (mask & DUCHEM_MASK) atomKK->k_duChem.sync_device();
  } else if (space == Host) {
    if (mask & X_MASK) atomKK->k_x.sync_host();
    if (mask & V_MASK) atomKK->k_v.sync_host();
    if (mask & F_MASK) atomKK->k_f.sync_host();
    if (mask & TAG_MASK) atomKK->k_tag.sync_host();
    if (mask & TYPE_MASK) atomKK->k_type.sync_host();
    if (mask & MASK_MASK) atomKK->k_mask.sync_host();
    if (mask & IMAGE_MASK) atomKK->k_image.sync_host();
    if (mask & DPDRHO_MASK) atomKK->k_rho.sync_host();
    if (mask & DPDTHETA_MASK) atomKK->k_dpdTheta.sync_host();
    if (mask & UCOND_MASK) atomKK->k_uCond.sync_host();
    if (mask & UMECH_MASK) atomKK->k_uMech.sync_host();
    if (mask & UCHEM_MASK) atomKK->k_uChem.sync_host();
    if (mask & UCG_MASK) atomKK->k_uCG.sync_host();
    if (mask & UCGNEW_MASK) atomKK->k_uCGnew.sync_host();
    if (mask & DUCHEM_MASK) atomKK->k_duChem.sync_host();
  } else if (space == HostKK) {
    if (mask & X_MASK) atomKK->k_x.sync_hostkk();
    if (mask & V_MASK) atomKK->k_v.sync_hostkk();
    if (mask & F_MASK) atomKK->k_f.sync_hostkk();
    if (mask & TAG_MASK) atomKK->k_tag.sync_host();
    if (mask & TYPE_MASK) atomKK->k_type.sync_host();
    if (mask & MASK_MASK) atomKK->k_mask.sync_host();
    if (mask & IMAGE_MASK) atomKK->k_image.sync_host();
    if (mask & DPDRHO_MASK) atomKK->k_rho.sync_hostkk();
    if (mask & DPDTHETA_MASK) atomKK->k_dpdTheta.sync_hostkk();
    if (mask & UCOND_MASK) atomKK->k_uCond.sync_hostkk();
    if (mask & UMECH_MASK) atomKK->k_uMech.sync_hostkk();
    if (mask & UCHEM_MASK) atomKK->k_uChem.sync_hostkk();
    if (mask & UCG_MASK) atomKK->k_uCG.sync_hostkk();
    if (mask & UCGNEW_MASK) atomKK->k_uCGnew.sync_hostkk();
    if (mask & DUCHEM_MASK) atomKK->k_duChem.sync_hostkk();
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecDPDKokkos::sync_pinned(ExecutionSpace space, uint64_t mask, int async_flag)
{
  if (space == Device) {
    if ((mask & X_MASK) && atomKK->k_x.need_sync_device())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d_3_lr>(atomKK->k_x,space,async_flag);
    if ((mask & V_MASK) && atomKK->k_v.need_sync_device())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d_3>(atomKK->k_v,space,async_flag);
    if ((mask & F_MASK) && atomKK->k_f.need_sync_device())
      perform_pinned_copy_transform<DAT::ttransform_kkacc_1d_3>(atomKK->k_f,space,async_flag);
    if ((mask & TAG_MASK) && atomKK->k_tag.need_sync_device())
      perform_pinned_copy<DAT::tdual_tagint_1d>(atomKK->k_tag,space,async_flag);
    if ((mask & TYPE_MASK) && atomKK->k_type.need_sync_device())
      perform_pinned_copy<DAT::tdual_int_1d>(atomKK->k_type,space,async_flag);
    if ((mask & MASK_MASK) && atomKK->k_mask.need_sync_device())
      perform_pinned_copy<DAT::tdual_int_1d>(atomKK->k_mask,space,async_flag);
    if ((mask & IMAGE_MASK) && atomKK->k_image.need_sync_device())
      perform_pinned_copy<DAT::tdual_imageint_1d>(atomKK->k_image,space,async_flag);
    if ((mask & DPDRHO_MASK) && atomKK->k_rho.need_sync_device())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d>(atomKK->k_rho,space,async_flag);
    if ((mask & DPDTHETA_MASK) && atomKK->k_dpdTheta.need_sync_device())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d>(atomKK->k_dpdTheta,space,async_flag);
    if ((mask & UCOND_MASK) && atomKK->k_uCond.need_sync_device())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d>(atomKK->k_uCond,space,async_flag);
    if ((mask & UMECH_MASK) && atomKK->k_uMech.need_sync_device())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d>(atomKK->k_uMech,space,async_flag);
    if ((mask & UCHEM_MASK) && atomKK->k_uChem.need_sync_device())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d>(atomKK->k_uChem,space,async_flag);
    if ((mask & UCG_MASK) && atomKK->k_uCG.need_sync_device())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d>(atomKK->k_uCG,space,async_flag);
    if ((mask & UCGNEW_MASK) && atomKK->k_uCGnew.need_sync_device())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d>(atomKK->k_uCGnew,space,async_flag);
    if ((mask & DUCHEM_MASK) && atomKK->k_duChem.need_sync_device())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d>(atomKK->k_duChem,space,async_flag);
  } else {
    if ((mask & X_MASK) && atomKK->k_x.need_sync_host())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d_3_lr>(atomKK->k_x,space,async_flag);
    if ((mask & V_MASK) && atomKK->k_v.need_sync_host())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d_3>(atomKK->k_v,space,async_flag);
    if ((mask & F_MASK) && atomKK->k_f.need_sync_host())
      perform_pinned_copy_transform<DAT::ttransform_kkacc_1d_3>(atomKK->k_f,space,async_flag);
    if ((mask & TAG_MASK) && atomKK->k_tag.need_sync_host())
      perform_pinned_copy<DAT::tdual_tagint_1d>(atomKK->k_tag,space,async_flag);
    if ((mask & TYPE_MASK) && atomKK->k_type.need_sync_host())
      perform_pinned_copy<DAT::tdual_int_1d>(atomKK->k_type,space,async_flag);
    if ((mask & MASK_MASK) && atomKK->k_mask.need_sync_host())
      perform_pinned_copy<DAT::tdual_int_1d>(atomKK->k_mask,space,async_flag);
    if ((mask & IMAGE_MASK) && atomKK->k_image.need_sync_host())
      perform_pinned_copy<DAT::tdual_imageint_1d>(atomKK->k_image,space,async_flag);
    if ((mask & DPDRHO_MASK) && atomKK->k_rho.need_sync_host())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d>(atomKK->k_rho,space,async_flag);
    if ((mask & DPDTHETA_MASK) && atomKK->k_dpdTheta.need_sync_host())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d>(atomKK->k_dpdTheta,space,async_flag);
    if ((mask & UCOND_MASK) && atomKK->k_uCond.need_sync_host())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d>(atomKK->k_uCond,space,async_flag);
    if ((mask & UMECH_MASK) && atomKK->k_uMech.need_sync_host())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d>(atomKK->k_uMech,space,async_flag);
    if ((mask & UCHEM_MASK) && atomKK->k_uChem.need_sync_host())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d>(atomKK->k_uChem,space,async_flag);
    if ((mask & UCG_MASK) && atomKK->k_uCG.need_sync_host())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d>(atomKK->k_uCG,space,async_flag);
    if ((mask & UCGNEW_MASK) && atomKK->k_uCGnew.need_sync_host())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d>(atomKK->k_uCGnew,space,async_flag);
    if ((mask & DUCHEM_MASK) && atomKK->k_duChem.need_sync_host())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d>(atomKK->k_duChem,space,async_flag);
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecDPDKokkos::modified(ExecutionSpace space, uint64_t mask)
{
  if (space == Device) {
    if (mask & X_MASK) atomKK->k_x.modify_device();
    if (mask & V_MASK) atomKK->k_v.modify_device();
    if (mask & F_MASK) atomKK->k_f.modify_device();
    if (mask & TAG_MASK) atomKK->k_tag.modify_device();
    if (mask & TYPE_MASK) atomKK->k_type.modify_device();
    if (mask & MASK_MASK) atomKK->k_mask.modify_device();
    if (mask & IMAGE_MASK) atomKK->k_image.modify_device();
    if (mask & DPDRHO_MASK) atomKK->k_rho.modify_device();
    if (mask & DPDTHETA_MASK) atomKK->k_dpdTheta.modify_device();
    if (mask & UCOND_MASK) atomKK->k_uCond.modify_device();
    if (mask & UMECH_MASK) atomKK->k_uMech.modify_device();
    if (mask & UCHEM_MASK) atomKK->k_uChem.modify_device();
    if (mask & UCG_MASK) atomKK->k_uCG.modify_device();
    if (mask & UCGNEW_MASK) atomKK->k_uCGnew.modify_device();
    if (mask & DUCHEM_MASK) atomKK->k_duChem.modify_device();
  } else if (space == Host) {
    if (mask & X_MASK) atomKK->k_x.modify_host();
    if (mask & V_MASK) atomKK->k_v.modify_host();
    if (mask & F_MASK) atomKK->k_f.modify_host();
    if (mask & TAG_MASK) atomKK->k_tag.modify_host();
    if (mask & TYPE_MASK) atomKK->k_type.modify_host();
    if (mask & MASK_MASK) atomKK->k_mask.modify_host();
    if (mask & IMAGE_MASK) atomKK->k_image.modify_host();
    if (mask & DPDRHO_MASK) atomKK->k_rho.modify_host();
    if (mask & DPDTHETA_MASK) atomKK->k_dpdTheta.modify_host();
    if (mask & UCOND_MASK) atomKK->k_uCond.modify_host();
    if (mask & UMECH_MASK) atomKK->k_uMech.modify_host();
    if (mask & UCHEM_MASK) atomKK->k_uChem.modify_host();
    if (mask & UCG_MASK) atomKK->k_uCG.modify_host();
    if (mask & UCGNEW_MASK) atomKK->k_uCGnew.modify_host();
    if (mask & DUCHEM_MASK) atomKK->k_duChem.modify_host();
  } else if (space == HostKK) {
    if (mask & X_MASK) atomKK->k_x.modify_hostkk();
    if (mask & V_MASK) atomKK->k_v.modify_hostkk();
    if (mask & F_MASK) atomKK->k_f.modify_hostkk();
    if (mask & TAG_MASK) atomKK->k_tag.modify_host();
    if (mask & TYPE_MASK) atomKK->k_type.modify_host();
    if (mask & MASK_MASK) atomKK->k_mask.modify_host();
    if (mask & IMAGE_MASK) atomKK->k_image.modify_host();
    if (mask & DPDRHO_MASK) atomKK->k_rho.modify_hostkk();
    if (mask & DPDTHETA_MASK) atomKK->k_dpdTheta.modify_hostkk();
    if (mask & UCOND_MASK) atomKK->k_uCond.modify_hostkk();
    if (mask & UMECH_MASK) atomKK->k_uMech.modify_hostkk();
    if (mask & UCHEM_MASK) atomKK->k_uChem.modify_hostkk();
    if (mask & UCG_MASK) atomKK->k_uCG.modify_hostkk();
    if (mask & UCGNEW_MASK) atomKK->k_uCGnew.modify_hostkk();
    if (mask & DUCHEM_MASK) atomKK->k_duChem.modify_hostkk();
  }
}
