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

/* ------------------------------------------------------------------------
   Contributing authors: Julien Tranchida (SNL)
                         Aidan Thompson (SNL)

   Please cite the related publication:
   Tranchida, J., Plimpton, S. J., Thibaudeau, P., & Thompson, A. P. (2018).
   Massively parallel symplectic algorithm for coupled magnetic spin dynamics
   and molecular dynamics. Journal of Computational Physics.
------------------------------------------------------------------------- */

#include "atom_vec_spin_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "memory_kokkos.h"
#include "modify.h"

using namespace LAMMPS_NS;

static constexpr int DELTA = 10;

/* ---------------------------------------------------------------------- */

AtomVecSpinKokkos::AtomVecSpinKokkos(LAMMPS *lmp) : AtomVec(lmp),
AtomVecKokkos(lmp), AtomVecSpin(lmp)
{
}

/* ----------------------------------------------------------------------
   process field strings to initialize data structs for all other methods
------------------------------------------------------------------------- */

void AtomVecSpinKokkos::init()
{
  AtomVecSpin::init();

  set_atom_masks();
}

/* ----------------------------------------------------------------------
   grow atom arrays
   n = 0 grows arrays by a chunk
   n > 0 allocates arrays to size n
------------------------------------------------------------------------- */

void AtomVecSpinKokkos::grow(int n)
{
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

  // allocating mech. quantities

  memoryKK->grow_kokkos(atomKK->k_x,atomKK->x,nmax,"atom:x");
  memoryKK->grow_kokkos(atomKK->k_v,atomKK->v,nmax,"atom:v");
  memoryKK->grow_kokkos(atomKK->k_f,atomKK->f,nmax,"atom:f");

  // allocating mag. quantities

  memoryKK->grow_kokkos(atomKK->k_sp,atomKK->sp,nmax,"atom:sp");
  memoryKK->grow_kokkos(atomKK->k_fm,atomKK->fm,nmax,"atom:fm");
  memoryKK->grow_kokkos(atomKK->k_fm_long,atomKK->fm_long,nmax,"atom:fm_long");

  grow_pointers();
  atomKK->sync(Host,ALL_MASK);

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->grow_arrays(nmax);
}

/* ----------------------------------------------------------------------
   reset local array ptrs
------------------------------------------------------------------------- */

void AtomVecSpinKokkos::grow_pointers()
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

  sp = atomKK->sp;
  d_sp = atomKK->k_sp.view_device();
  h_sp = atomKK->k_sp.view_hostkk();
  fm = atomKK->fm;
  d_fm = atomKK->k_fm.view_device();
  h_fm = atomKK->k_fm.view_hostkk();
  fm_long = atomKK->fm_long;
  d_fm_long = atomKK->k_fm_long.view_device();
  h_fm_long = atomKK->k_fm_long.view_hostkk();
}

/* ----------------------------------------------------------------------
   sort atom arrays on device
------------------------------------------------------------------------- */

void AtomVecSpinKokkos::sort_kokkos(Kokkos::BinSort<KeyViewType, BinOp> &Sorter)
{
  atomKK->sync(Device, TAG_MASK|TYPE_MASK|MASK_MASK|IMAGE_MASK|X_MASK|V_MASK|SP_MASK);

  Sorter.sort(LMPDeviceType(), d_tag);
  Sorter.sort(LMPDeviceType(), d_type);
  Sorter.sort(LMPDeviceType(), d_mask);
  Sorter.sort(LMPDeviceType(), d_image);
  Sorter.sort(LMPDeviceType(), d_x);
  Sorter.sort(LMPDeviceType(), d_v);
  Sorter.sort(LMPDeviceType(), d_sp);

  atomKK->modified(Device, TAG_MASK|TYPE_MASK|MASK_MASK|IMAGE_MASK|X_MASK|V_MASK|SP_MASK);
}

/* ----------------------------------------------------------------------
   clear extra forces starting at atom N
   nbytes = # of bytes to clear for a per-atom vector
   include f b/c this is invoked from within SPIN pair styles
------------------------------------------------------------------------- */

void AtomVecSpinKokkos::force_clear(int /*n*/, size_t nbytes)
{
  int nzero = (double)nbytes/sizeof(double);

  if (nzero) {
    atomKK->k_fm.clear_sync_state(); // will be cleared below
    atomKK->k_fm_long.clear_sync_state(); // will be cleared below

    // local variables for lambda capture

    auto l_fm = atomKK->k_fm.view_device();
    auto l_fm_long = atomKK->k_fm_long.view_device();

    Kokkos::parallel_for(nzero, LAMMPS_LAMBDA(int i) {
      l_fm(i,0) = 0.0;
      l_fm(i,1) = 0.0;
      l_fm(i,2) = 0.0;
      l_fm_long(i,0) = 0.0;
      l_fm_long(i,1) = 0.0;
      l_fm_long(i,2) = 0.0;
    });

    atomKK->modified(Device,FM_MASK|FML_MASK);
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecSpinKokkos::sync(ExecutionSpace space, uint64_t mask)
{
  if (space == Device) {
    if (mask & X_MASK) atomKK->k_x.sync_device();
    if (mask & V_MASK) atomKK->k_v.sync_device();
    if (mask & F_MASK) atomKK->k_f.sync_device();
    if (mask & TAG_MASK) atomKK->k_tag.sync_device();
    if (mask & TYPE_MASK) atomKK->k_type.sync_device();
    if (mask & MASK_MASK) atomKK->k_mask.sync_device();
    if (mask & IMAGE_MASK) atomKK->k_image.sync_device();
    if (mask & SP_MASK) atomKK->k_sp.sync_device();
    if (mask & FM_MASK) atomKK->k_fm.sync_device();
    if (mask & FML_MASK) atomKK->k_fm_long.sync_device();
  } else if (space == Host) {
    if (mask & X_MASK) atomKK->k_x.sync_host();
    if (mask & V_MASK) atomKK->k_v.sync_host();
    if (mask & F_MASK) atomKK->k_f.sync_host();
    if (mask & TAG_MASK) atomKK->k_tag.sync_host();
    if (mask & TYPE_MASK) atomKK->k_type.sync_host();
    if (mask & MASK_MASK) atomKK->k_mask.sync_host();
    if (mask & IMAGE_MASK) atomKK->k_image.sync_host();
    if (mask & SP_MASK) atomKK->k_sp.sync_host();
    if (mask & FM_MASK) atomKK->k_fm.sync_host();
    if (mask & FML_MASK) atomKK->k_fm_long.sync_host();
  } else if (space == HostKK) {
    if (mask & X_MASK) atomKK->k_x.sync_hostkk();
    if (mask & V_MASK) atomKK->k_v.sync_hostkk();
    if (mask & F_MASK) atomKK->k_f.sync_hostkk();
    if (mask & TAG_MASK) atomKK->k_tag.sync_host();
    if (mask & TYPE_MASK) atomKK->k_type.sync_host();
    if (mask & MASK_MASK) atomKK->k_mask.sync_host();
    if (mask & IMAGE_MASK) atomKK->k_image.sync_host();
    if (mask & SP_MASK) atomKK->k_sp.sync_hostkk();
    if (mask & FM_MASK) atomKK->k_fm.sync_hostkk();
    if (mask & FML_MASK) atomKK->k_fm_long.sync_hostkk();
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecSpinKokkos::modified(ExecutionSpace space, uint64_t mask)
{
  if (space == Device) {
    if (mask & X_MASK) atomKK->k_x.modify_device();
    if (mask & V_MASK) atomKK->k_v.modify_device();
    if (mask & F_MASK) atomKK->k_f.modify_device();
    if (mask & TAG_MASK) atomKK->k_tag.modify_device();
    if (mask & TYPE_MASK) atomKK->k_type.modify_device();
    if (mask & MASK_MASK) atomKK->k_mask.modify_device();
    if (mask & IMAGE_MASK) atomKK->k_image.modify_device();
    if (mask & SP_MASK) atomKK->k_sp.modify_device();
    if (mask & FM_MASK) atomKK->k_fm.modify_device();
    if (mask & FML_MASK) atomKK->k_fm_long.modify_device();
  } else if (space == Host) {
    if (mask & X_MASK) atomKK->k_x.modify_host();
    if (mask & V_MASK) atomKK->k_v.modify_host();
    if (mask & F_MASK) atomKK->k_f.modify_host();
    if (mask & TAG_MASK) atomKK->k_tag.modify_host();
    if (mask & TYPE_MASK) atomKK->k_type.modify_host();
    if (mask & MASK_MASK) atomKK->k_mask.modify_host();
    if (mask & IMAGE_MASK) atomKK->k_image.modify_host();
    if (mask & SP_MASK) atomKK->k_sp.modify_host();
    if (mask & FM_MASK) atomKK->k_fm.modify_host();
    if (mask & FML_MASK) atomKK->k_fm_long.modify_host();
  } else if (space == HostKK) {
    if (mask & X_MASK) atomKK->k_x.modify_hostkk();
    if (mask & V_MASK) atomKK->k_v.modify_hostkk();
    if (mask & F_MASK) atomKK->k_f.modify_hostkk();
    if (mask & TAG_MASK) atomKK->k_tag.modify_host();
    if (mask & TYPE_MASK) atomKK->k_type.modify_host();
    if (mask & MASK_MASK) atomKK->k_mask.modify_host();
    if (mask & IMAGE_MASK) atomKK->k_image.modify_host();
    if (mask & SP_MASK) atomKK->k_sp.modify_hostkk();
    if (mask & FM_MASK) atomKK->k_fm.modify_hostkk();
    if (mask & FML_MASK) atomKK->k_fm_long.modify_hostkk();
  }
}

void AtomVecSpinKokkos::sync_pinned(ExecutionSpace space, uint64_t mask, int async_flag)
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
    if ((mask & SP_MASK) && atomKK->k_sp.need_sync_device())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d_4>(atomKK->k_sp,space,async_flag);
    if ((mask & FM_MASK) && atomKK->k_sp.need_sync_device())
      perform_pinned_copy_transform<DAT::ttransform_kkacc_1d_3>(atomKK->k_fm,space,async_flag);
    if ((mask & FML_MASK) && atomKK->k_fm_long.need_sync_device())
      perform_pinned_copy_transform<DAT::ttransform_kkacc_1d_3>(atomKK->k_fm_long,space,async_flag);
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
    if ((mask & SP_MASK) && atomKK->k_sp.need_sync_host())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d_4>(atomKK->k_sp,space,async_flag);
    if ((mask & FM_MASK) && atomKK->k_fm.need_sync_host())
      perform_pinned_copy_transform<DAT::ttransform_kkacc_1d_3>(atomKK->k_fm,space,async_flag);
    if ((mask & FML_MASK) && atomKK->k_fm_long.need_sync_host())
      perform_pinned_copy_transform<DAT::ttransform_kkacc_1d_3>(atomKK->k_fm_long,space,async_flag);
  }
}
