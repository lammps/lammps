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

#include "atom_vec_bond_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "memory_kokkos.h"
#include "modify.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecBondKokkos::AtomVecBondKokkos(LAMMPS *lmp) : AtomVec(lmp),
AtomVecKokkos(lmp), AtomVecBond(lmp)
{
}

/* ----------------------------------------------------------------------
   process field strings to initialize data structs for all other methods
------------------------------------------------------------------------- */

void AtomVecBondKokkos::init()
{
  AtomVecBond::init();

  set_atom_masks();
}

/* ----------------------------------------------------------------------
   grow atom arrays
   n = 0 grows arrays by DELTA
   n > 0 allocates arrays to size n
------------------------------------------------------------------------- */

void AtomVecBondKokkos::grow(int n)
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

  memoryKK->grow_kokkos(atomKK->k_molecule,atomKK->molecule,nmax,"atom:molecule");
  memoryKK->grow_kokkos(atomKK->k_nspecial,atomKK->nspecial,nmax,3,"atom:nspecial");
  memoryKK->grow_kokkos(atomKK->k_special,atomKK->special,nmax,atomKK->maxspecial,"atom:special");
  memoryKK->grow_kokkos(atomKK->k_num_bond,atomKK->num_bond,nmax,"atom:num_bond");
  memoryKK->grow_kokkos(atomKK->k_bond_type,atomKK->bond_type,nmax,atomKK->bond_per_atom,"atom:bond_type");
  memoryKK->grow_kokkos(atomKK->k_bond_atom,atomKK->bond_atom,nmax,atomKK->bond_per_atom,"atom:bond_atom");

  grow_pointers();
  atomKK->sync(Host,ALL_MASK);

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atomKK->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->grow_arrays(nmax);
}

/* ----------------------------------------------------------------------
   reset local array ptrs
------------------------------------------------------------------------- */

void AtomVecBondKokkos::grow_pointers()
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

  molecule = atomKK->molecule;
  d_molecule = atomKK->k_molecule.view_device();
  h_molecule = atomKK->k_molecule.view_host();
  nspecial = atomKK->nspecial;
  d_nspecial = atomKK->k_nspecial.view_device();
  h_nspecial = atomKK->k_nspecial.view_hostkk();
  special = atomKK->special;
  d_special = atomKK->k_special.view_device();
  h_special = atomKK->k_special.view_hostkk();
  num_bond = atomKK->num_bond;
  d_num_bond = atomKK->k_num_bond.view_device();
  h_num_bond = atomKK->k_num_bond.view_host();
  bond_type = atomKK->bond_type;
  d_bond_type = atomKK->k_bond_type.view_device();
  h_bond_type = atomKK->k_bond_type.view_hostkk();
  bond_atom = atomKK->bond_atom;
  d_bond_atom = atomKK->k_bond_atom.view_device();
  h_bond_atom = atomKK->k_bond_atom.view_hostkk();
}

/* ----------------------------------------------------------------------
   sort atom arrays on device
------------------------------------------------------------------------- */

void AtomVecBondKokkos::sort_kokkos(Kokkos::BinSort<KeyViewType, BinOp> &Sorter)
{
  atomKK->sync(Device, ALL_MASK & ~F_MASK);

  Sorter.sort(LMPDeviceType(), d_tag);
  Sorter.sort(LMPDeviceType(), d_type);
  Sorter.sort(LMPDeviceType(), d_mask);
  Sorter.sort(LMPDeviceType(), d_image);
  Sorter.sort(LMPDeviceType(), d_x);
  Sorter.sort(LMPDeviceType(), d_v);
  Sorter.sort(LMPDeviceType(), d_molecule);
  Sorter.sort(LMPDeviceType(), d_num_bond);
  Sorter.sort(LMPDeviceType(), d_bond_type);
  Sorter.sort(LMPDeviceType(), d_bond_atom);
  Sorter.sort(LMPDeviceType(), d_nspecial);
  Sorter.sort(LMPDeviceType(), d_special);

  atomKK->modified(Device, ALL_MASK & ~F_MASK);
}

/* ---------------------------------------------------------------------- */

void AtomVecBondKokkos::sync(ExecutionSpace space, uint64_t mask)
{
  if (space == Device) {
    if (mask & X_MASK) atomKK->k_x.sync_device();
    if (mask & V_MASK) atomKK->k_v.sync_device();
    if (mask & F_MASK) atomKK->k_f.sync_device();
    if (mask & TAG_MASK) atomKK->k_tag.sync_device();
    if (mask & TYPE_MASK) atomKK->k_type.sync_device();
    if (mask & MASK_MASK) atomKK->k_mask.sync_device();
    if (mask & IMAGE_MASK) atomKK->k_image.sync_device();
    if (mask & MOLECULE_MASK) atomKK->k_molecule.sync_device();
    if (mask & SPECIAL_MASK) {
      atomKK->k_nspecial.sync_device();
      atomKK->k_special.sync_device();
    }
    if (mask & BOND_MASK) {
      atomKK->k_num_bond.sync_device();
      atomKK->k_bond_type.sync_device();
      atomKK->k_bond_atom.sync_device();
    }
  } else if (space == Host) {
    if (mask & X_MASK) atomKK->k_x.sync_host();
    if (mask & V_MASK) atomKK->k_v.sync_host();
    if (mask & F_MASK) atomKK->k_f.sync_host();
    if (mask & TAG_MASK) atomKK->k_tag.sync_host();
    if (mask & TYPE_MASK) atomKK->k_type.sync_host();
    if (mask & MASK_MASK) atomKK->k_mask.sync_host();
    if (mask & IMAGE_MASK) atomKK->k_image.sync_host();
    if (mask & MOLECULE_MASK) atomKK->k_molecule.sync_host();
    if (mask & SPECIAL_MASK) {
      atomKK->k_nspecial.sync_host();
      atomKK->k_special.sync_host();
    }
    if (mask & BOND_MASK) {
      atomKK->k_num_bond.sync_host();
      atomKK->k_bond_type.sync_host();
      atomKK->k_bond_atom.sync_host();
    }
  } else if (space == HostKK) {
    if (mask & X_MASK) atomKK->k_x.sync_hostkk();
    if (mask & V_MASK) atomKK->k_v.sync_hostkk();
    if (mask & F_MASK) atomKK->k_f.sync_hostkk();
    if (mask & TAG_MASK) atomKK->k_tag.sync_host();
    if (mask & TYPE_MASK) atomKK->k_type.sync_host();
    if (mask & MASK_MASK) atomKK->k_mask.sync_host();
    if (mask & IMAGE_MASK) atomKK->k_image.sync_host();
    if (mask & MOLECULE_MASK) atomKK->k_molecule.sync_host();
    if (mask & SPECIAL_MASK) {
      atomKK->k_nspecial.sync_hostkk();
      atomKK->k_special.sync_hostkk();
    }
    if (mask & BOND_MASK) {
      atomKK->k_num_bond.sync_host();
      atomKK->k_bond_type.sync_hostkk();
      atomKK->k_bond_atom.sync_hostkk();
    }
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecBondKokkos::sync_pinned(ExecutionSpace space, uint64_t mask, int async_flag)
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
    if ((mask & MOLECULE_MASK) && atomKK->k_molecule.need_sync_device())
      perform_pinned_copy<DAT::tdual_tagint_1d>(atomKK->k_molecule,space,async_flag);
    if (mask & SPECIAL_MASK) {
      if (atomKK->k_nspecial.need_sync_device())
        perform_pinned_copy_transform<DAT::ttransform_int_2d>(atomKK->k_nspecial,space,async_flag);
      if (atomKK->k_special.need_sync_device())
        perform_pinned_copy_transform<DAT::ttransform_tagint_2d>(atomKK->k_special,space,async_flag);
    }
    if (mask & BOND_MASK) {
      if (atomKK->k_num_bond.need_sync_device())
        perform_pinned_copy<DAT::tdual_int_1d>(atomKK->k_num_bond,space,async_flag);
      if (atomKK->k_bond_type.need_sync_device())
        perform_pinned_copy_transform<DAT::ttransform_int_2d>(atomKK->k_bond_type,space,async_flag);
      if (atomKK->k_bond_atom.need_sync_device())
        perform_pinned_copy_transform<DAT::ttransform_tagint_2d>(atomKK->k_bond_atom,space,async_flag);
    }
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
    if ((mask & MOLECULE_MASK) && atomKK->k_molecule.need_sync_host())
      perform_pinned_copy<DAT::tdual_tagint_1d>(atomKK->k_molecule,space,async_flag);
    if (mask & SPECIAL_MASK) {
      if (atomKK->k_nspecial.need_sync_host())
        perform_pinned_copy_transform<DAT::ttransform_int_2d>(atomKK->k_nspecial,space,async_flag);
      if (atomKK->k_special.need_sync_host())
        perform_pinned_copy_transform<DAT::ttransform_tagint_2d>(atomKK->k_special,space,async_flag);
    }
    if (mask & BOND_MASK) {
      if (atomKK->k_num_bond.need_sync_host())
        perform_pinned_copy<DAT::tdual_int_1d>(atomKK->k_num_bond,space,async_flag);
      if (atomKK->k_bond_type.need_sync_host())
        perform_pinned_copy_transform<DAT::ttransform_int_2d>(atomKK->k_bond_type,space,async_flag);
      if (atomKK->k_bond_atom.need_sync_host())
        perform_pinned_copy_transform<DAT::ttransform_tagint_2d>(atomKK->k_bond_atom,space,async_flag);
    }
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecBondKokkos::modified(ExecutionSpace space, uint64_t mask)
{
  if (space == Device) {
    if (mask & X_MASK) atomKK->k_x.modify_device();
    if (mask & V_MASK) atomKK->k_v.modify_device();
    if (mask & F_MASK) atomKK->k_f.modify_device();
    if (mask & TAG_MASK) atomKK->k_tag.modify_device();
    if (mask & TYPE_MASK) atomKK->k_type.modify_device();
    if (mask & MASK_MASK) atomKK->k_mask.modify_device();
    if (mask & IMAGE_MASK) atomKK->k_image.modify_device();
    if (mask & MOLECULE_MASK) atomKK->k_molecule.modify_device();
    if (mask & SPECIAL_MASK) {
      atomKK->k_nspecial.modify_device();
      atomKK->k_special.modify_device();
    }
    if (mask & BOND_MASK) {
      atomKK->k_num_bond.modify_device();
      atomKK->k_bond_type.modify_device();
      atomKK->k_bond_atom.modify_device();
    }
  } else if (space == Host) {
    if (mask & X_MASK) atomKK->k_x.modify_host();
    if (mask & V_MASK) atomKK->k_v.modify_host();
    if (mask & F_MASK) atomKK->k_f.modify_host();
    if (mask & TAG_MASK) atomKK->k_tag.modify_host();
    if (mask & TYPE_MASK) atomKK->k_type.modify_host();
    if (mask & MASK_MASK) atomKK->k_mask.modify_host();
    if (mask & IMAGE_MASK) atomKK->k_image.modify_host();
    if (mask & MOLECULE_MASK) atomKK->k_molecule.modify_host();
    if (mask & SPECIAL_MASK) {
      atomKK->k_nspecial.modify_host();
      atomKK->k_special.modify_host();
    }
    if (mask & BOND_MASK) {
      atomKK->k_num_bond.modify_host();
      atomKK->k_bond_type.modify_host();
      atomKK->k_bond_atom.modify_host();
    }
  } else if (space == HostKK) {
    if (mask & X_MASK) atomKK->k_x.modify_hostkk();
    if (mask & V_MASK) atomKK->k_v.modify_hostkk();
    if (mask & F_MASK) atomKK->k_f.modify_hostkk();
    if (mask & TAG_MASK) atomKK->k_tag.modify_host();
    if (mask & TYPE_MASK) atomKK->k_type.modify_host();
    if (mask & MASK_MASK) atomKK->k_mask.modify_host();
    if (mask & IMAGE_MASK) atomKK->k_image.modify_host();
    if (mask & MOLECULE_MASK) atomKK->k_molecule.modify_host();
    if (mask & SPECIAL_MASK) {
      atomKK->k_nspecial.modify_hostkk();
      atomKK->k_special.modify_hostkk();
    }
    if (mask & BOND_MASK) {
      atomKK->k_num_bond.modify_host();
      atomKK->k_bond_type.modify_hostkk();
      atomKK->k_bond_atom.modify_hostkk();
    }
  }
}
