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

#include "fix_property_atom_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "memory_kokkos.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixPropertyAtomKokkos::FixPropertyAtomKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixPropertyAtom(lmp, narg, arg)
{
  atomKK = (AtomKokkos *) atom;
  kokkosable = 1;

  dvector_flag = 0;
  ivector_flag = 0;
  iarray_flag = 0;
  darray_flag = 0;
  for (int nv = 0; nv < nvalue; nv++) {
    if (styles[nv] == IVEC) ivector_flag = 1;
    if (styles[nv] == DVEC) dvector_flag = 1;
    if (styles[nv] == IARRAY) iarray_flag = 1;
    if (styles[nv] == DARRAY) darray_flag = 1;
  }
}

/* ---------------------------------------------------------------------- */

void FixPropertyAtomKokkos::post_constructor()
{
  atomKK->update_property_atom();

  FixPropertyAtom::post_constructor();
}

/* ---------------------------------------------------------------------- */

FixPropertyAtomKokkos::~FixPropertyAtomKokkos()
{
  // deallocate per-atom vectors in Atom class
  // set ptrs to a null pointer, so they no longer exist for Atom class

  for (int nv = 0; nv < nvalue; nv++) {
    if (styles[nv] == MOLECULE) {
      atom->molecule_flag = 0;
      memoryKK->destroy_kokkos(atomKK->k_molecule,atom->molecule);
      atom->molecule = nullptr;
    } else if (styles[nv] == CHARGE) {
      atom->q_flag = 0;
      memoryKK->destroy_kokkos(atomKK->k_q,atom->q);
      atom->q = nullptr;
    } else if (styles[nv] == RMASS) {
      atom->rmass_flag = 0;
      memoryKK->destroy_kokkos(atomKK->k_rmass,atom->rmass);
      atom->rmass = nullptr;
    }
  }

  atomKK->update_property_atom();
}

/* ----------------------------------------------------------------------
   allocate atom-based arrays
   initialize new values to 0,
   since AtomVec class won't do it as atoms are added,
   e.g. in create_atom() or data_atom()
------------------------------------------------------------------------- */

void FixPropertyAtomKokkos::grow_arrays(int nmax)
{
  // ivector is one shared contiguous 2D view holding all custom int vectors,
  // so grow it once here rather than re-growing the whole view once per IVEC
  // value inside the loop below
  if (ivector_flag) {
    atomKK->sync(Device,IVECTOR_MASK);
    atomKK->modified(Device,IVECTOR_MASK);
    memoryKK->grow_kokkos(atomKK->k_ivector,atom->ivector,atomKK->k_ivector.extent(0),nmax,
                        "atom:ivector");
    atomKK->sync(Host,IVECTOR_MASK);
  }

  for (int nv = 0; nv < nvalue; nv++) {
    if (styles[nv] == MOLECULE) {
      atomKK->sync(Device,MOLECULE_MASK);
      atomKK->modified(Device,MOLECULE_MASK);
      memoryKK->grow_kokkos(atomKK->k_molecule,atom->molecule,nmax,"atom:molecule");
      atomKK->sync(Host,MOLECULE_MASK);
    } else if (styles[nv] == CHARGE) {
      atomKK->sync(Device,Q_MASK);
      atomKK->modified(Device,Q_MASK);
      memoryKK->grow_kokkos(atomKK->k_q,atom->q,nmax,"atom:q");
      atomKK->sync(Host,Q_MASK);
    } else if (styles[nv] == RMASS) {
      atomKK->sync(Device,RMASS_MASK);
      atomKK->modified(Device,RMASS_MASK);
      memoryKK->grow_kokkos(atomKK->k_rmass,atom->rmass,nmax,"atom:rmass");
      atomKK->sync(Host,RMASS_MASK);
    } else if (styles[nv] == TEMPERATURE) {
      memory->grow(atom->temperature, nmax, "atom:temperature");
      size_t nbytes = (nmax - nmax_old) * sizeof(double);
      memset(&atom->temperature[nmax_old], 0, nbytes);
    } else if (styles[nv] == HEATFLOW) {
      memory->grow(atom->heatflow, nmax, "atom:heatflow");
      size_t nbytes = (nmax - nmax_old) * sizeof(double);
      memset(&atom->heatflow[nmax_old], 0, nbytes);
    } else if (styles[nv] == IVEC) {
      // storage was grown once above (shared contiguous k_ivector); just zero
      // this vector's newly added tail
      size_t nbytes = (nmax-nmax_old) * sizeof(int);
      memset(&atom->ivector[index[nv]][nmax_old],0,nbytes);
    } else if (styles[nv] == DVEC) {
      atomKK->sync(Device,DVECTOR_MASK);
      atomKK->modified(Device,DVECTOR_MASK);
      memoryKK->grow_kokkos(atomKK->k_dvector,atom->dvector,atomKK->k_dvector.extent(0),nmax,
                          "atom:dvector");
      atomKK->sync(Host,DVECTOR_MASK);
    } else if (styles[nv] == IARRAY) {
      // ragged cols: grow this property's inner DualView in the view-of-views,
      // using the same device/host sync handshake as the DVEC branch so a grow
      // that follows a device-side modification does not lose data
      int idx = index[nv];
      auto& inner = atomKK->k_iarray.view_host()[idx].k_view;
      inner.sync_device();
      inner.modify_device();
      memoryKK->grow_kokkos(inner, atom->iarray[idx], nmax, cols[nv], "atom:iarray");
      inner.sync_host();
      size_t nbytes = (size_t) (nmax - nmax_old) * cols[nv] * sizeof(int);
      if (nbytes) {
        memset(&atom->iarray[idx][nmax_old][0], 0, nbytes);
        inner.modify_host();   // the memset wrote host; mark it so it reaches device
      }
      // re-sync the outer struct array: growing the inner view changed its
      // device pointer, which the device-side view-of-views must see
      atomKK->k_iarray.modify_host();
      atomKK->k_iarray.sync_device();
    } else if (styles[nv] == DARRAY) {
      int idx = index[nv];
      auto& inner = atomKK->k_darray.view_host()[idx].k_view;
      inner.sync_device();
      inner.modify_device();
      memoryKK->grow_kokkos(inner, atom->darray[idx], nmax, cols[nv], "atom:darray");
      inner.sync_host();
      size_t nbytes = (size_t) (nmax - nmax_old) * cols[nv] * sizeof(double);
      if (nbytes) {
        memset(&atom->darray[idx][nmax_old][0], 0, nbytes);
        inner.modify_host();
      }
      atomKK->k_darray.modify_host();
      atomKK->k_darray.sync_device();
    }
  }
  nmax_old = nmax;
}

/* ---------------------------------------------------------------------- */

void FixPropertyAtomKokkos::sync(ExecutionSpace space, uint64_t mask)
{
  if (space == Device) {
    if (molecule_flag && (mask & MOLECULE_MASK)) atomKK->k_molecule.sync_device();
    if (q_flag && (mask & Q_MASK)) atomKK->k_q.sync_device();
    if (rmass_flag && (mask & RMASS_MASK)) {atomKK->k_rmass.sync_device();}
    if (dvector_flag && (mask & DVECTOR_MASK)) atomKK->k_dvector.sync_device();
    if (ivector_flag && (mask & IVECTOR_MASK)) atomKK->k_ivector.sync_device();
    if (iarray_flag && (mask & IARRAY_MASK))
      for (int i = 0; i < atom->niarray; i++)
        atomKK->k_iarray.view_host()[i].k_view.sync_device();
    if (darray_flag && (mask & DARRAY_MASK))
      for (int i = 0; i < atom->ndarray; i++)
        atomKK->k_darray.view_host()[i].k_view.sync_device();
  } else if (space == Host) {
    if (molecule_flag && (mask & MOLECULE_MASK)) atomKK->k_molecule.sync_host();
    if (q_flag && (mask & Q_MASK)) atomKK->k_q.sync_host();
    if (rmass_flag && (mask & RMASS_MASK)) atomKK->k_rmass.sync_host();
    if (dvector_flag && (mask & DVECTOR_MASK)) atomKK->k_dvector.sync_host();
    if (ivector_flag && (mask & IVECTOR_MASK)) atomKK->k_ivector.sync_host();
    if (iarray_flag && (mask & IARRAY_MASK))
      for (int i = 0; i < atom->niarray; i++)
        atomKK->k_iarray.view_host()[i].k_view.sync_host();
    if (darray_flag && (mask & DARRAY_MASK))
      for (int i = 0; i < atom->ndarray; i++)
        atomKK->k_darray.view_host()[i].k_view.sync_host();
  } else if (space == HostKK) {
    if (molecule_flag && (mask & MOLECULE_MASK)) atomKK->k_molecule.sync_host();
    if (q_flag && (mask & Q_MASK)) atomKK->k_q.sync_hostkk();
    if (rmass_flag && (mask & RMASS_MASK)) atomKK->k_rmass.sync_hostkk();
    if (dvector_flag && (mask & DVECTOR_MASK)) atomKK->k_dvector.sync_hostkk();
    if (ivector_flag && (mask & IVECTOR_MASK)) atomKK->k_ivector.sync_host();
    if (iarray_flag && (mask & IARRAY_MASK))
      for (int i = 0; i < atom->niarray; i++)
        atomKK->k_iarray.view_host()[i].k_view.sync_host();
    if (darray_flag && (mask & DARRAY_MASK))
      for (int i = 0; i < atom->ndarray; i++)
        atomKK->k_darray.view_host()[i].k_view.sync_host();
  }
}

/* ---------------------------------------------------------------------- */

void FixPropertyAtomKokkos::sync_pinned(ExecutionSpace space, uint64_t mask, int async_flag)
{
  if (space == Device) {
    if ((mask & MOLECULE_MASK) && atomKK->k_molecule.need_sync_device())
      atomKK->avecKK->perform_pinned_copy<DAT::tdual_tagint_1d>(atomKK->k_molecule,space,async_flag);
    if ((mask & Q_MASK) && atomKK->k_q.need_sync_device())
      atomKK->avecKK->perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d>(atomKK->k_q,space,async_flag);
    if ((mask & RMASS_MASK) && atomKK->k_rmass.need_sync_device())
      atomKK->avecKK->perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d>(atomKK->k_rmass,space,async_flag);
    if ((mask & DVECTOR_MASK) && atomKK->k_dvector.need_sync_device())
      atomKK->avecKK->perform_pinned_copy_transform<DAT::ttransform_kkfloat_2d>(atomKK->k_dvector,space,async_flag);
    if ((mask & IVECTOR_MASK) && atomKK->k_ivector.need_sync_device())
      atomKK->avecKK->perform_pinned_copy<DAT::tdual_int_2d_lr>(atomKK->k_ivector,space,async_flag);
    if (mask & IARRAY_MASK)
      for (int i = 0; i < atom->niarray; i++)
        if (atomKK->k_iarray.view_host()[i].k_view.need_sync_device())
          atomKK->avecKK->perform_pinned_copy<DAT::tdual_int_2d_lr>(
              atomKK->k_iarray.view_host()[i].k_view, space, async_flag);
    if (mask & DARRAY_MASK)
      for (int i = 0; i < atom->ndarray; i++)
        if (atomKK->k_darray.view_host()[i].k_view.need_sync_device())
          atomKK->avecKK->perform_pinned_copy<DAT::tdual_double_2d_lr>(
              atomKK->k_darray.view_host()[i].k_view, space, async_flag);
  } else {
    if ((mask & MOLECULE_MASK) && atomKK->k_molecule.need_sync_host())
      atomKK->avecKK->perform_pinned_copy<DAT::tdual_tagint_1d>(atomKK->k_molecule,space,async_flag);
    if ((mask & Q_MASK) && atomKK->k_q.need_sync_host())
      atomKK->avecKK->perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d>(atomKK->k_q,space,async_flag);
    if ((mask & RMASS_MASK) && atomKK->k_rmass.need_sync_host())
      atomKK->avecKK->perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d>(atomKK->k_rmass,space,async_flag);
    if ((mask & DVECTOR_MASK) && atomKK->k_dvector.need_sync_host())
      atomKK->avecKK->perform_pinned_copy_transform<DAT::ttransform_kkfloat_2d>(atomKK->k_dvector,space,async_flag);
    if ((mask & IVECTOR_MASK) && atomKK->k_ivector.need_sync_host())
      atomKK->avecKK->perform_pinned_copy<DAT::tdual_int_2d_lr>(atomKK->k_ivector,space,async_flag);
    if (mask & IARRAY_MASK)
      for (int i = 0; i < atom->niarray; i++)
        if (atomKK->k_iarray.view_host()[i].k_view.need_sync_host())
          atomKK->avecKK->perform_pinned_copy<DAT::tdual_int_2d_lr>(
              atomKK->k_iarray.view_host()[i].k_view, space, async_flag);
    if (mask & DARRAY_MASK)
      for (int i = 0; i < atom->ndarray; i++)
        if (atomKK->k_darray.view_host()[i].k_view.need_sync_host())
          atomKK->avecKK->perform_pinned_copy<DAT::tdual_double_2d_lr>(
              atomKK->k_darray.view_host()[i].k_view, space, async_flag);
  }
}

/* ---------------------------------------------------------------------- */

void FixPropertyAtomKokkos::modified(ExecutionSpace space, uint64_t mask)
{
  if (space == Device) {
    if (molecule_flag && (mask & MOLECULE_MASK)) atomKK->k_molecule.modify_device();
    if (q_flag && (mask & Q_MASK)) atomKK->k_q.modify_device();
    if (rmass_flag && (mask & RMASS_MASK)) atomKK->k_rmass.modify_device();
    if (dvector_flag && (mask & DVECTOR_MASK)) atomKK->k_dvector.modify_device();
    if (ivector_flag && (mask & IVECTOR_MASK)) atomKK->k_ivector.modify_device();
    if (iarray_flag && (mask & IARRAY_MASK))
      for (int i = 0; i < atom->niarray; i++)
        atomKK->k_iarray.view_host()[i].k_view.modify_device();
    if (darray_flag && (mask & DARRAY_MASK))
      for (int i = 0; i < atom->ndarray; i++)
        atomKK->k_darray.view_host()[i].k_view.modify_device();
  } else if (space == Host) {
    if (molecule_flag && (mask & MOLECULE_MASK)) atomKK->k_molecule.modify_host();
    if (q_flag && (mask & Q_MASK)) atomKK->k_q.modify_host();
    if (rmass_flag && (mask & RMASS_MASK)) atomKK->k_rmass.modify_host();
    if (dvector_flag && (mask & DVECTOR_MASK)) atomKK->k_dvector.modify_host();
    if (ivector_flag && (mask & IVECTOR_MASK)) atomKK->k_ivector.modify_host();
    if (iarray_flag && (mask & IARRAY_MASK))
      for (int i = 0; i < atom->niarray; i++)
        atomKK->k_iarray.view_host()[i].k_view.modify_host();
    if (darray_flag && (mask & DARRAY_MASK))
      for (int i = 0; i < atom->ndarray; i++)
        atomKK->k_darray.view_host()[i].k_view.modify_host();
  } else if (space == HostKK) {
    if (molecule_flag && (mask & MOLECULE_MASK)) atomKK->k_molecule.modify_host();
    if (q_flag && (mask & Q_MASK)) atomKK->k_q.modify_hostkk();
    if (rmass_flag && (mask & RMASS_MASK)) atomKK->k_rmass.modify_hostkk();
    if (dvector_flag && (mask & DVECTOR_MASK)) atomKK->k_dvector.modify_hostkk();
    if (ivector_flag && (mask & IVECTOR_MASK)) atomKK->k_ivector.modify_host();
    if (iarray_flag && (mask & IARRAY_MASK))
      for (int i = 0; i < atom->niarray; i++)
        atomKK->k_iarray.view_host()[i].k_view.modify_host();
    if (darray_flag && (mask & DARRAY_MASK))
      for (int i = 0; i < atom->ndarray; i++)
        atomKK->k_darray.view_host()[i].k_view.modify_host();
  }
}
