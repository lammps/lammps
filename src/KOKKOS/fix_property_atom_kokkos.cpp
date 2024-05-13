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
  for (int nv = 0; nv < nvalue; nv++)
    if (styles[nv] == DVEC) dvector_flag = 1;
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
  for (int nv = 0; nv < nvalue; nv++) {
    if (styles[nv] == MOLECULE) {
      atomKK->sync(Device,MOLECULE_MASK);
      memoryKK->grow_kokkos(atomKK->k_molecule,atom->molecule,nmax,"atom:molecule");
      atomKK->modified(Device,MOLECULE_MASK);
    } else if (styles[nv] == CHARGE) {
      atomKK->sync(Device,Q_MASK);
      memoryKK->grow_kokkos(atomKK->k_q,atom->q,nmax,"atom:q");
      atomKK->modified(Device,Q_MASK);
    } else if (styles[nv] == RMASS) {
      atomKK->sync(Device,RMASS_MASK);
      memoryKK->grow_kokkos(atomKK->k_rmass,atom->rmass,nmax,"atom:rmass");
      atomKK->modified(Device,RMASS_MASK);
    } else if (styles[nv] == TEMPERATURE) {
      memory->grow(atom->temperature, nmax, "atom:temperature");
      size_t nbytes = (nmax - nmax_old) * sizeof(double);
      memset(&atom->temperature[nmax_old], 0, nbytes);
    } else if (styles[nv] == HEATFLOW) {
      memory->grow(atom->heatflow, nmax, "atom:heatflow");
      size_t nbytes = (nmax - nmax_old) * sizeof(double);
      memset(&atom->heatflow[nmax_old], 0, nbytes);
    } else if (styles[nv] == IVEC) {
      memory->grow(atom->ivector[index[nv]],nmax,"atom:ivector");
      size_t nbytes = (nmax-nmax_old) * sizeof(int);
      memset(&atom->ivector[index[nv]][nmax_old],0,nbytes);
    } else if (styles[nv] == DVEC) {
      atomKK->sync(Device,DVECTOR_MASK);
      memoryKK->grow_kokkos(atomKK->k_dvector,atom->dvector,atomKK->k_dvector.extent(0),nmax,
                          "atom:dvector");
      atomKK->modified(Device,DVECTOR_MASK);
    } else if (styles[nv] == IARRAY) {
      memory->grow(atom->iarray[index[nv]], nmax, cols[nv], "atom:iarray");
      size_t nbytes = (size_t) (nmax - nmax_old) * cols[nv] * sizeof(int);
      if (nbytes) memset(&atom->iarray[index[nv]][nmax_old][0], 0, nbytes);
    } else if (styles[nv] == DARRAY) {
      memory->grow(atom->darray[index[nv]], nmax, cols[nv], "atom:darray");
      size_t nbytes = (size_t) (nmax - nmax_old) * cols[nv] * sizeof(double);
      if (nbytes) memset(&atom->darray[index[nv]][nmax_old][0], 0, nbytes);
    }
  }
  nmax_old = nmax;
}

/* ---------------------------------------------------------------------- */

void FixPropertyAtomKokkos::sync(ExecutionSpace space, unsigned int mask)
{
  if (space == Device) {
    if (molecule_flag && (mask & MOLECULE_MASK)) atomKK->k_molecule.sync<LMPDeviceType>();
    if (q_flag && (mask & Q_MASK)) atomKK->k_q.sync<LMPDeviceType>();
    if (rmass_flag && (mask & RMASS_MASK)) {atomKK->k_rmass.sync<LMPDeviceType>();}
    if (dvector_flag && (mask & DVECTOR_MASK)) atomKK->k_dvector.sync<LMPDeviceType>();
  } else {
    if (molecule_flag && (mask & MOLECULE_MASK)) atomKK->k_molecule.sync<LMPHostType>();
    if (q_flag && (mask & Q_MASK)) atomKK->k_q.sync<LMPHostType>();
    if (rmass_flag && (mask & RMASS_MASK)) atomKK->k_rmass.sync<LMPHostType>();
    if (dvector_flag && (mask & DVECTOR_MASK)) atomKK->k_dvector.sync<LMPHostType>();
  }
}

/* ---------------------------------------------------------------------- */

void FixPropertyAtomKokkos::sync_overlapping_device(ExecutionSpace space, unsigned int mask)
{
  if (space == Device) {
    if ((mask & MOLECULE_MASK) && atomKK->k_molecule.need_sync<LMPDeviceType>())
      atomKK->avecKK->perform_async_copy<DAT::tdual_tagint_1d>(atomKK->k_molecule,space);
    if ((mask & Q_MASK) && atomKK->k_q.need_sync<LMPDeviceType>())
      atomKK->avecKK->perform_async_copy<DAT::tdual_float_1d>(atomKK->k_q,space);
    if ((mask & RMASS_MASK) && atomKK->k_rmass.need_sync<LMPDeviceType>())
      atomKK->avecKK->perform_async_copy<DAT::tdual_float_1d>(atomKK->k_rmass,space);
    if ((mask & DVECTOR_MASK) && atomKK->k_dvector.need_sync<LMPDeviceType>())
      atomKK->avecKK->perform_async_copy<DAT::tdual_float_2d>(atomKK->k_dvector,space);
  } else {
    if ((mask & MOLECULE_MASK) && atomKK->k_molecule.need_sync<LMPHostType>())
      atomKK->avecKK->perform_async_copy<DAT::tdual_tagint_1d>(atomKK->k_molecule,space);
    if ((mask & Q_MASK) && atomKK->k_q.need_sync<LMPHostType>())
      atomKK->avecKK->perform_async_copy<DAT::tdual_float_1d>(atomKK->k_q,space);
    if ((mask & RMASS_MASK) && atomKK->k_rmass.need_sync<LMPHostType>())
      atomKK->avecKK->perform_async_copy<DAT::tdual_float_1d>(atomKK->k_rmass,space);
    if ((mask & DVECTOR_MASK) && atomKK->k_dvector.need_sync<LMPHostType>())
      atomKK->avecKK->perform_async_copy<DAT::tdual_float_2d>(atomKK->k_dvector,space);
  }
}

/* ---------------------------------------------------------------------- */

void FixPropertyAtomKokkos::modified(ExecutionSpace space, unsigned int mask)
{
  if (space == Device) {
    if (molecule_flag && (mask & MOLECULE_MASK)) atomKK->k_molecule.modify<LMPDeviceType>();
    if (q_flag && (mask & Q_MASK)) atomKK->k_q.modify<LMPDeviceType>();
    if (rmass_flag && (mask & RMASS_MASK)) atomKK->k_rmass.modify<LMPDeviceType>();
    if (dvector_flag && (mask & DVECTOR_MASK)) atomKK->k_dvector.modify<LMPDeviceType>();
  } else {
    if (molecule_flag && (mask & MOLECULE_MASK)) atomKK->k_molecule.modify<LMPHostType>();
    if (q_flag && (mask & Q_MASK)) atomKK->k_q.modify<LMPHostType>();
    if (rmass_flag && (mask & RMASS_MASK)) atomKK->k_rmass.modify<LMPHostType>();
    if (dvector_flag && (mask & DVECTOR_MASK)) atomKK->k_dvector.modify<LMPHostType>();
  }
}
