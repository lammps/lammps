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

#include "compute_species_atom_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "pair_reaxff_kokkos.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ComputeSpeciesAtomKokkos<DeviceType>::ComputeSpeciesAtomKokkos(LAMMPS *lmp, int narg, char **arg) :
  ComputeSpeciesAtom(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  
  auto reaxff_kk = dynamic_cast<PairReaxFFKokkos<DeviceType> *>(force->pair_match("^reax/kk",0));
  if (reaxff_kk)
    d_tmpbo = reaxff_kk->k_tmpbo.template view<DeviceType>();
  else
    error->all(FLERR,"Cannot use compute species/atom/kk without Kokkos ReaxFF");

  if (nvalues == 1) memoryKK->create_kokkos(k_vector,vector_atom,nmax,"species/atom:vector");
  else memoryKK->create_kokkos(k_array,array_atom,nmax,nvalues,"species/atom:array");
  memoryKK->create_kokkos(k_properties,nvalues,"species/atom:k_properties");

  for (int i = 0; i < nvalues; i++)
    k_properties.h_view(i) = get_property_type(pack_choice[i]);

  k_properties.modify<LMPHostType>();
  k_properties.sync<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ComputeSpeciesAtomKokkos<DeviceType>::~ComputeSpeciesAtomKokkos()
{
  if (copymode) return;
  memoryKK->destroy_kokkos(k_vector);
  memoryKK->destroy_kokkos(k_array);
  memoryKK->destroy_kokkos(k_properties);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeSpeciesAtomKokkos<DeviceType>::compute_peratom()
{
  invoked_peratom = update->ntimestep;
  
  // Resize Kokkos views if necessary
  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    if (nvalues == 1) memoryKK->grow_kokkos(k_vector,vector_atom,nmax,"species/atom:vector");
    else memoryKK->grow_kokkos(k_array,array_atom,nmax,nvalues,"species/atom:array");
  }

  // Get necessary data from atom class
  atomKK->sync(execution_space, X_MASK | V_MASK | Q_MASK | MASK_MASK);
  auto x = atomKK->k_x.template view<DeviceType>();
  auto v = atomKK->k_v.template view<DeviceType>();
  auto q = atomKK->k_q.template view<DeviceType>();
  auto mask = atomKK->k_mask.template view<DeviceType>();

  // Get ReaxFF bond order data
  auto tmpbo = d_tmpbo;
  
  // Local copies for lambda capture
  int groupbit_local = groupbit;
  int nvalues_local = nvalues;
  auto d_properties = k_properties.template view<DeviceType>();

  // Mark that we'll be modifying the device side of our DualView
  if (nvalues == 1) {
    k_vector.template modify<DeviceType>();
    auto d_vec = k_vector.d_view;
    
    // Initialize the vector to zeros
    Kokkos::deep_copy(d_vec, 0.0);
    
    // Vector version (single property)
    unsigned int type = d_properties(0);
    
    Kokkos::parallel_for("ComputeSpeciesAtomKokkos:vector",
        Kokkos::RangePolicy<DeviceType>(0,atom->nlocal),
        KOKKOS_LAMBDA(const int i) {
          if (mask(i) & groupbit_local) {
            // Switch based on property type
            switch (type) {
              case 1: d_vec(i) = q(i);   break;
              case 2: d_vec(i) = x(i,0); break;
              case 3: d_vec(i) = x(i,1); break;
              case 4: d_vec(i) = x(i,2); break;
              case 5: d_vec(i) = v(i,0); break;
              case 6: d_vec(i) = v(i,1); break;
              case 7: d_vec(i) = v(i,2); break;
              default:
                // bond orders (abo01-abo24)
                if (type >= 8 && type <= 31) d_vec(i) = tmpbo(i, type-8);
                break;
            }
          }
        });
  } else {
    k_array.template modify<DeviceType>();
    auto d_arr = k_array.d_view;
    
    // Initialize the array to zeros
    Kokkos::deep_copy(d_arr, 0.0);
    
    // Array version (multiple properties)
    Kokkos::parallel_for("ComputeSpeciesAtomKokkos:array",
        Kokkos::RangePolicy<DeviceType>(0,atom->nlocal),
        KOKKOS_LAMBDA(const int i) {
          if (mask(i) & groupbit_local) {
            for (int n = 0; n < nvalues_local; n++) {
              unsigned int type = d_properties(n);
              
              // Switch based on property type
              switch (type) {
                case 1: d_arr(i,n) = q(i);   break;
                case 2: d_arr(i,n) = x(i,0); break;
                case 3: d_arr(i,n) = x(i,1); break;
                case 4: d_arr(i,n) = x(i,2); break;
                case 5: d_arr(i,n) = v(i,0); break;
                case 6: d_arr(i,n) = v(i,1); break;
                case 7: d_arr(i,n) = v(i,2); break;
                default:
                  // bond orders (abo01-abo24)
                  if (type >= 8 && type <= 31) d_arr(i,n) = tmpbo(i, type-8);
                  break;
              }
            }
          }
        });
  }
  
  // Sync data back to host views
  if (nvalues == 1) k_vector.template sync<LMPHostType>();
  else k_array.template sync<LMPHostType>();

}

// Explicit template instantiations
template class ComputeSpeciesAtomKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class ComputeSpeciesAtomKokkos<LMPHostType>;
#endif
