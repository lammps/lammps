/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Matt Bettencourt (NVIDIA)
------------------------------------------------------------------------- */

#ifndef LMP_MLIAP_DATA_KOKKOS_H
#define LMP_MLIAP_DATA_KOKKOS_H

#include "mliap_data.h"

#include "kokkos_type.h"
#include "memory_kokkos.h"
#include "pair_mliap_kokkos.h"
#include "pointers.h"

namespace LAMMPS_NS {
// clang-format off
enum {
  IATOMS_MASK      = 0x00000001,
  IELEMS_MASK      = 0x00000002,
  JATOMS_MASK      = 0x00000004,
  JELEMS_MASK      = 0x00000008,
  IJ_MASK          = 0x00000010,
  BETAS_MASK       = 0x00000020,
  DESCRIPTORS_MASK = 0x00000040,
  EATOMS_MASK      = 0x00000080,
  RIJ_MASK         = 0x00000100,
  GRADFORCE_MASK   = 0x00000200,
  GRADDESC_MASK    = 0x00000400,
  NUMNEIGHS_MASK   = 0x00000800,
  GAMMA_MASK_MASK  = 0x00001000,
  GAMMA_ROW_MASK   = 0x00002000,
  GAMMA_COL_MASK   = 0x00004000,
};
// clang-format on

template <class DeviceType> class MLIAPDataKokkos : public MLIAPData {
 public:
  MLIAPDataKokkos(class LAMMPS *, int, int *, class MLIAPModel *, class MLIAPDescriptor *,
                  class PairMLIAPKokkos<DeviceType> * = nullptr);
  ~MLIAPDataKokkos();
  ExecutionSpace execution_space;

  void generate_neighdata(class NeighList *, int = 0, int = 0);
  void grow_neigharrays();

  void modified(ExecutionSpace space, unsigned int mask, bool ignore_auto_sync = false);

  void sync(ExecutionSpace space, unsigned int mask, bool ignore_auto_sync = false);

  PairMLIAPKokkos<DeviceType> *k_pairmliap;

  DAT::tdual_int_1d k_iatoms;           // index of each atom
  DAT::tdual_int_1d k_ielems;           // element of each atom
  DAT::tdual_int_1d k_jatoms;           // index of each neighbor
  DAT::tdual_int_1d k_jelems;           // element of each neighbor
  DAT::tdual_int_1d k_ij;               // Start location for each particle
  DAT::tdual_float_2d k_betas;          // betas for all atoms in list
  DAT::tdual_float_2d k_descriptors;    // descriptors for all atoms in list
  DAT::tdual_float_1d k_eatoms;         // energies for all atoms in list
  DAT::tdual_float_2d k_rij;            // distance vector of each neighbor
  DAT::tdual_float_2d k_gradforce;
  DAT::tdual_float_3d k_graddesc;         // descriptor gradient w.r.t. each neighbor
  DAT::tdual_int_1d k_numneighs;          // neighbors count for each atom
  DAT::tdual_float_2d k_gamma;            // gamma element
  DAT::tdual_int_2d k_gamma_row_index;    // row (parameter) index
  DAT::tdual_int_2d k_gamma_col_index;    // column (descriptor) index

  int nij_total;

 protected:
  class LAMMPS *lmp;
};
}    // namespace LAMMPS_NS
#endif
