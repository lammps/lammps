/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Matt Bettencourt (NVIDIA)
------------------------------------------------------------------------- */

#ifndef SRC_KOKKOS_MLIAP_DATA_KOKKOS_H
#define SRC_KOKKOS_MLIAP_DATA_KOKKOS_H

#include "kokkos_type.h"
#include "mliap_data.h"
#include "pointers.h"
#include "memory_kokkos.h"
#include "pair_mliap_kokkos.h"

namespace LAMMPS_NS {
#define IATOMS_MASK           0x00000001
#define IELEMS_MASK           0x00000002
#define JATOMS_MASK           0x00000004
#define JELEMS_MASK           0x00000008
#define IJ_MASK               0x00000010
#define BETAS_MASK            0x00000020
#define DESCRIPTORS_MASK      0x00000040
#define EATOMS_MASK           0x00000080
#define RIJ_MASK              0x00000100
#define GRADFORCE_MASK        0x00000200
#define GRADDESC_MASK         0x00000400
#define NUMNEIGHS_MASK        0x00000800
#define GAMMA_MASK_MASK       0x00001000
#define GAMMA_ROW_MASK        0x00002000
#define GAMMA_COL_MASK        0x00004000

template<class DeviceType>
class MLIAPDataKokkos: public MLIAPData {
public:

  MLIAPDataKokkos(class LAMMPS *, int, int *, class MLIAPModel *, class MLIAPDescriptor *,
                  class PairMLIAPKokkos<DeviceType>* = nullptr);
  ~MLIAPDataKokkos();
  ExecutionSpace execution_space;

  void generate_neighdata(class NeighList *, int = 0, int = 0);
  void grow_neigharrays();

  void modified(ExecutionSpace space, unsigned int mask, bool ignore_auto_sync=false);

  void sync(ExecutionSpace space, unsigned int mask, bool ignore_auto_sync=false);

  PairMLIAPKokkos<DeviceType> *k_pairmliap;

  DAT::tdual_int_1d k_iatoms;// index of each atom
  DAT::tdual_int_1d k_ielems;                   // element of each atom
  DAT::tdual_int_1d k_jatoms;                   // index of each neighbor
  DAT::tdual_int_1d k_jelems;                   // element of each neighbor
  DAT::tdual_int_1d k_ij;                       // Start location for each particle
  DAT::tdual_float_2d k_betas;                  // betas for all atoms in list
  DAT::tdual_float_2d k_descriptors;            // descriptors for all atoms in list
  DAT::tdual_float_1d k_eatoms;                 // energies for all atoms in list
  DAT::tdual_float_2d k_rij;                    // distance vector of each neighbor
  DAT::tdual_float_2d k_gradforce;
  DAT::tdual_float_3d k_graddesc;               // descriptor gradient w.r.t. each neighbor
  DAT::tdual_int_1d k_numneighs;                // neighbors count for each atom
  DAT::tdual_float_2d k_gamma;                  // gamma element
  DAT::tdual_int_2d k_gamma_row_index;          // row (parameter) index
  DAT::tdual_int_2d k_gamma_col_index;          // column (descriptor) index

  int nij_total;
protected:
  class LAMMPS *lmp;
};

}

#endif /* SRC_KOKKOS_MLIAP_DATA_KOKKOS_H_ */
