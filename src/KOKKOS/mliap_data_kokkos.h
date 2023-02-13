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
  PAIR_I_MASK      = 0x00008000,
  ELEMS_MASK       = 0x00010000,
};
// clang-format on

template <class DeviceType> class MLIAPDataKokkos : public MLIAPData {
 public:
  MLIAPDataKokkos(class LAMMPS *, int, int *, class MLIAPModel *, class MLIAPDescriptor *,
                  class PairMLIAPKokkos<DeviceType> * = nullptr);
  ~MLIAPDataKokkos() override;
  ExecutionSpace execution_space;

  void generate_neighdata(class NeighList *, int = 0, int = 0) override;
  void grow_neigharrays() override;

  void modified(ExecutionSpace space, unsigned int mask, bool ignore_auto_sync = false);

  void sync(ExecutionSpace space, unsigned int mask, bool ignore_auto_sync = false);

  PairMLIAPKokkos<DeviceType> *k_pairmliap;

  DAT::tdual_int_1d k_iatoms;           // index of each atom
  DAT::tdual_int_1d k_ielems;           // element of each atom
  DAT::tdual_int_1d k_jatoms;           // index of each neighbor
  DAT::tdual_int_1d k_elems;            // element of each atom in or not in the neighborlist
  DAT::tdual_int_1d k_pair_i;           // index of each i atom for each ij pair
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

  // Just cached for python interface
  double *f_device;

 protected:
  class LAMMPS *lmp;
};

// Now we need a specific device version for communication with python
class MLIAPDataKokkosDevice {
public:

  MLIAPDataKokkosDevice(MLIAPDataKokkos<LMPDeviceType> &base) :
    size_array_rows(base.size_array_rows),
    size_array_cols(base.size_array_cols),
    natoms(base.natoms),
    yoffset(base.yoffset),
    zoffset(base.zoffset),
    ndims_force(base.ndims_force),
    ndims_virial(base.ndims_virial),
    size_gradforce(base.size_gradforce),
    f(base.f_device),
    gradforce(base.k_gradforce.d_view.data()),
    betas(base.k_betas.d_view.data()),
    descriptors(base.k_descriptors.d_view.data()),
    eatoms(base.k_eatoms.d_view.data()),
    energy(&base.energy),
    ndescriptors(base.ndescriptors),
    nparams(base.nparams),
    nelements(base.nelements),
    gamma_nnz(base.gamma_nnz),
    gamma(base.k_gamma.d_view.data()),
    gamma_row_index(base.k_gamma_row_index.d_view.data()),
    gamma_col_index(base.k_gamma_col_index.d_view.data()),
    egradient(nullptr),
    ntotal(base.ntotal),
    nlistatoms(base.nlistatoms),
    natomneigh(base.natomneigh),
    numneighs(base.numneighs),
    iatoms(base.k_iatoms.d_view.data()),
    pair_i(base.k_pair_i.d_view.data()),
    ielems(base.k_ielems.d_view.data()),
    nneigh_max(base.nneigh_max),
    npairs(base.npairs),
    jatoms(base.k_jatoms.d_view.data()),
    jelems(base.k_jelems.d_view.data()),
    elems(base.k_elems.d_view.data()),
    rij(base.k_rij.d_view.data()),
    graddesc(base.k_graddesc.d_view.data()),
    eflag(base.eflag),
    vflag(base.vflag),
    pairmliap(dynamic_cast<PairMLIAPKokkos<LMPDeviceType> *>(base.pairmliap)),
#if defined(KOKKOS_ENABLE_CUDA)
    dev(1)
#else
    dev(0)
#endif
    {  }
  int size_array_rows;
  int size_array_cols;
  int natoms;
  int yoffset;
  int zoffset;
  int ndims_force;
  int ndims_virial;
  int size_gradforce;

  //Write only
  double *f;
  double *gradforce;
  double *betas;
  double *descriptors;
  double *eatoms;
  double *energy;

  // sizing
  const int ndescriptors;
  const int nparams;
  const int nelements;

  //Ignored for now
  int gamma_nnz;
  double *gamma;
  int *gamma_row_index;
  int *gamma_col_index;
  double *egradient;

  // Neighborlist stuff
  const int ntotal;
  const int nlistatoms;
  const int natomneigh;
  int *numneighs;
  int *iatoms;
  int *pair_i;
  int *ielems;
  const int nneigh_max;
  const int npairs;
  int *jatoms;
  int *jelems;
  int *elems;
  double *rij;
  double *graddesc;
  int eflag;
  int vflag;

  class PairMLIAPKokkos<LMPDeviceType> *pairmliap;    // access to pair tally functions

  int dev;

#ifdef LMP_KOKKOS_GPU
  MLIAPDataKokkosDevice(MLIAPDataKokkos<LMPHostType> &base) : ndescriptors(-1),nparams(-1),nelements(-1),ntotal(-1),nlistatoms(-1),natomneigh(-1),
      nneigh_max(-1),npairs(-1)
  {
    // It cannot get here, but needed for compilation
  }
#endif
};


}    // namespace LAMMPS_NS
#endif
