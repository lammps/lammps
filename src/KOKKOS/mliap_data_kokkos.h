/* -*- c++ -*- ----------------------------------------------------------
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
   Contributing author: Matt Bettencourt (NVIDIA)
------------------------------------------------------------------------- */

#ifndef LMP_MLIAP_DATA_KOKKOS_H
#define LMP_MLIAP_DATA_KOKKOS_H

#include "mliap_data.h"

#include "kokkos_type.h"
#include "memory_kokkos.h"
#include "pair_mliap_kokkos.h"
#include "pointers.h"
#include <unordered_map>

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

template <class DeviceType> struct ExtraPropertiesKokkos {
  using execution_space = typename DeviceType::execution_space;
  using memory_space = typename DeviceType::memory_space;

  std::unordered_map<std::string, DAT::tdual_float_2d> data;
  std::unordered_map<std::string, int> dims;
  int nproperties; //May not be needed
  int nlistatoms;

  ExtraPropertiesKokkos() : nproperties(0), nlistatoms(0) {}
                      
  int get_dim(std::string name) {
    if (dims.find(name) == dims.end()) {
      return -1;
    }
    return dims[name];
  }

  LMP_FLOAT* get_device_pointer(std::string name) {
    if (data.find(name) == data.end()) {
      return nullptr;
    }
    return data[name].d_view.data();
  }

  Kokkos::View<LMP_FLOAT**, Kokkos::LayoutRight, Kokkos::HostSpace> get_host_view(std::string name) {
    if (data.find(name) == data.end()) {
      return DAT::tdual_float_2d().h_view;
    }
    return data[name].h_view;
  }

  void register_extra_property(std::string name, int dim) {
    data[name] = DAT::tdual_float_2d("extra_property_register", nlistatoms, dim);
    dims[name] = dim;
    nproperties++;
  }

  void grow(int new_nlistatoms) {
    //Check new size, return if smaller than or equal tonmax
    if (new_nlistatoms <= nlistatoms) return;
    //Iterate through each extra property and expand its size
    for (auto& pair : data) {
      DAT::tdual_float_2d new_view("extra_property_grow", new_nlistatoms, dims[pair.first]);
      auto newSubViewD = Kokkos::subview(new_view.d_view, std::make_pair(0, nlistatoms), std::make_pair(0,dims[pair.first]));
      auto newSubViewH = Kokkos::subview(new_view.h_view, std::make_pair(0, nlistatoms), std::make_pair(0,dims[pair.first]));
      Kokkos::deep_copy(newSubViewD, pair.second.d_view);
      Kokkos::deep_copy(newSubViewH, pair.second.h_view);
      pair.second = new_view;
    }
    nlistatoms = new_nlistatoms;
  }
    
  void modify_host() {
    for (auto& pair : data) {
      pair.second.template modify<Kokkos::HostSpace>();
    }
  }

  void modify_device() {
    for (auto& pair : data) {                                                                       
      pair.second.template modify<execution_space>();
    }
  }

  void sync_host() {
    for (auto& pair : data) {                                                                       
      pair.second.template sync<Kokkos::HostSpace>();
    }
  }

  void sync_device() {
    for (auto& pair : data) {                                                                       
      pair.second.template sync<execution_space>();
    }
  }
};

//Forward declaration due to circular reference to PairMLIAPKokkos.
template <class DeviceType>
class PairMLIAPKokkos;

template <class DeviceType> class MLIAPDataKokkos : public MLIAPData {
 public:
  MLIAPDataKokkos(class LAMMPS *, int, int *, class MLIAPModel *, class MLIAPDescriptor *,
                  class PairMLIAPKokkos<DeviceType> * = nullptr);
  ~MLIAPDataKokkos() override;
  ExecutionSpace execution_space;

  void generate_neighdata(class NeighList *, int = 0, int = 0) override;
  void grow_neigharrays() override;
  void register_extra_property(std::string, int);

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
  ExtraPropertiesKokkos<DeviceType> k_extra_properties;   // extra_properties data structure
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

  MLIAPDataKokkosDevice(MLIAPDataKokkos<LMPDeviceType> &base);

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
  const int nlocal;
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

  ExtraPropertiesKokkos<LMPDeviceType> extra_properties;
  //Wrapper functions for access to ExtraPropertiesKokkos struct of MLIAPDataKokkos
  int get_extra_property_dim(const char*);
  LMP_FLOAT* get_extra_property_device_pointer(const char*);

  //forward_exchange writes into ghosts
  template <typename CommType>
  void forward_exchange(CommType* copy_from, CommType* copy_to, const int vec_len){
    pairmliap->forward_comm(copy_from, copy_to, vec_len);
  }

  //reverse_exchange adds from ghosts and zeros out ghosts afterwards
  template <typename CommType>
  void reverse_exchange(CommType* copy_from, CommType* copy_to, const int vec_len){
    pairmliap->reverse_comm(copy_from, copy_to, vec_len);
  }


#ifdef LMP_KOKKOS_GPU
  MLIAPDataKokkosDevice(MLIAPDataKokkos<LMPHostType> &base) : ndescriptors(-1),nparams(-1),nelements(-1),ntotal(-1),nlistatoms(-1),nlocal(-1),natomneigh(-1),
      nneigh_max(-1),npairs(-1)
  {
    // It cannot get here, but needed for compilation
  }
#endif
};


}    // namespace LAMMPS_NS
#endif
