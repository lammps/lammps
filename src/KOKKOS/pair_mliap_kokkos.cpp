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

#include "pair_mliap_kokkos.h"
#include "memory_kokkos.h"
#include "atom_kokkos.h"
#include "atom_masks.h"
#include "mliap_data_kokkos.h"
#include "mliap_descriptor_so3_kokkos.h"
#include "mliap_model_linear_kokkos.h"
#include "error.h"
#include "neigh_request.h"
#include "lammps.h"
#include "kokkos.h"
#include "pointers.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairMLIAPKokkos<DeviceType>::PairMLIAPKokkos(class LAMMPS* l) : PairMLIAP(l)
{
  kokkosable = 1;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairMLIAPKokkos<DeviceType>::~PairMLIAPKokkos()
{
  memoryKK->destroy_kokkos(k_map, map);
  memoryKK->destroy_kokkos(k_cutsq, cutsq);
  memoryKK->destroy_kokkos(k_setflag, setflag);
  memoryKK->destroy_kokkos(k_eatom,eatom);
  memoryKK->destroy_kokkos(k_vatom,vatom);
  allocated = 0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairMLIAPKokkos<DeviceType>::compute(int eflag, int vflag)
{
  atomKK->sync(Host,F_MASK | ENERGY_MASK | VIRIAL_MASK);
  atomKK->sync(execution_space,X_MASK | TYPE_MASK );
  MLIAPDataKokkos<DeviceType> *k_data = (MLIAPDataKokkos<DeviceType>*)(data);

  int is_kokkos_model = (dynamic_cast<MLIAPModelKokkos<DeviceType>*>(model)) != nullptr;
  int is_kokkos_descriptor = (dynamic_cast<MLIAPDescriptorKokkos<DeviceType>*>(descriptor)) != nullptr;

  // consistency checks
  if (data->ndescriptors != model->ndescriptors)
    error->all(FLERR, "Incompatible model and descriptor descriptor count");

  if (data->nelements != model->nelements)
    error->all(FLERR, "Incompatible model and descriptor element count");

  ev_init(eflag, vflag);
  if (eflag_atom && k_eatom.h_view.extent(0) < maxeatom) {
     memoryKK->destroy_kokkos(k_eatom,eatom);
     memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"pair:eatom");
  }

  if (vflag_atom && k_vatom.h_view.extent(0) < maxeatom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxeatom,6,"pair:eatom");
  }

  data->generate_neighdata(list, eflag, vflag);

  // compute descriptors, if needed
  if (!is_kokkos_descriptor) {
    // put sync stuff here
  }
  if (model->nonlinearflag || eflag) descriptor->compute_descriptors(data);
  if (!is_kokkos_descriptor)
    k_data->modified(Host, DESCRIPTORS_MASK);
  if (is_kokkos_model)
    k_data->sync(execution_space, DESCRIPTORS_MASK);

  // compute E_i and beta_i = dE_i/dB_i for all i in list
  model->compute_gradients(data);

  // This is kokkos
  if (!is_kokkos_model)
    k_data->sync(execution_space, IATOMS_MASK | EATOMS_MASK);
  e_tally(data);

  // calculate force contributions beta_i*dB_i/dR_j
  if (is_kokkos_model)
    k_data->modified(execution_space, BETAS_MASK | EATOMS_MASK);
  if (!is_kokkos_descriptor)
    k_data->sync(Host, BETAS_MASK | EATOMS_MASK);
  descriptor->compute_forces(data);

  if (!is_kokkos_descriptor) {
    if (evflag) {
      atomKK->modified(Host,F_MASK | ENERGY_MASK | VIRIAL_MASK);
    } else {
      atomKK->modified(Host,F_MASK);
    }
  } else {
    if (evflag)
      atomKK->modified(execution_space,F_MASK | ENERGY_MASK | VIRIAL_MASK);
    else
      atomKK->modified(execution_space,F_MASK);
  }
  // calculate stress
  if (vflag_fdotr) {
    pair_virial_fdotr_compute(this);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairMLIAPKokkos<DeviceType>::allocate()
{
  int n = atom->ntypes;
  memoryKK->destroy_kokkos(k_map, map);
  memoryKK->destroy_kokkos(k_cutsq, cutsq);
  memoryKK->destroy_kokkos(k_setflag, setflag);
  memoryKK->create_kokkos(k_map, map, n+1, "pair_mliap:map");
  memoryKK->create_kokkos(k_cutsq, cutsq, n+1, n+1, "pair_mliap:cutsq");
  memoryKK->create_kokkos(k_setflag, setflag, n+1, n+1, "pair_mliap:setflag");

  auto h_cutsq=k_cutsq.template view<LMPHostType>();
  n = descriptor->nelements;
  for (int i=0;i<n;++i)
    for (int j=0;j<n;++j) {
      h_cutsq(i,j) = descriptor->cutsq[i][j];
    }
  k_cutsq.modify<LMPHostType>();
  k_cutsq.sync<LMPDeviceType>();

  // this is for the base class so it doesn't double delete
  allocated = 1;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairMLIAPKokkos<DeviceType>::settings(int narg, char ** arg)
{
  PairMLIAP::settings(narg, arg);
  int iarg=0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"model") == 0) {
      if (strcmp(arg[iarg+1],"linear") == 0) {
        if (iarg+3 > narg) error->all(FLERR,"Illegal pair_style mliap command");
        delete model;
        model = new MLIAPModelLinearKokkos<DeviceType>(lmp,arg[iarg+2]);
        iarg += 3;
      } else
        iarg += 2;
    } else if (strcmp(arg[iarg],"descriptor") == 0) {
      if (strcmp(arg[iarg+1],"so3") == 0) {
        if (iarg+3 > narg) error->all(FLERR,"Illegal pair_style mliap command");
        delete descriptor;
        descriptor = new MLIAPDescriptorSO3Kokkos<DeviceType>(lmp,arg[iarg+2]);
        iarg += 3;
      } else
        iarg ++;
    } else
      iarg++;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairMLIAPKokkos<DeviceType>::coeff(int narg, char **arg) {
  if (narg < 3) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) {
    PairMLIAP::allocate();
    allocate();
  }

  char* type1 = arg[0];
  char* type2 = arg[1];
  char** elemtypes = &arg[2];

  // insure I,J args are * *

  if (strcmp(type1,"*") != 0 || strcmp(type2,"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // read args that map atom types to elements
  // map[i] = which element the Ith atom type is, -1 if not mapped
  // map[0] is not used

  for (int i = 1; i <= atom->ntypes; i++) {
    char* elemname = elemtypes[i-1];
    int jelem;
    for (jelem = 0; jelem < descriptor->nelements; jelem++)
      if (strcmp(elemname,descriptor->elements[jelem]) == 0)
        break;

    if (jelem < descriptor->nelements)
      map[i] = jelem;
    else if (strcmp(elemname,"NULL") == 0) map[i] = -1;
    else error->all(FLERR,"Incorrect args for pair coefficients");
  }
  k_map.modify<LMPHostType>();
  k_map.sync<LMPDeviceType>();

  // clear setflag since coeff() called once with I,J = * *

  int n = atom->ntypes;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  // set setflag i,j for type pairs where both are mapped to elements

  int count = 0;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      if (map[i] >= 0 && map[j] >= 0) {
        setflag[i][j] = 1;
        count++;
      }
  k_setflag.modify<LMPHostType>();
  k_setflag.sync<LMPDeviceType>();

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");

  // set up model, descriptor, and mliap data structures
  model->init();
  descriptor->init();
  int gradgradflag = -1;
  delete data;
  data = new MLIAPDataKokkos<DeviceType>(lmp, gradgradflag, map, model, descriptor, this);
  data->init();
}

/* ----------------------------------------------------------------------
   add energies to eng_vdwl and per-atom energy
------------------------------------------------------------------------- */

template<class DeviceType>
void PairMLIAPKokkos<DeviceType>::e_tally(MLIAPData* data)
{
  if (eflag_global) eng_vdwl += data->energy;
  if (eflag_atom) {
    MLIAPDataKokkos<DeviceType> *k_data = static_cast<MLIAPDataKokkos<DeviceType>*>(data);
    k_data->sync(execution_space, IATOMS_MASK);
    auto d_iatoms = k_data->k_iatoms.template view<DeviceType>();
    auto d_eatoms = k_data->k_eatoms.template view<DeviceType>();
    auto d_eatom = k_eatom.template view<DeviceType>();
    Kokkos::parallel_for(data->nlistatoms, KOKKOS_LAMBDA (int ii) {
      d_eatom(d_iatoms(ii)) = d_eatoms(ii);
    });
    k_eatom.modify<DeviceType>();
    // This sync has to be here for the hybrid pair type
    k_eatom.sync<LMPHostType>();
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairMLIAPKokkos<DeviceType>::init_style()
{

  PairMLIAP::init_style();
  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same<DeviceType,LMPHostType>::value &&
                           !std::is_same<DeviceType,LMPDeviceType>::value);
  request->set_kokkos_device(std::is_same<DeviceType,LMPDeviceType>::value);
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class PairMLIAPKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairMLIAPKokkos<LMPHostType>;
#endif
}
