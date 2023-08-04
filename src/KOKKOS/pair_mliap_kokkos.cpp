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

#include "pair_mliap_kokkos.h"
#include "memory_kokkos.h"
#include "atom_kokkos.h"
#include "atom_masks.h"
#include "mliap_data_kokkos.h"
#include "mliap_descriptor_so3_kokkos.h"
#include "mliap_model_linear_kokkos.h"
#ifdef MLIAP_PYTHON
#include "mliap_model_python_kokkos.h"
#include "mliap_unified_kokkos.h"
#endif
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
  datamask_modify = 0;
  is_child=true;
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
  delete model;
  delete descriptor;
  model=nullptr;
  descriptor=nullptr;
  allocated = 0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairMLIAPKokkos<DeviceType>::compute(int eflag, int vflag)
{
  atomKK->sync(execution_space,X_MASK | TYPE_MASK );
  MLIAPDataKokkos<DeviceType> *k_data = (MLIAPDataKokkos<DeviceType>*)(data);

  int is_kokkos_model = (dynamic_cast<MLIAPModelKokkos<DeviceType>*>(model)) != nullptr;
  int is_kokkos_descriptor = (dynamic_cast<MLIAPDescriptorKokkos<DeviceType>*>(descriptor)) != nullptr;
  auto model_space = is_kokkos_model ? execution_space : Host;
  auto descriptor_space = is_kokkos_descriptor? execution_space : Host;
  // consistency checks
  if (data->ndescriptors != model->ndescriptors)
    error->all(FLERR, "Incompatible model and descriptor descriptor count");

  if (data->nelements != model->nelements)
    error->all(FLERR, "Incompatible model and descriptor element count");

  ev_init(eflag, vflag);
  if (eflag_atom && (int)k_eatom.h_view.extent(0) < maxeatom) {
     memoryKK->destroy_kokkos(k_eatom,eatom);
     memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"pair:eatom");
  }

  if (vflag_atom && (int)k_vatom.h_view.extent(0) < maxeatom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxeatom,6,"pair:eatom");
  }

  data->generate_neighdata(list, eflag, vflag);

  // compute descriptors, if needed
  if (model->nonlinearflag || eflag)  {
    k_data->sync(descriptor_space, NUMNEIGHS_MASK | IATOMS_MASK | IELEMS_MASK | ELEMS_MASK | JATOMS_MASK | PAIR_I_MASK | JELEMS_MASK | RIJ_MASK );
    descriptor->compute_descriptors(data);
    if (!is_kokkos_descriptor)
      k_data->modified(descriptor_space, DESCRIPTORS_MASK);
  }

  // compute E_i and beta_i = dE_i/dB_i for all i in list
  k_data->sync(model_space, IELEMS_MASK | DESCRIPTORS_MASK);
  model->compute_gradients(data);
  k_data->modified(model_space, BETAS_MASK);
  if (eflag_atom) {
    k_data->modified(model_space, EATOMS_MASK);
  }

  // calculate force contributions beta_i*dB_i/dR_j
  atomKK->sync(descriptor_space,F_MASK);
  k_data->sync(descriptor_space, NUMNEIGHS_MASK | IATOMS_MASK | IELEMS_MASK | ELEMS_MASK | BETAS_MASK | JATOMS_MASK | PAIR_I_MASK | JELEMS_MASK | RIJ_MASK );

  descriptor->compute_forces(data);

  e_tally(data);

  if (evflag) {
    atomKK->modified(descriptor_space,F_MASK | ENERGY_MASK | VIRIAL_MASK);
    atomKK->sync(execution_space,F_MASK | ENERGY_MASK | VIRIAL_MASK);
  } else {
    atomKK->modified(descriptor_space,F_MASK);
    atomKK->sync(execution_space,F_MASK);
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

  // this is for the base class so it doesn't double delete
  allocated = 1;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairMLIAPKokkos<DeviceType>::settings(int narg, char ** arg)
{
  std::vector<char*> new_args;
  int iarg=0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"model") == 0) {
      if (strcmp(arg[iarg+1],"linear") == 0) {
        if (iarg+3 > narg) error->all(FLERR,"Illegal pair_style mliap command");
        delete model;
        model = new MLIAPModelLinearKokkos<DeviceType>(lmp,arg[iarg+2]);
        iarg += 3;
      } else if (strcmp(arg[iarg+1],"mliappy") == 0) {
#ifdef MLIAP_PYTHON
        if (iarg+3 > narg) utils::missing_cmd_args(FLERR, "pair_style mliap mliappy", error);
        delete model;
        model = new MLIAPModelPythonKokkos<DeviceType>(lmp,arg[iarg+2]);
        iarg += 3;
#else
        error->all(FLERR,"Using pair_style mliap model mliappy requires ML-IAP with python support");
#endif
      } else {
        new_args.push_back(arg[iarg++]);
        new_args.push_back(arg[iarg++]);
      }
    } else if (strcmp(arg[iarg],"descriptor") == 0) {
      if (strcmp(arg[iarg+1],"so3") == 0) {
        if (iarg+3 > narg) error->all(FLERR,"Illegal pair_style mliap command");
        delete descriptor;
        descriptor = new MLIAPDescriptorSO3Kokkos<DeviceType>(lmp,arg[iarg+2]);
        iarg += 3;
      } else
        new_args.push_back(arg[iarg++]);
    } else if (strcmp(arg[iarg], "unified") == 0) {
#ifdef MLIAP_PYTHON
      if (model != nullptr) error->all(FLERR,"Illegal multiple pair_style mliap model definitions");
      if (descriptor != nullptr) error->all(FLERR,"Illegal multiple pair_style mliap descriptor definitions");
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "pair_style mliap unified", error);
      MLIAPBuildUnifiedKokkos_t<DeviceType> build = build_unified(arg[iarg+1], dynamic_cast<MLIAPDataKokkos<DeviceType>*>(data), lmp);
      if (iarg+3 > narg) {
        ghostneigh = 0;
      } else {
        ghostneigh = utils::logical(FLERR, arg[iarg+2], false, lmp);
      }

      iarg += 3;
      model = build.model;
      descriptor = build.descriptor;
#else
      error->all(FLERR,"Using pair_style mliap unified requires ML-IAP with python support");
#endif
    } else
      new_args.push_back(arg[iarg++]);
  }
  PairMLIAP::settings(new_args.size(), new_args.data());

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

  // ensure I,J args are * *

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

  auto h_cutsq=k_cutsq.template view<LMPHostType>();
  for (int itype=1; itype <= atom->ntypes; ++itype)
    for (int jtype=1; jtype <= atom->ntypes; ++jtype)
      h_cutsq(itype,jtype) = descriptor->cutsq[map[itype]][map[jtype]];
  k_cutsq.modify<LMPHostType>();
  k_cutsq.sync<DeviceType>();
  constexpr int gradgradflag = -1;
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
    k_data->sync(execution_space, IATOMS_MASK | EATOMS_MASK, true);
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
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class PairMLIAPKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairMLIAPKokkos<LMPHostType>;
#endif
}
