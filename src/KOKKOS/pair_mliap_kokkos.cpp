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
  reverse_comm_device = 1;
  comm_type=COMM_TYPE::UNSET;

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairMLIAPKokkos<DeviceType>::~PairMLIAPKokkos()
{
  memoryKK->destroy_kokkos(k_map, map);
  memoryKK->destroy_kokkos(k_cutsq, cutsq);
  memoryKK->destroy_kokkos(k_setflag, setflag);
  memoryKK->destroy_kokkos(k_eatom, eatom);
  memoryKK->destroy_kokkos(k_vatom, vatom);
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
  if (narg < 3) error->all(FLERR,"Incorrect args for pair coefficients" + utils::errorurl(21));
  if (!allocated) {
    PairMLIAP::allocate();
    allocate();
  }

  char* type1 = arg[0];
  char* type2 = arg[1];
  char** elemtypes = &arg[2];

  // ensure I,J args are * *

  if (strcmp(type1,"*") != 0 || strcmp(type2,"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients" + utils::errorurl(21));

  // read args that map atom types to elements
  // map[i] = which element the Ith atom type is, -1 if not mapped
  // map[0] is not used

  for (int i = 1; i <= atom->ntypes; i++) {
    char* elemname = elemtypes[i-1];
    int jelem;
    for (jelem = 0; jelem < descriptor->nelements; jelem++)
      if (strcmp(elemname,descriptor->elements[jelem]) == 0)
        break;

    //printf(">>> nelements: %d\n", descriptor->nelements);
    if (jelem < descriptor->nelements)
      map[i] = jelem;
    else if (strcmp(elemname,"NULL") == 0) map[i] = -1;
    else error->all(FLERR,"Incorrect args for pair coefficients" + utils::errorurl(21));
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

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients" + utils::errorurl(21));

  // set up model, descriptor, and mliap data structures
  model->init();
  descriptor->init();

  auto h_cutsq=k_cutsq.template view<LMPHostType>();
  for (int itype=1; itype <= atom->ntypes; ++itype)
    for (int jtype=1; jtype <= atom->ntypes; ++jtype)
      // do not set cuts for NULL atoms
      if (map[itype] >= 0 && map[jtype] >= 0) {
        h_cutsq(itype,jtype) = descriptor->cutsq[map[itype]][map[jtype]];
      }
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
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0,data->nlistatoms), KOKKOS_LAMBDA (int ii) {
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

template<class DeviceType>
template<class CommType>
int PairMLIAPKokkos<DeviceType>::forward_comm(CommType* copy_from_, CommType* copy_to_, const int vl)
{
  static_assert( std::is_same_v<CommType,float>
              || std::is_same_v<CommType,double>,
                 "Unsupported CommType");
  if constexpr ( std::is_same_v<CommType,float> ) {
    comm_type = COMM_TYPE::FLOAT;
  } else if constexpr ( std::is_same_v<CommType,double> ) {
    comm_type = COMM_TYPE::DOUBLE;
  }
  copy_to = copy_to_;
  copy_from = copy_from_;
  comm_forward = vec_len=vl;

  Kokkos::parallel_for((atom->nlocal+atom->nghost)*vl, KOKKOS_LAMBDA (int i) {
    copy_to_[i] = copy_from_[i];
  });
  //call comm
  comm->forward_comm(this);

  return 0;
}

/* ---------------------------------------------------------------------- */
template<class DeviceType>
template<class CommType>
int PairMLIAPKokkos<DeviceType>::reverse_comm(CommType* copy_from_, CommType* copy_to_, const int vl)
{
  static_assert( std::is_same_v<CommType,float>
              || std::is_same_v<CommType,double>,
                 "Unsupported CommType");
  if constexpr ( std::is_same_v<CommType,float> ) {
    comm_type = COMM_TYPE::FLOAT;
  } else if constexpr ( std::is_same_v<CommType,double> ) {
    comm_type = COMM_TYPE::DOUBLE;
  }
  copy_to = copy_to_;
  copy_from = copy_from_;
  comm_reverse = vec_len = vl;

  Kokkos::parallel_for((atom->nlocal+atom->nghost)*vl, KOKKOS_LAMBDA (int i) {
    copy_to_[i] = copy_from_[i]; // Copy inputs
  });

  comm->reverse_comm(this);

  Kokkos::parallel_for(
    Kokkos::RangePolicy<>( (atom->nlocal)*vl,(atom->nlocal + atom->nghost)*vl),
    KOKKOS_LAMBDA (int i) {
    copy_to_[i] = 0; //Zero out ghosts
  });

  return 0;
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
int PairMLIAPKokkos<DeviceType>::pack_forward_comm_kokkos(
    int nv, DAT::tdual_int_1d idx_v, DAT::tdual_xfloat_1d &fill, int int2,
    int *intp) {
  switch( comm_type ) {
    case COMM_TYPE::FLOAT:
      return pack_forward_comm_kokkos(nv,idx_v,fill,int2,intp,std::get<float*>(copy_to));
    case COMM_TYPE::DOUBLE:
      return pack_forward_comm_kokkos(nv,idx_v,fill,int2,intp,std::get<double*>(copy_to));
    case COMM_TYPE::UNSET:
    default:
      error->all(FLERR,"comm_type was never set");
      return -1;
  }
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
template <typename CommType>
int PairMLIAPKokkos<DeviceType>::pack_forward_comm_kokkos(
    int nv, DAT::tdual_int_1d idx_v, DAT::tdual_xfloat_1d &fill, int int2,
    int *intp, CommType *copy_to) {
  auto idx=idx_v.view<DeviceType>();
  auto val=fill.view<DeviceType>();
  int nf=vec_len;
  auto to=copy_to;
  Kokkos::parallel_for(nv, KOKKOS_LAMBDA (int i) {
    int gstart=idx(i)*nf;
    int start=i*nf;
    for (int j=0;j<nf;++j)
      val(start++) = static_cast<double>(to[gstart++]);
  });
  return nv*nf;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int PairMLIAPKokkos<DeviceType>::pack_forward_comm(int nv, int* idx_v, double *fill,
                                                   int int2, int *intp)
{
  static bool first=true;
  if (first) {
    error->warning(FLERR,"PackForwardComm has only been tested on Kokkos devices");
    first=false;
  }
  switch( comm_type ) {
    case COMM_TYPE::FLOAT:
      return pack_forward_comm(nv,idx_v,fill,int2,intp,std::get<float*>(copy_to));
    case COMM_TYPE::DOUBLE:
      return pack_forward_comm(nv,idx_v,fill,int2,intp,std::get<double*>(copy_to));
    case COMM_TYPE::UNSET:
    default:
      error->all(FLERR,"comm_type was never set");
      return -1;
  }
}

template<class DeviceType>
template <typename CommType>
int PairMLIAPKokkos<DeviceType>::pack_forward_comm(int nv, int* idx_v, double *fill,
                                                   int int2, int *intp, CommType *copy_to)
{
  for (int i=0;i<nv;++i) {
    int gstart=idx_v[i]*vec_len;
    int start=i*vec_len;
    for (int j=0;j<vec_len;++j)
      fill[start++] = static_cast<double>(copy_to[gstart++]);
  }
  return nv*vec_len;
}
/* ---------------------------------------------------------------------- */

template <class DeviceType>
void PairMLIAPKokkos<DeviceType>::unpack_forward_comm_kokkos(
    int nv, int first_up, DAT::tdual_xfloat_1d &fill) {
  switch( comm_type ) {
    case COMM_TYPE::FLOAT:
      return unpack_forward_comm_kokkos(nv,first_up,fill,std::get<float*>(copy_to));
    case COMM_TYPE::DOUBLE:
      return unpack_forward_comm_kokkos(nv,first_up,fill,std::get<double*>(copy_to));
    case COMM_TYPE::UNSET:
    default:
      error->all(FLERR,"comm_type was never set");
  }
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
template <typename CommType>
void PairMLIAPKokkos<DeviceType>::unpack_forward_comm_kokkos(
    int nv, int first_up, DAT::tdual_xfloat_1d &fill, CommType *copy_to) {
  auto val=fill.view<DeviceType>();
  int nf=vec_len;

  Kokkos::parallel_for(nv, KOKKOS_LAMBDA (int i) {
    int gstart=(first_up+i)*nf;
    int start=i*nf;
    for (int j=0;j<nf;++j) {
      copy_to[gstart+j] = static_cast<CommType>(val(start+j));
    }
  });
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairMLIAPKokkos<DeviceType>::unpack_forward_comm(int nv, int first_up, double *fill)
{
  switch( comm_type ) {
    case COMM_TYPE::FLOAT:
      return unpack_forward_comm(nv,first_up,fill,std::get<float*>(copy_to));
    case COMM_TYPE::DOUBLE:
      return unpack_forward_comm(nv,first_up,fill,std::get<double*>(copy_to));
    case COMM_TYPE::UNSET:
    default:
      error->all(FLERR,"comm_type was never set");
  }
}

template <class DeviceType>
template <typename CommType>
void PairMLIAPKokkos<DeviceType>::unpack_forward_comm(
    int nv, int first_up, double *fill, CommType *copy_to) {
  for (int i=0; i<nv; ++i) {
    int gstart=(first_up+i)*vec_len;
    int start=i*vec_len;
    for (int j=0;j<vec_len;++j) {
      copy_to[gstart+j] = static_cast<CommType>(fill[start+j]);
    }
  }
}
/* ---------------------------------------------------------------------- */

template<class DeviceType>
int PairMLIAPKokkos<DeviceType>::pack_reverse_comm_kokkos(int nv, int first_up, DAT::tdual_xfloat_1d &fill)
{
  switch( comm_type ) {
    case COMM_TYPE::FLOAT:
      return pack_reverse_comm_kokkos(nv,first_up,fill,std::get<float*>(copy_to));
    case COMM_TYPE::DOUBLE:
      return pack_reverse_comm_kokkos(nv,first_up,fill,std::get<double*>(copy_to));
    case COMM_TYPE::UNSET:
    default:
      error->all(FLERR,"comm_type was never set");
  }
  return -1;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<typename CommType>
int PairMLIAPKokkos<DeviceType>::pack_reverse_comm_kokkos(int nv, int first_up, DAT::tdual_xfloat_1d &fill, CommType *copy_to)
{
  int nf=vec_len;
  auto val=fill.view<DeviceType>();
  Kokkos::parallel_for(nv, KOKKOS_LAMBDA (int i) {
    int gstart=(first_up+i)*nf;
    int start=i*nf;
    for (int j=0;j<nf;++j) {
      val(start++) = static_cast<double>(copy_to[gstart++]);
    }
  });
  return nv*nf;
}
/* ---------------------------------------------------------------------- */

template<class DeviceType>
int PairMLIAPKokkos<DeviceType>::pack_reverse_comm(int nv, int first_up, double *fill)
{
  switch( comm_type ) {
    case COMM_TYPE::FLOAT:
      return pack_reverse_comm(nv,first_up,fill,std::get<float*>(copy_to));
    case COMM_TYPE::DOUBLE:
      return pack_reverse_comm(nv,first_up,fill,std::get<double*>(copy_to));
    case COMM_TYPE::UNSET:
    default:
      error->all(FLERR,"comm_type was never set");
  }
  return -1;
}

/* ---------------------------------------------------------------------- */
template<class DeviceType>
template<typename CommType>
int PairMLIAPKokkos<DeviceType>::pack_reverse_comm(int nv, int first_up, double *fill, CommType *copy_to)
{
  for (int i=0;i<nv;++i) {
    int gstart=(first_up+i)*vec_len;
    int start=i*vec_len;
    for (int j=0;j<vec_len;++j) {
      fill[start++] = static_cast<double>(copy_to[gstart++]);
    }
  }
  return nv*vec_len;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairMLIAPKokkos<DeviceType>::unpack_reverse_comm_kokkos(int nv, DAT::tdual_int_1d idx_v, DAT::tdual_xfloat_1d &fill)
{
  switch( comm_type ) {
    case COMM_TYPE::FLOAT:
      return unpack_reverse_comm_kokkos(nv,idx_v,fill,std::get<float*>(copy_to));
    case COMM_TYPE::DOUBLE:
      return unpack_reverse_comm_kokkos(nv,idx_v,fill,std::get<double*>(copy_to));
    case COMM_TYPE::UNSET:
    default:
      error->all(FLERR,"comm_type was never set");
      return;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<typename CommType>
void PairMLIAPKokkos<DeviceType>::unpack_reverse_comm_kokkos(int nv, DAT::tdual_int_1d idx_v, DAT::tdual_xfloat_1d &fill, CommType *copy_to)
{
  int nf=vec_len;
  auto val=fill.view<DeviceType>();
  auto idx=idx_v.view<DeviceType>();
  auto to=copy_to;
  Kokkos::parallel_for(nv, KOKKOS_LAMBDA (int i) {
    int gstart=idx(i)*nf;
    int start=i*nf;
    for (int j=0;j<nf;++j)
      to[gstart++] += static_cast<CommType>(val(start++));
  });
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairMLIAPKokkos<DeviceType>::unpack_reverse_comm(int nv, int *idx, double *fill)
{
  switch( comm_type ) {
    case COMM_TYPE::FLOAT:
      return unpack_reverse_comm(nv,idx,fill,std::get<float*>(copy_to));
    case COMM_TYPE::DOUBLE:
      return unpack_reverse_comm(nv,idx,fill,std::get<double*>(copy_to));
    case COMM_TYPE::UNSET:
    default:
      error->all(FLERR,"comm_type was never set");
      return;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<typename CommType>
void PairMLIAPKokkos<DeviceType>::unpack_reverse_comm(int nv, int *idx, double *fill, CommType *copy_to)
{
  for (int i=0;i<nv;++i) {
    int gstart=idx[i]*vec_len;
    int start=i*vec_len;
    for (int j=0;j<vec_len;++j)
      copy_to[gstart++] += static_cast<CommType>(fill[start++]);
  }
}
namespace LAMMPS_NS {
template class PairMLIAPKokkos<LMPDeviceType>;
template int PairMLIAPKokkos<LMPDeviceType>::forward_comm<float>(float*,float*,const int);
template int PairMLIAPKokkos<LMPDeviceType>::forward_comm<double>(double*,double*,const int);
template int PairMLIAPKokkos<LMPDeviceType>::reverse_comm<float>(float*,float*,const int);
template int PairMLIAPKokkos<LMPDeviceType>::reverse_comm<double>(double*,double*,const int);
#ifdef LMP_KOKKOS_GPU
template class PairMLIAPKokkos<LMPHostType>;
template int PairMLIAPKokkos<LMPHostType>::forward_comm<float>(float*,float*,const int);
template int PairMLIAPKokkos<LMPHostType>::forward_comm<double>(double*,double*,const int);
template int PairMLIAPKokkos<LMPHostType>::reverse_comm<float>(float*,float*,const int);
template int PairMLIAPKokkos<LMPHostType>::reverse_comm<double>(double*,double*,const int);
#endif
}
