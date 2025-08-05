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

#include "compute_mliap_property_atom_kokkos.h"
#include "compute.h"
#include "atom.h"
#include "force.h"
#include "memory.h"
#include "error.h"
#include "pair_mliap.h"
#include "kokkos_type.h"
#include <cstring>
#include <cstdio>
using namespace LAMMPS_NS;

template <class DeviceType>
ComputeMLIAPPropertyAtomKokkos<DeviceType>::ComputeMLIAPPropertyAtomKokkos(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg)
{
  kokkosable = 1;
  peratom_flag = 1;
  //Check for args here
  int curArg = 3;
  int nameFlag = 0, dimFlag = 0;
  while (curArg < narg - 1) //All args come in pairs of two argName argValue
  //If this changes you need to do bounds checking in each argName
  {
    if(strcmp(arg[curArg], "name") == 0) { //Parse name
      property_name = arg[curArg+1];
      nameFlag = 1; 
    }
    if(strcmp(arg[curArg], "dim") == 0) { //Parse dim
      size_peratom_cols = atoi(arg[curArg+1]);
      dimFlag = 1;
      if (size_peratom_cols <= 0) {
        error->all(FLERR, "Compute mliap/property/atom/kk requires dim > 0.");
      }
    }
    curArg += 2; //Must put inside ifs if variable length args
  }

  //If not all essential properties were set. Error
  if (nameFlag == 0 || dimFlag == 0)
  {
    error->all(FLERR, "Compute mliap/property/atom/kk missing args (Expecting name and dim).");
  }

  //Initalize Member Variables
  k_data = nullptr;
  nlocal = 0;
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
ComputeMLIAPPropertyAtomKokkos<DeviceType>::~ComputeMLIAPPropertyAtomKokkos()
{
  if (array_atom != nullptr) {
    delete[] array_atom;
  }
  //data is managed by PairMLIAP so we should only set it to nullptr
  k_data = nullptr;
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void ComputeMLIAPPropertyAtomKokkos<DeviceType>::init()
{
  // Check if a pair style has been defined
  if (force->pair == nullptr)
    error->all(FLERR,"Compute mliap/property/atom/kk requires a pair style be defined");

  // Check if it is safe downcast to PairMLIAP pair style
  PairMLIAPKokkos<DeviceType>* castedPair = dynamic_cast<PairMLIAPKokkos<DeviceType>*>(force->pair);
  if (castedPair == nullptr)
    error->all(FLERR, "Compute mliap/property/atom/kk requires pair mliap to be active");  

  //Get the pointer
  k_data = castedPair->get_k_data();
  if (k_data == nullptr)
    error->all(FLERR, "Compute mliap/property/atom/kk requires a MLIAPDataKokkos");

  //Register the extra_property compute with the data class
  k_data->register_extra_property(property_name, size_peratom_cols);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void ComputeMLIAPPropertyAtomKokkos<DeviceType>::compute_peratom()
{
  k_data->k_extra_properties.sync_host();

  nlocal = atom->nlocal;
  if (array_atom != nullptr) {
    delete[] array_atom;
  }

  auto extra_property_view = k_data->k_extra_properties.get_host_view(property_name);

  array_atom = new LMP_FLOAT*[nlocal];  // manually allocated

  for (int i = 0; i < nlocal; ++i) {
    array_atom[i] = &extra_property_view(i, 0);  // each row starts at column 0
  }
}

template class ComputeMLIAPPropertyAtomKokkos<LMPDeviceType>;
