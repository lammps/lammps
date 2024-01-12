// clang-format off
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
#ifdef MLIAP_PYTHON

#include "mliap_model_python_kokkos.h"

#include "mliap_data_kokkos.h"
#include "comm.h"
#include "error.h"
#include "utils.h"
#include "mliap_model_python_couple_kokkos.h"
#include "lmppython.h"
#include "python_compat.h"
#include <type_traits>

using namespace LAMMPS_NS;

template<class DeviceType>
MLIAPModelPythonKokkos<DeviceType>::~MLIAPModelPythonKokkos() {
  auto nontemplated_this = static_cast<MLIAPModelPythonKokkosDevice*>((void*)this);
  if (model_loaded)
    MLIAPPYKokkos_unload_model(nontemplated_this);
  model_loaded=false;
}

template<class DeviceType>
MLIAPModelPythonKokkos<DeviceType>::MLIAPModelPythonKokkos(LAMMPS *lmp, char *coefffilename) :
  MLIAPModelPython(lmp,coefffilename,true),
  MLIAPModelKokkos<DeviceType>(lmp, this)
{
  if  (!std::is_same_v<DeviceType,LMPDeviceType> )
    MLIAPModelKokkos<DeviceType>::error->all(FLERR, "MLIAP Kokkos version of the python interface is ONLY available on device");

  model_loaded = 0;
  MLIAPModelKokkos<DeviceType>::python->init();
  PyGILState_STATE gstate = PyGILState_Ensure();

  PyObject *pyMain = PyImport_AddModule("__main__");

  if (!pyMain) {
    PyGILState_Release(gstate);
    MLIAPModelKokkos<DeviceType>::error->all(FLERR, "Could not initialize embedded Python");
  }

  PyObject *coupling_module = PyImport_ImportModule("mliap_model_python_couple_kokkos");

  if (!coupling_module) {
    PyErr_Print();
    PyErr_Clear();
    PyGILState_Release(gstate);
    MLIAPModelKokkos<DeviceType>::error->all(FLERR, "Loading MLIAPPYKokkos coupling module failure.");
  }
  // Recipe from lammps/src/pair_python.cpp :
  // add current directory to PYTHONPATH
  PyObject *py_path = PySys_GetObject((char *) "path");
  PyList_Append(py_path, PY_STRING_FROM_STRING("."));

  // if LAMMPS_POTENTIALS environment variable is set, add it to PYTHONPATH as well
  const char *potentials_path = getenv("LAMMPS_POTENTIALS");
  if (potentials_path != nullptr) {
    PyList_Append(py_path, PY_STRING_FROM_STRING(potentials_path));
  }
  PyGILState_Release(gstate);
  if (coefffilename) read_coeffs(coefffilename);

  if (coefffilename) MLIAPModelKokkos<DeviceType>::set_k_coeffelem();

  nonlinearflag = 1;
}
/* ---------------------------------------------------------------------- */

template<class DeviceType>
void MLIAPModelPythonKokkos<DeviceType>::read_coeffs(char *fname)
{
  PyGILState_STATE gstate = PyGILState_Ensure();
  auto nontemplated_this = static_cast<MLIAPModelPythonKokkosDevice*>((void*)this);
  model_loaded = MLIAPPYKokkos_load_model(nontemplated_this, fname);
  if (PyErr_Occurred()) {
    PyErr_Print();
    PyErr_Clear();
    PyGILState_Release(gstate);
    MLIAPModelKokkos<DeviceType>::error->all(FLERR, "Loading python model failure.");
  }
  PyGILState_Release(gstate);

  if (model_loaded) {
    this->connect_param_counts();
  } else {
    if (MLIAPModelKokkos<DeviceType>::comm->me == 0)
      utils::logmesg(MLIAPModelKokkos<DeviceType>::lmp, "Loading python model deferred.\n");
  }
}

/* ---------------------------------------------------------------------- */
// Finalize loading of the model.
template<class DeviceType>
void MLIAPModelPythonKokkos<DeviceType>::connect_param_counts()
{
  PyGILState_STATE gstate = PyGILState_Ensure();
  auto nontemplated_this = static_cast<MLIAPModelPythonKokkosDevice*>((void*)this);
  nelements = MLIAPPYKokkos_nelements(nontemplated_this);
  nparams = MLIAPPYKokkos_nparams(nontemplated_this);
  ndescriptors = MLIAPPYKokkos_ndescriptors(nontemplated_this);

  if (PyErr_Occurred()) {
    PyErr_Print();
    PyErr_Clear();
    PyGILState_Release(gstate);
    MLIAPModelKokkos<DeviceType>::error->all(FLERR, "Loading python model failure.");
  }
  PyGILState_Release(gstate);
  model_loaded = 1;
  utils::logmesg(MLIAPModelKokkos<DeviceType>::lmp, "Loading python model complete.\n");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void MLIAPModelPythonKokkos<DeviceType>::compute_gradients(class MLIAPData *data)
{
  if (!model_loaded) { MLIAPModelKokkos<DeviceType>::error->all(FLERR, "Model not loaded."); }

  PyGILState_STATE gstate = PyGILState_Ensure();

  auto nontemplated_this = static_cast<MLIAPModelPythonKokkosDevice*>((void*)this);
  auto *kokkos_data = dynamic_cast<MLIAPDataKokkos<DeviceType>*>(data);
  MLIAPDataKokkosDevice raw_data(*kokkos_data);
  MLIAPPYKokkos_compute_gradients(nontemplated_this, &raw_data);
  if (PyErr_Occurred()) {
    PyErr_Print();
    PyErr_Clear();
    PyGILState_Release(gstate);
    MLIAPModelKokkos<DeviceType>::error->all(FLERR, "Running python model failure.");
  }
  PyGILState_Release(gstate);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void MLIAPModelPythonKokkos<DeviceType>::compute_gradgrads(class MLIAPData *data)
{
  MLIAPModelPython::compute_gradgrads(data);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void MLIAPModelPythonKokkos<DeviceType>::compute_force_gradients(class MLIAPData *data)
{
  MLIAPModelPython::compute_force_gradients(data);
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class MLIAPModelPythonKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class MLIAPModelPythonKokkos<LMPHostType>;
#endif
}

#endif
