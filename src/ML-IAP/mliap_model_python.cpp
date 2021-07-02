/* ----------------------------------------------------------------------
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
   Contributing author: Nicholas Lubbers (LANL)
------------------------------------------------------------------------- */

#ifdef MLIAP_PYTHON

#include "mliap_model_python.h"

#include "error.h"
#include "lmppython.h"
#include "mliap_data.h"
#include "mliap_model_python_couple.h"
#include "pair_mliap.h"
#include "python_compat.h"
#include "utils.h"

#include <Python.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

MLIAPModelPython::MLIAPModelPython(LAMMPS *lmp, char *coefffilename) :
    MLIAPModel(lmp, coefffilename)
{
  model_loaded = 0;
  python->init();
  PyGILState_STATE gstate = PyGILState_Ensure();

  PyObject *pyMain = PyImport_AddModule("__main__");

  if (!pyMain) {
    PyGILState_Release(gstate);
    error->all(FLERR, "Could not initialize embedded Python");
  }

  PyObject *coupling_module = PyImport_ImportModule("mliap_model_python_couple");

  if (!coupling_module) {
    PyErr_Print();
    PyErr_Clear();
    PyGILState_Release(gstate);
    error->all(FLERR, "Loading MLIAPPY coupling module failure.");
  }
  // Recipe from lammps/src/pair_python.cpp :
  // add current directory to PYTHONPATH
  PyObject *py_path = PySys_GetObject((char *) "path");
  PyList_Append(py_path, PY_STRING_FROM_STRING("."));

  // if LAMMPS_POTENTIALS environment variable is set, add it to PYTHONPATH as well
  const char *potentials_path = getenv("LAMMPS_POTENTIALS");
  if (potentials_path != NULL) { PyList_Append(py_path, PY_STRING_FROM_STRING(potentials_path)); }
  PyGILState_Release(gstate);

  if (coefffilename) read_coeffs(coefffilename);

  nonlinearflag = 1;
}

/* ---------------------------------------------------------------------- */

MLIAPModelPython::~MLIAPModelPython()
{
  MLIAPPY_unload_model(this);
}

/* ----------------------------------------------------------------------
   get number of parameters
   ---------------------------------------------------------------------- */

int MLIAPModelPython::get_nparams()
{
  return nparams;
}

void MLIAPModelPython::read_coeffs(char *fname)
{
  PyGILState_STATE gstate = PyGILState_Ensure();

  int loaded = MLIAPPY_load_model(this, fname);
  if (PyErr_Occurred()) {
    PyErr_Print();
    PyErr_Clear();
    PyGILState_Release(gstate);
    error->all(FLERR, "Loading python model failure.");
  }
  PyGILState_Release(gstate);

  if (loaded) {
    this->connect_param_counts();
  } else {
    utils::logmesg(lmp, "Loading python model deferred.\n");
  }
}

// Finalize loading of the model.
void MLIAPModelPython::connect_param_counts()
{
  PyGILState_STATE gstate = PyGILState_Ensure();
  nelements = MLIAPPY_nelements(this);
  nparams = MLIAPPY_nparams(this);
  ndescriptors = MLIAPPY_ndescriptors(this);

  if (PyErr_Occurred()) {
    PyErr_Print();
    PyErr_Clear();
    PyGILState_Release(gstate);
    error->all(FLERR, "Loading python model failure.");
  }
  PyGILState_Release(gstate);
  model_loaded = 1;
  utils::logmesg(lmp, "Loading python model complete.\n");
}

/* ----------------------------------------------------------------------
   Calculate model gradients w.r.t descriptors
   for each atom beta_i = dE(B_i)/dB_i
   ---------------------------------------------------------------------- */

void MLIAPModelPython::compute_gradients(MLIAPData *data)
{
  if (not model_loaded) { error->all(FLERR, "Model not loaded."); }

  PyGILState_STATE gstate = PyGILState_Ensure();
  MLIAPPY_compute_gradients(this, data);
  if (PyErr_Occurred()) {
    PyErr_Print();
    PyErr_Clear();
    PyGILState_Release(gstate);
    error->all(FLERR, "Running python model failure.");
  }
  PyGILState_Release(gstate);
}

/* ----------------------------------------------------------------------
   Calculate model double gradients w.r.t descriptors and parameters
   for each atom energy gamma_lk = d2E(B)/dB_k/dsigma_l,
   where sigma_l is a parameter, B_k a descriptor,
   and atom subscript i is omitted

   gamma is in CSR format:
      nnz = number of non-zero values
      gamma_row_index[inz] = l indices, 0 <= l < nparams
      gamma_col_indexiinz] = k indices, 0 <= k < ndescriptors
      gamma[i][inz] = non-zero values, 0 <= inz < nnz

   egradient is derivative of energy w.r.t. parameters
   ---------------------------------------------------------------------- */

void MLIAPModelPython::compute_gradgrads(class MLIAPData *)
{
  error->all(FLERR, "compute_gradgrads not implemented");
}

/* ----------------------------------------------------------------------
   calculate gradients of forces w.r.t. parameters
   egradient is derivative of energy w.r.t. parameters
   ---------------------------------------------------------------------- */

void MLIAPModelPython::compute_force_gradients(class MLIAPData *)
{
  error->all(FLERR, "compute_force_gradients not implemented");
}

/* ----------------------------------------------------------------------
   count the number of non-zero entries in gamma matrix
   ---------------------------------------------------------------------- */

int MLIAPModelPython::get_gamma_nnz(class MLIAPData *)
{
  // todo: get_gamma_nnz
  return 0;
}

double MLIAPModelPython::memory_usage()
{
  // todo: get approximate memory usage in coupling code.
  return 0;
}

#endif
