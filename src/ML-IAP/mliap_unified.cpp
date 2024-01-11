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

/* ----------------------------------------------------------------------
   Contributing author: Steven Ray Anaya (LANL)
------------------------------------------------------------------------- */

#ifdef MLIAP_PYTHON

#include "mliap_unified.h"
#include <Python.h>

#include "error.h"
#include "lmppython.h"
#include "memory.h"
#include "mliap_data.h"
#include "mliap_unified_couple.h"
#include "pair_mliap.h"
#include "python_compat.h"
#include "utils.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

MLIAPDummyDescriptor::MLIAPDummyDescriptor(LAMMPS *_lmp) : Pointers(_lmp), MLIAPDescriptor(_lmp) {}

MLIAPDummyDescriptor::~MLIAPDummyDescriptor()
{
  // manually decrement borrowed reference from Python
  Py_DECREF(unified_interface);
}

/* ----------------------------------------------------------------------
   invoke compute_descriptors from Cython interface
   ---------------------------------------------------------------------- */

void MLIAPDummyDescriptor::compute_descriptors(class MLIAPData *data)
{
  PyGILState_STATE gstate = PyGILState_Ensure();
  compute_descriptors_python(unified_interface, data);
  if (PyErr_Occurred()) {
    PyErr_Print();
    PyErr_Clear();
    PyGILState_Release(gstate);
    lmp->error->all(FLERR, "Running mliappy unified compute_descriptors failure.");
  }
  PyGILState_Release(gstate);
}

/* ----------------------------------------------------------------------
   invoke compute_forces from Cython interface
   ---------------------------------------------------------------------- */

void MLIAPDummyDescriptor::compute_forces(class MLIAPData *data)
{
  PyGILState_STATE gstate = PyGILState_Ensure();
  compute_forces_python(unified_interface, data);
  if (PyErr_Occurred()) {
    PyErr_Print();
    PyErr_Clear();
    PyGILState_Release(gstate);
    lmp->error->all(FLERR, "Running mliappy unified compute_forces failure.");
  }
  PyGILState_Release(gstate);
}

// not implemented
void MLIAPDummyDescriptor::compute_force_gradients(class MLIAPData *)
{
  error->all(FLERR, "compute_force_gradients not implemented");
}

// not implemented
void MLIAPDummyDescriptor::compute_descriptor_gradients(class MLIAPData *)
{
  error->all(FLERR, "compute_descriptor_gradients not implemented");
}

void MLIAPDummyDescriptor::init()
{
  memory->create(radelem, nelements, "mliap_dummy_descriptor:radelem");
  for (int ielem = 0; ielem < nelements; ielem++) { radelem[ielem] = 1; }

  double cut;
  cutmax = 0.0;
  memory->create(cutsq, nelements, nelements, "mliap/descriptor/dummy:cutsq");
  memory->create(cutghost, nelements, nelements, "mliap/descriptor/dummy:cutghost");
  for (int ielem = 0; ielem < nelements; ielem++) {
    // rcutfac set from python, is global cutoff for all elements
    cut = 2.0 * radelem[ielem] * rcutfac;
    if (cut > cutmax) cutmax = cut;
    cutsq[ielem][ielem] = cut * cut;
    cutghost[ielem][ielem] = cut * cut;
    for (int jelem = ielem + 1; jelem < nelements; jelem++) {
      cut = (radelem[ielem] + radelem[jelem]) * rcutfac;
      cutsq[ielem][jelem] = cutsq[jelem][ielem] = cut * cut;
      cutghost[ielem][jelem] = cutghost[jelem][ielem] = cut * cut;
    }
  }
}

void MLIAPDummyDescriptor::set_elements(char **elems, int nelems)
{
  nelements = nelems;
  elements = new char *[nelems];
  for (int i = 0; i < nelems; i++) { elements[i] = utils::strdup(elems[i]); }
}

/* ---------------------------------------------------------------------- */

MLIAPDummyModel::MLIAPDummyModel(LAMMPS *lmp, char *coefffilename) : MLIAPModel(lmp, coefffilename)
{
  nonlinearflag = 1;
}

MLIAPDummyModel::~MLIAPDummyModel()
{
  // manually decrement borrowed reference from Python
  Py_DECREF(unified_interface);
}

int MLIAPDummyModel::get_nparams()
{
  return nparams;
}

int MLIAPDummyModel::get_gamma_nnz(class MLIAPData *)
{
  // TODO: get_gamma_nnz
  return 0;
}

/* ----------------------------------------------------------------------
   invoke compute_gradients from Cython interface
   ---------------------------------------------------------------------- */

void MLIAPDummyModel::compute_gradients(class MLIAPData *data)
{
  PyGILState_STATE gstate = PyGILState_Ensure();
  compute_gradients_python(unified_interface, data);
  if (PyErr_Occurred()) {
    PyErr_Print();
    PyErr_Clear();
    PyGILState_Release(gstate);
    lmp->error->all(FLERR, "Running mliappy unified compute_gradients failure.");
  }
  PyGILState_Release(gstate);
}

// not implemented
void MLIAPDummyModel::compute_gradgrads(class MLIAPData *)
{
  error->all(FLERR, "compute_gradgrads not implemented");
}

// not implemented
void MLIAPDummyModel::compute_force_gradients(class MLIAPData *)
{
  error->all(FLERR, "compute_force_gradients not implemented");
}

/* ----------------------------------------------------------------------
   memory usage unclear due to Cython/Python implementation
   ---------------------------------------------------------------------- */

double MLIAPDummyModel::memory_usage()
{
  // TODO: implement memory usage in Cython(?)
  return 0;
}

// not implemented
void MLIAPDummyModel::read_coeffs(char *)
{
  error->all(FLERR, "read_coeffs not implemented");
}

/* ----------------------------------------------------------------------
   build the unified interface object, connect to dummy model and descriptor
   ---------------------------------------------------------------------- */

MLIAPBuildUnified_t LAMMPS_NS::build_unified(char *unified_fname, MLIAPData *data, LAMMPS *lmp,
                                             char *coefffilename)
{
  lmp->python->init();
  PyGILState_STATE gstate = PyGILState_Ensure();

  PyObject *pyMain = PyImport_AddModule("__main__");

  if (!pyMain) {
    PyGILState_Release(gstate);
    lmp->error->all(FLERR, "Could not initialize embedded Python");
  }

  PyImport_ImportModule("mliap_unified_couple");

  if (PyErr_Occurred()) {
    PyErr_Print();
    PyErr_Clear();
    PyGILState_Release(gstate);
    lmp->error->all(FLERR, "Loading mliappy unified module failure.");
  }

  // Connect dummy model, dummy descriptor, data to Python unified
  MLIAPDummyModel *model = new MLIAPDummyModel(lmp, coefffilename);
  MLIAPDummyDescriptor *descriptor = new MLIAPDummyDescriptor(lmp);

  PyObject *unified_interface = mliap_unified_connect(unified_fname, model, descriptor);
  if (PyErr_Occurred()) {
    PyErr_Print();
    PyErr_Clear();
    PyGILState_Release(gstate);
    lmp->error->all(FLERR, "Running mliappy unified module failure.");
  }

  // Borrowed references must be manually incremented
  model->unified_interface = unified_interface;
  Py_INCREF(unified_interface);
  descriptor->unified_interface = unified_interface;
  Py_INCREF(unified_interface);

  PyGILState_Release(gstate);

  MLIAPBuildUnified_t build = {data, descriptor, model};
  return build;
}

/* ---------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   set energy for ij atom pairs
   ---------------------------------------------------------------------- */

void LAMMPS_NS::update_pair_energy(MLIAPData *data, double *eij)
{
  double e_total = 0.0;
  const auto nlistatoms = data->nlistatoms;
  const auto nlocal = data->nlocal;
  for (int ii = 0; ii < nlistatoms; ii++) data->eatoms[ii] = 0;

  for (int ii = 0; ii < data->npairs; ii++) {
    int i = data->pair_i[ii];
    double e = 0.5 * eij[ii];

    // must not count any contribution where i is not a local atom
    if (i < nlocal) {
      data->eatoms[i] += e;
      e_total += e;
    }
  }
  data->energy = e_total;
}

/* ----------------------------------------------------------------------
   set forces for ij atom pairs
   ---------------------------------------------------------------------- */

void LAMMPS_NS::update_pair_forces(MLIAPData *data, double *fij)
{
  //Bugfix: need to account for Null atoms in local atoms
  //const auto nlistatoms = data->nlistatoms;
  const auto nlocal = data->nlocal;
  double **f = data->f;
  for (int ii = 0; ii < data->npairs; ii++) {
    int ii3 = ii * 3;
    int i = data->pair_i[ii];
    int j = data->jatoms[ii];

    // must not count any contribution where i is not a local atom
    if (i < nlocal) {
      f[i][0] += fij[ii3];
      f[i][1] += fij[ii3 + 1];
      f[i][2] += fij[ii3 + 2];
      f[j][0] -= fij[ii3];
      f[j][1] -= fij[ii3 + 1];
      f[j][2] -= fij[ii3 + 2];

      if (data->vflag) data->pairmliap->v_tally(i, j, &fij[ii3], data->rij[ii]);
    }
  }
}

#endif
