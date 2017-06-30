/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Richard Berger (Temple U)
------------------------------------------------------------------------- */

#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fix_python_integrate.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "python.h"
#include "error.h"
#include "python_compat.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixPythonIntegrate::FixPythonIntegrate(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  dynamic_group_allow = 1;
  time_integrate = 1;

  python->init();

  py_integrator = NULL;

  PyGILState_STATE gstate = PyGILState_Ensure();

  // add current directory to PYTHONPATH
  PyObject * py_path = PySys_GetObject("path");
  PyList_Append(py_path, PY_STRING_FROM_STRING("."));


  // create integrator instance
  char * full_cls_name = arg[2];
  char * lastpos = strrchr(full_cls_name, '.');

  if (lastpos == NULL) {
    error->all(FLERR,"Python pair style requires fully qualified class name");
  }

  size_t module_name_length = strlen(full_cls_name) - strlen(lastpos);
  size_t cls_name_length = strlen(lastpos)-1;

  char * module_name = new char[module_name_length+1];
  char * cls_name = new char[cls_name_length+1];
  strncpy(module_name, full_cls_name, module_name_length);
  module_name[module_name_length] = 0;

  strcpy(cls_name, lastpos+1);

  PyObject * pModule = PyImport_ImportModule(module_name);
  if (!pModule) {
    PyErr_Print();
    PyErr_Clear();
    PyGILState_Release(gstate);
    error->all(FLERR,"Loading python integrator module failure");
  }

  // create LAMMPS atom type to potential file type mapping in python class
  // by calling 'lammps_pair_style.map_coeff(name,type)'

  PyObject *py_integrator_type = PyObject_GetAttrString(pModule, cls_name);
  if (!py_integrator_type) {
    PyErr_Print();
    PyErr_Clear();
    PyGILState_Release(gstate);
    error->all(FLERR,"Could not find integrator class in module'");
  }

  delete [] module_name;
  delete [] cls_name;

  PyObject * py_integrator_obj = PyObject_CallObject(py_integrator_type, NULL);
  if (!py_integrator_obj) {
    PyErr_Print();
    PyErr_Clear();
    PyGILState_Release(gstate);
    error->all(FLERR,"Could not instantiate instance of integrator class'");
  }

  // check object interface
  py_integrator = (void *) py_integrator_obj;

  PyGILState_Release(gstate);
}

/* ---------------------------------------------------------------------- */

FixPythonIntegrate::~FixPythonIntegrate()
{
  PyGILState_STATE gstate = PyGILState_Ensure();
  if(py_integrator) Py_DECREF((PyObject*) py_integrator);
  PyGILState_Release(gstate);
}

/* ---------------------------------------------------------------------- */

int FixPythonIntegrate::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= FINAL_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPythonIntegrate::init()
{
  PyGILState_STATE gstate = PyGILState_Ensure();
  PyObject *py_integrator_obj = (PyObject *) py_integrator;
  PyObject *py_init = PyObject_GetAttrString(py_integrator_obj,"init");
  if (!py_init) {
    PyErr_Print();
    PyErr_Clear();
    PyGILState_Release(gstate);
    error->all(FLERR,"Could not find 'init' method'");
  }
  PyObject * ptr = PY_VOID_POINTER(lmp);
  PyObject * arglist = Py_BuildValue("(O)", ptr);
  PyObject * result = PyEval_CallObject(py_init, arglist);
  Py_DECREF(arglist);
  PyGILState_Release(gstate);
}

/* ---------------------------------------------------------------------- */

void FixPythonIntegrate::initial_integrate(int vflag)
{
  PyGILState_STATE gstate = PyGILState_Ensure();
  PyObject *py_integrator_obj = (PyObject *) py_integrator;
  PyObject *py_initial_integrate = PyObject_GetAttrString(py_integrator_obj,"initial_integrate");
  if (!py_initial_integrate) {
    PyErr_Print();
    PyErr_Clear();
    PyGILState_Release(gstate);
    error->all(FLERR,"Could not find 'initial_integrate' method'");
  }
  PyObject * ptr = PY_VOID_POINTER(lmp);
  PyObject * arglist = Py_BuildValue("(Oi)", ptr, vflag);
  PyObject * result = PyEval_CallObject(py_initial_integrate, arglist);
  Py_DECREF(arglist);
  PyGILState_Release(gstate);
}

/* ---------------------------------------------------------------------- */

void FixPythonIntegrate::final_integrate()
{
  PyGILState_STATE gstate = PyGILState_Ensure();
  PyObject *py_integrator_obj = (PyObject *) py_integrator;
  PyObject *py_final_integrate = PyObject_GetAttrString(py_integrator_obj,"final_integrate");
  if (!py_final_integrate) {
    PyErr_Print();
    PyErr_Clear();
    PyGILState_Release(gstate);
    error->all(FLERR,"Could not find 'final_integrate' method'");
  }
  PyObject * ptr = PY_VOID_POINTER(lmp);
  PyObject * arglist = Py_BuildValue("(O)", ptr);
  PyObject * result = PyEval_CallObject(py_final_integrate, arglist);
  Py_DECREF(arglist);
  PyGILState_Release(gstate);
}

/* ---------------------------------------------------------------------- */

void FixPythonIntegrate::initial_integrate_respa(int vflag, int ilevel, int iloop)
{
  PyGILState_STATE gstate = PyGILState_Ensure();
  PyObject *py_integrator_obj = (PyObject *) py_integrator;
  PyObject *py_initial_integrate_respa = PyObject_GetAttrString(py_integrator_obj,"initial_integrate_respa");
  if (!py_initial_integrate_respa) {
    PyErr_Print();
    PyErr_Clear();
    PyGILState_Release(gstate);
    error->all(FLERR,"Could not find 'initial_integrate_respa' method'");
  }
  PyObject * ptr = PY_VOID_POINTER(lmp);
  PyObject * arglist = Py_BuildValue("(Oiii)", ptr, vflag, ilevel, iloop);
  PyObject * result = PyEval_CallObject(py_initial_integrate_respa, arglist);
  Py_DECREF(arglist);
  PyGILState_Release(gstate);
}

/* ---------------------------------------------------------------------- */

void FixPythonIntegrate::final_integrate_respa(int ilevel, int iloop)
{
  PyGILState_STATE gstate = PyGILState_Ensure();
  PyObject *py_integrator_obj = (PyObject *) py_integrator;
  PyObject *py_final_integrate_respa = PyObject_GetAttrString(py_integrator_obj,"final_integrate_respa");
  if (!py_final_integrate_respa) {
    PyErr_Print();
    PyErr_Clear();
    PyGILState_Release(gstate);
    error->all(FLERR,"Could not find 'final_integrate_respa' method'");
  }
  PyObject * ptr = PY_VOID_POINTER(lmp);
  PyObject * arglist = Py_BuildValue("(Oii)", ptr, ilevel, iloop);
  PyObject * result = PyEval_CallObject(py_final_integrate_respa, arglist);
  Py_DECREF(arglist);
  PyGILState_Release(gstate);
}

/* ---------------------------------------------------------------------- */

void FixPythonIntegrate::reset_dt()
{
  PyGILState_STATE gstate = PyGILState_Ensure();
  PyObject *py_integrator_obj = (PyObject *) py_integrator;
  PyObject *py_reset_dt = PyObject_GetAttrString(py_integrator_obj,"reset_dt");
  if (!py_reset_dt) {
    PyErr_Print();
    PyErr_Clear();
    PyGILState_Release(gstate);
    error->all(FLERR,"Could not find 'reset_dt' method'");
  }
  PyObject * ptr = PY_VOID_POINTER(lmp);
  PyObject * arglist = Py_BuildValue("(O)", ptr);
  PyObject * result = PyEval_CallObject(py_reset_dt, arglist);
  Py_DECREF(arglist);
  PyGILState_Release(gstate);
}
