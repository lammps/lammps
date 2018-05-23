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
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "fix_python_move.h"
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

FixPythonMove::FixPythonMove(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  dynamic_group_allow = 1;
  time_integrate = 1;

  python->init();

  py_move = NULL;

  PyGILState_STATE gstate = PyGILState_Ensure();

  // add current directory to PYTHONPATH
  PyObject * py_path = PySys_GetObject((char *)"path");
  PyList_Append(py_path, PY_STRING_FROM_STRING("."));


  // create integrator instance
  char * full_cls_name = arg[3];
  char * lastpos = strrchr(full_cls_name, '.');

  if (lastpos == NULL) {
    error->all(FLERR,"Fix python/integrate requires fully qualified class name");
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

  PyObject *py_move_type = PyObject_GetAttrString(pModule, cls_name);
  if (!py_move_type) {
    PyErr_Print();
    PyErr_Clear();
    PyGILState_Release(gstate);
    error->all(FLERR,"Could not find integrator class in module'");
  }

  delete [] module_name;
  delete [] cls_name;

  PyObject * ptr = PY_VOID_POINTER(lmp);
  PyObject * arglist = Py_BuildValue("(O)", ptr);
  PyObject * py_move_obj = PyObject_CallObject(py_move_type, arglist);
  Py_DECREF(arglist);

  if (!py_move_obj) {
    PyErr_Print();
    PyErr_Clear();
    PyGILState_Release(gstate);
    error->all(FLERR,"Could not instantiate instance of integrator class'");
  }

  // check object interface
  py_move = (void *) py_move_obj;

  PyGILState_Release(gstate);
}

/* ---------------------------------------------------------------------- */

FixPythonMove::~FixPythonMove()
{
  PyGILState_STATE gstate = PyGILState_Ensure();
  if(py_move) Py_DECREF((PyObject*) py_move);
  PyGILState_Release(gstate);
}

/* ---------------------------------------------------------------------- */

int FixPythonMove::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= FINAL_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPythonMove::init()
{
  PyGILState_STATE gstate = PyGILState_Ensure();
  PyObject *py_move_obj = (PyObject *) py_move;
  PyObject *py_init = PyObject_GetAttrString(py_move_obj,(char *)"init");
  if (!py_init) {
    PyErr_Print();
    PyErr_Clear();
    PyGILState_Release(gstate);
    error->all(FLERR,"Could not find 'init' method'");
  }
  PyObject * result = PyEval_CallObject(py_init, NULL);
  PyGILState_Release(gstate);
}

/* ---------------------------------------------------------------------- */

void FixPythonMove::initial_integrate(int vflag)
{
  PyGILState_STATE gstate = PyGILState_Ensure();
  PyObject *py_move_obj = (PyObject *) py_move;
  PyObject *py_initial_integrate = PyObject_GetAttrString(py_move_obj,"initial_integrate");
  if (!py_initial_integrate) {
    PyErr_Print();
    PyErr_Clear();
    PyGILState_Release(gstate);
    error->all(FLERR,"Could not find 'initial_integrate' method'");
  }
  PyObject * arglist = Py_BuildValue("(i)", vflag);
  PyObject * result = PyEval_CallObject(py_initial_integrate, arglist);
  Py_DECREF(arglist);
  PyGILState_Release(gstate);
}

/* ---------------------------------------------------------------------- */

void FixPythonMove::final_integrate()
{
  PyGILState_STATE gstate = PyGILState_Ensure();
  PyObject *py_move_obj = (PyObject *) py_move;
  PyObject *py_final_integrate = PyObject_GetAttrString(py_move_obj,"final_integrate");
  if (!py_final_integrate) {
    PyErr_Print();
    PyErr_Clear();
    PyGILState_Release(gstate);
    error->all(FLERR,"Could not find 'final_integrate' method'");
  }
  PyObject * result = PyEval_CallObject(py_final_integrate, NULL);
  PyGILState_Release(gstate);
}

/* ---------------------------------------------------------------------- */

void FixPythonMove::initial_integrate_respa(int vflag, int ilevel, int iloop)
{
  PyGILState_STATE gstate = PyGILState_Ensure();
  PyObject *py_move_obj = (PyObject *) py_move;
  PyObject *py_initial_integrate_respa = PyObject_GetAttrString(py_move_obj,"initial_integrate_respa");
  if (!py_initial_integrate_respa) {
    PyErr_Print();
    PyErr_Clear();
    PyGILState_Release(gstate);
    error->all(FLERR,"Could not find 'initial_integrate_respa' method'");
  }
  PyObject * arglist = Py_BuildValue("(iii)", vflag, ilevel, iloop);
  PyObject * result = PyEval_CallObject(py_initial_integrate_respa, arglist);
  Py_DECREF(arglist);
  PyGILState_Release(gstate);
}

/* ---------------------------------------------------------------------- */

void FixPythonMove::final_integrate_respa(int ilevel, int iloop)
{
  PyGILState_STATE gstate = PyGILState_Ensure();
  PyObject *py_move_obj = (PyObject *) py_move;
  PyObject *py_final_integrate_respa = PyObject_GetAttrString(py_move_obj,"final_integrate_respa");
  if (!py_final_integrate_respa) {
    PyErr_Print();
    PyErr_Clear();
    PyGILState_Release(gstate);
    error->all(FLERR,"Could not find 'final_integrate_respa' method'");
  }
  PyObject * arglist = Py_BuildValue("(ii)", ilevel, iloop);
  PyObject * result = PyEval_CallObject(py_final_integrate_respa, arglist);
  Py_DECREF(arglist);
  PyGILState_Release(gstate);
}

/* ---------------------------------------------------------------------- */

void FixPythonMove::reset_dt()
{
  PyGILState_STATE gstate = PyGILState_Ensure();
  PyObject *py_move_obj = (PyObject *) py_move;
  PyObject *py_reset_dt = PyObject_GetAttrString(py_move_obj,"reset_dt");
  if (!py_reset_dt) {
    PyErr_Print();
    PyErr_Clear();
    PyGILState_Release(gstate);
    error->all(FLERR,"Could not find 'reset_dt' method'");
  }
  PyObject * result = PyEval_CallObject(py_reset_dt, NULL);
  PyGILState_Release(gstate);
}
