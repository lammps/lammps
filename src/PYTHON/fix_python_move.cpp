// clang-format off
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
   Contributing author: Richard Berger (Temple U)
------------------------------------------------------------------------- */

#include "fix_python_move.h"

#include "error.h"
#include "lmppython.h"
#include "python_compat.h"
#include "python_utils.h"

#include <Python.h>   // IWYU pragma: export

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixPythonMove::FixPythonMove(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  dynamic_group_allow = 1;
  time_integrate = 1;

  python->init();

  py_move = nullptr;

  PyUtils::GIL lock;

  // add current directory to PYTHONPATH
  PyObject *py_path = PySys_GetObject((char *)"path");
  PyList_Append(py_path, PY_STRING_FROM_STRING("."));


  // create integrator instance
  std::string full_cls_name = arg[3];
  std::string module_name = "__main__";
  std::string cls_name = full_cls_name;
  size_t lastpos = full_cls_name.rfind(".");

  if (lastpos != std::string::npos) {
    module_name = full_cls_name.substr(0, lastpos);
    cls_name = full_cls_name.substr(lastpos+1);
  }

  PyObject *pModule = PyImport_ImportModule(module_name.c_str());
  if (!pModule) {
    PyUtils::Print_Errors();
    error->all(FLERR,"Loading python integrator module failure");
  }

  // create LAMMPS atom type to potential file type mapping in python class
  // by calling 'lammps_pair_style.map_coeff(name,type)'

  PyObject *py_move_type = PyObject_GetAttrString(pModule, cls_name.c_str());
  if (!py_move_type) {
    PyUtils::Print_Errors();
    error->all(FLERR,"Could not find integrator class {} in module {}", cls_name, module_name);
  }

  PyObject *ptr = PY_VOID_POINTER(lmp);
  PyObject *py_move_obj = PyObject_CallFunction(py_move_type, (char *)"O", ptr);
  Py_CLEAR(ptr);

  if (!py_move_obj) {
    PyUtils::Print_Errors();
    error->all(FLERR,"Could not instantiate instance of integrator class'");
  }

  // check object interface
  py_move = (void *) py_move_obj;
}

/* ---------------------------------------------------------------------- */

FixPythonMove::~FixPythonMove()
{
  PyUtils::GIL lock;
  Py_CLEAR(py_move);
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
  PyUtils::GIL lock;
  PyObject * result = PyObject_CallMethod((PyObject *)py_move, (char *)"init", nullptr);

  if (!result) {
    PyUtils::Print_Errors();
    error->all(FLERR,"Fix python/move init() method failed");
  }
  Py_CLEAR(result);
}

/* ---------------------------------------------------------------------- */

void FixPythonMove::initial_integrate(int vflag)
{
  PyUtils::GIL lock;
  PyObject * result = PyObject_CallMethod((PyObject*)py_move, (char *)"initial_integrate", (char *)"i", vflag);

  if (!result) {
    PyUtils::Print_Errors();
    error->all(FLERR,"Fix python/move initial_integrate() method failed");
  }
  Py_CLEAR(result);
}

/* ---------------------------------------------------------------------- */

void FixPythonMove::final_integrate()
{
  PyUtils::GIL lock;
  PyObject * result = PyObject_CallMethod((PyObject*)py_move, (char *)"final_integrate", nullptr);

  if (!result) {
    PyUtils::Print_Errors();
    error->all(FLERR,"Fix python/move final_integrate() method failed");
  }
  Py_CLEAR(result);
}

/* ---------------------------------------------------------------------- */

void FixPythonMove::initial_integrate_respa(int vflag, int ilevel, int iloop)
{
  PyUtils::GIL lock;
  PyObject * result = PyObject_CallMethod((PyObject*)py_move, (char *)"initial_integrate_respa", (char *)"iii", vflag, ilevel, iloop);

  if (!result) {
    PyUtils::Print_Errors();
    error->all(FLERR,"Fix python/move initial_integrate_respa() method failed");
  }
  Py_CLEAR(result);
}

/* ---------------------------------------------------------------------- */

void FixPythonMove::final_integrate_respa(int ilevel, int iloop)
{
  PyUtils::GIL lock;
  PyObject * result = PyObject_CallMethod((PyObject*)py_move, (char *)"final_integrate_respa", (char *)"ii", ilevel, iloop);

  if (!result) {
    PyUtils::Print_Errors();
    error->all(FLERR,"Fix python/move final_integrate_respa() method failed");
  }
  Py_CLEAR(result);
}

/* ---------------------------------------------------------------------- */

void FixPythonMove::reset_dt()
{
  PyUtils::GIL lock;
  PyObject * result = PyObject_CallMethod((PyObject*)py_move, (char *)"reset_dt", nullptr);

  if (!result) {
    PyUtils::Print_Errors();
    error->all(FLERR,"Fix python/move reset_dt() method failed");
  }
  Py_CLEAR(result);
}
