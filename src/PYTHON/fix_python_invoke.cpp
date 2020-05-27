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

#include "fix_python_invoke.h"
#include <Python.h>   // IWYU pragma: keep
#include <cstring>
#include "force.h"
#include "update.h"
#include "error.h"
#include "lmppython.h"
#include "python_compat.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixPythonInvoke::FixPythonInvoke(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 6) error->all(FLERR,"Illegal fix python/invoke command");

  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery <= 0) error->all(FLERR,"Illegal fix python/invoke command");

  // ensure Python interpreter is initialized
  python->init();

  if (strcmp(arg[4],"post_force") == 0) {
    selected_callback = POST_FORCE;
  } else if (strcmp(arg[4],"end_of_step") == 0) {
    selected_callback = END_OF_STEP;
  } else {
    error->all(FLERR,"Unsupported callback name for fix python/invoke");
  }

  // get Python function
  PyGILState_STATE gstate = PyGILState_Ensure();

  PyObject * pyMain = PyImport_AddModule("__main__");

  if (!pyMain) {
    PyGILState_Release(gstate);
    error->all(FLERR,"Could not initialize embedded Python");
  }

  char * fname = arg[5];
  pFunc = PyObject_GetAttrString(pyMain, fname);

  if (!pFunc) {
    PyGILState_Release(gstate);
    error->all(FLERR,"Could not find Python function");
  }

  PyGILState_Release(gstate);
}

/* ---------------------------------------------------------------------- */

int FixPythonInvoke::setmask()
{
  return selected_callback;
}

/* ---------------------------------------------------------------------- */

void FixPythonInvoke::end_of_step()
{
  PyGILState_STATE gstate = PyGILState_Ensure();

  PyObject * ptr = PY_VOID_POINTER(lmp);
  PyObject * arglist = Py_BuildValue("(O)", ptr);

  PyObject * result = PyEval_CallObject((PyObject*)pFunc, arglist);
  Py_DECREF(arglist);

  PyGILState_Release(gstate);
}

/* ---------------------------------------------------------------------- */

void FixPythonInvoke::post_force(int vflag)
{
  if (update->ntimestep % nevery != 0) return;

  PyGILState_STATE gstate = PyGILState_Ensure();

  PyObject * ptr = PY_VOID_POINTER(lmp);
  PyObject * arglist = Py_BuildValue("(Oi)", ptr, vflag);

  PyObject * result = PyEval_CallObject((PyObject*)pFunc, arglist);
  Py_DECREF(arglist);

  PyGILState_Release(gstate);
}
