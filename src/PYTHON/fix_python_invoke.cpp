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

/* ----------------------------------------------------------------------
   Contributing author: Richard Berger (Temple U)
------------------------------------------------------------------------- */

#include "fix_python_invoke.h"

#include "error.h"
#include "lmppython.h"
#include "python_compat.h"
#include "python_utils.h"
#include "modify.h"
#include "update.h"

#include <cstring>
#include <Python.h>   // IWYU pragma: export

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixPythonInvoke::FixPythonInvoke(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 6) error->all(FLERR,"Illegal fix python/invoke command");

  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
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
  PyUtils::GIL lock;

  PyObject *pyMain = PyImport_AddModule("__main__");

  if (!pyMain) {
    PyUtils::Print_Errors();
    error->all(FLERR,"Could not initialize embedded Python");
  }

  char *fname = arg[5];
  pFunc = PyObject_GetAttrString(pyMain, fname);

  if (!pFunc) {
    PyUtils::Print_Errors();
    error->all(FLERR,"Could not find Python function");
  }

  lmpPtr = PyCapsule_New((void *)lmp, nullptr, nullptr);

  // nvalid = next step on which end_of_step or post_force does something
  // add nextvalid() to all computes that store invocation times
  // since we don't know a priori which are invoked by python code
  nvalid = nextvalid();
  modify->addstep_compute_all(nvalid);
}

/* ---------------------------------------------------------------------- */

FixPythonInvoke::~FixPythonInvoke()
{
  PyUtils::GIL lock;
  Py_CLEAR(lmpPtr);
}

/* ---------------------------------------------------------------------- */

int FixPythonInvoke::setmask()
{
  return selected_callback;
}

/* ---------------------------------------------------------------------- */

void FixPythonInvoke::init()
{
  // need to reset nvalid if nvalid < ntimestep b/c minimize was performed

  if (nvalid < update->ntimestep) {
    nvalid = nextvalid();
    modify->addstep_compute_all(nvalid);
  }
}

/* ---------------------------------------------------------------------- */

void FixPythonInvoke::end_of_step()
{
  // python code may invoke computes so wrap with clear/add
  modify->clearstep_compute();

  PyUtils::GIL lock;

  PyObject * result = PyObject_CallFunction((PyObject*)pFunc, (char *)"O", (PyObject*)lmpPtr);

  if (!result) {
    PyUtils::Print_Errors();
    error->all(FLERR,"Fix python/invoke end_of_step() method failed");
  }

  Py_CLEAR(result);

  nvalid = nextvalid();
  modify->addstep_compute(nvalid);
}

/* ---------------------------------------------------------------------- */

void FixPythonInvoke::setup(int vflag)
{
  if (selected_callback == POST_FORCE) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixPythonInvoke::post_force(int vflag)
{
  if (update->ntimestep % nevery != 0) return;

  // python code may invoke computes so wrap with clear/add
  modify->clearstep_compute();

  PyUtils::GIL lock;
  char fmt[] = "Oi";

  PyObject * result = PyObject_CallFunction((PyObject*)pFunc, fmt, (PyObject*)lmpPtr, vflag);

  if (!result) {
    PyUtils::Print_Errors();
    error->all(FLERR,"Fix python/invoke post_force() method failed");
  }

  Py_CLEAR(result);

  nvalid = nextvalid();
  modify->addstep_compute(nvalid);
}

/* ---------------------------------------------------------------------- */

bigint FixPythonInvoke::nextvalid()
{
  return (update->ntimestep/nevery + 1)*nevery;
}
