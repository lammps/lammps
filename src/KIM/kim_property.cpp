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
   Contributing author: Yaser Afshar (UMN)
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the Free
   Software Foundation; either version 2 of the License, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
   more details.

   You should have received a copy of the GNU General Public License along with
   this program; if not, see <https://www.gnu.org/licenses>.

   Linking LAMMPS statically or dynamically with other modules is making a
   combined work based on LAMMPS. Thus, the terms and conditions of the GNU
   General Public License cover the whole combination.

   In addition, as a special exception, the copyright holders of LAMMPS give
   you permission to combine LAMMPS with free software programs or libraries
   that are released under the GNU LGPL and with code included in the standard
   release of the "kim-api" under the CDDL (or modified versions of such code,
   with unchanged license). You may copy and distribute such a system following
   the terms of the GNU GPL for LAMMPS and the licenses of the other code
   concerned, provided that you include the source code of that other code
   when and as the GNU GPL requires distribution of source code.

   Note that people who make modified versions of LAMMPS are not obligated to
   grant this special exception for their modified versions; it is their choice
   whether to do so. The GNU General Public License gives permission to release
   a modified version without this exception; this exception also makes it
   possible to release a modified version which carries forward this exception.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Designed for use with the kim-api-2.1.0 (and newer) package
------------------------------------------------------------------------- */

#if LMP_PYTHON
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#endif

#include "kim_property.h"

#include "comm.h"
#include "input.h"
#include "variable.h"
#include "utils.h"
#include "error.h"

#include <string>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

kimProperty::kimProperty(LAMMPS *lmp) : Pointers(lmp)
{
#if LMP_PYTHON
#if PY_MAJOR_VERSION != 3
  error->all(FLERR, "Invalid Python version.\n"
                    "The kim-property Python package requires Python "
                    "3 >= 3.6 support.");
#endif
  // one-time initialization of Python interpreter
  if (!Py_IsInitialized()) {
    Py_Initialize();
    PyEval_InitThreads();
  }
#else
  error->all(FLERR, "Error Python support missing! Compile with PYTHON "
                    "package installed!");
#endif // LMP_PYTHON
}

void kimProperty::command(int narg, char **arg)
{
#if LMP_PYTHON
#if PY_MAJOR_VERSION == 3
  if (narg < 2)
    error->all(FLERR, "Invalid `kim_property` command.");

  if (!(strcmp(arg[0], "create") == 0) &&
      !(strcmp(arg[0], "destroy") == 0) &&
      !(strcmp(arg[0], "modify") == 0) &&
      !(strcmp(arg[0], "remove") == 0) &&
      !(strcmp(arg[0], "dump") == 0)) {
    std::string msg("Error incorrect arguments in kim_property command.\n");
    msg += "`kim_property create/destroy/modify/remove/dump` ";
    msg += "is mandatory.";
    error->all(FLERR, msg.c_str());
  }

  if (comm->me == 0) {
    std::string msg;
    msg = "#=== kim-property ===========================================\n";
    input->write_echo(msg.c_str());
  }

  // Get the kim_str ptr to the data associated with a kim_property_str
  // variable
  char *kim_str =
    input->variable->retrieve(const_cast<char *>("kim_property_str"));

  char **kim_str_cmd = new char *[3];
  kim_str_cmd[0] = const_cast<char *>("kim_property_str");

  PyGILState_STATE gstate = PyGILState_Ensure();

  // kim_property module
  PyObject *kim_property = NULL;

  // import kim_property
  {
    PyObject *obj = PyUnicode_FromString("kim_property");
    if (!obj) {
      PyGILState_Release(gstate);
      error->all(FLERR, "Error creating a `PyObject`!");
    }

    kim_property = PyImport_Import(obj);
    if (!kim_property) {
      PyGILState_Release(gstate);
      error->all(FLERR, "Error unable to import Python `kim_property` module!"
                        "\nkim-property Python package can be installed "
                        "with pip:\n`pip install kim-property`\n"
                        "For further information, please refer to the "
                        "to the installation instructions at :\n"
                        "https://github.com/openkim/kim-property");
    }

    // Decrementing of the reference count
    Py_XDECREF(obj);
  }

  // kim_property create 1 cohesive-potential-energy-cubic-crystal
  if (strcmp(arg[0], "create") == 0) {
    if (narg != 3) {
      PyGILState_Release(gstate);
      error->all(FLERR, "Error invalid `kim_property create` command.");
    }

    int const ID = utils::inumeric(FLERR, arg[1], true, lmp);

    Py_ssize_t const nSize = (kim_str ? 3 : 2);

    // Python function
    // This is the equivalent of the Python expression:
    // kim_property.kim_property_create.
    PyObject *pFunc =
      PyObject_GetAttrString(kim_property, "kim_property_create");
    if (!pFunc) {
      PyGILState_Release(gstate);
      error->all(FLERR, "Error unable to get an attribute named "
                        "`kim_property_create` from a kim_property object!");
    }

    // Decrementing of the reference count
    Py_XDECREF(kim_property);

    // create Python tuple of input arguments
    PyObject *pArgs = PyTuple_New(nSize);
    if (!pArgs) {
      PyGILState_Release(gstate);
      error->all(FLERR, "Error could not create Python function arguments.");
    }

    // Python object to set the tuple
    // Property ID
    PyObject *pValue = PyLong_FromLong(ID);
    PyTuple_SetItem(pArgs, 0, pValue);

    // Property name
    pValue = PyUnicode_FromString(arg[2]);
    PyTuple_SetItem(pArgs, 1, pValue);

    if (nSize == 3) {
      pValue = PyUnicode_FromString(kim_str);
      PyTuple_SetItem(pArgs, 2, pValue);
    }

    // call the Python kim_property_create function
    // error check with one() since only some procs may fail
    pValue = PyObject_CallObject(pFunc, pArgs);
    if (!pValue) {
      PyErr_Print();
      PyGILState_Release(gstate);
      error->one(FLERR, "Error Python `kim_property_create` function "
                        "evaluation failed!");
    }

    // Python function returned a string value
    const char *pystr = PyUnicode_AsUTF8(pValue);

    kim_str_cmd[2] = const_cast<char *>(pystr);

    if (kim_str) {
      input->variable->set_string(kim_str_cmd[0], kim_str_cmd[2]);
    } else {
      kim_str_cmd[1] = const_cast<char *>("string");
      input->variable->set(3, kim_str_cmd);
    }

    Py_XDECREF(pArgs);
    Py_XDECREF(pFunc);
    Py_XDECREF(pValue);
  }
  else if (strcmp(arg[0], "destroy") == 0) {
    if (narg != 2) {
      PyGILState_Release(gstate);
      error->all(FLERR, "Error invalid `kim_property destroy` command.");
    }

    if (!kim_str) {
      PyGILState_Release(gstate);
      return;
    }

    int const ID = utils::inumeric(FLERR, arg[1], true, lmp);

    // Python function
    // This is the equivalent of the Python expression kim_property.kim_property_destroy
    PyObject *pFunc =
      PyObject_GetAttrString(kim_property, "kim_property_destroy");
    if (!pFunc) {
      PyGILState_Release(gstate);
      error->all(FLERR, "Error unable to get an attribute named "
                        "`kim_property_destroy` from a kim_property object!");
    }

    // Decrementing of the reference count
    Py_XDECREF(kim_property);

    // create Python tuple of input arguments
    PyObject *pArgs = PyTuple_New(2);
    if (!pArgs) {
      PyGILState_Release(gstate);
      error->all(FLERR, "Error could not create Python function arguments.");
    }

    // Python object to set the tuple
    PyObject *pValue = PyUnicode_FromString(kim_str);
    PyTuple_SetItem(pArgs, 0, pValue);

    pValue = PyLong_FromLong(ID);
    PyTuple_SetItem(pArgs, 1, pValue);

    // call the Python kim_property_destroy function
    // error check with one() since only some procs may fail
    pValue = PyObject_CallObject(pFunc, pArgs);
    if (!pValue) {
      PyErr_Print();
      PyGILState_Release(gstate);
      error->one(FLERR, "Error Python `kim_property_destroy` function "
                        "evaluation failed!");
    }

    // Python function returned a string value
    const char *pystr = PyUnicode_AsUTF8(pValue);

    kim_str_cmd[2] = const_cast<char *>(pystr);

    input->variable->set_string(kim_str_cmd[0], kim_str_cmd[2]);

    Py_XDECREF(pArgs);
    Py_XDECREF(pFunc);
    Py_XDECREF(pValue);
  }
  else if (strcmp(arg[0], "modify") == 0) {
    if (narg < 6) {
      PyGILState_Release(gstate);
      error->all(FLERR, "Error invalid `kim_property modify` command.");
    }

    if (!kim_str) {
      PyGILState_Release(gstate);
      error->all(FLERR, "Error There is no property instance to modify "
                        "the content.");
    }

    int const ID = utils::inumeric(FLERR, arg[1], true, lmp);

    // Python function
    // This is the equivalent of the Python expression
    // kim_property.kim_property_modify
    PyObject *pFunc =
      PyObject_GetAttrString(kim_property, "kim_property_modify");
    if (!pFunc) {
      PyGILState_Release(gstate);
      error->all(FLERR, "Error unable to get an attribute named "
                        "`kim_property_modify` from a kim_property object!");
    }

    // Decrementing of the reference count
    Py_XDECREF(kim_property);

    // create Python tuple of input arguments
    PyObject *pArgs = PyTuple_New(static_cast<Py_ssize_t>(narg));
    if (!pArgs) {
      PyGILState_Release(gstate);
      error->all(FLERR, "Error could not create Python function arguments.");
    }

    // Python object to set the tuple
    PyObject *pValue = PyUnicode_FromString(kim_str);
    PyTuple_SetItem(pArgs, 0, pValue);

    pValue = PyLong_FromLong(ID);
    PyTuple_SetItem(pArgs, 1, pValue);

    for (Py_ssize_t i = 2; i < static_cast<Py_ssize_t>(narg); ++i) {
      pValue = PyUnicode_FromString(arg[i]);
      PyTuple_SetItem(pArgs, i, pValue);
    }

    // call the Python kim_property_modify function
    // error check with one() since only some procs may fail
    pValue = PyObject_CallObject(pFunc, pArgs);
    if (!pValue) {
      PyErr_Print();
      PyGILState_Release(gstate);
      error->one(FLERR, "Error Python `kim_property_modify` function "
                        "evaluation failed!");
    }

    // Python function returned a string value
    const char *pystr = PyUnicode_AsUTF8(pValue);

    kim_str_cmd[2] = const_cast<char *>(pystr);

    input->variable->set_string(kim_str_cmd[0], kim_str_cmd[2]);

    Py_XDECREF(pArgs);
    Py_XDECREF(pFunc);
    Py_XDECREF(pValue);
  }
  else if (strcmp(arg[0], "remove") == 0) {
    if (narg < 4) {
      PyGILState_Release(gstate);
      error->all(FLERR, "Error invalid `kim_property remove` command.");
    }

    if (!kim_str) {
      PyGILState_Release(gstate);
      error->all(FLERR, "Error There is no property instance to remove "
                        "the content.");
    }

    int const ID = utils::inumeric(FLERR, arg[1], true, lmp);

    // Python function
    // This is the equivalent of the Python expression
    // kim_property.kim_property_remove
    PyObject *pFunc =
      PyObject_GetAttrString(kim_property, "kim_property_remove");
    if (!pFunc) {
      PyGILState_Release(gstate);
      error->all(FLERR, "Error unable to get an attribute named "
                        "`kim_property_remove` from a kim_property object!");
    }

    // Decrementing of the reference count
    Py_XDECREF(kim_property);

    // create Python tuple of input arguments
    PyObject *pArgs = PyTuple_New(static_cast<Py_ssize_t>(narg));
    if (!pArgs) {
      PyGILState_Release(gstate);
      error->all(FLERR, "Error could not create Python function arguments.");
    }

    // Python object to set the tuple
    PyObject *pValue = PyUnicode_FromString(kim_str);
    PyTuple_SetItem(pArgs, 0, pValue);

    pValue = PyLong_FromLong(ID);
    PyTuple_SetItem(pArgs, 1, pValue);

    for (Py_ssize_t i = 2; i < static_cast<Py_ssize_t>(narg); ++i) {
      pValue = PyUnicode_FromString(arg[i]);
      PyTuple_SetItem(pArgs, i, pValue);
    }

    // call the Python kim_property_remove function
    // error check with one() since only some procs may fail
    pValue = PyObject_CallObject(pFunc, pArgs);
    if (!pValue) {
      PyErr_Print();
      PyGILState_Release(gstate);
      error->one(FLERR, "Error Python `kim_property_remove` function "
                        "evaluation failed!");
    }

    // Python function returned a string value
    const char *pystr = PyUnicode_AsUTF8(pValue);

    kim_str_cmd[2] = const_cast<char *>(pystr);

    input->variable->set_string(kim_str_cmd[0], kim_str_cmd[2]);

    Py_XDECREF(pArgs);
    Py_XDECREF(pFunc);
    Py_XDECREF(pValue);
  }
  else if (strcmp(arg[0], "dump") == 0) {
    if (narg != 2) {
      PyGILState_Release(gstate);
      error->all(FLERR, "Error invalid `kim_property dump` command.");
    }

    if (!kim_str) {
      PyGILState_Release(gstate);
      error->all(FLERR, "Error There is no property instance to dump "
                        "the content.");
    }

    // Python function
    // This is the equivalent of the Python expression
    // kim_property.kim_property_dump
    PyObject *pFunc =
      PyObject_GetAttrString(kim_property, "kim_property_dump");
    if (!pFunc) {
      PyGILState_Release(gstate);
      error->all(FLERR, "Error unable to get an attribute named "
                        "`kim_property_dump` from a kim_property object!");
    }

    // Decrementing of the reference count
    Py_XDECREF(kim_property);

    // create Python tuple of input arguments
    PyObject *pArgs = PyTuple_New(2);
    if (!pArgs) {
      PyGILState_Release(gstate);
      error->all(FLERR, "Could not create Python function arguments.");
    }

    // Python object to set the tuple
    PyObject *pValue = PyUnicode_FromString(kim_str);
    PyTuple_SetItem(pArgs, 0, pValue);

    pValue = PyUnicode_FromString(arg[1]);
    PyTuple_SetItem(pArgs, 1, pValue);

    if (comm->me == 0) {
      // call the Python kim_property_dump function
      // error check with one() since only some procs may fail
      pValue = PyObject_CallObject(pFunc, pArgs);
      if (!pValue) {
        PyErr_Print();
        PyGILState_Release(gstate);
        error->one(FLERR, "Error Python `kim_property_dump` function "
                          "evaluation failed!");
      }
    }

    // Destroy the variable
    kim_str_cmd[1] = const_cast<char *>("delete");
    input->variable->set(2, kim_str_cmd);

    Py_XDECREF(pArgs);
    Py_XDECREF(pFunc);
    Py_XDECREF(pValue);
  }

  PyGILState_Release(gstate);

  delete[] kim_str_cmd;
#endif // PY_MAJOR_VERSION
#endif // LMP_PYTHON
}
