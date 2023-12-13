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
   Contributing authors: Richard Berger and Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "python_impl.h"

#include "error.h"
#include "input.h"
#include "memory.h"
#include "python_compat.h"
#include "python_utils.h"
#include "variable.h"

#include <Python.h>    // IWYU pragma: export
#include <cstring>

#ifdef MLIAP_PYTHON
#include "mliap_model_python.h"
#if defined(__PYX_EXTERN_C) && !defined(CYTHON_EXTERN_C)
#undef __PYX_EXTERN_C
#endif
#include "mliap_unified.h"
// The above should somehow really be included in the next file.
// We could get around this with cython --capi-reexport-cincludes
// However, that exposes -too many- headers.
#include "mliap_model_python_couple.h"
#if defined(__PYX_EXTERN_C) && !defined(CYTHON_EXTERN_C)
#undef __PYX_EXTERN_C
#endif
#include "mliap_unified_couple.h"
#ifdef LMP_KOKKOS
#include "mliap_model_python_kokkos.h"
#if defined(__PYX_EXTERN_C) && !defined(CYTHON_EXTERN_C)
#undef __PYX_EXTERN_C
#endif
#include "mliap_unified_kokkos.h"
// The above should somehow really be included in the next file.
// We could get around this with cython --capi-reexport-cincludes
// However, that exposes -too many- headers.
#include "mliap_model_python_couple_kokkos.h"
#if defined(__PYX_EXTERN_C) && !defined(CYTHON_EXTERN_C)
#undef __PYX_EXTERN_C
#endif
#include "mliap_unified_couple_kokkos.h"


#endif
#endif

using namespace LAMMPS_NS;

enum { NONE, INT, DOUBLE, STRING, PTR };

/* ---------------------------------------------------------------------- */

PythonImpl::PythonImpl(LAMMPS *lmp) : Pointers(lmp)
{
  // pfuncs stores interface info for each Python function

  nfunc = 0;
  pfuncs = nullptr;

#if PY_MAJOR_VERSION >= 3 && !defined(Py_LIMITED_API)
  // check for PYTHONUNBUFFERED environment variable
  const char *PYTHONUNBUFFERED = getenv("PYTHONUNBUFFERED");
  // Force the stdout and stderr streams to be unbuffered.
  bool unbuffered = PYTHONUNBUFFERED != nullptr && strcmp(PYTHONUNBUFFERED, "1") == 0;

#if PY_VERSION_HEX >= 0x030800f0
  PyConfig config;
  PyConfig_InitPythonConfig(&config);
  config.buffered_stdio = !unbuffered;
#else
  // Python Global configuration variable
  Py_UnbufferedStdioFlag = unbuffered;
#endif
#endif

#ifdef MLIAP_PYTHON
  // Inform python intialization scheme of the mliappy module.
  // This -must- happen before python is initialized.
  int err = PyImport_AppendInittab("mliap_model_python_couple", PyInit_mliap_model_python_couple);
  if (err) error->all(FLERR, "Could not register MLIAPPY embedded python module.");

  err = PyImport_AppendInittab("mliap_unified_couple", PyInit_mliap_unified_couple);
  if (err) error->all(FLERR, "Could not register MLIAPPY unified embedded python module.");
#ifdef LMP_KOKKOS
  // Inform python intialization scheme of the mliappy module.
  // This -must- happen before python is initialized.
  err = PyImport_AppendInittab("mliap_model_python_couple_kokkos", PyInit_mliap_model_python_couple_kokkos);
  if (err) error->all(FLERR, "Could not register MLIAPPY embedded python module.");

  err = PyImport_AppendInittab("mliap_unified_couple_kokkos", PyInit_mliap_unified_couple_kokkos);
  if (err) error->all(FLERR, "Could not register MLIAPPY unified embedded python module.");

#endif
#endif

#if PY_VERSION_HEX >= 0x030800f0 && !defined(Py_LIMITED_API)
  Py_InitializeFromConfig(&config);
  PyConfig_Clear(&config);
#else
  Py_Initialize();
#endif

  // only needed for Python 2.x and Python 3 < 3.7
  // With Python 3.7 this function is now called by Py_Initialize()
  // Deprecated since version 3.9, will be removed in version 3.11
#if PY_VERSION_HEX < 0x030700f0
  if (!PyEval_ThreadsInitialized()) { PyEval_InitThreads(); }
#endif

  PyUtils::GIL lock;

  PyObject *pModule = PyImport_AddModule("__main__");
  if (!pModule) error->all(FLERR, "Could not initialize embedded Python");

  pyMain = (void *) pModule;
}

/* ---------------------------------------------------------------------- */

PythonImpl::~PythonImpl()
{
  if (pyMain) {
    // clean up
    PyUtils::GIL lock;

    for (int i = 0; i < nfunc; i++) {
      delete[] pfuncs[i].name;
      deallocate(i);
      Py_CLEAR(pfuncs[i].pFunc);
    }
  }

  memory->sfree(pfuncs);
}

/* ---------------------------------------------------------------------- */

void PythonImpl::command(int narg, char **arg)
{
  if (narg < 2) utils::missing_cmd_args(FLERR, "python", error);

  // if invoke is only keyword, invoke the previously defined function

  if (narg == 2 && strcmp(arg[1], "invoke") == 0) {
    int ifunc = find(arg[0]);
    if (ifunc < 0) error->all(FLERR, "Python invoke of unknown function: {}", arg[0]);

    char *str = nullptr;
    if (pfuncs[ifunc].noutput) {
      str = input->variable->pythonstyle(pfuncs[ifunc].ovarname, pfuncs[ifunc].name);
      if (!str)
        error->all(FLERR,
                   "Python variable {} does not match variable {} "
                   "registered with Python function {}",
                   arg[0], pfuncs[ifunc].ovarname, pfuncs[ifunc].name);
    }

    invoke_function(ifunc, str);
    return;
  }

  // if source is only keyword, execute the python code in file

  if ((narg > 1) && (strcmp(arg[0], "source") == 0)) {
    int err = -1;

    if ((narg > 2) && (strcmp(arg[1], "here") == 0)) {
      err = execute_string(arg[2]);
    } else {
      if (platform::file_is_readable(arg[1]))
        err = execute_file(arg[1]);
      else
        error->all(FLERR, "Could not open python source file {} for processing", arg[1]);
    }
    if (err) error->all(FLERR, "Failure in python source command");

    return;
  }

  // parse optional args, invoke is not allowed in this mode

  int ninput = 0;
  int noutput = 0;
  char **istr = nullptr;
  char *ostr = nullptr;
  char *format = nullptr;
  int length_longstr = 0;
  char *pyfile = nullptr;
  char *herestr = nullptr;
  int existflag = 0;

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "input") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "python input", error);
      ninput = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (ninput < 0) error->all(FLERR, "Invalid number of python input arguments: {}", ninput);
      iarg += 2;
      delete[] istr;
      istr = new char *[ninput];
      if (iarg + ninput > narg) utils::missing_cmd_args(FLERR, "python input", error);
      for (int i = 0; i < ninput; i++) istr[i] = arg[iarg + i];
      iarg += ninput;
    } else if (strcmp(arg[iarg], "return") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "python return", error);
      noutput = 1;
      ostr = arg[iarg + 1];
      iarg += 2;
    } else if (strcmp(arg[iarg], "format") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "python format", error);
      format = utils::strdup(arg[iarg + 1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "length") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "python length", error);
      length_longstr = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (length_longstr <= 0) error->all(FLERR, "Invalid python return value length");
      iarg += 2;
    } else if (strcmp(arg[iarg], "file") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "python file", error);
      delete[] pyfile;
      pyfile = utils::strdup(arg[iarg + 1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "here") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "python here", error);
      herestr = arg[iarg + 1];
      iarg += 2;
    } else if (strcmp(arg[iarg], "exists") == 0) {
      existflag = 1;
      iarg++;
    } else
      error->all(FLERR, "Unknown python command keyword: {}", arg[iarg]);
  }

  if (pyfile && herestr)
    error->all(FLERR, "Must not use python 'file' and 'here' keywords at the same time");
  if (pyfile && existflag)
    error->all(FLERR, "Must not use python 'file' and 'exists' keywords at the same time");
  if (herestr && existflag)
    error->all(FLERR, "Must not use python 'here' and 'exists' keywords at the same time");

  // create or overwrite entry in pfuncs vector with name = arg[0]

  int ifunc = create_entry(arg[0], ninput, noutput, length_longstr, istr, ostr, format);

  PyUtils::GIL lock;

  // send Python code to Python interpreter
  // file: read the file via PyRun_SimpleFile()
  // here: process the here string directly
  // exist: do nothing, assume code has already been run

  if (pyfile) {
    FILE *fp = fopen(pyfile, "r");

    if (fp == nullptr) {
      PyUtils::Print_Errors();
      error->all(FLERR, "Could not open Python file: {}", pyfile);
    }

    int err = PyRun_SimpleFile(fp, pyfile);
    if (err) {
      PyUtils::Print_Errors();
      error->all(FLERR, "Could not process Python file: {}", pyfile);
    }
    fclose(fp);

  } else if (herestr) {
    int err = PyRun_SimpleString(herestr);
    if (err) {
      PyUtils::Print_Errors();
      error->all(FLERR, "Could not process Python string: {}", herestr);
    }
  }

  // pFunc = function object for requested function

  auto pModule = (PyObject *) pyMain;
  PyObject *pFunc = PyObject_GetAttrString(pModule, pfuncs[ifunc].name);

  if (!pFunc) {
    PyUtils::Print_Errors();
    error->all(FLERR, "Could not find Python function {}", pfuncs[ifunc].name);
  }

  if (!PyCallable_Check(pFunc)) {
    PyUtils::Print_Errors();
    error->all(FLERR, "Python function {} is not callable", pfuncs[ifunc].name);
  }

  pfuncs[ifunc].pFunc = (void *) pFunc;

  // clean-up input storage

  delete[] istr;
  delete[] format;
  delete[] pyfile;
}

/* ------------------------------------------------------------------ */

void PythonImpl::invoke_function(int ifunc, char *result)
{
  PyUtils::GIL lock;
  PyObject *pValue;
  char *str;

  auto pFunc = (PyObject *) pfuncs[ifunc].pFunc;

  // create Python tuple of input arguments

  int ninput = pfuncs[ifunc].ninput;
  PyObject *pArgs = PyTuple_New(ninput);

  if (!pArgs)
    error->all(FLERR, "Could not prepare arguments for Python function {}", pfuncs[ifunc].name);

  for (int i = 0; i < ninput; i++) {
    int itype = pfuncs[ifunc].itype[i];
    if (itype == INT) {
      if (pfuncs[ifunc].ivarflag[i]) {
        str = input->variable->retrieve(pfuncs[ifunc].svalue[i]);
        if (!str)
          error->all(FLERR, "Could not evaluate Python function {} input variable: {}",
                     pfuncs[ifunc].name, pfuncs[ifunc].svalue[i]);
        pValue = PY_INT_FROM_LONG(atoi(str));
      } else {
        pValue = PY_INT_FROM_LONG(pfuncs[ifunc].ivalue[i]);
      }
    } else if (itype == DOUBLE) {
      if (pfuncs[ifunc].ivarflag[i]) {
        str = input->variable->retrieve(pfuncs[ifunc].svalue[i]);
        if (!str)
          error->all(FLERR, "Could not evaluate Python function {} input variable: {}",
                     pfuncs[ifunc].name, pfuncs[ifunc].svalue[i]);
        pValue = PyFloat_FromDouble(atof(str));
      } else {
        pValue = PyFloat_FromDouble(pfuncs[ifunc].dvalue[i]);
      }
    } else if (itype == STRING) {
      if (pfuncs[ifunc].ivarflag[i]) {
        str = input->variable->retrieve(pfuncs[ifunc].svalue[i]);
        if (!str)
          error->all(FLERR, "Could not evaluate Python function {} input variable: {}",
                     pfuncs[ifunc].name, pfuncs[ifunc].svalue[i]);
        pValue = PY_STRING_FROM_STRING(str);
      } else {
        pValue = PY_STRING_FROM_STRING(pfuncs[ifunc].svalue[i]);
      }
    } else if (itype == PTR) {
      pValue = PY_VOID_POINTER(lmp);
    } else {
      error->all(FLERR, "Unsupported variable type: {}", itype);
    }
    PyTuple_SetItem(pArgs, i, pValue);
  }

  // call the Python function
  // error check with one() since only some procs may fail

  pValue = PyObject_CallObject(pFunc, pArgs);
  Py_CLEAR(pArgs);

  if (!pValue) {
    PyUtils::Print_Errors();
    error->one(FLERR, "Python evaluation of function {} failed", pfuncs[ifunc].name);
  }

  // function returned a value
  // assign it to result string stored by python-style variable
  // or if user specified a length, assign it to longstr

  if (pfuncs[ifunc].noutput) {
    int otype = pfuncs[ifunc].otype;
    if (otype == INT) {
      auto value = fmt::format("{}", PY_INT_AS_LONG(pValue));
      strncpy(result, value.c_str(), Variable::VALUELENGTH - 1);
    } else if (otype == DOUBLE) {
      auto value = fmt::format("{:.15g}", PyFloat_AsDouble(pValue));
      strncpy(result, value.c_str(), Variable::VALUELENGTH - 1);
    } else if (otype == STRING) {
      const char *pystr = PY_STRING_AS_STRING(pValue);
      if (pfuncs[ifunc].longstr)
        strncpy(pfuncs[ifunc].longstr, pystr, pfuncs[ifunc].length_longstr);
      else
        strncpy(result, pystr, Variable::VALUELENGTH - 1);
    }
  }
  Py_CLEAR(pValue);
}

/* ------------------------------------------------------------------ */

int PythonImpl::find(const char *name)
{
  for (int i = 0; i < nfunc; i++)
    if (strcmp(name, pfuncs[i].name) == 0) return i;
  return -1;
}

/* ------------------------------------------------------------------ */

int PythonImpl::variable_match(const char *name, const char *varname, int numeric)
{
  int ifunc = find(name);
  if (ifunc < 0) return -1;
  if (pfuncs[ifunc].noutput == 0) return -2;
  if (strcmp(pfuncs[ifunc].ovarname, varname) != 0) return -3;
  if (numeric && pfuncs[ifunc].otype == STRING) return -4;
  return ifunc;
}

/* ------------------------------------------------------------------ */

char *PythonImpl::long_string(int ifunc)
{
  return pfuncs[ifunc].longstr;
}

/* ------------------------------------------------------------------ */

int PythonImpl::create_entry(char *name, int ninput, int noutput, int length_longstr, char **istr,
                             char *ostr, char *format)
{
  // ifunc = index to entry by name in pfuncs vector, can be old or new
  // free old vectors if overwriting old pfunc

  int ifunc = find(name);

  if (ifunc < 0) {
    ifunc = nfunc;
    nfunc++;
    pfuncs = (PyFunc *) memory->srealloc(pfuncs, nfunc * sizeof(struct PyFunc), "python:pfuncs");
    pfuncs[ifunc].name = utils::strdup(name);
  } else
    deallocate(ifunc);

  pfuncs[ifunc].ninput = ninput;
  pfuncs[ifunc].noutput = noutput;

  if (!format && ninput + noutput)
    error->all(FLERR, "Missing python format keyword");
  else if (format && ((int) strlen(format) != ninput + noutput))
    error->all(FLERR, "Input/output arguments ({}) and format characters ({}) are inconsistent",
               (ninput + noutput), strlen(format));

  // process inputs as values or variables

  pfuncs[ifunc].itype = new int[ninput];
  pfuncs[ifunc].ivarflag = new int[ninput];
  pfuncs[ifunc].ivalue = new int[ninput];
  pfuncs[ifunc].dvalue = new double[ninput];
  pfuncs[ifunc].svalue = new char *[ninput];

  for (int i = 0; i < ninput; i++) {
    pfuncs[ifunc].svalue[i] = nullptr;
    char type = format[i];
    if (type == 'i') {
      pfuncs[ifunc].itype[i] = INT;
      if (utils::strmatch(istr[i], "^v_")) {
        pfuncs[ifunc].ivarflag[i] = 1;
        pfuncs[ifunc].svalue[i] = utils::strdup(istr[i] + 2);
      } else {
        pfuncs[ifunc].ivarflag[i] = 0;
        pfuncs[ifunc].ivalue[i] = utils::inumeric(FLERR, istr[i], false, lmp);
      }
    } else if (type == 'f') {
      pfuncs[ifunc].itype[i] = DOUBLE;
      if (utils::strmatch(istr[i], "^v_")) {
        pfuncs[ifunc].ivarflag[i] = 1;
        pfuncs[ifunc].svalue[i] = utils::strdup(istr[i] + 2);
      } else {
        pfuncs[ifunc].ivarflag[i] = 0;
        pfuncs[ifunc].dvalue[i] = utils::numeric(FLERR, istr[i], false, lmp);
      }
    } else if (type == 's') {
      pfuncs[ifunc].itype[i] = STRING;
      if (utils::strmatch(istr[i], "^v_")) {
        pfuncs[ifunc].ivarflag[i] = 1;
        pfuncs[ifunc].svalue[i] = utils::strdup(istr[i] + 2);
      } else {
        pfuncs[ifunc].ivarflag[i] = 0;
        pfuncs[ifunc].svalue[i] = utils::strdup(istr[i]);
      }
    } else if (type == 'p') {
      pfuncs[ifunc].ivarflag[i] = 0;
      pfuncs[ifunc].itype[i] = PTR;
      if (strcmp(istr[i], "SELF") != 0) error->all(FLERR, "Invalid python command");

    } else
      error->all(FLERR, "Invalid python format character: {}", type);
  }

  // process output as value or variable

  pfuncs[ifunc].ovarname = nullptr;
  pfuncs[ifunc].longstr = nullptr;
  if (!noutput) return ifunc;

  char type = format[ninput];
  if (type == 'i')
    pfuncs[ifunc].otype = INT;
  else if (type == 'f')
    pfuncs[ifunc].otype = DOUBLE;
  else if (type == 's')
    pfuncs[ifunc].otype = STRING;
  else
    error->all(FLERR, "Invalid python return format character: {}", type);

  if (length_longstr) {
    if (pfuncs[ifunc].otype != STRING)
      error->all(FLERR, "Python command length keyword cannot be used unless output is a string");
    pfuncs[ifunc].length_longstr = length_longstr;
    pfuncs[ifunc].longstr = new char[length_longstr + 1];
    pfuncs[ifunc].longstr[length_longstr] = '\0';
  }

  if (strstr(ostr, "v_") != ostr) error->all(FLERR, "Invalid python command");
  pfuncs[ifunc].ovarname = utils::strdup(ostr + 2);

  return ifunc;
}

/* ---------------------------------------------------------------------- */

int PythonImpl::execute_string(char *cmd)
{
  PyUtils::GIL lock;
  int err = PyRun_SimpleString(cmd);
  if (err) PyUtils::Print_Errors();
  return err;
}

/* ---------------------------------------------------------------------- */

int PythonImpl::execute_file(char *fname)
{
  FILE *fp = fopen(fname, "r");
  if (fp == nullptr) return -1;

  PyUtils::GIL lock;
  int err = PyRun_SimpleFile(fp, fname);
  if (err) PyUtils::Print_Errors();

  if (fp) fclose(fp);
  return err;
}

/* ------------------------------------------------------------------ */

void PythonImpl::deallocate(int i)
{
  delete[] pfuncs[i].itype;
  delete[] pfuncs[i].ivarflag;
  delete[] pfuncs[i].ivalue;
  delete[] pfuncs[i].dvalue;
  for (int j = 0; j < pfuncs[i].ninput; j++) delete[] pfuncs[i].svalue[j];
  delete[] pfuncs[i].svalue;
  delete[] pfuncs[i].ovarname;
  delete[] pfuncs[i].longstr;
}

/* ------------------------------------------------------------------ */

bool PythonImpl::has_minimum_version(int major, int minor)
{
  return (PY_MAJOR_VERSION == major && PY_MINOR_VERSION >= minor) || (PY_MAJOR_VERSION > major);
}

/* ------------------------------------------------------------------ */

void PythonImpl::finalize()
{
  if (Py_IsInitialized()) Py_Finalize();
}
