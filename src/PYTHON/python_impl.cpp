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
   Contributing author: Richard Berger and Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include <cstdlib>
#include <cstring>
#include <Python.h>
#include "python_impl.h"
#include "force.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "python_compat.h"

using namespace LAMMPS_NS;

enum{NONE,INT,DOUBLE,STRING,PTR};

#define VALUELENGTH 64               // also in variable.cpp


/* ---------------------------------------------------------------------- */

PythonImpl::PythonImpl(LAMMPS *lmp) : Pointers(lmp)
{
  ninput = noutput = 0;
  istr = NULL;
  ostr = NULL;
  format = NULL;
  length_longstr = 0;

  // pfuncs stores interface info for each Python function

  nfunc = 0;
  pfuncs = NULL;

  // one-time initialization of Python interpreter
  // pyMain stores pointer to main module
  external_interpreter = Py_IsInitialized();

  Py_Initialize();
  PyEval_InitThreads();

  PyGILState_STATE gstate = PyGILState_Ensure();

  PyObject *pModule = PyImport_AddModule("__main__");
  if (!pModule) error->all(FLERR,"Could not initialize embedded Python");

  pyMain = (void *) pModule;
  PyGILState_Release(gstate);
}

/* ---------------------------------------------------------------------- */

PythonImpl::~PythonImpl()
{
  if(pyMain) {
    // clean up
    PyGILState_STATE gstate = PyGILState_Ensure();

    for (int i = 0; i < nfunc; i++) {
      delete [] pfuncs[i].name;
      deallocate(i);
      PyObject *pFunc = (PyObject *) pfuncs[i].pFunc;
      Py_XDECREF(pFunc);
    }

    // shutdown Python interpreter

    if (!external_interpreter) {
      Py_Finalize();
    }
    else {
      PyGILState_Release(gstate);
    }
  }

  memory->sfree(pfuncs);
}

/* ---------------------------------------------------------------------- */

void PythonImpl::command(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Invalid python command");

  // if invoke is only keyword, invoke the previously defined function

  if (narg == 2 && strcmp(arg[1],"invoke") == 0) {
    int ifunc = find(arg[0]);
    if (ifunc < 0) error->all(FLERR,"Python invoke of undefined function");

    char *str = NULL;
    if (noutput) {
      str = input->variable->pythonstyle(pfuncs[ifunc].ovarname,
                                         pfuncs[ifunc].name);
      if (!str)
        error->all(FLERR,"Python variable does not match Python function");
    }

    invoke_function(ifunc,str);
    return;
  }

  // if source is only keyword, execute the python code

  if (narg == 3 && strcmp(arg[1],"source") == 0) {
    int err;

    FILE *fp = fopen(arg[2],"r");
    if (fp == NULL)
      err = execute_string(arg[2]);
    else
      err = execute_file(arg[2]);

    if (fp) fclose(fp);
    if (err) error->all(FLERR,"Could not process Python source command");

    return;
  }

  // parse optional args, invoke is not allowed in this mode

  ninput = noutput = 0;
  istr = NULL;
  ostr = NULL;
  format = NULL;
  length_longstr = 0;
  char *pyfile = NULL;
  char *herestr = NULL;
  int existflag = 0;

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"input") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Invalid python command");
      ninput = force->inumeric(FLERR,arg[iarg+1]);
      if (ninput < 0) error->all(FLERR,"Invalid python command");
      iarg += 2;
      istr = new char*[ninput];
      if (iarg+ninput > narg) error->all(FLERR,"Invalid python command");
      for (int i = 0; i < ninput; i++) istr[i] = arg[iarg+i];
      iarg += ninput;
    } else if (strcmp(arg[iarg],"return") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Invalid python command");
      noutput = 1;
      ostr = arg[iarg+1];
      iarg += 2;
    } else if (strcmp(arg[iarg],"format") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Invalid python command");
      int n = strlen(arg[iarg+1]) + 1;
      format = new char[n];
      strcpy(format,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"length") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Invalid python command");
      length_longstr = force->inumeric(FLERR,arg[iarg+1]);
      if (length_longstr <= 0) error->all(FLERR,"Invalid python command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"file") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Invalid python command");
      delete[] pyfile;
      int n = strlen(arg[iarg+1]) + 1;
      pyfile = new char[n];
      strcpy(pyfile,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"here") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Invalid python command");
      herestr = arg[iarg+1];
      iarg += 2;
    } else if (strcmp(arg[iarg],"exists") == 0) {
      existflag = 1;
      iarg++;
    } else error->all(FLERR,"Invalid python command");
  }

  if (pyfile && herestr) error->all(FLERR,"Invalid python command");
  if (pyfile && existflag) error->all(FLERR,"Invalid python command");
  if (herestr && existflag) error->all(FLERR,"Invalid python command");

  // create or overwrite entry in pfuncs vector with name = arg[0]

  int ifunc = create_entry(arg[0]);

  PyGILState_STATE gstate = PyGILState_Ensure();

  // send Python code to Python interpreter
  // file: read the file via PyRun_SimpleFile()
  // here: process the here string directly
  // exist: do nothing, assume code has already been run

  if (pyfile) {
    FILE *fp = fopen(pyfile,"r");

    if (fp == NULL) {
      PyGILState_Release(gstate);
      error->all(FLERR,"Could not open Python file");
    }

    int err = PyRun_SimpleFile(fp,pyfile);

    if (err) {
      PyGILState_Release(gstate);
      error->all(FLERR,"Could not process Python file");
    }

    fclose(fp);
  } else if (herestr) {
    int err = PyRun_SimpleString(herestr);

    if (err) {
      PyGILState_Release(gstate);
      error->all(FLERR,"Could not process Python string");
    }
  }

  // pFunc = function object for requested function

  PyObject *pModule = (PyObject *) pyMain;
  PyObject *pFunc = PyObject_GetAttrString(pModule,pfuncs[ifunc].name);

  if (!pFunc) {
    PyGILState_Release(gstate);
    error->all(FLERR,"Could not find Python function");
  }

  if (!PyCallable_Check(pFunc)) {
    PyGILState_Release(gstate);
    error->all(FLERR,"Python function is not callable");
  }

  pfuncs[ifunc].pFunc = (void *) pFunc;

  // clean-up input storage

  delete [] istr;
  delete [] format;
  delete [] pyfile;
  PyGILState_Release(gstate);
}

/* ------------------------------------------------------------------ */

void PythonImpl::invoke_function(int ifunc, char *result)
{
  PyGILState_STATE gstate = PyGILState_Ensure();
  PyObject *pValue;
  char *str;

  PyObject *pFunc = (PyObject *) pfuncs[ifunc].pFunc;

  // create Python tuple of input arguments

  int ninput = pfuncs[ifunc].ninput;
  PyObject *pArgs = PyTuple_New(ninput);

  if (!pArgs) {
    PyGILState_Release(gstate);
    error->all(FLERR,"Could not create Python function arguments");
  }

  for (int i = 0; i < ninput; i++) {
    int itype = pfuncs[ifunc].itype[i];
    if (itype == INT) {
      if (pfuncs[ifunc].ivarflag[i]) {
        str = input->variable->retrieve(pfuncs[ifunc].svalue[i]);

        if (!str) {
          PyGILState_Release(gstate);
          error->all(FLERR,"Could not evaluate Python function input variable");
        }

        pValue = PY_INT_FROM_LONG(atoi(str));
      } else {
        pValue = PY_INT_FROM_LONG(pfuncs[ifunc].ivalue[i]);
      }
    } else if (itype == DOUBLE) {
      if (pfuncs[ifunc].ivarflag[i]) {
        str = input->variable->retrieve(pfuncs[ifunc].svalue[i]);

        if (!str) {
          PyGILState_Release(gstate);
          error->all(FLERR,"Could not evaluate Python function input variable");
        }

        pValue = PyFloat_FromDouble(atof(str));
      } else {
        pValue = PyFloat_FromDouble(pfuncs[ifunc].dvalue[i]);
      }
    } else if (itype == STRING) {
      if (pfuncs[ifunc].ivarflag[i]) {
        str = input->variable->retrieve(pfuncs[ifunc].svalue[i]);
        if (!str) {
          PyGILState_Release(gstate);
          error->all(FLERR,"Could not evaluate Python function input variable");
        }

        pValue = PY_STRING_FROM_STRING(str);
      } else {
        pValue = PY_STRING_FROM_STRING(pfuncs[ifunc].svalue[i]);
      }
    } else if (itype == PTR) {
      pValue = PY_VOID_POINTER(lmp);
    } else {
      PyGILState_Release(gstate);
      error->all(FLERR,"Unsupported variable type");
    }
    PyTuple_SetItem(pArgs,i,pValue);
  }

  // call the Python function
  // error check with one() since only some procs may fail

  pValue = PyObject_CallObject(pFunc,pArgs);

  if (!pValue) {
    PyErr_Print();
    PyGILState_Release(gstate);
    error->one(FLERR,"Python function evaluation failed");
  }

  Py_DECREF(pArgs);

  // function returned a value
  // assign it to result string stored by python-style variable
  // or if user specified a length, assign it to longstr

  if (pfuncs[ifunc].noutput) {
    int otype = pfuncs[ifunc].otype;
    if (otype == INT) {
      sprintf(result,"%ld",PY_INT_AS_LONG(pValue));
    } else if (otype == DOUBLE) {
      sprintf(result,"%.15g",PyFloat_AsDouble(pValue));
    } else if (otype == STRING) {
      const char *pystr = PY_STRING_AS_STRING(pValue);
      if (pfuncs[ifunc].longstr)
        strncpy(pfuncs[ifunc].longstr,pystr,pfuncs[ifunc].length_longstr);
      else strncpy(result,pystr,VALUELENGTH-1);
    }
    Py_DECREF(pValue);
  }

  PyGILState_Release(gstate);
}

/* ------------------------------------------------------------------ */

int PythonImpl::find(char *name)
{
  for (int i = 0; i < nfunc; i++)
    if (strcmp(name,pfuncs[i].name) == 0) return i;
  return -1;
}

/* ------------------------------------------------------------------ */

int PythonImpl::variable_match(char *name, char *varname, int numeric)
{
  int ifunc = find(name);
  if (ifunc < 0) return -1;
  if (pfuncs[ifunc].noutput == 0) return -1;
  if (strcmp(pfuncs[ifunc].ovarname,varname) != 0) return -1;
  if (numeric && pfuncs[ifunc].otype == STRING) return -1;
  return ifunc;
}

/* ------------------------------------------------------------------ */

char *PythonImpl::long_string(int ifunc)
{
  return pfuncs[ifunc].longstr;
}

/* ------------------------------------------------------------------ */

int PythonImpl::create_entry(char *name)
{
  // ifunc = index to entry by name in pfuncs vector, can be old or new
  // free old vectors if overwriting old pfunc

  int ifunc = find(name);

  if (ifunc < 0) {
    ifunc = nfunc;
    nfunc++;
    pfuncs = (PyFunc *)
      memory->srealloc(pfuncs,nfunc*sizeof(struct PyFunc),"python:pfuncs");
    int n = strlen(name) + 1;
    pfuncs[ifunc].name = new char[n];
    strcpy(pfuncs[ifunc].name,name);
  } else deallocate(ifunc);

  pfuncs[ifunc].ninput = ninput;
  pfuncs[ifunc].noutput = noutput;

  if (!format && ninput+noutput)
    error->all(FLERR,"Invalid python command");
  else if (format && strlen(format) != ninput+noutput)
    error->all(FLERR,"Invalid python command");

  // process inputs as values or variables

  pfuncs[ifunc].itype = new int[ninput];
  pfuncs[ifunc].ivarflag = new int[ninput];
  pfuncs[ifunc].ivalue = new int[ninput];
  pfuncs[ifunc].dvalue = new double[ninput];
  pfuncs[ifunc].svalue = new char*[ninput];

  for (int i = 0; i < ninput; i++) {
    pfuncs[ifunc].svalue[i] = NULL;
    char type = format[i];
    if (type == 'i') {
      pfuncs[ifunc].itype[i] = INT;
      if (strstr(istr[i],"v_") == istr[i]) {
        pfuncs[ifunc].ivarflag[i] = 1;
        int n = strlen(&istr[i][2]) + 1;
        pfuncs[ifunc].svalue[i] = new char[n];
        strcpy(pfuncs[ifunc].svalue[i],&istr[i][2]);
      } else {
        pfuncs[ifunc].ivarflag[i] = 0;
        pfuncs[ifunc].ivalue[i] = force->inumeric(FLERR,istr[i]);
      }
    } else if (type == 'f') {
      pfuncs[ifunc].itype[i] = DOUBLE;
      if (strstr(istr[i],"v_") == istr[i]) {
        pfuncs[ifunc].ivarflag[i] = 1;
        int n = strlen(&istr[i][2]) + 1;
        pfuncs[ifunc].svalue[i] = new char[n];
        strcpy(pfuncs[ifunc].svalue[i],&istr[i][2]);
      } else {
        pfuncs[ifunc].ivarflag[i] = 0;
        pfuncs[ifunc].dvalue[i] = force->numeric(FLERR,istr[i]);
      }
    } else if (type == 's') {
      pfuncs[ifunc].itype[i] = STRING;
      if (strstr(istr[i],"v_") == istr[i]) {
        pfuncs[ifunc].ivarflag[i] = 1;
        int n = strlen(&istr[i][2]) + 1;
        pfuncs[ifunc].svalue[i] = new char[n];
        strcpy(pfuncs[ifunc].svalue[i],&istr[i][2]);
      } else {
        pfuncs[ifunc].ivarflag[i] = 0;
        int n = strlen(istr[i]) + 1;
        pfuncs[ifunc].svalue[i] = new char[n];
        strcpy(pfuncs[ifunc].svalue[i],istr[i]);
      }
    } else if (type == 'p') {
      pfuncs[ifunc].ivarflag[i] = 0;
      pfuncs[ifunc].itype[i] = PTR;
      if (strcmp(istr[i],"SELF") != 0)
        error->all(FLERR,"Invalid python command");

    } else error->all(FLERR,"Invalid python command");
  }

  // process output as value or variable

  pfuncs[ifunc].ovarname = NULL;
  pfuncs[ifunc].longstr = NULL;
  if (!noutput) return ifunc;

  char type = format[ninput];
  if (type == 'i') pfuncs[ifunc].otype = INT;
  else if (type == 'f') pfuncs[ifunc].otype = DOUBLE;
  else if (type == 's') pfuncs[ifunc].otype = STRING;
  else error->all(FLERR,"Invalid python command");

  if (length_longstr) {
    if (pfuncs[ifunc].otype != STRING)
      error->all(FLERR,"Python command length keyword "
                 "cannot be used unless output is a string");
    pfuncs[ifunc].length_longstr = length_longstr;
    pfuncs[ifunc].longstr = new char[length_longstr+1];
  }

  if (strstr(ostr,"v_") != ostr) error->all(FLERR,"Invalid python command");
  int n = strlen(&ostr[2]) + 1;
  pfuncs[ifunc].ovarname = new char[n];
  strcpy(pfuncs[ifunc].ovarname,&ostr[2]);

  return ifunc;
}

/* ---------------------------------------------------------------------- */

int PythonImpl::execute_string(char *cmd)
{
  PyGILState_STATE gstate = PyGILState_Ensure();
  int err = PyRun_SimpleString(cmd);
  PyGILState_Release(gstate);

  return err;
}

/* ---------------------------------------------------------------------- */

int PythonImpl::execute_file(char *fname)
{
  FILE *fp = fopen(fname,"r");
  if (fp == NULL) return -1;

  PyGILState_STATE gstate = PyGILState_Ensure();
  int err = PyRun_SimpleFile(fp,fname);
  PyGILState_Release(gstate);

  if (fp) fclose(fp);
  return err;
}

/* ------------------------------------------------------------------ */

void PythonImpl::deallocate(int i)
{
  delete [] pfuncs[i].itype;
  delete [] pfuncs[i].ivarflag;
  delete [] pfuncs[i].ivalue;
  delete [] pfuncs[i].dvalue;
  for (int j = 0; j < pfuncs[i].ninput; j++)
    delete [] pfuncs[i].svalue[j];
  delete [] pfuncs[i].svalue;
  delete [] pfuncs[i].ovarname;
  delete [] pfuncs[i].longstr;
}
