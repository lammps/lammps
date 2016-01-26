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

#include <Python.h>
#include "python.h"
#include "force.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{NONE,INT,DOUBLE,STRING,PTR};

#define VALUELENGTH 64               // also in variable.cpp

/* ---------------------------------------------------------------------- */

Python::Python(LAMMPS *lmp) : Pointers(lmp)
{
  python_exists = 1;

  pyMain = NULL;

  // pfuncs stores interface info for each Python function

  nfunc = 0;
  pfuncs = NULL;
}

/* ---------------------------------------------------------------------- */

Python::~Python()
{
  // clean up

  for (int i = 0; i < nfunc; i++) {
    delete [] pfuncs[i].name;
    deallocate(i);
    PyObject *pFunc = (PyObject *) pfuncs[i].pFunc;
    Py_XDECREF(pFunc);
  }

  // shutdown Python interpreter

  if (pyMain) Py_Finalize();

  memory->sfree(pfuncs);
}

/* ---------------------------------------------------------------------- */

void Python::command(int narg, char **arg)
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

  // parse optional args, invoke is not allowed in this mode

  ninput = noutput = 0;
  istr = NULL;
  ostr = NULL;
  format = NULL;
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

  // one-time intitialization of Python interpreter
  // Py_SetArgv() enables finding of *.py module files in current dir
  //   only needed for module load, not for direct file read into __main__
  // pymain stores pointer to main module

  if (pyMain == NULL) {
    if (Py_IsInitialized())
      error->all(FLERR,"Cannot embed Python when also "
                 "extending Python with LAMMPS");
    Py_Initialize();

    //char *arg = (char *) "./lmp";
    //PySys_SetArgv(1,&arg);

    //PyObject *pName = PyString_FromString("__main__");
    //if (!pName) errorX->all(FLERR,"Bad pName");
    //PyObject *pModule = PyImport_Import(pName);
    //Py_DECREF(pName);

    PyObject *pModule = PyImport_AddModule("__main__");
    if (!pModule) error->all(FLERR,"Could not initialize embedded Python");
    pyMain = (void *) pModule;
  }

  // send Python code to Python interpreter
  // file: read the file via PyRun_SimpleFile()
  // here: process the here string directly
  // exist: do nothing, assume code has already been run

  if (pyfile) {
    FILE *fp = fopen(pyfile,"r");
    if (fp == NULL) error->all(FLERR,"Could not open Python file");
    int err = PyRun_SimpleFile(fp,pyfile);
    if (err) error->all(FLERR,"Could not process Python file");
    fclose(fp);
  } else if (herestr) {
    int err = PyRun_SimpleString(herestr);
    if (err) error->all(FLERR,"Could not process Python string");
  }

  // pFunc = function object for requested function

  PyObject *pModule = (PyObject *) pyMain;
  PyObject *pFunc = PyObject_GetAttrString(pModule,pfuncs[ifunc].name);
  if (!pFunc) error->all(FLERR,"Could not find Python function");
  if (!PyCallable_Check(pFunc))
    error->all(FLERR,"Python function is not callable");
  pfuncs[ifunc].pFunc = (void *) pFunc;

  // clean-up input storage

  delete [] istr;
  delete [] format;
  delete [] pyfile;
}

/* ------------------------------------------------------------------ */

void Python::invoke_function(int ifunc, char *result)
{
  PyObject *pValue;
  char *str;

  PyObject *pFunc = (PyObject *) pfuncs[ifunc].pFunc;

  // create Python tuple of input arguments

  int ninput = pfuncs[ifunc].ninput;
  PyObject *pArgs = PyTuple_New(ninput);
  if (!pArgs) error->all(FLERR,"Could not create Python function arguments");

  for (int i = 0; i < ninput; i++) {
    int itype = pfuncs[ifunc].itype[i];
    if (itype == INT) {
      if (pfuncs[ifunc].ivarflag[i]) {
        str = input->variable->retrieve(pfuncs[ifunc].svalue[i]);
        if (!str)
          error->all(FLERR,"Could not evaluate Python function input variable");
        pValue = PyInt_FromLong(atoi(str));
      } else pValue = PyInt_FromLong(pfuncs[ifunc].ivalue[i]);
    } else if (itype == DOUBLE) {
      if (pfuncs[ifunc].ivarflag[i]) {
        str = input->variable->retrieve(pfuncs[ifunc].svalue[i]);
        if (!str)
          error->all(FLERR,"Could not evaluate Python function input variable");
        pValue = PyFloat_FromDouble(atof(str));
      } else pValue = PyFloat_FromDouble(pfuncs[ifunc].dvalue[i]);
    } else if (itype == STRING) {
      if (pfuncs[ifunc].ivarflag[i]) {
        str = input->variable->retrieve(pfuncs[ifunc].svalue[i]);
        if (!str)
          error->all(FLERR,"Could not evaluate Python function input variable");
        pValue = PyString_FromString(str);
      } else pValue = PyString_FromString(pfuncs[ifunc].svalue[i]);
    } else if (itype == PTR) {
      pValue = PyCObject_FromVoidPtr((void *) lmp,NULL);
    }
    PyTuple_SetItem(pArgs,i,pValue);
  }

  // call the Python function
  // error check with one() since only some procs may fail

  pValue = PyObject_CallObject(pFunc,pArgs);
  if (!pValue) error->one(FLERR,"Python function evaluation failed");
  Py_DECREF(pArgs);

  // function returned a value
  // assign it to result string stored by python-style variable

  if (pfuncs[ifunc].noutput) {
    int otype = pfuncs[ifunc].otype;
    if (otype == INT) {
      sprintf(result,"%ld",PyInt_AsLong(pValue));
    } else if (otype == DOUBLE) {
      sprintf(result,"%.15g",PyFloat_AsDouble(pValue));
    } else if (otype == STRING) {
      char *pystr = PyString_AsString(pValue);
      strncpy(result,pystr,VALUELENGTH-1);
    }
    Py_DECREF(pValue);
  }
}

/* ------------------------------------------------------------------ */

int Python::find(char *name)
{
  for (int i = 0; i < nfunc; i++)
    if (strcmp(name,pfuncs[i].name) == 0) return i;
  return -1;
}

/* ------------------------------------------------------------------ */

int Python::variable_match(char *name, char *varname, int numeric)
{
  int ifunc = find(name);
  if (ifunc < 0) return -1;
  if (pfuncs[ifunc].noutput == 0) return -1;
  if (strcmp(pfuncs[ifunc].ovarname,varname) != 0) return -1;
  if (numeric && pfuncs[ifunc].otype == STRING) return -1;
  return ifunc;
}

/* ------------------------------------------------------------------ */

int Python::create_entry(char *name)
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
  if (!noutput) return ifunc;

  char type = format[ninput];
  if (type == 'i') pfuncs[ifunc].otype = INT;
  else if (type == 'f') pfuncs[ifunc].otype = DOUBLE;
  else if (type == 's') pfuncs[ifunc].otype = STRING;
  else error->all(FLERR,"Invalid python command");

  if (strstr(ostr,"v_") != ostr) error->all(FLERR,"Invalid python command");
  int n = strlen(&ostr[2]) + 1;
  pfuncs[ifunc].ovarname = new char[n];
  strcpy(pfuncs[ifunc].ovarname,&ostr[2]);

  return ifunc;
}

/* ------------------------------------------------------------------ */

void Python::deallocate(int i)
{
  delete [] pfuncs[i].itype;
  delete [] pfuncs[i].ivarflag;
  delete [] pfuncs[i].ivalue;
  delete [] pfuncs[i].dvalue;
  for (int j = 0; j < pfuncs[i].ninput; j++)
    delete [] pfuncs[i].svalue[j];
  delete [] pfuncs[i].svalue;
  delete [] pfuncs[i].ovarname;
}
