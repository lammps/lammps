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
   Contributing authors: Axel Kohlmeyer and Richard Berger (Temple U)
------------------------------------------------------------------------- */

#include "pair_python.h"

#include "atom.h"
#include "error.h"
#include "force.h"
#include "lmppython.h"
#include "memory.h"
#include "neigh_list.h"
#include "python_compat.h"
#include "python_utils.h"
#include "update.h"

#include <cstring>
#include <Python.h>  // IWYU pragma: export

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairPython::PairPython(LAMMPS *lmp) : Pair(lmp) {
  respa_enable = 0;
  single_enable = 1;
  writedata = 0;
  restartinfo = 0;
  one_coeff = 1;
  reinitflag = 0;
  cut_global = 0.0;
  centroidstressflag = CENTROID_SAME;

  py_potential = nullptr;
  skip_types = nullptr;

  python->init();

  // add current directory to PYTHONPATH

  PyUtils::GIL lock;
  PyObject *py_path = PySys_GetObject((char *)"path");
  PyList_Append(py_path, PY_STRING_FROM_STRING("."));

  // if LAMMPS_POTENTIALS environment variable is set,
  // add it to PYTHONPATH as well

  const char *potentials_path = getenv("LAMMPS_POTENTIALS");
  if (potentials_path != nullptr) {
    PyList_Append(py_path, PY_STRING_FROM_STRING(potentials_path));
  }
}

/* ---------------------------------------------------------------------- */

PairPython::~PairPython()
{
  PyUtils::GIL lock;
  Py_CLEAR(py_potential);
  delete[] skip_types;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
  }
}

/* ---------------------------------------------------------------------- */

void PairPython::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,factor_lj;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // prepare access to compute_force and compute_energy functions

  PyUtils::GIL lock;
  PyObject *py_pair_instance = (PyObject *) py_potential;
  PyObject *py_compute_force = PyObject_GetAttrString(py_pair_instance,"compute_force");
  if (!py_compute_force) {
    PyUtils::Print_Errors();
    error->all(FLERR,"Could not find 'compute_force' method'");
  }
  if (!PyCallable_Check(py_compute_force)) {
    PyUtils::Print_Errors();
    error->all(FLERR,"Python 'compute_force' is not callable");
  }

  PyObject *py_compute_energy = PyObject_GetAttrString(py_pair_instance,"compute_energy");
  if (!py_compute_energy) {
    PyUtils::Print_Errors();
    error->all(FLERR,"Could not find 'compute_energy' method'");
  }
  if (!PyCallable_Check(py_compute_energy)) {
    PyUtils::Print_Errors();
    error->all(FLERR,"Python 'compute_energy' is not callable");
  }

  PyObject *py_compute_args = PyTuple_New(3);
  if (!py_compute_args) {
    PyUtils::Print_Errors();
    error->all(FLERR,"Could not create tuple for 'compute' function arguments");
  }

  PyObject *py_rsq, *py_itype, *py_jtype, *py_value;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    py_itype = PY_INT_FROM_LONG(itype);
    PyTuple_SetItem(py_compute_args,1,py_itype);

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      // with hybrid/overlay we might get called for skipped types
      if (skip_types[itype] || skip_types[jtype]) continue;

      py_jtype = PY_INT_FROM_LONG(jtype);
      PyTuple_SetItem(py_compute_args,2,py_jtype);

      if (rsq < cutsq[itype][jtype]) {
        py_rsq = PyFloat_FromDouble(rsq);
        PyTuple_SetItem(py_compute_args,0,py_rsq);
        py_value = PyObject_CallObject(py_compute_force,py_compute_args);
        if (!py_value) {
          PyUtils::Print_Errors();
          error->all(FLERR,"Calling 'compute_force' function failed");
        }
        fpair = factor_lj*PyFloat_AsDouble(py_value);
        Py_CLEAR(py_value);

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          py_value = PyObject_CallObject(py_compute_energy,py_compute_args);
          if (!py_value) {
            PyUtils::Print_Errors();
            error->all(FLERR,"Calling 'compute_energy' function failed");
          }
          evdwl = factor_lj*PyFloat_AsDouble(py_value);
          Py_CLEAR(py_value);
        } else evdwl = 0.0;

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);
      }
    }
  }
  Py_CLEAR(py_compute_args);

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairPython::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairPython::settings(int narg, char **arg)
{
  if (narg != 1)
    error->all(FLERR,"Illegal pair_style command");

  cut_global = utils::numeric(FLERR,arg[0],false,lmp);
}

/* ----------------------------------------------------------------------
   set coeffs for all type pairs
------------------------------------------------------------------------- */

void PairPython::coeff(int narg, char **arg)
{
  const int ntypes = atom->ntypes;

  if (narg != 3+ntypes)
    error->all(FLERR,"Incorrect args for pair coefficients");

  if (!allocated) allocate();

  // make sure I,J args are * *

  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients");


  // check if python potential class type exists
  // load module if necessary
  std::string full_cls_name = arg[2];
  std::string module_name = "__main__";
  std::string cls_name = full_cls_name;

  size_t lastpos = full_cls_name.rfind(".");

  if (lastpos != std::string::npos) {
    module_name = full_cls_name.substr(0, lastpos);
    cls_name = full_cls_name.substr(lastpos+1);
  }

  PyUtils::GIL lock;

  PyObject * pModule = PyImport_ImportModule(module_name.c_str());

  if (!pModule) {
    PyUtils::Print_Errors();
    error->all(FLERR,"Loading python pair style module failure");
  }

  // create LAMMPS atom type to potential file type mapping in python class
  // by calling 'lammps_pair_style.map_coeff(name,type)'

  PyObject *py_pair_type = PyObject_GetAttrString(pModule, cls_name.c_str());
  if (!py_pair_type) {
    PyUtils::Print_Errors();
    error->all(FLERR, "Could not find pair style class {} in module {}", cls_name, module_name);
  }

  PyObject * py_pair_instance = PyObject_CallObject(py_pair_type, nullptr);
  if (!py_pair_instance) {
    PyUtils::Print_Errors();
    error->all(FLERR,"Could not instantiate instance of pair style class'");
  }

  py_potential = (void *) py_pair_instance;

  PyObject *py_value = PyObject_CallMethod(py_pair_instance, (char *)"check_units", (char *)"s", update->unit_style);
  if (!py_value) {
    PyUtils::Print_Errors();
    error->all(FLERR,"Calling 'check_units' function failed");
  }
  Py_CLEAR(py_value);


  delete[] skip_types;
  skip_types = new int[ntypes+1];
  skip_types[0] = 1;
  for (int i = 1; i <= ntypes ; i++) {
    if (strcmp(arg[2+i],"NULL") == 0) {
      skip_types[i] = 1;
      continue;
    } else skip_types[i] = 0;
    const int type = i;
    const char * name = arg[2+i];
    py_value = PyObject_CallMethod(py_pair_instance, (char *)"map_coeff", (char *)"si", name, type);
    if (!py_value) {
      PyUtils::Print_Errors();
      error->all(FLERR,"Calling 'map_coeff' function failed");
    }
    Py_CLEAR(py_value);

    for (int j = i; j <= ntypes ; j++) {
      setflag[i][j] = 1;
      cutsq[i][j] = cut_global*cut_global;
    }
  }
}

/* ---------------------------------------------------------------------- */

double PairPython::init_one(int, int)
{
  return cut_global;
}

/* ---------------------------------------------------------------------- */

double PairPython::single(int /* i */, int /* j */, int itype, int jtype,
                         double rsq, double /* factor_coul */,
                         double factor_lj, double &fforce)
{
  // with hybrid/overlay we might get called for skipped types
  if (skip_types[itype] || skip_types[jtype]) {
    fforce = 0.0;
    return 0.0;
  }

  // prepare access to compute_force and compute_energy functions

  PyUtils::GIL lock;
  PyObject *py_compute_force = (PyObject *) get_member_function("compute_force");
  PyObject *py_compute_energy = (PyObject *) get_member_function("compute_energy");
  PyObject *py_compute_args = Py_BuildValue("(dii)", rsq, itype, jtype);

  if (!py_compute_args) {
    PyUtils::Print_Errors();
    error->all(FLERR,"Could not create tuple for 'compute' function arguments");
  }

  PyObject * py_value = PyObject_CallObject(py_compute_force, py_compute_args);
  if (!py_value) {
    PyUtils::Print_Errors();
    error->all(FLERR,"Calling 'compute_force' function failed");
  }
  fforce = factor_lj*PyFloat_AsDouble(py_value);
  Py_CLEAR(py_value);

  py_value = PyObject_CallObject(py_compute_energy, py_compute_args);
  if (!py_value) {
    PyUtils::Print_Errors();
    error->all(FLERR,"Calling 'compute_energy' function failed");
  }
  double evdwl = factor_lj*PyFloat_AsDouble(py_value);

  Py_CLEAR(py_value);
  Py_CLEAR(py_compute_args);

  return evdwl;
}


/* ---------------------------------------------------------------------- */

void * PairPython::get_member_function(const char * name)
{
  PyUtils::GIL lock;
  PyObject *py_pair_instance = (PyObject *) py_potential;
  PyObject * py_mfunc = PyObject_GetAttrString(py_pair_instance, name);
  if (!py_mfunc) {
    PyUtils::Print_Errors();
    error->all(FLERR,"Could not find '{}' method'", name);
  }
  if (!PyCallable_Check(py_mfunc)) {
    PyUtils::Print_Errors();
    error->all(FLERR,"Python '{}' is not callable", name);
  }
  return py_mfunc;
}
