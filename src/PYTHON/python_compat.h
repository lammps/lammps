/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_PYTHON_COMPAT_H
#define LMP_PYTHON_COMPAT_H

#include <Python.h>

// Wrap API changes between Python 2 and 3 using macros
#if PY_MAJOR_VERSION == 2
#if defined(_MSC_VER) || defined(__MINGW32__)
#define PY_INT_FROM_LONG(X) PyLong_FromLongLong(X)
#else
#define PY_INT_FROM_LONG(X) PyInt_FromLong(X)
#endif
#if defined(_MSC_VER) || defined(__MINGW32__)
#define PY_INT_AS_LONG(X) PyLong_AsLongLong(X)
#else
#define PY_INT_AS_LONG(X) PyInt_AsLong(X)
#endif
#define PY_STRING_FROM_STRING(X) PyString_FromString(X)
#define PY_VOID_POINTER(X) PyCObject_FromVoidPtr((void *) X, nullptr)
#define PY_STRING_AS_STRING(X) PyString_AsString(X)

#elif PY_MAJOR_VERSION == 3
#if defined(_MSC_VER) || defined(__MINGW32__)
#define PY_INT_FROM_LONG(X) PyLong_FromLongLong(X)
#else
#define PY_INT_FROM_LONG(X) PyLong_FromLong(X)
#endif
#if defined(_MSC_VER) || defined(__MINGW32__)
#define PY_INT_AS_LONG(X) PyLong_AsLongLong(X)
#else
#define PY_INT_AS_LONG(X) PyLong_AsLong(X)
#endif
#define PY_STRING_FROM_STRING(X) PyUnicode_FromString(X)
#define PY_VOID_POINTER(X) PyCapsule_New((void *) X, nullptr, nullptr)
#define PY_STRING_AS_STRING(X) PyUnicode_AsUTF8(X)
#endif

#endif
