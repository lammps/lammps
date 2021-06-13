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

#ifndef LMP_PYTHON_UTILS_H
#define LMP_PYTHON_UTILS_H

#include <Python.h>

namespace LAMMPS_NS {

namespace PyUtils {

  class GIL {
    PyGILState_STATE gstate;

   public:
    GIL() : gstate(PyGILState_Ensure()) {}
    ~GIL() { PyGILState_Release(gstate); }
  };

  static void Print_Errors()
  {
    PyErr_Print();
    PyErr_Clear();
  }

}    // namespace PyUtils

}    // namespace LAMMPS_NS

#endif
