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
   Contributing author: Aidan Thompson (SNL)
------------------------------------------------------------------------- */

#include "mliap_descriptor.h"

#include "memory.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

MLIAPDescriptor::MLIAPDescriptor(LAMMPS *lmp) :
    Pointers(lmp), ndescriptors(0), nelements(0), elements(nullptr), cutsq(nullptr),
    radelem(nullptr), wjelem(nullptr)
{
  cutmax = 0.0;
}

/* ---------------------------------------------------------------------- */

MLIAPDescriptor::~MLIAPDescriptor()
{
  for (int i = 0; i < nelements; i++) delete[] elements[i];
  delete[] elements;
  memory->destroy(cutsq);
  memory->destroy(radelem);
  memory->destroy(wjelem);
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double MLIAPDescriptor::memory_usage()
{
  double bytes = (double) nelements * sizeof(double);          // radelem
  bytes += (double) nelements * sizeof(double);                // welem
  bytes += (double) nelements * nelements * sizeof(double);    // cutsq

  return bytes;
}
