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

// collection of smart pointers for specific purposes

#include "safe_pointers.h"

#include "platform.h"

using namespace LAMMPS_NS;

SafeFilePtr::~SafeFilePtr()
{
  if (fp) {
    if (use_pclose)
      platform::pclose(fp);
    else
      fclose(fp);
  }
}

SafeFilePtr &SafeFilePtr::operator=(FILE *_fp)
{
  if (fp && (fp != _fp)) {
    if (use_pclose)
      platform::pclose(fp);
    else
      fclose(fp);
  }
  fp = _fp;
  if (_fp == nullptr) use_pclose = false;
  return *this;
}
