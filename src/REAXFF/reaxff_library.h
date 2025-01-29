/* -*- c -*- ------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LAMMPS_REAXFF_LIBRARY_H
#define LAMMPS_REAXFF_LIBRARY_H

#include "library.h"

// Ifdefs to allow this file to be included in C and C++ programs
#ifdef __cplusplus
extern "C" {
#endif

void lammps_set_reaxff_atm_parameter(void *handle, int type, int parameter_index, double value);
void lammps_set_reaxff_bnd_parameter(void *handle, int type1, int type2, int parameter_index, double value);
void lammps_set_reaxff_ofd_parameter(void *handle, int type1, int type2, int parameter_index, double value);
void lammps_set_reaxff_ang_parameter(void *handle, int type1, int type2, int type3, int parameter_index, double value);
void lammps_set_reaxff_tor_parameter(void *handle, int type1, int type2, int type3, int type4, int parameter_index, double value);
void lammps_set_reaxff_hbd_parameter(void *handle, int type1, int type2, int type3, int parameter_index, double value);

#ifdef __cplusplus
}
#endif

#endif /* LAMMPS_REAXFF_LIBRARY_H */

