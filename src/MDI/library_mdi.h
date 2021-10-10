/* -*- c -*- ------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LAMMPS_LIBRARY_MDI_H
#define LAMMPS_LIBRARY_MDI_H

/* C style library calls to LAMMPS when a LAMMPS shared library is
 *  used as a plugin through MolSSI Driver Interface (MDI). */

#include <mdi.h>

extern "C" {
int MDI_Plugin_init_lammps();
int lammps_execute_mdi_command(const char *, MDI_Comm, void *);
}
#endif
