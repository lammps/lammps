/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_PLUGIN_H
#define LMP_PLUGIN_H

#include "lammpsplugin.h"

namespace LAMMPS_NS
{
  class LAMMPS;

  void plugin_load(const char *, LAMMPS *);
  void plugin_register(lammpsplugin_t *, void *);

  void plugin_unload(const char *, const char *, LAMMPS *);
  void plugin_erase(const char *, const char *);

  int plugin_get_num_plugins();
  int plugin_find(const char *, const char *);
  const lammpsplugin_t *plugin_get_info(int);
}

#endif
