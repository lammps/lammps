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

#ifdef COMMAND_CLASS

CommandStyle(plugin,Plugin)

#else

#ifndef LMP_PLUGIN_H
#define LMP_PLUGIN_H

#include "lammpsplugin.h"
#include "pointers.h"

namespace LAMMPS_NS
{

  class Plugin : protected Pointers {
  public:
    Plugin(class LAMMPS *);
    void command(int, char **);
  };

  void plugin_load(const char *, LAMMPS *);
  void plugin_register(lammpsplugin_t *, void *);

  void plugin_unload(const char *, const char *, LAMMPS *);
  void plugin_erase(const char *, const char *);
  void plugin_clear(LAMMPS *);

  int plugin_get_num_plugins();
  int plugin_find(const char *, const char *);
  const lammpsplugin_t *plugin_get_info(int);
}

#endif
#endif
