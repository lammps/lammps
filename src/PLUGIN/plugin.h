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

#ifdef COMMAND_CLASS
// clang-format off
CommandStyle(plugin,Plugin);
// clang-format on
#else

#ifndef LMP_PLUGIN_H
#define LMP_PLUGIN_H

#include "command.h"
#include "lammpsplugin.h"    // IWYU pragma: export

namespace LAMMPS_NS {

class Plugin : public Command {
 public:
  Plugin(class LAMMPS *);
  void command(int, char **) override;
};

void plugin_auto_load(LAMMPS *);
int plugin_load(const char *, LAMMPS *);
void plugin_register(lammpsplugin_t *, void *);

void plugin_unload(const char *, const char *, LAMMPS *);
void plugin_erase(const char *, const char *);
void plugin_clear(LAMMPS *);

int plugin_get_num_plugins();
int plugin_find(const char *, const char *);
const lammpsplugin_t *plugin_get_info(int);
}    // namespace LAMMPS_NS

#endif
#endif
