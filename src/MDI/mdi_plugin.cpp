/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Taylor Barnes (MolSSI)
   MolSSI Driver Interface (MDI) support for LAMMPS
------------------------------------------------------------------------- */

#include "mdi_plugin.h"

#include "error.h"
#include "input.h"
#include "memory.h"

#include <cstring>

#include <mdi.h>

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   trigger LAMMPS to load an MDI plugin engine
   after loading the plugin library, it executes a LAMMPS command such as "run"
   the command will use other LAMMPS commands, such as fix mdi/aimd
     which act as an MDI driver, issuing MDI commands to the engine
   when MDI_Launch_plugin() exits, the engine is shut down and
     this class is destroyed
---------------------------------------------------------------------- */

MDIPlugin::MDIPlugin(LAMMPS *_lmp, int narg, char **arg) : Pointers(_lmp)
{
  if (narg < 1) error->all(FLERR, "Illegal mdi plugin command");

  char *plugin_name = arg[0];

  char *mdi_arg = nullptr;
  char *infile_arg = nullptr;
  char *extra_arg = nullptr;
  lammps_command = nullptr;

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "mdi") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal mdi plugin command");
      mdi_arg = arg[iarg + 1];
      iarg += 2;
    } else if (strcmp(arg[iarg], "infile") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal mdi plugin command");
      infile_arg = arg[iarg + 1];

      iarg += 2;
    } else if (strcmp(arg[iarg], "extra") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal mdi plugin command");
      extra_arg = arg[iarg + 1];

      // do variable substitution in multiple word extra_arg

      int ncopy = strlen(extra_arg) + 1;
      char *copy = (char *) memory->smalloc(ncopy,"mdi_plugin:copy");
      strncpy(copy, extra_arg, ncopy);
      char *work = (char *) memory->smalloc(ncopy,"mdi_plugin:work");
      int nwork = ncopy;
      input->substitute(copy, work, ncopy, nwork, 0);
      memory->sfree(work);
      extra_arg = copy;

      iarg += 2;
    } else if (strcmp(arg[iarg], "command") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal mdi plugin command");
      lammps_command = arg[iarg + 1];

      // do variable substitution in multiple word lammps_command

      int ncopy = strlen(lammps_command) + 1;
      char *copy = (char *) memory->smalloc(ncopy,"mdi_plugin:work");
      strncpy(copy, lammps_command, ncopy);
      char *work = (char *) memory->smalloc(ncopy,"mdi_plugin:work");
      int nwork = ncopy;
      input->substitute(copy, work, ncopy, nwork, 0);
      memory->sfree(work);
      lammps_command = copy;

      iarg += 2;
    } else
      error->all(FLERR, "Illegal mdi plugin command");
  }

  // error checks

  if (!mdi_arg || !lammps_command)
    error->all(FLERR, "MDI plugin must specify mdi and command keywords");

  // build full plugin_args string for args to plugin library

  int n = strlen(mdi_arg) + 16;
  if (infile_arg) n += strlen(infile_arg);
  if (extra_arg) n += strlen(extra_arg);
  auto plugin_args = new char[n];
  plugin_args[0] = 0;
  strcat(plugin_args, "-mdi \"");
  strcat(plugin_args, mdi_arg);
  strcat(plugin_args, "\"");
  if (infile_arg) {
    strcat(plugin_args, " -in ");
    strcat(plugin_args, infile_arg);
  }
  if (extra_arg) {
    strcat(plugin_args, " ");
    strcat(plugin_args, extra_arg);
  }

  // launch the MDI plugin library
  // path for lib was specified in -mdi command-line arg when LAMMPS started
  // this calls back to plugin_wrapper(), which issues MDI EXIT at end & returns
  // plugin_wrapper() must be a static method

  MDI_Launch_plugin(plugin_name, plugin_args, &world, plugin_wrapper, (void *) this);

  delete[] plugin_args;
  memory->sfree(extra_arg);
  memory->sfree(lammps_command);
}

/* ----------------------------------------------------------------------
   wrapper on entire interaction of LAMMPS as a driver with the plugin engine
   invoked as a callback by MDI once plugin library engine is launched
   this is a static method in mdi_plugin.h
---------------------------------------------------------------------- */

int MDIPlugin::plugin_wrapper(void * /*pmpicomm*/, MDI_Comm mdicomm, void *vptr)
{
  auto ptr = (MDIPlugin *) vptr;
  LAMMPS *lammps = ptr->lmp;
  char *lammps_command = ptr->lammps_command;

  // invoke the specified LAMMPS command
  // that operation will issue MDI commands to the plugin engine

  lammps->input->one(lammps_command);

  // send MDI exit to plugin, which unloads the plugin

  MDI_Send_command("EXIT", mdicomm);

  return 0;
}
