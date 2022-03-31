/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

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
#include "fix_mdi_aimd.h"
#include "input.h"
#include "modify.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   mdi command: plugin
   may later have other MDI command variants
---------------------------------------------------------------------- */

void MDIPlugin::command(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal mdi/plugin command");

  char *plugin_name = arg[0];
  char *plugin_args = nullptr;
  plugin_command = nullptr;

  printf("NARG %d\n",narg);

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"args") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal mdi/plugin command");
      plugin_args = arg[iarg+1];
      iarg += 2;
    } else if (strcmp(arg[iarg],"command") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal mdi/plugin command");
      plugin_command = arg[iarg+1];
      iarg += 2;
    } else error->all(FLERR,"Illegal mdi/plugin command");
  }

  // error if no command was specified

  if (!plugin_command) error->all(FLERR,"MDI/plugin must specify command");

  // find FixMDIAimd instance so can reset its mdicomm

  fixptr = modify->get_fix_by_style("mdi/aimd")[0];

  // launch the MDI plugin library
  // path for lib was specified in -mdi command-line arg when LAMMPS started
  // this calls back to plugin_wrapper, which must issue MDI EXIT at end

  printf("PRE-LAUNCH\n");
  printf("NAME %s\n",plugin_name);
  printf("ARGS %s\n",plugin_args);

  MDI_Launch_plugin(plugin_name,plugin_args,world,plugin_wrapper,(void *)this);
}

/* ----------------------------------------------------------------------
   callback function from MDI_Launch_plugin()
   this function must wrap entire interaction of LAMMPS as a driver
     with the plugin
---------------------------------------------------------------------- */

int MDIPlugin::plugin_wrapper(void *pmpicomm, MDI_Comm mdicomm, 
                              void *ptr)
{
  printf("INSIDE CALLBACK\n");

  MPI_Comm mpicomm = *(MPI_Comm *) pmpicomm;
  MDIPlugin *thisptr = (MDIPlugin *) ptr;
  LAMMPS *lammps = thisptr->lmp;
 
  // set FixMDIAimd mdicomm to this mdicomm

  FixMDIAimd *aimdptr = (FixMDIAimd *) (thisptr->fixptr);
  aimdptr->mdicomm = mdicomm;

  // invoke the specified LAMMPS command
  // that operation will issue MDI commands to the plugin engine

  printf("PRE RUN command\n");

  lammps->input->one(thisptr->plugin_command);

  // send MDI exit to plugin, which unloads the plugin

  MDI_Send_command("EXIT",mdicomm);

  return 0;
}
