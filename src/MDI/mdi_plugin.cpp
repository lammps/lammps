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

#include <cstring>
#include "error.h"
#include "fix_mdi_aimd.h"
#include "input.h"
#include "modify.h"

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

MDIPlugin::MDIPlugin(LAMMPS *lmp, int narg, char **arg) : Pointers(lmp)
{
  if (narg < 1) error->all(FLERR,"Illegal mdi plugin command");

  char *plugin_name = arg[0];

  char *mdi_arg = nullptr;
  char *infile_arg = nullptr;
  char *extra_arg = nullptr;
  lammps_command = nullptr;

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"mdi") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal mdi plugin command");
      mdi_arg = arg[iarg+1];
      iarg += 2;
    } else if (strcmp(arg[iarg],"infile") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal mdi plugin command");
      infile_arg = arg[iarg+1];
      iarg += 2;
    } else if (strcmp(arg[iarg],"extra") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal mdi plugin command");
      extra_arg = arg[iarg+1];
      iarg += 2;
    } else if (strcmp(arg[iarg],"command") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal mdi plugin command");
      int n = strlen(arg[iarg+1]) + 1;
      lammps_command = new char[n];
      strcpy(lammps_command,arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal mdi plugin command");
  }

  // error checks

  if (!mdi_arg || !infile_arg || !lammps_command)
    error->all(FLERR,"MDI plugin must specify mdi, infile, command keywords");

  // build full plugin_args string for args to plugin library

  int n = strlen(mdi_arg) + strlen(infile_arg) + strlen(extra_arg) + 16;
  char *plugin_args = new char[n];
  plugin_args[0] = 0;
  strcat(plugin_args,"-mdi \"");
  strcat(plugin_args,mdi_arg);
  strcat(plugin_args,"\" -in ");
  strcat(plugin_args,infile_arg);
  if (extra_arg) {
    strcat(plugin_args," ");
    strcat(plugin_args,extra_arg);
  }

  // find FixMDIAimd instance so can reset its mdicomm
  // NOTE: this is a kludge - need better way to handle this

  fixptr = modify->get_fix_by_style("mdi/aimd")[0];

  // launch the MDI plugin library
  // path for lib was specified in -mdi command-line arg when LAMMPS started
  // this calls back to plugin_wrapper, which must issue MDI EXIT at end

  MDI_Launch_plugin(plugin_name,plugin_args,&world,plugin_wrapper,(void *)this);

  delete [] plugin_args;
}

/* ----------------------------------------------------------------------
   callback function from MDI_Launch_plugin()
   this function must wrap entire interaction of LAMMPS as a driver
     with the plugin
---------------------------------------------------------------------- */

int MDIPlugin::plugin_wrapper(void *pmpicomm, MDI_Comm mdicomm,
                              void *vptr)
{
  MPI_Comm mpicomm = *(MPI_Comm *) pmpicomm;
  MDIPlugin *ptr = (MDIPlugin *) vptr;
  LAMMPS *lammps = ptr->lmp;
  char *lammps_command = ptr->lammps_command;

  // set FixMDIAimd mdicomm to driver's mdicomm passed to this callback

  FixMDIAimd *aimdptr = (FixMDIAimd *) (ptr->fixptr);
  aimdptr->mdicomm = mdicomm;

  // invoke the specified LAMMPS command
  // that operation will issue MDI commands to the plugin engine

  lammps->input->one(lammps_command);
  delete [] lammps_command;

  // send MDI exit to plugin, which unloads the plugin

  MDI_Send_command("EXIT",mdicomm);

  return 0;
}
