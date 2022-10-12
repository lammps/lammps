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

#include "mdi_command.h"

#include "error.h"
#include "mdi_engine.h"
#include "mdi_plugin.h"
#include "memory.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   mdi command: engine or plugin or connect or exit
   engine is used when LAMMPS is an MDI engine, to start listening for requests
   plugin is used when LAMMPS is an MDI driver to load a plugin library
   connect and exit are used when LAMMPS is an MDI driver to
     (a) connect = setup comm with a stand-alone MDI engine
     (b) exit = terminate comm with a stand-alone MDI engine
---------------------------------------------------------------------- */

void MDICommand::command(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR, "Illegal mdi command");

  if (strcmp(arg[0], "engine") == 0) {
    MDIEngine(lmp, narg - 1, &arg[1]);

  } else if (strcmp(arg[0], "plugin") == 0) {
    MDIPlugin(lmp, narg - 1, &arg[1]);

  } else if (strcmp(arg[0], "connect") == 0) {

    if (lmp->mdicomm != nullptr)
      error->all(FLERR, "MDI cannot connect to already connected engine");

    MDI_Comm mdicomm;
    MDI_Get_communicator(&mdicomm, 0);

    if (mdicomm == MDI_COMM_NULL) {
      MDI_Accept_communicator(&mdicomm);
      if (mdicomm == MDI_COMM_NULL)
        error->all(FLERR, "MDI unable to connect to stand-alone engine");
    } else
      error->all(FLERR, "Cannot use mdi connect with plugin engine");

    int nbytes = sizeof(MDI_Comm);
    char *ptrcomm = (char *) memory->smalloc(nbytes, "mdi:mdicomm");
    memcpy(ptrcomm, &mdicomm, nbytes);

    lmp->mdicomm = (void *) ptrcomm;

  } else if (strcmp(arg[0], "exit") == 0) {

    if (lmp->mdicomm == nullptr) error->all(FLERR, "MDI cannot send exit to unconnected engine");

    MDI_Comm mdicomm;
    int nbytes = sizeof(MDI_Comm);
    char *ptrcomm = (char *) lmp->mdicomm;
    memcpy(&mdicomm, ptrcomm, nbytes);

    int ierr = MDI_Send_command("EXIT", mdicomm);
    if (ierr) error->all(FLERR, "MDI: EXIT command");

    memory->sfree(ptrcomm);
    lmp->mdicomm = nullptr;

  } else
    error->all(FLERR, "Illegal mdi command");
}
