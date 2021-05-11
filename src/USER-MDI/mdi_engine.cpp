/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
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

#include "mdi_engine.h"

#include "atom.h"
#include "error.h"
#include "fix_mdi_engine.h"
#include "force.h"
#include "mdi.h"
#include "min.h"
#include "minimize.h"
#include "modify.h"
#include "output.h"
#include "timer.h"
#include "update.h"
#include "verlet.h"

#include <string.h>
#include <limits>

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   trigger LAMMPS to start acting as an MDI engine
   endlessly loop over receiving commands from driver and responding
   much of the logic for this is in FixMDIEngine
   when EXIT command is received, mdi_engine command exits
---------------------------------------------------------------------- */

void MDIEngine::command(int narg, char **arg)
{
  // list of nodes and commands that a MDI-compliant MD code should support

  // default node and its commands

  MDI_Register_Node("@DEFAULT");
  MDI_Register_Command("@DEFAULT", "<@");
  MDI_Register_Command("@DEFAULT", "<CELL");
  MDI_Register_Command("@DEFAULT", "<CHARGES");
  MDI_Register_Command("@DEFAULT", "<COORDS");
  MDI_Register_Command("@DEFAULT", "<LABELS");
  MDI_Register_Command("@DEFAULT", "<NATOMS");
  MDI_Register_Command("@DEFAULT", "<MASSES");
  MDI_Register_Command("@DEFAULT", ">COORDS");
  MDI_Register_Command("@DEFAULT", "@INIT_MD");
  MDI_Register_Command("@DEFAULT", "@INIT_OPTG");
  MDI_Register_Command("@DEFAULT", "EXIT");

  // node for setting up and running a dynamics simulation

  MDI_Register_Node("@INIT_MD");
  MDI_Register_Command("@INIT_MD", "<@");
  MDI_Register_Command("@INIT_MD", "<CELL");
  MDI_Register_Command("@INIT_MD", "<CHARGES");
  MDI_Register_Command("@INIT_MD", "<COORDS");
  MDI_Register_Command("@INIT_MD", "<ENERGY");
  MDI_Register_Command("@INIT_MD", "<FORCES");
  MDI_Register_Command("@INIT_MD", "<KE");
  MDI_Register_Command("@INIT_MD", "<LABELS");
  MDI_Register_Command("@INIT_MD", "<MASSES");
  MDI_Register_Command("@INIT_MD", "<NATOMS");
  MDI_Register_Command("@INIT_MD", "<PE");
  MDI_Register_Command("@INIT_MD", ">COORDS");
  MDI_Register_Command("@INIT_MD", ">FORCES");
  MDI_Register_Command("@INIT_MD", "@");
  MDI_Register_Command("@INIT_MD", "@COORDS");
  MDI_Register_Command("@INIT_MD", "@DEFAULT");
  MDI_Register_Command("@INIT_MD", "@FORCES");
  MDI_Register_Command("@INIT_MD", "@PRE-FORCES");
  MDI_Register_Command("@INIT_MD", "EXIT");

  // node for setting up and running a minimization

  MDI_Register_Node("@INIT_OPTG");
  MDI_Register_Command("@INIT_OPTG", "<@");
  MDI_Register_Command("@INIT_OPTG", "<CELL");
  MDI_Register_Command("@INIT_OPTG", "<CHARGES");
  MDI_Register_Command("@INIT_OPTG", "<COORDS");
  MDI_Register_Command("@INIT_OPTG", "<ENERGY");
  MDI_Register_Command("@INIT_OPTG", "<FORCES");
  MDI_Register_Command("@INIT_OPTG", "<KE");
  MDI_Register_Command("@INIT_OPTG", "<LABELS");
  MDI_Register_Command("@INIT_OPTG", "<MASSES");
  MDI_Register_Command("@INIT_OPTG", "<NATOMS");
  MDI_Register_Command("@INIT_OPTG", "<PE");
  MDI_Register_Command("@INIT_OPTG", ">COORDS");
  MDI_Register_Command("@INIT_OPTG", ">FORCES");
  MDI_Register_Command("@INIT_OPTG", "@");
  MDI_Register_Command("@INIT_OPTG", "@COORDS");
  MDI_Register_Command("@INIT_OPTG", "@DEFAULT");
  MDI_Register_Command("@INIT_OPTG", "@FORCES");
  MDI_Register_Command("@INIT_OPTG", "EXIT");

  // node at POST_FORCE location in timestep

  MDI_Register_Node("@FORCES");
  MDI_Register_Callback("@FORCES", ">FORCES");
  MDI_Register_Command("@FORCES", "<@");
  MDI_Register_Command("@FORCES", "<CELL");
  MDI_Register_Command("@FORCES", "<CHARGES");
  MDI_Register_Command("@FORCES", "<COORDS");
  MDI_Register_Command("@FORCES", "<ENERGY");
  MDI_Register_Command("@FORCES", "<FORCES");
  MDI_Register_Command("@FORCES", "<KE");
  MDI_Register_Command("@FORCES", "<LABELS");
  MDI_Register_Command("@FORCES", "<MASSES");
  MDI_Register_Command("@FORCES", "<NATOMS");
  MDI_Register_Command("@FORCES", "<PE");
  MDI_Register_Command("@FORCES", ">COORDS");
  MDI_Register_Command("@FORCES", ">FORCES");
  MDI_Register_Command("@FORCES", "@");
  MDI_Register_Command("@FORCES", "@COORDS");
  MDI_Register_Command("@FORCES", "@DEFAULT");
  MDI_Register_Command("@FORCES", "@FORCES");
  MDI_Register_Command("@FORCES", "@PRE-FORCES");
  MDI_Register_Command("@FORCES", "EXIT");

  // node at POST_INTEGRATE location in timestep

  MDI_Register_Node("@COORDS");
  MDI_Register_Command("@COORDS", "<@");
  MDI_Register_Command("@COORDS", "<CELL");
  MDI_Register_Command("@COORDS", "<CHARGES");
  MDI_Register_Command("@COORDS", "<COORDS");
  MDI_Register_Command("@COORDS", "<ENERGY");
  MDI_Register_Command("@COORDS", "<FORCES");
  MDI_Register_Command("@COORDS", "<KE");
  MDI_Register_Command("@COORDS", "<LABELS");
  MDI_Register_Command("@COORDS", "<MASSES");
  MDI_Register_Command("@COORDS", "<NATOMS");
  MDI_Register_Command("@COORDS", "<PE");
  MDI_Register_Command("@COORDS", ">COORDS");
  MDI_Register_Command("@COORDS", ">FORCES");
  MDI_Register_Command("@COORDS", "@");
  MDI_Register_Command("@COORDS", "@COORDS");
  MDI_Register_Command("@COORDS", "@DEFAULT");
  MDI_Register_Command("@COORDS", "@FORCES");
  MDI_Register_Command("@COORDS", "@PRE-FORCES");
  MDI_Register_Command("@COORDS", "EXIT");

  // if the mdi_engine fix is not already present, add it now

  int ifix = modify->find_fix_by_style("mdi/engine");
  bool added_mdi_engine_fix = false;
  if (ifix < 0) {
    modify->add_fix("MDI_ENGINE_INTERNAL all mdi/engine");
    added_mdi_engine_fix = true;
  }

  // identify the mdi_engine fix

  ifix = modify->find_fix_by_style("mdi/engine");
  mdi_fix = static_cast<FixMDIEngine*>(modify->fix[ifix]);

  // check that LAMMPS is setup as a compatible MDI engine

  if (narg > 0) error->all(FLERR,"Illegal mdi_engine command");

  if (atom->tag_enable == 0)
    error->all(FLERR,"Cannot use mdi_engine without atom IDs");

  if (atom->tag_consecutive() == 0)
    error->all(FLERR,"mdi_engine requires consecutive atom IDs");

  // endless engine loop, responding to driver commands

  char *command;

  while (1) {

    // mdi_engine command only recognizes three nodes
    // DEFAULT, INIT_MD, INIT_OPTG

    command = mdi_fix->engine_mode("@DEFAULT");
    
    // MDI commands for dynamics or minimization

    if (strcmp(command,"@INIT_MD") == 0 ) {
      command = mdi_md();
      if (strcmp(command,"EXIT")) break;
      
    } else if (strcmp(command,"@INIT_OPTG") == 0 ) {
      command = mdi_optg();
      if (strcmp(command,"EXIT")) break;

    } else if (strcmp(command,"EXIT") == 0) {
      break;
  
    } else
      error->all(FLERR,
                 fmt::format("MDI node exited with "
                             "invalid command: {}",command));
  }

  // remove mdi/engine fix that mdi_engine instantiated

  if (added_mdi_engine_fix) modify->delete_fix("MDI_ENGINE_INTERNAL");
}

/* ----------------------------------------------------------------------
   run an MD simulation under control of driver
---------------------------------------------------------------------- */

char *MDIEngine::mdi_md()
{
  // initialize an MD simulation

  update->whichflag = 1;
  timer->init_timeout();
  update->nsteps = 1;
  update->ntimestep = 0;
  update->firststep = update->ntimestep;
  update->laststep = update->ntimestep + update->nsteps;
  update->beginstep = update->firststep;
  update->endstep = update->laststep;

  lmp->init();

  // engine is now at @INIT_MD node

  char *command = NULL;
  command = mdi_fix->engine_mode("@INIT_MD");

  if (strcmp(command,"@DEFAULT") == 0 || strcmp(command,"EXIT") == 0)
    return command;

  // setup the MD simulation

  update->integrate->setup(1);

  command = mdi_fix->engine_mode("@FORCES");

  if (strcmp(command,"@DEFAULT") == 0 || strcmp(command,"EXIT") == 0)
    return command;

  // run MD one step at a time

  while (1) {
    update->whichflag = 1;
    timer->init_timeout();
    update->nsteps += 1;
    update->laststep += 1;
    update->endstep = update->laststep;
    output->next = update->ntimestep + 1;

    // single MD timestep

    update->integrate->run(1);

    // done with MD if driver sends @DEFAULT or EXIT

    command = mdi_fix->command;

    if (strcmp(command,"@DEFAULT") == 0 || strcmp(command,"EXIT") == 0)
      return command;
  }

  return NULL;
}

/* ----------------------------------------------------------------------
   perform minimization under control of driver
---------------------------------------------------------------------- */

char *MDIEngine::mdi_optg()
{
  // initialize an energy minization

  Minimize *minimizer = new Minimize(lmp);

  // setup the minimizer in a way that ensures optimization 
  // will continue until MDI driver exits

  update->etol = std::numeric_limits<double>::min();
  update->ftol = std::numeric_limits<double>::min();
  update->nsteps = std::numeric_limits<int>::max();
  update->max_eval = std::numeric_limits<int>::max();

  update->whichflag = 2;
  update->beginstep = update->firststep = update->ntimestep;
  update->endstep = update->laststep = update->firststep + update->nsteps;

  lmp->init();

  // engine is now at @INIT_OPTG node

  char *command = NULL;
  command = mdi_fix->engine_mode("@INIT_OPTG");

  if (strcmp(command,"@DEFAULT") == 0 || strcmp(command,"EXIT") == 0)
    return command;

  // setup the minimization

  update->minimize->setup();

  // get new command

  command = mdi_fix->command;

  if (strcmp(command,"@DEFAULT") == 0 || strcmp(command,"EXIT") == 0)
    return command;

  // Start a minimization, which is configured to run (essentially)
  //       infinite steps.  When the driver sends the EXIT command,
  //       the minimizer's energy and force tolerances are set to
  //       extremely large values, causing the minimization to end.

  update->minimize->iterate(update->nsteps);

  // return if driver sends @DEFAULT or EXIT

  command = mdi_fix->command;

  if (strcmp(command,"@DEFAULT") == 0 || strcmp(command,"EXIT") == 0)
    return command;

  error->all(FLERR,
             fmt::format("MDI reached end of OPTG simulation "
                         "with invalid command: {}",command));
  return NULL;
}
