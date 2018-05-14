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

/* ---------------------------------------------------------------------- */

CommandMDIEngine::CommandMDIEngine(LAMMPS *lmp) : Command(lmp) {
  return;
}

CommandMDIEngine::~CommandMDIEngine() {
  return;
}

/* ---------------------------------------------------------------------- */

void CommandMDIEngine::command(int narg, char **arg)
{

  // register the default node
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

  // register the MD initialization node
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

  // register the OPTG initialization node
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

  // register the pre-forces node
  MDI_Register_Node("@PRE-FORCES");
  MDI_Register_Command("@PRE-FORCES", "<@");
  MDI_Register_Command("@PRE-FORCES", "<CELL");
  MDI_Register_Command("@PRE-FORCES", "<CHARGES");
  MDI_Register_Command("@PRE-FORCES", "<COORDS");
  MDI_Register_Command("@PRE-FORCES", "<ENERGY");
  MDI_Register_Command("@PRE-FORCES", "<FORCES");
  MDI_Register_Command("@PRE-FORCES", "<KE");
  MDI_Register_Command("@PRE-FORCES", "<LABELS");
  MDI_Register_Command("@PRE-FORCES", "<MASSES");
  MDI_Register_Command("@PRE-FORCES", "<NATOMS");
  MDI_Register_Command("@PRE-FORCES", "<PE");
  MDI_Register_Command("@PRE-FORCES", ">COORDS");
  MDI_Register_Command("@PRE-FORCES", ">FORCES");
  MDI_Register_Command("@PRE-FORCES", "@");
  MDI_Register_Command("@PRE-FORCES", "@COORDS");
  MDI_Register_Command("@PRE-FORCES", "@DEFAULT");
  MDI_Register_Command("@PRE-FORCES", "@FORCES");
  MDI_Register_Command("@PRE-FORCES", "@PRE-FORCES");
  MDI_Register_Command("@PRE-FORCES", "EXIT");

  // register the forces node
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

  // register the coordinates node
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


  // identify the mdi_engine fix
  int ifix = modify->find_fix_by_style("mdi/engine");
  if (ifix < 0) error->all(FLERR,"The mdi_engine command requires the mdi/engine fix");
  mdi_fix = static_cast<FixMDIEngine*>(modify->fix[ifix]);

  /* format for MDI Engine command:
   * mdi_engine
   */
  if (narg > 0) error->all(FLERR,"Illegal MDI command");

  if (atom->tag_enable == 0)
    error->all(FLERR,"Cannot use MDI command without atom IDs");

  if (atom->tag_consecutive() == 0)
    error->all(FLERR,"MDI command requires consecutive atom IDs");

  // begin engine_mode
  char *command = NULL;
  while ( true ) {
    // listen for MDI commands at the default command
    // the response to most MDI commands is handled here
    command = mdi_fix->engine_mode("@DEFAULT");

    // MDI commands that involve large-scale program flow are handled here
    if (strcmp(command,"@INIT_MD") == 0 ) {
      // enter MD control loop
      int received_exit = mdi_md();
      if ( received_exit == 1 ) {
	return;
      }
    }
    if (strcmp(command,"@INIT_OPTG") == 0 ) {
      // enter minimizer control loop
      int received_exit = mdi_optg();
      if ( received_exit == 1 ) {
	return;
      }
    }
    else if (strcmp(command,"EXIT") == 0 ) {
      return;
    }
    else {
      error->all(FLERR,fmt::format("MDI node exited with invalid command: {}",command));
    }
  }

  return;
}



int CommandMDIEngine::mdi_md()
{
  // initialize an MD simulation
  update->whichflag = 1; // 1 for dynamics
  timer->init_timeout();
  update->nsteps = 1;
  update->ntimestep = 0;
  update->firststep = update->ntimestep;
  update->laststep = update->ntimestep + update->nsteps;
  update->beginstep = update->firststep;
  update->endstep = update->laststep;
  lmp->init();

  // the MD simulation is now at the @INIT_MD node
  char *command = NULL;
  command = mdi_fix->engine_mode("@INIT_MD");

  if (strcmp(command,"@DEFAULT") == 0 ) {
    // return, and flag for @DEFAULT node
    return 0;
  }
  else if (strcmp(command,"EXIT") == 0 ) {
    // return, and flag for global exit
    return 1;
  }

  // continue the MD simulation
  update->integrate->setup(1);

  // the MD simulation is now at the @FORCES node
  command = mdi_fix->engine_mode("@FORCES");

  if (strcmp(command,"@DEFAULT") == 0 ) {
    // return, and flag for @DEFAULT node
    return 0;
  }
  else if (strcmp(command,"EXIT") == 0 ) {
    // return, and flag for global exit
    return 1;
  }

  // do MD iterations until told to exit
  while ( true ) {

    // run an MD timestep
    update->whichflag = 1; // 1 for dynamics
    timer->init_timeout();
    update->nsteps += 1;
    update->laststep += 1;
    update->endstep = update->laststep;
    output->next = update->ntimestep + 1;
    update->integrate->run(1);

    // get the most recent command the MDI engine received
    command = mdi_fix->command;

    if (strcmp(command,"@DEFAULT") == 0 ) {
      // return, and flag for @DEFAULT node
      return 0;
    }
    else if (strcmp(command,"EXIT") == 0 ) {
      // return, and flag for global exit
      return 1;
    }

  }

}



int CommandMDIEngine::mdi_optg()
{
  char *command = NULL;

  // create instance of the Minimizer class
  Minimize *minimizer = new Minimize(lmp);

  // initialize the minimizer in a way that ensures optimization will continue until MDI exits
  update->etol = std::numeric_limits<double>::min();
  update->ftol = std::numeric_limits<double>::min();
  update->nsteps = std::numeric_limits<int>::max();
  update->max_eval = std::numeric_limits<int>::max();

  update->whichflag = 2; // 2 for minimization
  update->beginstep = update->firststep = update->ntimestep;
  update->endstep = update->laststep = update->firststep + update->nsteps;
  lmp->init();

  command = mdi_fix->engine_mode("@INIT_OPTG");

  if (strcmp(command,"@DEFAULT") == 0 ) {
    // return, and flag for @DEFAULT node
    return 0;
  }
  else if (strcmp(command,"EXIT") == 0 ) {
    // return, and flag for global exit
    return 1;
  }

  update->minimize->setup();
  command = mdi_fix->command;

  if (strcmp(command,"@DEFAULT") == 0 ) {
    // return, and flag for @DEFAULT node
    return 0;
  }
  else if (strcmp(command,"EXIT") == 0 ) {
    // return, and flag for global exit
    return 1;
  }

  update->minimize->iterate(update->nsteps);
  command = mdi_fix->command;

  if (strcmp(command,"@DEFAULT") == 0 ) {
    // return, and flag for @DEFAULT node
    return 0;
  }
  else if (strcmp(command,"EXIT") == 0 ) {
    // return, and flag for global exit
    return 1;
  }

  error->all(FLERR,fmt::format("MDI reached end of OPTG simulation with invalid command: {}",command));
  return 0;
}
