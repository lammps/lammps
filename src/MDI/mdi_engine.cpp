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

#include "mdi_engine.h"

#include <limits>
#include <cstring>
#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "fix_mdi_engine.h"
#include "force.h"
#include "group.h"
#include "input.h"
//#include "irregular.h"
#include "integrate.h"
#include "library.h"
#include "memory.h"
#include "min.h"
#include "modify.h"
#include "neighbor.h"
#include "output.h"
#include "thermo.h"
#include "timer.h"
#include "update.h"

#include <mdi.h>

using namespace LAMMPS_NS;

enum{NATIVE,REAL,METAL};       // LAMMPS units which MDI supports
enum{DEFAULT,MD,OPT,SYS};      // top-level MDI engine mode

// per-atom data which engine commands access

enum{TYPE,CHARGE,MASS,COORD,VELOCITY,FORCE};  

// stages of CREATE_ATOM commands

enum{CREATE_ATOM,CREATE_ID,CREATE_TYPE,CREATE_X,CREATE_V,CREATE_IMAGE,CREATE_GO};

/* ----------------------------------------------------------------------
   mdi command: engine
   NOTE: may later have other MDI command variants?
---------------------------------------------------------------------- */

void MDIEngine::command(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal mdi command");

  if (strcmp(arg[0],"engine") == 0) mdi_engine(narg-1,&arg[1]);
  else error->all(FLERR,"Illegal mdi command");
}

/* ----------------------------------------------------------------------
   trigger LAMMPS to start acting as an MDI engine
   either in standalone mode or plugin mode
     MDI_Init() for standalone mode is in main.cpp
     MDI_Init() for plugin mode is in library_mdi.cpp::MDI_Plugin_init_lammps()
   endlessly loop over receiving commands from driver and responding
   when EXIT command is received, mdi engine command exits
---------------------------------------------------------------------- */

void MDIEngine::mdi_engine(int narg, char **arg)
{
  // process args

  enable_fix = 0;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"nodes") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal mdi engine command");
      if (strcmp(arg[iarg+1],"yes") == 0) enable_fix = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) enable_fix = 0;
      iarg += 2;
    } else error->all(FLERR,"Illegal mdi engine command");
  }

  // check requirements for LAMMPS to work with MDI as an engine

  if (atom->tag_enable == 0) 
    error->all(FLERR,"Cannot use MDI engine without atom IDs");

  if (atom->natoms && atom->tag_consecutive() == 0) 
    error->all(FLERR,"MDI engine requires consecutive atom IDs");

  // confirm LAMMPS is being run as an engine
  
  int role;
  MDI_Get_role(&role);
  if (role != MDI_ENGINE) 
    error->all(FLERR,"Must invoke LAMMPS as an MDI engine to use mdi engine");

  // root = 1 for proc 0, otherwise 0

  root = (comm->me == 0) ? 1 : 0;

  // MDI setup

  mdicmd = new char[MDI_COMMAND_LENGTH];
  node_engine = new char[MDI_COMMAND_LENGTH];
  strncpy(node_engine,"@DEFAULT",MDI_COMMAND_LENGTH);
  node_driver = new char[MDI_COMMAND_LENGTH];
  strncpy(node_driver,"\0",MDI_COMMAND_LENGTH);

  // create computes for KE. PE, pressure
  // pressure compute only calculates virial, no kinetic term
 
  id_ke = utils::strdup(std::string("MDI_ENGINE") + "_ke");
  modify->add_compute(fmt::format("{} all ke", id_ke));

  id_pe = utils::strdup(std::string("MDI_ENGINE") + "_pe");
  modify->add_compute(fmt::format("{} all pe", id_pe));

  id_press = utils::strdup(std::string("MDI_ENGINE") + "_press");
  modify->add_compute(fmt::format("{} all pressure NULL virial", id_press));

  int icompute_ke = modify->find_compute(id_ke);
  int icompute_pe = modify->find_compute(id_pe);
  int icompute_press = modify->find_compute(id_press);

  ke = modify->compute[icompute_ke];
  pe = modify->compute[icompute_pe];
  press = modify->compute[icompute_press];

  // set unit conversion factors

  if (strcmp(update->unit_style, "real") == 0) lmpunits = REAL;
  else if (strcmp(update->unit_style, "metal") == 0) lmpunits = METAL;
  else lmpunits = NATIVE;

  unit_conversions();

  // irregular class and data structs used by MDI

  buf1 = buf1all = nullptr;
  buf3 = buf3all = nullptr;
  ibuf1 = ibuf1all = nullptr;
  maxbuf = 0;

  // internal state of engine

  nbytes = -1;
  need_evaluation = 1;
  create_atoms_flag = 0;

  // define MDI commands that LAMMPS engine recognizes

  mdi_commands();

  // one-time operation to establish a connection with the driver
  
  MDI_Accept_communicator(&mdicomm);
  if (mdicomm <= 0) error->all(FLERR,"Unable to connect to MDI driver");

  // endless engine loop, responding to driver commands

  mode = DEFAULT;
  node_match = true;
  exit_command = false;

  while (1) {

    // top-level mdi engine only recognizes three nodes
    // DEFAULT, INIT_MD, INIT_OPTG, INIT_SYS

    engine_node("@DEFAULT");

    // MDI commands for dynamics or minimization

    if (strcmp(mdicmd,"@INIT_MD") == 0) {
      mdi_md();
      if (exit_command) break;

    } else if (strcmp(mdicmd,"@INIT_OPTG") == 0) {
      mdi_optg();
      if (exit_command) break;

    } else if (strcmp(mdicmd,"@INIT_SYS") == 0) {
      mdi_sys();
      if (exit_command) break;

    } else if (exit_command) {
      break;

    } else
      error->all(FLERR,
                 fmt::format("MDI engine exited with invalid command: {}",
                             mdicmd));
  }

  // clean up
  
  delete [] mdicmd;
  delete [] node_engine;
  delete [] node_driver;

  modify->delete_compute(id_ke);
  modify->delete_compute(id_pe);
  modify->delete_compute(id_press);

  delete [] id_ke;
  delete [] id_pe;
  delete [] id_press;

  // delete irregular;

  memory->destroy(ibuf1);
  memory->destroy(buf1);
  memory->destroy(buf3);

  memory->destroy(ibuf1all);
  memory->destroy(buf1all);
  memory->destroy(buf3all);
}

/* ----------------------------------------------------------------------
   engine is now at this MDI node
   loop over received commands so long as driver is also at this node
   return when not the case or EXIT command received
---------------------------------------------------------------------- */

void MDIEngine::engine_node(const char *node)
{
  int ierr;

  // do not process commands if engine and driver request are not the same

  strncpy(node_engine,node,MDI_COMMAND_LENGTH);

  if (strcmp(node_driver,"\0") != 0 && strcmp(node_driver,node_engine) != 0)
    node_match = false;

  // respond to commands from the driver

  while (!exit_command && node_match) {

    // read the next command from the driver
    // all procs call this, but only proc 0 receives the command

    ierr = MDI_Recv_command(mdicmd,mdicomm);
    if (ierr) error->all(FLERR,"MDI: Unable to receive command from driver");

    // broadcast command to the other MPI tasks

    MPI_Bcast(mdicmd,MDI_COMMAND_LENGTH,MPI_CHAR,0,world);

    // execute the command

    execute_command(mdicmd,mdicomm);

    // check if driver request is now different than engine node

    if (strcmp(node_driver,"\0") != 0 && strcmp(node_driver,node_engine) != 0)
      node_match = false;
  }

  // node exit was triggered so reset node_match

  node_match = true;
}

/* ----------------------------------------------------------------------
   process a single driver command
   called by engine_node() in loop
   also called by MDI itself via lib::lammps_execute_mdi_command()
     when LAMMPS is running as a plugin
---------------------------------------------------------------------- */

int MDIEngine::execute_command(const char *command, MDI_Comm mdicomm)
{
  int ierr;

  // confirm this command is supported at this node
  // otherwise is error

  int command_exists;
  if (root) {
    ierr = MDI_Check_command_exists(node_engine,command,MDI_COMM_NULL,
                                    &command_exists);
    if (ierr) 
      error->one(FLERR, 
                 "MDI: Unable to check whether current command is supported");
  }

  MPI_Bcast(&command_exists,1,MPI_INT,0,world);
  if (!command_exists)
    error->all(FLERR,"MDI: Received a command unsupported by engine node");

  // ---------------------------------------
  // respond to MDI standard commands
  // receives first, sends second, node commands third
  // ---------------------------------------

  if (strcmp(command,">CELL") == 0) {
    receive_cell();

  } else if (strcmp(command,">CELL_DISPL") == 0) {
    receive_cell_displ();

  } else if (strcmp(command,">CHARGES") == 0) {
    receive_charges();

  } else if (strcmp(command,">COORDS") == 0) {
    receive_coords();
    need_evaluation = 1;

  } else if (strcmp(command,">NATOMS") == 0) {
    receive_natoms();

  } else if (strcmp(command,">TYPES") == 0) {
    receive_types();

  } else if (strcmp(command,">VELOCITIES") == 0) {
    receive_velocities();

  // -----------------------------------------------

  } else if (strcmp(command,"<@") == 0) {
    ierr = MDI_Send(node_engine,MDI_NAME_LENGTH,MDI_CHAR,mdicomm);
    if (ierr) error->all(FLERR,"MDI: <@ data");

  } else if (strcmp(command,"<CELL") == 0) {
    send_cell();

  } else if (strcmp(command,"<CELL_DISPL") == 0) {
    send_cell_displ();

  } else if (strcmp(command,"<CHARGES") == 0) {
    send_double1(CHARGE);

  } else if (strcmp(command,"<COORDS") == 0) {
    send_double3(COORD);

  } else if (strcmp(command,"<ENERGY") == 0) {
    if (need_evaluation) {
      evaluate();
      need_evaluation = 0;
    }
    send_total_energy();

  } else if (strcmp(command,"<FORCES") == 0) {
    if (need_evaluation) {
      evaluate();
      need_evaluation = 0;
    }
    send_double3(FORCE);

  } else if (strcmp(command,"<LABELS") == 0) {
    send_labels();

  } else if (strcmp(command,"<MASSES") == 0) {
    send_double1(MASS);

  } else if (strcmp(command,"<NATOMS") == 0) {
    send_natoms();

  } else if (strcmp(command,"<PE") == 0) {
    if (need_evaluation) {
      evaluate();
      need_evaluation = 0;
    }
    send_pe();

  } else if (strcmp(command,"<STRESS") == 0) {
    if (need_evaluation) {
      evaluate();
      need_evaluation = 0;
    }
    send_stress();

  } else if (strcmp(command,"<TYPES") == 0) {
    send_int1(TYPE);

  } else if (strcmp(command,"<VELOCITIES") == 0) {
    send_double3(VELOCITY);

  // MDI node commands

  } else if (strcmp(command,"@INIT_MD") == 0) {
    if (mode != DEFAULT) 
      error->all(FLERR,"MDI: MDI engine is already performing a simulation");
    mode = MD;
    node_match = false;

  } else if (strcmp(command,"@INIT_OPTG") == 0) {
    if (mode != DEFAULT) 
      error->all(FLERR,"MDI: MDI engine is already performing a simulation");
    mode = OPT;
    node_match = false;

  } else if (strcmp(command,"@") == 0) {
    strncpy(node_driver,"\0",MDI_COMMAND_LENGTH);
    node_match = false;

  } else if (strcmp(command,"@DEFAULT") == 0) {
    mode = DEFAULT;
    strncpy(node_driver,"@DEFAULT",MDI_COMMAND_LENGTH);
    node_match = false;

    // if minimization in progress, force it to quit

    if (mode == OPT) {
      update->etol = std::numeric_limits<double>::max();
      update->ftol = std::numeric_limits<double>::max();
      update->max_eval = 0;
    }

  } else if (strcmp(command,"@COORDS") == 0) {
    strncpy(node_driver,"@COORDS",MDI_COMMAND_LENGTH);
    node_match = false;

  } else if (strcmp(command,"@FORCES") == 0) {
    strncpy(node_driver,"@FORCES",MDI_COMMAND_LENGTH);
    node_match = false;

  } else if (strcmp(command,"@ENDSTEP") == 0) {
    strncpy(node_driver,"@ENDSTEP",MDI_COMMAND_LENGTH);
    node_match = false;

  // exit command

  } else if (strcmp(command,"EXIT") == 0) {
    exit_command = true;

    // if minimization in progress, force it to quit

    if (mode == OPT) {
      update->etol = std::numeric_limits<double>::max();
      update->ftol = std::numeric_limits<double>::max();
      update->max_eval = 0;
    }

  // -------------------------------------------------------
  // custom LAMMPS commands
  // -------------------------------------------------------

  } else if (strcmp(command,"@INIT_SYS") == 0) {
    if (mode != DEFAULT) 
      error->all(FLERR,"MDI: MDI engine is already performing a simulation");
    mode = SYS;
    node_match = false;

  } else if (strcmp(command,"NBYTES") == 0) {
    nbytes_command();
  } else if (strcmp(command,"COMMAND") == 0) {
    single_command();
  } else if (strcmp(command,"COMMANDS") == 0) {
    many_commands();
  } else if (strcmp(command,"INFILE") == 0) {
    infile();
  } else if (strcmp(command,"<KE") == 0) {
    send_ke();

  // -------------------------------------------------------
  // unknown command
  // -------------------------------------------------------

  } else {
    error->all(FLERR,"MDI: Unknown command received from driver");
  }

  return 0;
}

/* ----------------------------------------------------------------------
   define which MDI commands the LAMMPS engine recognizes at each node
   both standard MDI commands and custom LAMMPS commands
   max length for a command is currently 11 chars
---------------------------------------------------------------------- */

void MDIEngine::mdi_commands()
{
  // default node, MDI standard commands

  MDI_Register_node("@DEFAULT");
  MDI_Register_command("@DEFAULT", "<@");
  MDI_Register_command("@DEFAULT", "<CELL");
  MDI_Register_command("@DEFAULT", "<CELL_DISPL");
  MDI_Register_command("@DEFAULT", "<CHARGES");
  MDI_Register_command("@DEFAULT", "<COORDS");
  MDI_Register_command("@DEFAULT", "<ENERGY");
  MDI_Register_command("@DEFAULT", "<FORCES");
  MDI_Register_command("@DEFAULT", "<LABELS");
  MDI_Register_command("@DEFAULT", "<MASSES");
  MDI_Register_command("@DEFAULT", "<NATOMS");
  MDI_Register_command("@DEFAULT", "<PE");
  MDI_Register_command("@DEFAULT", "<STRESS");
  MDI_Register_command("@DEFAULT", "<TYPES");
  MDI_Register_command("@DEFAULT", "<VELOCITIES");
  MDI_Register_command("@DEFAULT", ">CELL");
  MDI_Register_command("@DEFAULT", ">CELL_DISPL");
  MDI_Register_command("@DEFAULT", ">CHARGES");
  MDI_Register_command("@DEFAULT", ">COORDS");
  MDI_Register_command("@DEFAULT", ">NATOMS");
  MDI_Register_command("@DEFAULT", ">TYPES");
  MDI_Register_command("@DEFAULT", ">VELOCITIES");
  MDI_Register_command("@DEFAULT", "@INIT_MD");
  MDI_Register_command("@DEFAULT", "@INIT_OPTG");
  MDI_Register_command("@DEFAULT", "EXIT");

  // default node, custom commands added by LAMMPS

  MDI_Register_command("@DEFAULT", "NBYTES");
  MDI_Register_command("@DEFAULT", "COMMAND");
  MDI_Register_command("@DEFAULT", "COMMANDS");
  MDI_Register_command("@DEFAULT", "INFILE");
  MDI_Register_command("@DEFAULT", "<KE");
  MDI_Register_command("@DEFAULT", "@INIT_SYS");

  // custom nodes added by LAMMPS

  MDI_Register_node("@INIT_SYS");
  MDI_Register_command("@INIT_SYS", ">CELL");
  MDI_Register_command("@INIT_SYS", ">CELL_DISPL");
  MDI_Register_command("@INIT_SYS", ">CHARGES");
  MDI_Register_command("@INIT_SYS", ">COORDS");
  MDI_Register_command("@INIT_SYS", ">NATOMS");
  MDI_Register_command("@INIT_SYS", ">TYPES");
  MDI_Register_command("@INIT_SYS", ">VELOCITIES");
  MDI_Register_command("@INIT_SYS", "EXIT");

  // node for setting up and running a dynamics simulation

  MDI_Register_node("@INIT_MD");
  MDI_Register_command("@INIT_MD", "<@");
  MDI_Register_command("@INIT_MD", "<NATOMS");
  MDI_Register_command("@INIT_MD", "@");
  MDI_Register_command("@INIT_MD", "@DEFAULT");
  MDI_Register_command("@INIT_MD", "@COORDS");
  MDI_Register_command("@INIT_MD", "@FORCES");
  MDI_Register_command("@INIT_MD", "@ENDSTEP");
  MDI_Register_command("@INIT_MD", "EXIT");

  // node for setting up and running a minimization

  MDI_Register_node("@INIT_OPTG");
  MDI_Register_command("@INIT_OPTG", "<@");
  MDI_Register_command("@INIT_OPTG", "<NATOMS");
  MDI_Register_command("@INIT_OPTG", "@");
  MDI_Register_command("@INIT_OPTG", "@DEFAULT");
  MDI_Register_command("@INIT_OPTG", "@COORDS");
  MDI_Register_command("@INIT_OPTG", "@FORCES");
  MDI_Register_command("@INIT_OPTG", "@ENDSTEP");
  MDI_Register_command("@INIT_OPTG", "EXIT");

  // node at POST_INTEGRATE location in timestep
  // only if fix MDI/ENGINE is instantiated

  if (enable_fix) {
    MDI_Register_node("@COORDS");
    MDI_Register_command("@COORDS", "<@");
    MDI_Register_command("@COORDS", "<COORDS");
    MDI_Register_command("@COORDS", ">COORDS");
    MDI_Register_command("@COORDS", "@");
    MDI_Register_command("@COORDS", "@DEFAULT");
    MDI_Register_command("@COORDS", "@COORDS");
    MDI_Register_command("@COORDS", "@FORCES");
    MDI_Register_command("@COORDS", "@ENDSTEP");
    MDI_Register_command("@COORDS", "EXIT");
  }

  // node at POST_FORCE location in timestep
  // only if fix MDI/ENGINE is instantiated

  if (enable_fix) {
    MDI_Register_node("@FORCES");
    MDI_Register_command("@FORCES", "<@"); 
    MDI_Register_command("@FORCES", "<ENERGY");
    MDI_Register_command("@FORCES", "<FORCES");
    MDI_Register_command("@FORCES", "<PE");
    MDI_Register_command("@FORCES", "<STRESS");
    MDI_Register_callback("@FORCES", ">FORCES");
    MDI_Register_callback("@FORCES", ">+FORCES");
    MDI_Register_command("@FORCES", "<KE");
    MDI_Register_command("@FORCES", "@");
    MDI_Register_command("@FORCES", "@DEFAULT");
    MDI_Register_command("@FORCES", "@COORDS");
    MDI_Register_command("@FORCES", "@FORCES");
    MDI_Register_command("@FORCES", "@ENDSTEP");
    MDI_Register_command("@FORCES", "EXIT");
  }

  // node at END_OF_STEP location in timestep
  // only if fix MDI/ENGINE is instantiated

  if (enable_fix) {
    MDI_Register_node("@ENDSTEP");
    MDI_Register_command("@ENDSTEP", "<@");
    MDI_Register_command("@FORCES", "<ENERGY");
    MDI_Register_command("@FORCES", "<FORCES");
    MDI_Register_command("@FORCES", "<PE");
    MDI_Register_command("@FORCES", "<STRESS");
    MDI_Register_command("@ENDSTEP", "@");
    MDI_Register_command("@ENDSTEP", "@DEFAULT");
    MDI_Register_command("@ENDSTEP", "@COORDS");
    MDI_Register_command("@ENDSTEP", "@FORCES");
    MDI_Register_command("@ENDSTEP", "@ENDSTEP");
    MDI_Register_command("@ENDSTEP", "EXIT");
  }
}

/* ----------------------------------------------------------------------
   run MD simulation under control of driver one step at a time
   use of fix MDI/ENGINE allows MDI comm within timesteps
---------------------------------------------------------------------- */

void MDIEngine::mdi_md()
{
  // add an instance of fix MDI/ENGINE
  // delete the instance before this method returns

  if (!enable_fix) 
    error->all(FLERR,"MDI engine command did not enable @INIT_MD support");

  modify->add_fix("MDI_ENGINE_INTERNAL all MDI/ENGINE");
  mdi_fix = (FixMDIEngine *) modify->get_fix_by_id("MDI_ENGINE_INTERNAL");
  mdi_fix->mdi_engine = this;

  // initialize a new MD simulation

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
  // receive any commands driver wishes to send

  engine_node("@INIT_MD");

  if (strcmp(mdicmd,"@DEFAULT") == 0 || strcmp(mdicmd,"EXIT") == 0) {
    modify->delete_fix("MDI_ENGINE_INTERNAL");
    return;
  }

  // setup the MD simulation

  update->integrate->setup(1);

  // run MD one step at a time until driver sends @DEFAULT or EXIT
  // driver can communicate with LAMMPS within each timestep 
  // by sending a node command which matches a method in FixMDIEngine

  while (true) {
    update->whichflag = 1;
    timer->init_timeout();
    update->nsteps += 1;
    update->laststep += 1;
    update->endstep = update->laststep;
    output->next = update->ntimestep + 1;

    // single MD timestep

    update->integrate->run(1);

    // driver triggers end of MD loop by sending @DEFAULT or EXIT

    if (strcmp(mdicmd,"@DEFAULT") == 0 || strcmp(mdicmd,"EXIT") == 0) {
      modify->delete_fix("MDI_ENGINE_INTERNAL");
      return;
    }
  }
}

/* ----------------------------------------------------------------------
   perform minimization under control of driver one iteration at a time
   use of fix MDI/ENGINE allows MDI comm within iteration
---------------------------------------------------------------------- */

void MDIEngine::mdi_optg()
{
  // add an instance of fix MDI/ENGINE
  // delete the instance before this method returns

  if (!enable_fix) 
    error->all(FLERR,"MDI engine command did not enable @INIT_OPTG support");

  modify->add_fix("MDI_ENGINE_INTERNAL all MDI/ENGINE");
  mdi_fix = (FixMDIEngine *) modify->get_fix_by_id("MDI_ENGINE_INTERNAL");
  mdi_fix->mdi_engine = this;

  // set tolerances to epsilon and iteration limits huge
  // allows MDI driver to force an exit

  update->etol = std::numeric_limits<double>::min();
  update->ftol = std::numeric_limits<double>::min();
  update->nsteps = std::numeric_limits<int>::max();
  update->max_eval = std::numeric_limits<int>::max();

  update->whichflag = 2;
  update->beginstep = update->firststep = update->ntimestep;
  update->endstep = update->laststep = update->firststep + update->nsteps;

  lmp->init();

  // engine is now at @INIT_OPTG node
  // receive any commands driver wishes to send

  engine_node("@INIT_OPTG");

  if (strcmp(mdicmd,"@DEFAULT") == 0 || strcmp(mdicmd,"EXIT") == 0) {
    modify->delete_fix("MDI_ENGINE_INTERNAL");
    return;
  }

  // setup the minimization

  update->minimize->setup();

  if (strcmp(mdicmd,"@DEFAULT") == 0 || strcmp(mdicmd,"EXIT") == 0) return;

  // start minimization
  // when the driver sends @DEFAULT or EXIT minimizer tolerances are 
  // set to large values to force it to exit

  update->minimize->iterate(update->nsteps);

  // return if driver sends @DEFAULT or EXIT

  if (strcmp(mdicmd,"@DEFAULT") == 0 || strcmp(mdicmd,"EXIT") == 0) {
    modify->delete_fix("MDI_ENGINE_INTERNAL");
    return;
  }
}

/* ----------------------------------------------------------------------
   initialize a new simulation for single-point calc or dynamics or min
---------------------------------------------------------------------- */

void MDIEngine::mdi_sys()
{
  // engine is now at @INIT_SYS node
  // receive commands driver sends to define the system
  // GO command will force return from INIT_SYS

  engine_node("@INIT_SYS");

  if (strcmp(mdicmd,"@DEFAULT") == 0 || strcmp(mdicmd,"EXIT") == 0) return;

  // clear system via delete_atoms command
  // lib->reset_box()
  // lib->create_atoms()
  // then init LAMMPS for new system for MD or min

  //lmp->input->command();

  // initialize a new simulation

  update->whichflag = 1;
  timer->init_timeout();
  update->nsteps = 1;
  update->ntimestep = 0;
  update->firststep = update->ntimestep;
  update->laststep = update->ntimestep + update->nsteps;
  update->beginstep = update->firststep;
  update->endstep = update->laststep;

  lmp->init();

  // setup the MD simulation

  update->integrate->setup(1);

  // run MD one step at a time until driver sends @DEFAULT or EXIT
  // driver can communicate with LAMMPS within each timestep 
  // by sending a node command which matches a method in FixMDIEngine

  while (true) {
    update->whichflag = 1;
    timer->init_timeout();
    update->nsteps += 1;
    update->laststep += 1;
    update->endstep = update->laststep;
    output->next = update->ntimestep + 1;

    // single MD timestep

    update->integrate->run(1);

    // driver triggers end of MD loop by sending @DEFAULT or EXIT

    if (strcmp(mdicmd,"@DEFAULT") == 0 || strcmp(mdicmd,"EXIT") == 0) {
      modify->delete_fix("MDI_ENGINE_INTERNAL");
      return;
    }
  }

  // check that return is NOT GO

  engine_node("@INIT_SYS");

  if (strcmp(mdicmd,"@DEFAULT") == 0 || strcmp(mdicmd,"EXIT") == 0) return;
}

/* ----------------------------------------------------------------------
   evaluate() = compute forces, energy, pressure of current system
   usage modes:
   (1) called many times by a driver for system that is continuously evolving
   (2) called once per system by a driver that is passing in many systems
   distinguishes between:
     (a) first-time call
     (b) system needs reneighboring
     (c) system does not need reneighboring
  for (b) and (c), timestep is advanced
---------------------------------------------------------------------- */

void MDIEngine::evaluate()
{
  // NOTE: ago is not a good test
  //       caller should decide which verion of evaluate() is needed?
  // currently cannot call it interspersed with "delete atoms" calls
  //   b/c whichflag stays set
  // separate issue, need to unset whichflag when done

  if (neighbor->ago < 0) {

    update->whichflag = 1;
    lmp->init(); 
    update->integrate->setup(1);
    update->whichflag = 0;

  } else {

    // insure potential energy and virial are tallied on this step

    update->ntimestep++;
    pe->addstep(update->ntimestep);
    press->addstep(update->ntimestep);

    int nflag = neighbor->decide();
    if (nflag == 0) {
      comm->forward_comm();
      update->integrate->setup_minimal(0);
      modify->clearstep_compute();
      output->thermo->compute(1);
      modify->addstep_compute(update->ntimestep+1);
    } else {
      update->integrate->setup_minimal(1);
      modify->clearstep_compute();
      output->thermo->compute(1);
      modify->addstep_compute(update->ntimestep+1);
    }
  }
}

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// responses to ">" MDI driver commands that send data
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   >CELL command
   reset the simulation box edge vectors
   in conjunction with >CELL_DISPL this can adjust box arbitrarily
---------------------------------------------------------------------- */

void MDIEngine::receive_cell()
{
  if (mode == DEFAULT) receive_cell_default();
  else if (mode == SYS) receive_cell_sys();
}

void MDIEngine::receive_cell_default()
{
  double cell[9];

  int ierr = MDI_Recv(cell,9,MDI_DOUBLE,mdicomm);
  if (ierr) error->all(FLERR, "MDI: >CELL data");
  MPI_Bcast(cell,9,MPI_DOUBLE,0,world);

  for (int icell = 0; icell < 9; icell++) 
    cell[icell] *= mdi2lmp_length;

  // error check that edge vectors match LAMMPS triclinic requirement

  if (cell[1] != 0.0 || cell[2] != 0.0 || cell[5] != 0.0)
    error->all(FLERR,"MDI: Received cell edges are not LAMMPS compatible");

  // convert atoms to lamda coords before changing box

  domain->x2lamda(atom->nlocal);

  // convert celldata to new boxlo, boxhi, and tilt factors

  domain->boxhi[0] = cell[0] + domain->boxlo[0];
  domain->boxhi[1] = cell[4] + domain->boxlo[1];
  domain->boxhi[2] = cell[8] + domain->boxlo[2];

  domain->xy = cell[3];
  domain->xz = cell[6];
  domain->yz = cell[7];

  // reset all Domain variables that depend on box size/shape
  // convert atoms coords back to new box coords

  domain->set_global_box();
  domain->set_local_box();
  domain->lamda2x(atom->nlocal);
}

void MDIEngine::receive_cell_sys()
{
  int ierr = MDI_Recv(sys_cell,9,MDI_DOUBLE,mdicomm);
  if (ierr) error->all(FLERR, "MDI: >CELL data");
  MPI_Bcast(sys_cell,9,MPI_DOUBLE,0,world);

  for (int icell = 0; icell < 9; icell++) 
    sys_cell[icell] *= mdi2lmp_length;

  // error check that edge vectors match LAMMPS triclinic requirement

  if (sys_cell[1] != 0.0 || sys_cell[2] != 0.0 || sys_cell[5] != 0.0)
    error->all(FLERR,"MDI: Received cell edges are not LAMMPS compatible");
}

/* ----------------------------------------------------------------------
   >CELL_DISPL command
   reset simulation box origin = lower-left corner
---------------------------------------------------------------------- */

void MDIEngine::receive_cell_displ()
{
  if (mode == DEFAULT) receive_cell_displ_default();
  else if (mode == SYS) receive_cell_displ_sys();
}

void MDIEngine::receive_cell_displ_default()
{
  double cell_displ[3];
  int ierr = MDI_Recv(cell_displ,3,MDI_DOUBLE,mdicomm);
  if (ierr) 
    error->all(FLERR,"MDI: >CELL_DISPLS data");
  MPI_Bcast(cell_displ,3,MPI_DOUBLE,0,world);

  for (int icell = 0; icell < 3; icell++)
    cell_displ[icell] *= mdi2lmp_length;

  // convert atoms to lamda coords before changing box

  domain->x2lamda(atom->nlocal);

  // convert cell_displ to new boxlo and boxhi
  
  double old_boxlo[3];
  old_boxlo[0] = domain->boxlo[0];
  old_boxlo[1] = domain->boxlo[1];
  old_boxlo[2] = domain->boxlo[2];

  domain->boxlo[0] = cell_displ[0];
  domain->boxlo[1] = cell_displ[1];
  domain->boxlo[2] = cell_displ[2];

  domain->boxhi[0] += domain->boxlo[0] - old_boxlo[0];
  domain->boxhi[1] += domain->boxlo[1] - old_boxlo[1];
  domain->boxhi[2] += domain->boxlo[2] - old_boxlo[2];

  // reset all Domain variables that depend on box origin
  // convert atoms coords back to new box coords

  domain->set_global_box();
  domain->set_local_box();
  domain->lamda2x(atom->nlocal);
}

void MDIEngine::receive_cell_displ_sys()
{
  int ierr = MDI_Recv(sys_cell_displ,3,MDI_DOUBLE,mdicomm);
  if (ierr) error->all(FLERR,"MDI: >CELL_DISPLS data");
  MPI_Bcast(sys_cell_displ,3,MPI_DOUBLE,0,world);

  for (int icell = 0; icell < 3; icell++)
    sys_cell_displ[icell] *= mdi2lmp_length;
}

/* ----------------------------------------------------------------------
   >CHARGES command
---------------------------------------------------------------------- */

void MDIEngine::receive_charges()
{
  if (mode == DEFAULT) receive_double1(CHARGE);
  else if (mode == SYS) receive_charges_sys();
}

void MDIEngine::receive_charges_sys()
{
  int ierr = MDI_Recv(sys_charges,sys_natoms,MDI_DOUBLE,mdicomm);
  if (ierr) error->all(FLERR,"MDI: >CHARGES data");
  MPI_Bcast(sys_charges,sys_natoms,MPI_DOUBLE,0,world);
}

/* ----------------------------------------------------------------------
   >COORDS command
---------------------------------------------------------------------- */

void MDIEngine::receive_coords()
{
  if (mode == DEFAULT) receive_double3(COORD,0);
  else if (mode == SYS) receive_coords_sys();
}

void MDIEngine::receive_coords_sys()
{
  int ierr = MDI_Recv(sys_coords,3*sys_natoms,MDI_DOUBLE,mdicomm);
  if (ierr) error->all(FLERR,"MDI: >COORDS data");
  MPI_Bcast(sys_coords,3*sys_natoms,MPI_DOUBLE,0,world);
}

/* ----------------------------------------------------------------------
   >NATOMS command
   natoms cannot exceed 32-bit int for use with MDI
---------------------------------------------------------------------- */

void MDIEngine::receive_natoms()
{
  if (mode == DEFAULT) receive_natoms_default();
  else if (mode == SYS) receive_natoms_sys();
}

void MDIEngine::receive_natoms_default()
{
  int natoms;
  int ierr = MDI_Recv(&natoms,1,MDI_INT,mdicomm);
  if (ierr) error->all(FLERR,"MDI: >NATOMS data");
  MPI_Bcast(&natoms,1,MPI_INT,0,world);
  if (natoms < 0) error->all(FLERR,"MDI received natoms < 0");
  atom->natoms = natoms;
}

void MDIEngine::receive_natoms_sys()
{
  int ierr = MDI_Recv(&sys_natoms,1,MDI_INT,mdicomm);
  if (ierr) error->all(FLERR,"MDI: >NATOMS data");
  MPI_Bcast(&sys_natoms,1,MPI_INT,0,world);
  if (sys_natoms < 0) error->all(FLERR,"MDI received natoms < 0");
}

/* ----------------------------------------------------------------------
   >TYPES command
---------------------------------------------------------------------- */

void MDIEngine::receive_types()
{
  if (mode == DEFAULT) receive_int1(TYPE);
  else if (mode == SYS) receive_types_sys();
}

void MDIEngine::receive_types_sys()
{
  int ierr = MDI_Recv(sys_types,sys_natoms,MDI_INT,mdicomm);
  if (ierr) error->all(FLERR,"MDI: >TYPES data");
  MPI_Bcast(sys_types,sys_natoms,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   >VELOCITIES command
---------------------------------------------------------------------- */

void MDIEngine::receive_velocities()
{
  if (mode == DEFAULT) receive_double3(VELOCITY,0);
  else if (mode == SYS) receive_velocities_sys();
}

void MDIEngine::receive_velocities_sys()
{
  int ierr = MDI_Recv(sys_velocities,3*sys_natoms,MDI_DOUBLE,mdicomm);
  if (ierr) error->all(FLERR,"MDI: >VELOCITIES data");
  MPI_Bcast(sys_velocities,3*sys_natoms,MPI_DOUBLE,0,world);
}

/* ----------------------------------------------------------------------
   receive vector of 1 double for all atoms
   atoms are ordered by atomID, 1 to Natoms
   assumes all atoms already exist
   used by >CHARGES command
---------------------------------------------------------------------- */

void MDIEngine::receive_double1(int which)
{
  reallocate();

  int ierr = MDI_Recv(buf1,atom->natoms,MDI_DOUBLE,mdicomm);
  if (ierr) error->all(FLERR,"MDI: >double1 data");
  MPI_Bcast(buf1,atom->natoms,MPI_DOUBLE,0,world);

  // extract owned atom values
  // use atomID to index into ordered buf

  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  int ilocal;

  if (which == CHARGE) {
    double *q = atom->q;
    for (int i = 0; i < nlocal; i++) {
      ilocal = static_cast<int> (tag[i]) - 1;
      q[i] = buf1[ilocal];
    }
  }
}

/* ----------------------------------------------------------------------
   receive vector of 1 int for all atoms
   atoms are ordered by atomID, 1 to Natoms
   assumes all atoms already exist
   used by >TYPES command
---------------------------------------------------------------------- */

void MDIEngine::receive_int1(int which)
{
  reallocate();

  int ierr = MDI_Recv(ibuf1,atom->natoms,MDI_INT,mdicomm);
  if (ierr) error->all(FLERR,"MDI: >int1 data");
  MPI_Bcast(ibuf1,atom->natoms,MPI_INT,0,world);

  // extract owned atom values
  // use atomID to index into ordered buf

  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  int ilocal;

  if (which == TYPE) {
    int *type = atom->type;
    for (int i = 0; i < nlocal; i++) {
      ilocal = static_cast<int> (tag[i]) - 1;
      type[i] = ibuf1[ilocal];
    }
  }
}

/* ----------------------------------------------------------------------
   receive vector of 3 doubles for all atoms
   atoms are ordered by atomID, 1 to Natoms
   assumes all atoms already exist
   used by >COORDS, >FORCES, >VELOCITIES commands
   for COORD, assumes atom displacement is small
---------------------------------------------------------------------- */

void MDIEngine::receive_double3(int which, int addflag)
{
  reallocate();

  int ierr = MDI_Recv(buf3,3*atom->natoms,MDI_DOUBLE,mdicomm);
  if (ierr) error->all(FLERR,"MDI: >double3 data");
  MPI_Bcast(buf3,3*atom->natoms,MPI_DOUBLE,0,world);

  // extract owned atom values
  // use atomID to index into ordered buf

  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  int ilocal;

  if (which == COORD) {
    double **x = atom->x;
    for (int i = 0; i < nlocal; i++) {
      ilocal = static_cast<int> (tag[i]) - 1;
      x[i][0] = buf3[3*ilocal+0] * mdi2lmp_length;
      x[i][1] = buf3[3*ilocal+1] * mdi2lmp_length;
      x[i][2] = buf3[3*ilocal+2] * mdi2lmp_length;
    }
  } else if (which == FORCE) {
    if (!addflag) {
      double **f = atom->f;
      for (int i = 0; i < nlocal; i++) {
        ilocal = static_cast<int> (tag[i]) - 1;
        f[i][0] = buf3[3*ilocal+0] * mdi2lmp_force;
        f[i][1] = buf3[3*ilocal+1] * mdi2lmp_force;
        f[i][2] = buf3[3*ilocal+2] * mdi2lmp_force;
      }
    } else {
      double **f = atom->f;
      for (int i = 0; i < nlocal; i++) {
        ilocal = static_cast<int> (tag[i]) - 1;
        f[i][0] += buf3[3*ilocal+0] * mdi2lmp_force;
        f[i][1] += buf3[3*ilocal+1] * mdi2lmp_force;
        f[i][2] += buf3[3*ilocal+2] * mdi2lmp_force;
      }
    }
  } else if (which == VELOCITY) {
    double **v = atom->v;
    for (int i = 0; i < nlocal; i++) {
      ilocal = static_cast<int> (tag[i]) - 1;
      v[i][0] = buf3[3*ilocal+0] * mdi2lmp_velocity;
      v[i][1] = buf3[3*ilocal+1] * mdi2lmp_velocity;
      v[i][2] = buf3[3*ilocal+2] * mdi2lmp_velocity;
    }
  }
}

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// responses to "<" MDI driver commands that request data
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   <NATOMS command
   natoms cannot exceed 32-bit int for use with MDI
---------------------------------------------------------------------- */

void MDIEngine::send_natoms()
{
  int natoms = static_cast<int> (atom->natoms);
  int ierr = MDI_Send(&natoms,1,MDI_INT,mdicomm);
  if (ierr != 0) error->all(FLERR,"MDI: <NATOMS data");
}













/* ----------------------------------------------------------------------
   send vector of 1 double for all atoms
   atoms are ordered by atomID, 1 to Natoms
   used by <CHARGE, <MASSES commands
---------------------------------------------------------------------- */

void MDIEngine::send_double1(int which)
{
  reallocate();
  memset(buf1,0,atom->natoms*sizeof(double));

  // use atomID to index into ordered buf1

  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  int ilocal;

  if (which == CHARGE) {
    double *q = atom->q;
    for (int i = 0; i < nlocal; i++) {
      ilocal = static_cast<int> (tag[i]) - 1;
      buf1[ilocal] = q[i];
    }
  } else if (which == MASS) {
    double *mass = atom->mass;
    double *rmass = atom->rmass;
    if (rmass) {
      for (int i = 0; i < nlocal; i++) {
        ilocal = static_cast<int> (tag[i]) - 1;
        buf1[ilocal] = rmass[i];
      }
    } else {
      int *type = atom->type;
      for (int i = 0; i < nlocal; i++) {
        ilocal = static_cast<int> (tag[i]) - 1;
        buf1[ilocal] = mass[type[i]];
      }
    }
  }

  MPI_Reduce(buf1,buf1all,atom->natoms,MPI_DOUBLE,MPI_SUM,0,world);

  int ierr = MDI_Send(buf1all,atom->natoms,MDI_DOUBLE,mdicomm);
  if (ierr) error->all(FLERR,"MDI: <double1 data");
}

/* ----------------------------------------------------------------------
   send vector of 1 int for all atoms
   atoms are ordered by atomID, 1 to Natoms
   use by <TYPES command
---------------------------------------------------------------------- */

void MDIEngine::send_int1(int which)
{
  reallocate();
  memset(ibuf1,0,atom->natoms*sizeof(int));

  // use atomID to index into ordered buf1

  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  int ilocal;

  if (which == TYPE) {
    int *type = atom->type;
    for (int i = 0; i < nlocal; i++) {
      ilocal = static_cast<int> (tag[i]) - 1;
      ibuf1[ilocal] = type[i];
    }
  }

  MPI_Reduce(ibuf1,ibuf1all,atom->natoms,MPI_INT,MPI_SUM,0,world);

  int ierr = MDI_Send(ibuf1all,atom->natoms,MDI_INT,mdicomm);
  if (ierr) error->all(FLERR,"MDI: <int1 data");
}

/* ----------------------------------------------------------------------
   <COORDS, <FORCES, <VELOCITIES commands
   send vector of 3 doubles for all atoms
   atoms are ordered by atomID, 1 to Natoms
---------------------------------------------------------------------- */

void MDIEngine::send_double3(int which)
{
  reallocate();
  memset(buf3,0,3*atom->natoms*sizeof(double));

  // use atomID to index into ordered buf3

  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  int ilocal;

  if (which == COORD) {
    double **x = atom->x;
    for (int i = 0; i < nlocal; i++) {
      ilocal = static_cast<int> (tag[i]) - 1;
      buf3[3*ilocal+0] = x[i][0] * lmp2mdi_length;
      buf3[3*ilocal+1] = x[i][1] * lmp2mdi_length;
      buf3[3*ilocal+2] = x[i][2] * lmp2mdi_length;
    }
  } else if (which == FORCE) {
    double **f = atom->f;
    for (int i = 0; i < nlocal; i++) {
      ilocal = static_cast<int> (tag[i]) - 1;
      buf3[3*ilocal+0] = f[i][0] * lmp2mdi_force;
      buf3[3*ilocal+1] = f[i][1] * lmp2mdi_force;
      buf3[3*ilocal+2] = f[i][2] * lmp2mdi_force;
    }
  } else if (which == VELOCITY) {
    double **v = atom->v;
    for (int i = 0; i < nlocal; i++) {
      ilocal = static_cast<int> (tag[i]) - 1;
      buf3[3*ilocal+0] = v[i][0] * lmp2mdi_velocity;
      buf3[3*ilocal+1] = v[i][1] * lmp2mdi_velocity;
      buf3[3*ilocal+2] = v[i][2] * lmp2mdi_velocity;
    }
  }

  MPI_Reduce(buf3,buf3all,3*atom->natoms,MPI_DOUBLE,MPI_SUM,0,world);

  int ierr = MDI_Send(buf3all,3*atom->natoms,MDI_DOUBLE,mdicomm);
  if (ierr) error->all(FLERR,"MDI: <double3 data");
}

/* ----------------------------------------------------------------------
   <LABELS command
   convert atom type to string for each atom
   atoms are ordered by atomID, 1 to Natoms
---------------------------------------------------------------------- */

void MDIEngine::send_labels()
{
  char *labels = new char[atom->natoms * MDI_LABEL_LENGTH];
  memset(labels,' ',atom->natoms * MDI_LABEL_LENGTH);

  // NOTE: this loop will not work in parallel

  for (int iatom = 0; iatom < atom->natoms; iatom++) {
    std::string label = std::to_string(atom->type[iatom]);
    int label_len = std::min(int(label.length()), MDI_LABEL_LENGTH);
    strncpy(&labels[iatom * MDI_LABEL_LENGTH], label.c_str(), label_len);
  }

  int ierr = MDI_Send(labels,atom->natoms*MDI_LABEL_LENGTH,MDI_CHAR,mdicomm);
  if (ierr) error->all(FLERR,"MDI: <LABELS data");

  delete [] labels;
}

/* ----------------------------------------------------------------------
   <ENERGY command
   send total energy = PE + KE
---------------------------------------------------------------------- */

void MDIEngine::send_total_energy()
{
  double potential_energy = pe->compute_scalar();
  double kinetic_energy = ke->compute_scalar();
  double total_energy = potential_energy + kinetic_energy;
  total_energy *= lmp2mdi_energy;

  int ierr = MDI_Send(&total_energy,1,MDI_DOUBLE,mdicomm);
  if (ierr) error->all(FLERR,"MDI: <ENERGY data");
}

/* ----------------------------------------------------------------------
   <PE command
   send potential energy
---------------------------------------------------------------------- */

void MDIEngine::send_pe()
{
  double potential_energy = pe->compute_scalar();
  potential_energy *= lmp2mdi_energy;

  int ierr = MDI_Send(&potential_energy,1,MDI_DOUBLE,mdicomm);
  if (ierr) error->all(FLERR,"MDI: <PE data");
}

/* ----------------------------------------------------------------------
   <KE command
   send kinetic energy
---------------------------------------------------------------------- */

void MDIEngine::send_ke()
{
  double kinetic_energy = ke->compute_scalar();
  kinetic_energy *= lmp2mdi_energy;

  int ierr = MDI_Send(&kinetic_energy,1,MDI_DOUBLE,mdicomm);
  if (ierr) error->all(FLERR,"MDI: <KE data");
}

/* ----------------------------------------------------------------------
   <CELL command
   send simulation box edge vectors
---------------------------------------------------------------------- */

void MDIEngine::send_cell()
{
  double celldata[9];

  celldata[0] = domain->boxhi[0] - domain->boxlo[0];
  celldata[1] = 0.0;
  celldata[2] = 0.0;
  celldata[3] = domain->xy;
  celldata[4] = domain->boxhi[1] - domain->boxlo[1];
  celldata[5] = 0.0;
  celldata[6] = domain->xz;
  celldata[7] = domain->yz;
  celldata[8] = domain->boxhi[2] - domain->boxlo[2];

  for (int icell = 0; icell < 9; icell++) 
    celldata[icell] *= lmp2mdi_length;

  int ierr = MDI_Send(celldata,9,MDI_DOUBLE,mdicomm);
  if (ierr) error->all(FLERR,"MDI: <CELL data");
}

/* ----------------------------------------------------------------------
   <CELL_DISPL command
   send simulation box origin = lower-left corner
---------------------------------------------------------------------- */

void MDIEngine::send_cell_displ()
{
  double celldata[3];

  celldata[0] = domain->boxlo[0];
  celldata[1] = domain->boxlo[1];
  celldata[2] = domain->boxlo[2];

  for (int icell = 0; icell < 3; icell++)
    celldata[icell] *= lmp2mdi_length;

  int ierr = MDI_Send(celldata,3,MDI_DOUBLE,mdicomm);
  if (ierr) error->all(FLERR,"MDI: <CELL_DISPLS data");
}

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// responses to custom LAMMPS MDI commands
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   NBYTES command
   store received value in nbytes
   for use by a subsequent command, e.g. ones that send strings
---------------------------------------------------------------------- */

void MDIEngine::nbytes_command()
{
  int ierr = MDI_Recv(&nbytes,1,MDI_INT,mdicomm);
  if (ierr) error->all(FLERR,"MDI: NBYTES data");
  MPI_Bcast(&nbytes,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   COMMAND command
   store received value as string of length nbytes
   invoke as a LAMMPS command
---------------------------------------------------------------------- */

void MDIEngine::single_command()
{
  if (nbytes < 0) error->all(FLERR,"MDI: COMMAND nbytes has not been set");

  char *cmd = new char[nbytes+1];
  int ierr = MDI_Recv(cmd,nbytes+1,MDI_CHAR,mdicomm);
  if (ierr) error->all(FLERR,"MDI: COMMAND data");
  MPI_Bcast(cmd,nbytes+1,MPI_CHAR,0,world);
  cmd[nbytes] = '\0';

  lammps_command(lmp,cmd);

  delete [] cmd;
}

/* ----------------------------------------------------------------------
   COMMANDS command
   store received value as multi-line string of length nbytes
   invoke as multiple LAMMPS commands
---------------------------------------------------------------------- */

void MDIEngine::many_commands()
{
  if (nbytes < 0) error->all(FLERR,"MDI: COMMANDS nbytes has not been set");

  char *cmds = new char[nbytes+1];
  int ierr = MDI_Recv(cmds, nbytes+1, MDI_CHAR, mdicomm);
  if (ierr) error->all(FLERR,"MDI: COMMANDS data");
  MPI_Bcast(cmds,nbytes+1,MPI_CHAR,0,world);
  cmds[nbytes] = '\0';
  
  lammps_commands_string(lmp,cmds);
  
  delete [] cmds;
}

/* ----------------------------------------------------------------------
   INFILE command
   store received value as infile of length length_param
   invoke as a LAMMPS input script
---------------------------------------------------------------------- */

void MDIEngine::infile()
{
  if (nbytes < 0) error->all(FLERR,"MDI: INFILE nbytes has not been set");

  char *infile = new char[nbytes+1];
  int ierr = MDI_Recv(infile,nbytes+1,MDI_CHAR,mdicomm);
  if (ierr) error->all(FLERR,"MDI: INFILE data");
  MPI_Bcast(infile,nbytes+1,MPI_CHAR,0,world);
  infile[nbytes] = '\0';

  lammps_file(lmp,infile);

  delete [] infile;
}

/* ----------------------------------------------------------------------
   RESET_BOX command
   wrapper on library reset_box() method
   9 values = boxlo, boxhi, xy, yz, xz
   requires no atoms exist
   allows caller to define a new simulation box
---------------------------------------------------------------------- */

void MDIEngine::reset_box()
{
  int ierr;
  double values[9];

  if (atom->natoms > 0) 
    error->all(FLERR,"MDI RESET_BOX cannot be used when atoms exist");

  ierr = MDI_Recv(values,9,MDI_DOUBLE,mdicomm);
  MPI_Bcast(values,9,MPI_DOUBLE,0,world);

  lammps_reset_box(lmp,&values[0],&values[3],values[6],values[7],values[8]);
}

/* ----------------------------------------------------------------------
   CREATE_ATOM command
   wrapper on library create_atoms() method
   requires simulation box be defined
   allows caller to define a new set of atoms
     with their IDs, types, coords, velocities, image flags
   called in stages via flag
     since MDI plugin mode only allows 1 MDI Send/Recv per MDI command
   assumes current atom->natoms set by >NATOMS command is correct
---------------------------------------------------------------------- */

void MDIEngine::create_atoms(int flag)
{
  int ierr;

  // NOTE: error check on imageint = INT

  if (flag == CREATE_ATOM) {

    if (create_atoms_flag) 
      error->all(FLERR,"MDI CREATE_ATOM already in progress");

    ierr = MDI_Recv(&create_natoms,1,MDI_INT,mdicomm);
    MPI_Bcast(&create_natoms,1,MPI_INT,0,world);

    create_atoms_flag = 1;
    create_id = nullptr;
    create_type = nullptr;
    create_x = nullptr;
    create_v = nullptr;
    create_image = nullptr;

  } else if (flag == CREATE_ID) {

    if (!create_atoms_flag) error->all(FLERR,"MDI CREATE_ATOM not in progress");
    if (create_id) error->all(FLERR,"MDI CREATE_ATOM already in progress");
    
    int natom = create_natoms;
    memory->create(create_id,natom,"mdi:create_id");
    ierr = MDI_Recv(create_id,natom,MDI_INT,mdicomm);
    MPI_Bcast(create_id,natom,MPI_INT,0,world);

  } else if (flag == CREATE_TYPE) {

    if (!create_atoms_flag) error->all(FLERR,"MDI CREATE_ATOM not in progress");
    if (create_type) error->all(FLERR,"MDI CREATE_ATOM already in progress");

    int natom = create_natoms;
    if (create_type) error->all(FLERR,"MDI CREATE_ATOM already in progress");
    memory->create(create_type,natom,"mdi:create_type");
    ierr = MDI_Recv(create_type,natom,MDI_INT,mdicomm);
    MPI_Bcast(create_type,natom,MPI_INT,0,world);

  } else if (flag == CREATE_X) {

    if (!create_atoms_flag) error->all(FLERR,"MDI CREATE_ATOM not in progress");
    if (create_x) error->all(FLERR,"MDI CREATE_ATOM already in progress");

    int natom = create_natoms;
    if (create_x) error->all(FLERR,"MDI CREATE_ATOM already in progress");
    memory->create(create_x,3*natom,"mdi:create_x");
    ierr = MDI_Recv(create_x,3*natom,MDI_DOUBLE,mdicomm);
    MPI_Bcast(create_x,3*natom,MPI_DOUBLE,0,world);

  } else if (flag == CREATE_V) {

    if (!create_atoms_flag) error->all(FLERR,"MDI CREATE_ATOM not in progress");
    if (create_v) error->all(FLERR,"MDI CREATE_ATOM already in progress");

    int natom = create_natoms;
    if (create_v) error->all(FLERR,"MDI CREATE_ATOM already in progress");
    memory->create(create_v,3*natom,"mdi:create_x");
    ierr = MDI_Recv(create_v,3*natom,MDI_DOUBLE,mdicomm);
    MPI_Bcast(create_v,3*natom,MPI_DOUBLE,0,world);

  } else if (flag == CREATE_IMAGE) {

    if (!create_atoms_flag) error->all(FLERR,"MDI CREATE_ATOM not in progress");
    if (create_image) error->all(FLERR,"MDI CREATE_ATOM already in progress");

    int natom = create_natoms;
    if (create_image) error->all(FLERR,"MDI CREATE_ATOM already in progress");
    memory->create(create_image,natom,"mdi:create_image");
    ierr = MDI_Recv(create_image,natom,MDI_INT,mdicomm);
    MPI_Bcast(create_image,natom,MPI_INT,0,world);

  } else if (flag == CREATE_GO) {

    if (!create_atoms_flag) error->all(FLERR,"MDI CREATE_ATOM not in progress");
    if (!create_type || !create_x)
      error->all(FLERR,"MDI: CREATE_ATOM requires types and coords");
 
    int ncreate = lammps_create_atoms(lmp,create_natoms,create_id,create_type,
                                      create_x,create_v,create_image,1);
    
    if (ncreate != create_natoms) 
      error->all(FLERR, "MDI: CREATE ATOM created atoms != sent atoms");

    // clean up create_atoms state

    create_atoms_flag = 0;
    memory->destroy(create_id);
    memory->destroy(create_type);
    memory->destroy(create_x);
    memory->destroy(create_v);
    memory->destroy(create_image);
  }
}

/* ----------------------------------------------------------------------
   <STRESS command
   send 6-component stress tensor (no kinetic energy term)
---------------------------------------------------------------------- */

void MDIEngine::send_stress()
{
  double vtensor[6];
  press->compute_vector();
  for (int i = 0; i < 6; i++)
    vtensor[i] = press->vector[i] * lmp2mdi_pressure;

  int ierr = MDI_Send(vtensor,6,MDI_DOUBLE,mdicomm);
  if (ierr) error->all(FLERR,"MDI: <STRESS data");
}

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// additional methods
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   reallocate storage for all atoms if necessary
------------------------------------------------------------------------- */

void MDIEngine::reallocate()
{
  if (atom->natoms <= maxbuf) return;

  if (3*atom->natoms > MAXSMALLINT)
    error->all(FLERR,"Natoms too large to use with mdi engine");

  maxbuf = atom->natoms;

  memory->destroy(ibuf1);
  memory->destroy(buf1);
  memory->destroy(buf3);

  memory->destroy(ibuf1all);
  memory->destroy(buf1all);
  memory->destroy(buf3all);

  memory->create(ibuf1,maxbuf,"mdi:ibuf1");
  memory->create(buf1,maxbuf,"mdi:buf1");
  memory->create(buf3,3*maxbuf,"mdi:buf3");

  memory->create(ibuf1all,maxbuf,"mdi:ibuf1all");
  memory->create(buf1all,maxbuf,"mdi:buf1all");
  memory->create(buf3all,3*maxbuf,"mdi:buf3all");
}

/* ----------------------------------------------------------------------
   MDI to/from LAMMPS conversion factors
------------------------------------------------------------------------- */

void MDIEngine::unit_conversions()
{
  double angstrom_to_bohr,kelvin_to_hartree,ev_to_hartree,second_to_aut;

  MDI_Conversion_factor("angstrom","bohr",&angstrom_to_bohr);
  MDI_Conversion_factor("kelvin_energy","hartree",&kelvin_to_hartree);
  MDI_Conversion_factor("electron_volt","hartree",&ev_to_hartree);
  MDI_Conversion_Factor("second","atomic_unit_of_time",&second_to_aut);

  // length units

  mdi2lmp_length = 1.0;
  lmp2mdi_length = 1.0;

  if (lmpunits == REAL || lmpunits == METAL) {
    lmp2mdi_length = angstrom_to_bohr;
    mdi2lmp_length = 1.0 / angstrom_to_bohr;
  }

  // energy units

  mdi2lmp_energy = 1.0;
  lmp2mdi_energy = 1.0;

  if (lmpunits == REAL) {
    lmp2mdi_energy = kelvin_to_hartree / force->boltz;
    mdi2lmp_energy = force->boltz / kelvin_to_hartree;
  } else if (lmpunits == METAL) {
    lmp2mdi_energy = ev_to_hartree;
    mdi2lmp_energy = 1.0 / ev_to_hartree;
  }

  // force units = energy/length

  mdi2lmp_force = 1.0;
  lmp2mdi_force = 1.0;

  if (lmpunits == REAL) {
    lmp2mdi_force = (kelvin_to_hartree / force->boltz) / angstrom_to_bohr;
    mdi2lmp_force = 1.0 / lmp2mdi_force;
  } else if (lmpunits == METAL) {
    lmp2mdi_force = ev_to_hartree / angstrom_to_bohr;
    mdi2lmp_force = angstrom_to_bohr / ev_to_hartree;
  }

  // pressure or stress units = force/area = energy/volume

  mdi2lmp_pressure = 1.0;
  lmp2mdi_pressure = 1.0;

  if (lmpunits == REAL) {
    lmp2mdi_pressure = (kelvin_to_hartree / force->boltz) / 
      (angstrom_to_bohr * angstrom_to_bohr * angstrom_to_bohr) / force->nktv2p;
    mdi2lmp_pressure = 1.0 / lmp2mdi_pressure;
  } else if (lmpunits == METAL) {
    lmp2mdi_pressure = ev_to_hartree /
      (angstrom_to_bohr * angstrom_to_bohr * angstrom_to_bohr) / force->nktv2p;
    mdi2lmp_pressure = 1.0 / lmp2mdi_pressure;
  }

  // velocity units = distance/time

  mdi2lmp_velocity = 1.0;
  lmp2mdi_velocity = 1.0;

  if (lmpunits == REAL) {
    lmp2mdi_velocity = angstrom_to_bohr / (1.0e-15 * second_to_aut);
    mdi2lmp_velocity = 1.0 / lmp2mdi_velocity;
  } else if (lmpunits == METAL) {
    lmp2mdi_velocity = angstrom_to_bohr / (1.0e-12 * second_to_aut);
    mdi2lmp_velocity = 1.0 / lmp2mdi_velocity;
  }
}
