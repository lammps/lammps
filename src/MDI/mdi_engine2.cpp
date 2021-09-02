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

#include "mdi_engine2.h"

#include <limits>
#include <cstring>
#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "fix_mdi_engine2.h"
#include "force.h"
#include "group.h"
#include "input.h"
#include "irregular.h"
#include "library.h"
#include "memory.h"
#include "min.h"
#include "minimize.h"
#include "modify.h"
#include "neighbor.h"
#include "output.h"
#include "timer.h"
#include "update.h"
#include "verlet.h"

#include <mdi.h>

using namespace LAMMPS_NS;

enum{NATIVE,REAL,METAL};   // LAMMPS units which MDI supports
enum{DEFAULT,MD,OPT};      // top-level MDI engine mode
enum{TYPE,CHARGE,MASS,COORD,VELOCITY,FORCE};

/* ----------------------------------------------------------------------
   mdi command: quit or engine
---------------------------------------------------------------------- */

void MDIEngine2::command(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR, "Illegal mdi command");

  if (strcmp(arg[0],"engine") == 0) mdi_engine(narg-1,&arg[1]);
  else error->all(FLERR, "Illegal mdi command");
}

/* ----------------------------------------------------------------------
   trigger LAMMPS to start acting as an MDI engine
   either in standalone mode or plugin mode
     MDI_Init() for standalone mode is in main.cpp
     MDI_Init() for plugin mode is in library_mdi.cpp::MDI_Plugin_init_lammps()
   endlessly loop over receiving commands from driver and responding
   when EXIT command is received, mdi engine command exits
---------------------------------------------------------------------- */

void MDIEngine2::mdi_engine(int narg, char **arg)
{
  // args to choose what MDI nodes LAMMPS will interact with
  // also whether a fix mdi/engine is needed
  // NOTE: add args

  if (narg > 0) error->all(FLERR, "Illegal mdi engine command");

  // check requirements for LAMMPS to work with MDI as an engine

  if (atom->tag_enable == 0) 
    error->all(FLERR, "Cannot use mdi engined without atom IDs");

  if (atom->natoms && atom->tag_consecutive() == 0) 
    error->all(FLERR, "mdi engine requires consecutive atom IDs");

  if (strcmp(update->unit_style, "real") == 0) lmpunits = REAL;
  else if (strcmp(update->unit_style, "metal") == 0) lmpunits = METAL;
  else lmpunits = NATIVE;

  // confirm LAMMPS is being run as an engine
  
  int role;
  MDI_Get_role(&role);
  if (role != MDI_ENGINE) 
    error->all(FLERR,
               "Must invoke LAMMPS as an MDI engine to use mdi/engine command");

  // if the mdi/engine fix is not already present, add it now
  // NOTE: make this optional above

  //int ifix = modify->find_fix_by_style("mdi/engine");
  //bool added_mdi_engine_fix = false;
  //if (ifix < 0) {
  //  modify->add_fix("MDI_ENGINE_INTERNAL all mdi/engine");
  //  added_mdi_engine_fix = true;
  // }

  // identify the mdi_engine fix

  //ifix = modify->find_fix_by_style("mdi/engine");
  //mdi_fix = static_cast<FixMDIEngine2 *>(modify->fix[ifix]);
  //mdi_fix->mdi_engine = this;

  // root = 1 for proc 0, otherwise 0

  root = (comm->me == 0) ? 1 : 0;

  // MDI setup

  mode = DEFAULT;
  node_match = true;
  exit_command = false;

  mdicmd = new char[MDI_COMMAND_LENGTH];
  node_engine = new char[MDI_COMMAND_LENGTH];
  strncpy(node_engine,"@DEFAULT",MDI_COMMAND_LENGTH);
  node_driver = new char[MDI_COMMAND_LENGTH];
  strncpy(node_driver,"\0",MDI_COMMAND_LENGTH);

  // create computes for KE. PE, pressure
 
  id_ke = utils::strdup(std::string("MDI_ENGINE") + "_ke");
  modify->add_compute(fmt::format("{} all ke", id_ke));

  id_pe = utils::strdup(std::string("MDI_ENGINE") + "_pe");
  modify->add_compute(fmt::format("{} all pe", id_pe));

  id_press = utils::strdup(std::string("MDI_ENGINE") + "_press");
  modify->add_compute(fmt::format("{} all pressure thermo_temp", id_press));

  // set unit conversion factors

  unit_conversions();

  // store pointers to the new computes

  int icompute_ke = modify->find_compute(id_ke);
  int icompute_pe = modify->find_compute(id_pe);
  int icompute_press = modify->find_compute(id_press);

  ke = modify->compute[icompute_ke];
  pe = modify->compute[icompute_pe];
  press = modify->compute[icompute_press];

  // irregular class and data structs used by MDI

  irregular = new Irregular(lmp);

  buf1 = nullptr;
  buf3 = nullptr;
  ibuf1 = nullptr;
  maxbuf = 0;

  // define MDI commands that LAMMPS engine recognizes

  mdi_commands();

  // one-time operation to establish a connection with the driver
  
  MDI_Accept_communicator(&mdicomm);
  if (mdicomm <= 0) error->all(FLERR, "Unable to connect to MDI driver");

  // endless engine loop, responding to driver commands

  while (1) {

    // top-level mdi engine only recognizes three nodes
    // DEFAULT, INIT_MD, INIT_OPTG

    engine_node("@DEFAULT");

    // MDI commands for dynamics or minimization

    if (strcmp(mdicmd,"@INIT_MD") == 0) {
      mdi_md();
      if (strcmp(mdicmd,"EXIT")) break;

    } else if (strcmp(mdicmd,"@INIT_OPTG") == 0) {
      mdi_optg();
      if (strcmp(mdicmd,"EXIT")) break;

    } else if (strcmp(mdicmd,"EXIT") == 0) {
      break;

    } else
      error->all(FLERR,fmt::format("MDI node exited with invalid command: {}"
                                   ,mdicmd));
  }

  // clean up
  
  delete [] mdicmd;
  delete [] node_engine;
  delete [] node_driver;

  modify->delete_compute(id_pe);
  modify->delete_compute(id_ke);
  modify->delete_compute(id_press);
  delete irregular;

  memory->destroy(buf3);
  memory->destroy(buf1);
  memory->destroy(ibuf1);

  // remove mdi/engine fix that mdi/engine instantiated
  // NOTE: make this optional

  //if (added_mdi_engine_fix) modify->delete_fix("MDI_ENGINE_INTERNAL");
}

/* ---------------------------------------------------------------------- */

void MDIEngine2::engine_node(const char *node)
{
  int ierr;

  // do not process commands if engine and driver are not at same node

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

    // check if driver node is now something other than engine node

    if (strcmp(node_driver,"\0") != 0 && strcmp(node_driver,node_engine) != 0)
      node_match = false;
  }

  // node exit was triggered so reset node exit flag

  node_match = true;
}

/* ----------------------------------------------------------------------
   process a single driver command
   called by engine_node() in loop
   also called by MDI itself via lib::lammps_execute_mdi_command()
     when LAMMPS is running as a plugin
---------------------------------------------------------------------- */

int MDIEngine2::execute_command(const char *command, MDI_Comm mdicomm)
{
  int ierr;

  // confirm this command is supported at this node

  int command_exists = 1;
  if (root) {
    ierr = MDI_Check_command_exists(node_engine,command,MDI_COMM_NULL,
                                    &command_exists);
    if (ierr) 
      error->one(FLERR, 
                 "MDI: Unable to check whether current command is supported");
  }

  MPI_Bcast(&command_exists,1,MPI_INT,0,world);
  if (!command_exists)
    error->all(FLERR, "MDI: Received a command unsupported by engine node");

  // respond to each possible driver command

  if (strcmp(command, ">NATOMS") == 0) {
    // NOTE: needs to be 64-bit value or copied into 64-bit value
    ierr = MDI_Recv((char *) &atom->natoms,1,MDI_INT,mdicomm);
    if (ierr) 
      error->all(FLERR, "MDI: Unable to receive number of atoms from driver");
    MPI_Bcast(&atom->natoms, 1, MPI_INT, 0, world);

  } else if (strcmp(command, "<NATOMS") == 0) {
    // NOTE: does MDI expect 32-bit value?
    int64_t mdi_natoms = atom->natoms;
    ierr = MDI_Send((char *) &mdi_natoms, 1, MDI_INT64_T, mdicomm);
    if (ierr != 0) 
      error->all(FLERR, "MDI: Unable to send number of atoms to driver");

  } else if (strcmp(command, "<NTYPES") == 0) {
    ierr = MDI_Send((char *) &atom->ntypes, 1, MDI_INT, mdicomm);
    if (ierr != 0) 
      error->all(FLERR, "MDI: Unable to send number of atom types to driver");

  } else if (strcmp(command, "<TYPES") == 0) {
    send_int1(TYPE);

  } else if (strcmp(command, "<LABELS") == 0) {
    send_labels();

  } else if (strcmp(command, "<MASSES") == 0) {
    send_double1(MASS);

  } else if (strcmp(command, "<CELL") == 0) {
    send_cell();

  } else if (strcmp(command, ">CELL") == 0) {
    receive_cell();

  } else if (strcmp(command, "<CELL_DISPL") == 0) {
    send_celldispl();

  } else if (strcmp(command, ">CELL_DISPL") == 0) {
    receive_celldispl();

  } else if (strcmp(command, ">COORDS") == 0) {
    receive_double3(COORD,0);

  } else if (strcmp(command, "<COORDS") == 0) {
    send_double3(COORD);

  } else if (strcmp(command, "<CHARGES") == 0) {
    send_double1(CHARGE);

  } else if (strcmp(command, "<ENERGY") == 0) {
    send_energy();

  } else if (strcmp(command, "<FORCES") == 0) {
    send_double3(FORCE);

  } else if (strcmp(command, ">FORCES") == 0) {
    receive_double3(FORCE,0);

  } else if (strcmp(command, ">+FORCES") == 0) {
    receive_double3(FORCE,1);

  } else if (strcmp(command, "<KE") == 0) {
    send_ke();

  } else if (strcmp(command, "<PE") == 0) {
    send_pe();

  // node commands

  } else if (strcmp(command,"@INIT_MD") == 0) {
    if (mode != DEFAULT) 
      error->all(FLERR, "MDI: MDI engine is already performing a simulation");
    mode = MD;
    node_match = false;

  } else if (strcmp(command,"@INIT_OPTG") == 0) {
    if (mode != DEFAULT) 
      error->all(FLERR, "MDI: MDI engine is already performing a simulation");
    mode = OPT;
    node_match = false;

  } else if (strcmp(command,"@") == 0) {
    strncpy(node_driver,"\0",MDI_COMMAND_LENGTH);
    node_match = false;

  } else if (strcmp(command,"<@") == 0) {
    ierr = MDI_Send(node_engine,MDI_NAME_LENGTH,MDI_CHAR,mdicomm);
    if (ierr) error->all(FLERR,"MDI: Unable to send node to driver");

  } else if (strcmp(command,"@DEFAULT") == 0) {
    mode = DEFAULT;

    strncpy(node_driver,"@DEFAULT",MDI_COMMAND_LENGTH);
    node_match = false;

    // are we in the middle of a geometry optimization?
    if (mode == OPT) {
      // ensure that the energy and force tolerances are met
      update->etol = std::numeric_limits<double>::max();
      update->ftol = std::numeric_limits<double>::max();

      // set the maximum number of force evaluations to 0
      update->max_eval = 0;
    }

  } else if (strcmp(command,"@COORDS") == 0) {
    strncpy(node_driver,"@COORDS",MDI_COMMAND_LENGTH);
    node_match = false;

  } else if (strcmp(command,"@FORCES") == 0) {
    strncpy(node_driver,"@FORCES",MDI_COMMAND_LENGTH);
    node_match = false;

  // exit command

  } else if (strcmp(command, "EXIT") == 0) {
    exit_command = true;

    // are we in the middle of a geometry optimization?
    // if so, ensure energy and force tolerances are met
    // set the maximum number of force evaluations to 0

    if (mode == OPT) {
      update->etol = std::numeric_limits<double>::max();
      update->ftol = std::numeric_limits<double>::max();
      update->max_eval = 0;
    }

  // -------------------------------------------------------
  // custom LAMMPS commands
  // -------------------------------------------------------

  } else if (strcmp(command, "LENGTH") == 0) {
    length_command();
  } else if (strcmp(command, "COMMAND") == 0) {
    single_command();
  } else if (strcmp(command, "COMMANDS") == 0) {
    many_commands();
  } else if (strcmp(command, "INFILE") == 0) {
    infile();
  } else if (strcmp(command, "EVAL") == 0) {
    evaluate();
  } else if (strcmp(command, "RESET_BOX") == 0) {
    reset_box();
  } else if (strcmp(command, "CREATE_ATOM") == 0) {
    create_atoms();
  } else if (strcmp(command, "<PRESSURE") == 0) {
    send_pressure();
  } else if (strcmp(command, "<PTENSOR") == 0) {
    send_ptensor();

  // -------------------------------------------------------
  // unknown command
  // -------------------------------------------------------

  } else {
    error->all(FLERR, "MDI: Unknown command from driver");
  }

  return 0;
}

/* ----------------------------------------------------------------------
   define which MDI commands the LAMMPS engine recognizes at which node
   standard MDI commands and custom LAMMPS commands
---------------------------------------------------------------------- */

void MDIEngine2::mdi_commands()
{
  // ------------------------------------
  // commands and nodes that an MDI-compliant MD code supports
  // ------------------------------------

  // default node and its commands

  MDI_Register_node("@DEFAULT");
  MDI_Register_command("@DEFAULT", "<@");
  MDI_Register_command("@DEFAULT", "<CELL");
  MDI_Register_command("@DEFAULT", "<CELL_DISPL");
  MDI_Register_command("@DEFAULT", "<CHARGES");
  MDI_Register_command("@DEFAULT", "<COORDS");
  MDI_Register_command("@DEFAULT", "<LABELS");
  MDI_Register_command("@DEFAULT", "<MASSES");
  MDI_Register_command("@DEFAULT", "<NATOMS");
  MDI_Register_command("@DEFAULT", "<TYPES");
  MDI_Register_command("@DEFAULT", ">CELL");
  MDI_Register_command("@DEFAULT", ">CELL_DISPL");
  MDI_Register_command("@DEFAULT", ">COORDS");
  MDI_Register_command("@DEFAULT", "@INIT_MD");
  MDI_Register_command("@DEFAULT", "@INIT_OPTG");
  MDI_Register_command("@DEFAULT", "EXIT");

  // node for setting up and running a dynamics simulation

  MDI_Register_node("@INIT_MD");
  MDI_Register_command("@INIT_MD", "<@");
  MDI_Register_command("@INIT_MD", "<CELL");
  MDI_Register_command("@INIT_MD", "<CELL_DISPL");
  MDI_Register_command("@INIT_MD", "<CHARGES");
  MDI_Register_command("@INIT_MD", "<COORDS");
  MDI_Register_command("@INIT_MD", "<ENERGY");
  MDI_Register_command("@INIT_MD", "<FORCES");
  MDI_Register_command("@INIT_MD", "<KE");
  MDI_Register_command("@INIT_MD", "<LABELS");
  MDI_Register_command("@INIT_MD", "<MASSES");
  MDI_Register_command("@INIT_MD", "<NATOMS");
  MDI_Register_command("@INIT_MD", "<PE");
  MDI_Register_command("@INIT_MD", "<TYPES");
  MDI_Register_command("@INIT_MD", ">CELL");
  MDI_Register_command("@INIT_MD", ">CELL_DISPL");
  MDI_Register_command("@INIT_MD", ">COORDS");
  MDI_Register_command("@INIT_MD", ">FORCES");
  MDI_Register_command("@INIT_MD", ">+FORCES");
  MDI_Register_command("@INIT_MD", "@");
  MDI_Register_command("@INIT_MD", "@COORDS");
  MDI_Register_command("@INIT_MD", "@DEFAULT");
  MDI_Register_command("@INIT_MD", "@FORCES");
  MDI_Register_command("@INIT_MD", "EXIT");

  // node for setting up and running a minimization

  MDI_Register_node("@INIT_OPTG");
  MDI_Register_command("@INIT_OPTG", "<@");
  MDI_Register_command("@INIT_OPTG", "<CELL");
  MDI_Register_command("@INIT_OPTG", "<CELL_DISPL");
  MDI_Register_command("@INIT_OPTG", "<CHARGES");
  MDI_Register_command("@INIT_OPTG", "<COORDS");
  MDI_Register_command("@INIT_OPTG", "<ENERGY");
  MDI_Register_command("@INIT_OPTG", "<FORCES");
  MDI_Register_command("@INIT_OPTG", "<KE");
  MDI_Register_command("@INIT_OPTG", "<LABELS");
  MDI_Register_command("@INIT_OPTG", "<MASSES");
  MDI_Register_command("@INIT_OPTG", "<NATOMS");
  MDI_Register_command("@INIT_OPTG", "<PE");
  MDI_Register_command("@INIT_OPTG", "<TYPES");
  MDI_Register_command("@INIT_OPTG", ">CELL");
  MDI_Register_command("@INIT_OPTG", ">CELL_DISPL");
  MDI_Register_command("@INIT_OPTG", ">COORDS");
  MDI_Register_command("@INIT_OPTG", ">FORCES");
  MDI_Register_command("@INIT_OPTG", ">+FORCES");
  MDI_Register_command("@INIT_OPTG", "@");
  MDI_Register_command("@INIT_OPTG", "@COORDS");
  MDI_Register_command("@INIT_OPTG", "@DEFAULT");
  MDI_Register_command("@INIT_OPTG", "@FORCES");
  MDI_Register_command("@INIT_OPTG", "EXIT");

  // node at POST_FORCE location in timestep

  MDI_Register_node("@FORCES");
  MDI_Register_callback("@FORCES", ">FORCES");
  MDI_Register_callback("@FORCES", ">+FORCES");
  MDI_Register_command("@FORCES", "<@");
  MDI_Register_command("@FORCES", "<CELL");
  MDI_Register_command("@FORCES", "<CELL_DISPL");
  MDI_Register_command("@FORCES", "<CHARGES");
  MDI_Register_command("@FORCES", "<COORDS");
  MDI_Register_command("@FORCES", "<ENERGY");
  MDI_Register_command("@FORCES", "<FORCES");
  MDI_Register_command("@FORCES", "<KE");
  MDI_Register_command("@FORCES", "<LABELS");
  MDI_Register_command("@FORCES", "<MASSES");
  MDI_Register_command("@FORCES", "<NATOMS");
  MDI_Register_command("@FORCES", "<PE");
  MDI_Register_command("@FORCES", "<TYPES");
  MDI_Register_command("@FORCES", ">CELL");
  MDI_Register_command("@FORCES", ">CELL_DISPL");
  MDI_Register_command("@FORCES", ">COORDS");
  MDI_Register_command("@FORCES", ">FORCES");
  MDI_Register_command("@FORCES", ">+FORCES");
  MDI_Register_command("@FORCES", "@");
  MDI_Register_command("@FORCES", "@COORDS");
  MDI_Register_command("@FORCES", "@DEFAULT");
  MDI_Register_command("@FORCES", "@FORCES");
  MDI_Register_command("@FORCES", "EXIT");

  // node at POST_INTEGRATE location in timestep

  MDI_Register_node("@COORDS");
  MDI_Register_command("@COORDS", "<@");
  MDI_Register_command("@COORDS", "<CELL");
  MDI_Register_command("@COORDS", "<CELL_DISPL");
  MDI_Register_command("@COORDS", "<CHARGES");
  MDI_Register_command("@COORDS", "<COORDS");
  MDI_Register_command("@COORDS", "<ENERGY");
  MDI_Register_command("@COORDS", "<FORCES");
  MDI_Register_command("@COORDS", "<KE");
  MDI_Register_command("@COORDS", "<LABELS");
  MDI_Register_command("@COORDS", "<MASSES");
  MDI_Register_command("@COORDS", "<NATOMS");
  MDI_Register_command("@COORDS", "<PE");
  MDI_Register_command("@COORDS", "<TYPES");
  MDI_Register_command("@COORDS", ">CELL");
  MDI_Register_command("@COORDS", ">CELL_DISPL");
  MDI_Register_command("@COORDS", ">COORDS");
  MDI_Register_command("@COORDS", ">FORCES");
  MDI_Register_command("@COORDS", ">+FORCES");
  MDI_Register_command("@COORDS", "@");
  MDI_Register_command("@COORDS", "@COORDS");
  MDI_Register_command("@COORDS", "@DEFAULT");
  MDI_Register_command("@COORDS", "@FORCES");
  MDI_Register_command("@COORDS", "EXIT");

  // ------------------------------------
  // custom commands and nodes which LAMMPS adds support for
  // max length for a command is current 11 chars in MDI
  // ------------------------------------

  MDI_Register_command("@DEFAULT", "LENGTH");
  MDI_Register_command("@DEFAULT", "COMMAND");
  MDI_Register_command("@DEFAULT", "COMMANDS");
  MDI_Register_command("@DEFAULT", "INFILE");
  MDI_Register_command("@DEFAULT", "EVAL");
  MDI_Register_command("@DEFAULT", "RESET_BOX");
  MDI_Register_command("@DEFAULT", "CREATE_ATOM");
  MDI_Register_command("@DEFAULT", "<KE");
  MDI_Register_command("@DEFAULT", "<PE");
  MDI_Register_command("@DEFAULT", "<PRESSURE");
  MDI_Register_command("@DEFAULT", "<PTENSOR");
  MDI_Register_command("@DEFAULT", "<FORCES");
}

/* ----------------------------------------------------------------------
   run an MD simulation under control of driver
---------------------------------------------------------------------- */

void MDIEngine2::mdi_md()
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

  engine_node("@INIT_MD");
  if (strcmp(mdicmd,"@DEFAULT") == 0 || strcmp(mdicmd,"EXIT") == 0) return;

  // setup the MD simulation

  update->integrate->setup(1);
  engine_node("@FORCES");
  if (strcmp(mdicmd,"@DEFAULT") == 0 || strcmp(mdicmd,"EXIT") == 0) return;

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

    if (strcmp(mdicmd,"@DEFAULT") == 0 || strcmp(mdicmd,"EXIT") == 0) return;
  }
}

/* ----------------------------------------------------------------------
   perform minimization under control of driver
---------------------------------------------------------------------- */

void MDIEngine2::mdi_optg()
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

  engine_node("@INIT_OPTG");

  if (strcmp(mdicmd,"@DEFAULT") == 0 || strcmp(mdicmd,"EXIT") == 0) return;

  // setup the minimization

  update->minimize->setup();

  if (strcmp(mdicmd,"@DEFAULT") == 0 || strcmp(mdicmd,"EXIT") == 0) return;

  // Start a minimization, which is configured to run (essentially)
  //       infinite steps.  When the driver sends the EXIT command,
  //       the minimizer's energy and force tolerances are set to
  //       extremely large values, causing the minimization to end.

  update->minimize->iterate(update->nsteps);

  // return if driver sends @DEFAULT or EXIT

  if (strcmp(mdicmd,"@DEFAULT") == 0 || strcmp(mdicmd,"EXIT") == 0) return;

  error->all(FLERR,
             fmt::format("MDI reached end of OPTG simulation "
                         "with invalid command: {}",mdicmd));
}

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// responses to standard MDI driver commands
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   receive vector of 1 double for all atoms
   atoms are ordered by atomID, 1 to Natoms
   assumes all atoms already exist
---------------------------------------------------------------------- */

void MDIEngine2::receive_double1(int which)
{
  reallocate();

  int ierr = MDI_Recv((char *) buf1,atom->natoms,MDI_DOUBLE,mdicomm);
  if (ierr) error->all(FLERR,"MDI: Unable to receive double1 data from driver");
  MPI_Bcast(buf1,atom->natoms,MPI_DOUBLE,0,world);

  // extract onwed atom value
  // use atomID to index into ordered buf3

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
---------------------------------------------------------------------- */

void MDIEngine2::receive_int1(int which)
{
  reallocate();

  int ierr = MDI_Recv((char *) ibuf1,atom->natoms,MDI_INT,mdicomm);
  if (ierr) error->all(FLERR,"MDI: Unable to receive double1 data from driver");
  MPI_Bcast(ibuf1,atom->natoms,MPI_INT,0,world);

  // extract onwed atom value
  // use atomID to index into ordered buf3

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
   for which = COORD, assumes atom displacement is small
---------------------------------------------------------------------- */

void MDIEngine2::receive_double3(int which, int addflag)
{
  reallocate();

  int ierr = MDI_Recv((char *) buf3,3*atom->natoms,MDI_DOUBLE,mdicomm);
  if (ierr) error->all(FLERR,"MDI: Unable to receive double3 data from driver");
  MPI_Bcast(buf3,3*atom->natoms,MPI_DOUBLE,0,world);

  // extract owned atom values
  // use atomID to index into ordered buf3

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

  /*
  // ensure atoms are in current box & update box via shrink-wrap
  // has to be be done before invoking Irregular::migrate_atoms()
  //   since it requires atoms be inside simulation box

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  domain->reset_box();
  if (domain->triclinic) domain->lamda2x(atom->nlocal);

  // move atoms to new processors via irregular() only needed if
  // migrate_check() says an atom moves too far

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  if (irregular->migrate_check()) irregular->migrate_atoms();
  if (domain->triclinic) domain->lamda2x(atom->nlocal);
  */
}

/* ----------------------------------------------------------------------
   send vector of 1 double for all atoms
   atoms are ordered by atomID, 1 to Natoms
---------------------------------------------------------------------- */

void MDIEngine2::send_double1(int which)
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

  int ierr = MDI_Send((char *) buf1all,atom->natoms,MDI_DOUBLE,mdicomm);
  if (ierr) error->all(FLERR, "MDI: Unable to send double1 data to driver");
}

/* ----------------------------------------------------------------------
   send vector of 1 int for all atoms
   atoms are ordered by atomID, 1 to Natoms
---------------------------------------------------------------------- */

void MDIEngine2::send_int1(int which)
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

  int ierr = MDI_Send((char *) ibuf1all,atom->natoms,MDI_INT,mdicomm);
  if (ierr) error->all(FLERR, "MDI: Unable to send int1 data to driver");
}

/* ----------------------------------------------------------------------
   send vector of 3 doubles for all atoms
   atoms are ordered by atomID, 1 to Natoms
---------------------------------------------------------------------- */

void MDIEngine2::send_double3(int which)
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
    double **f = atom->f;
    for (int i = 0; i < nlocal; i++) {
      ilocal = static_cast<int> (tag[i]) - 1;
      buf3[3*ilocal+0] = f[i][0] * lmp2mdi_velocity;
      buf3[3*ilocal+1] = f[i][1] * lmp2mdi_velocity;
      buf3[3*ilocal+2] = f[i][2] * lmp2mdi_velocity;
    }
  }

  MPI_Reduce(buf3,buf3all,3*atom->natoms,MPI_DOUBLE,MPI_SUM,0,world);

  int ierr = MDI_Send((char *) buf3all,3*atom->natoms,MDI_DOUBLE,mdicomm);
  if (ierr) error->all(FLERR, "MDI: Unable to send double3 data to driver");
}

/* ---------------------------------------------------------------------- */

void MDIEngine2::send_labels()
{
  char *labels = new char[atom->natoms * MDI_LABEL_LENGTH];
  memset(labels, ' ', atom->natoms * MDI_LABEL_LENGTH);

  for (int iatom = 0; iatom < atom->natoms; iatom++) {
    std::string label = std::to_string(atom->type[iatom]);
    int label_len = std::min(int(label.length()), MDI_LABEL_LENGTH);
    strncpy(&labels[iatom * MDI_LABEL_LENGTH], label.c_str(), label_len);
  }

  int ierr = MDI_Send(labels, atom->natoms * MDI_LABEL_LENGTH, MDI_CHAR, mdicomm);
  if (ierr) error->all(FLERR, "MDI: Unable to send atom labels to driver");

  delete[] labels;
}

/* ---------------------------------------------------------------------- */

void MDIEngine2::send_energy()
{
  double convert = 1.0;

  double potential_energy = pe->compute_scalar();
  double kinetic_energy = ke->compute_scalar();
  double total_energy = potential_energy + kinetic_energy;
  total_energy *= lmp2mdi_energy;

  int ierr = MDI_Send((char *) &total_energy,1,MDI_DOUBLE,mdicomm);
  if (ierr) error->all(FLERR, "MDI: Unable to send total energy to driver");
}

/* ---------------------------------------------------------------------- */

void MDIEngine2::send_pe()
{
  double potential_energy = pe->compute_scalar();
  potential_energy *= lmp2mdi_energy;

  int ierr = MDI_Send((char *) &potential_energy,1,MDI_DOUBLE,mdicomm);
  if (ierr) error->all(FLERR, "MDI: Unable to send potential energy to driver");
}

/* ---------------------------------------------------------------------- */

void MDIEngine2::send_ke()
{
  double kinetic_energy = ke->compute_scalar();
  kinetic_energy *= lmp2mdi_energy;

  int ierr = MDI_Send((char *) &kinetic_energy,1,MDI_DOUBLE,mdicomm);
  if (ierr) error->all(FLERR,"MDI: Unable to send kinetic energy to driver");
}

/* ---------------------------------------------------------------------- */

void MDIEngine2::send_cell()
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

  int ierr = MDI_Send((char *) celldata,9,MDI_DOUBLE,mdicomm);
  if (ierr) error->all(FLERR, "MDI: Unable to send cell dimensions to driver");
}

/* ---------------------------------------------------------------------- */

void MDIEngine2::receive_cell()
{
  double celldata[9];

  int ierr = MDI_Recv((char *) celldata,9,MDI_DOUBLE,mdicomm);
  if (ierr) error->all(FLERR, "MDI: Unable to send cell dimensions to driver");
  MPI_Bcast(&celldata[0],9,MPI_DOUBLE,0,world);

  for (int icell = 0; icell < 9; icell++) 
    celldata[icell] /= mdi2lmp_length;

  // require the new cell vector be orthogonal (for now)

  double small = std::numeric_limits<double>::min();
  if (abs(celldata[1]) > small or abs(celldata[2]) > small or 
      abs(celldata[3]) > small or
      abs(celldata[5]) > small or abs(celldata[6]) > small or 
      abs(celldata[7]) > small) {
    error->all(FLERR,
               "MDI: LAMMPS currently only supports the >CELL command "
               "for orthogonal cell vectors");
  }

  domain->boxhi[0] = celldata[0] + domain->boxlo[0];
  domain->boxhi[1] = celldata[4] + domain->boxlo[1];
  domain->boxhi[2] = celldata[8] + domain->boxlo[2];
  domain->xy = 0.0;
  domain->xz = 0.0;
  domain->yz = 0.0;
}

/* ---------------------------------------------------------------------- */

void MDIEngine2::send_celldispl()
{
  double celldata[3];

  celldata[0] = domain->boxlo[0];
  celldata[1] = domain->boxlo[1];
  celldata[2] = domain->boxlo[2];

  for (int icell = 0; icell < 3; icell++)
    celldata[icell] *= lmp2mdi_length;

  int ierr = MDI_Send((char *) celldata,3,MDI_DOUBLE,mdicomm);
  if (ierr) error->all(FLERR, "MDI: Unable to send cell displacement to driver");
}

/* ---------------------------------------------------------------------- */

void MDIEngine2::receive_celldispl()
{
  double celldata[3];
  int ierr = MDI_Recv((char *) celldata,3,MDI_DOUBLE,mdicomm);
  if (ierr) 
    error->all(FLERR, "MDI: Unable to receive cell displacement from driver");
  MPI_Bcast(&celldata[0],3,MPI_DOUBLE,0,world);

  for (int icell = 0; icell < 3; icell++)
    celldata[icell] *= mdi2lmp_length;

  double old_boxlo[3];
  old_boxlo[0] = domain->boxlo[0];
  old_boxlo[1] = domain->boxlo[1];
  old_boxlo[2] = domain->boxlo[2];

  // adjust the values of boxlo and boxhi for the new cell displacement vector

  domain->boxlo[0] = celldata[0];
  domain->boxlo[1] = celldata[1];
  domain->boxlo[2] = celldata[2];
  domain->boxhi[0] += domain->boxlo[0] - old_boxlo[0];
  domain->boxhi[1] += domain->boxlo[1] - old_boxlo[1];
  domain->boxhi[2] += domain->boxlo[2] - old_boxlo[2];
}

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// responses to custom LAMMPS MDI commands
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   LENGTH command
   store received value in length_param
   for use by other commands, e.g. that send strings
---------------------------------------------------------------------- */

void MDIEngine2::length_command()
{
  int ierr = MDI_Recv(&length_param,1,MDI_INT,mdicomm);
  if (ierr) error->all(FLERR, "MDI LENGTH error");
  MPI_Bcast(&length_param,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   COMMAND command
   store received value as string of length length_param
   invoke as a LAMMPS command
---------------------------------------------------------------------- */

void MDIEngine2::single_command()
{
  char *cmd = new char[length_param+1];
  int ierr = MDI_Recv(cmd,length_param+1,MDI_CHAR,mdicomm);
  if (ierr) error->all(FLERR,"MDI COMMAND error");
  MPI_Bcast(cmd,length_param+1,MPI_CHAR,0,world);
  cmd[length_param+1] = '\0';

  lammps_command(lmp,cmd);

  delete [] cmd;
}

/* ----------------------------------------------------------------------
   COMMANDS command
   store received value as multi-line string of length length_param
   invoke as multiple LAMMPS commands
---------------------------------------------------------------------- */

void MDIEngine2::many_commands()
{
  char *cmds = new char[length_param+1];
  int ierr = MDI_Recv(cmds, length_param+1, MDI_CHAR, mdicomm);
  if (ierr) error->all(FLERR, "MDI COMMANDS error");
  MPI_Bcast(cmds,length_param+1,MPI_CHAR,0,world);
  cmds[length_param+1] = '\0';
  
  lammps_commands_string(lmp,cmds);
  
  delete [] cmds;
}

/* ----------------------------------------------------------------------
   INFILE command
   store received value as infile of length length_param
   invoke as a LAMMPS input script
---------------------------------------------------------------------- */

void MDIEngine2::infile()
{
  char *infile = new char[length_param+1];
  int ierr = MDI_Recv(infile,length_param+1,MDI_CHAR,mdicomm);
  if (ierr) error->all(FLERR, "MDI INFILE error");
  MPI_Bcast(infile,length_param+1,MPI_CHAR,0,world);
  infile[length_param+1] = '\0';

  lammps_file(lmp,infile);

  delete [] infile;
}

/* ---------------------------------------------------------------------- */

void MDIEngine2::evaluate()
{
  if (neighbor->ago < 0) {
    update->whichflag = 1;
    lmp->init();
    update->integrate->setup(0);

  } else {
    pe->addstep(update->ntimestep);
    press->addstep(update->ntimestep);

    int nflag = neighbor->decide();
    if (nflag == 0) {
      comm->forward_comm();
      update->integrate->setup_minimal(0);
    } else {
      update->integrate->setup_minimal(1);
    }
  }
}

/* ---------------------------------------------------------------------- */

void MDIEngine2::reset_box()
{
  int ierr;
  double boxlo[3],boxhi[3],tilts[3];
  ierr = MDI_Recv(boxlo, 3, MDI_DOUBLE, mdicomm);
  ierr = MDI_Recv(boxhi, 3, MDI_DOUBLE, mdicomm);
  ierr = MDI_Recv(tilts, 3, MDI_DOUBLE, mdicomm);
  // ierr check?
  MPI_Bcast(boxlo, 3, MPI_DOUBLE, 0, world);
  MPI_Bcast(boxhi, 3, MPI_DOUBLE, 0, world);
  MPI_Bcast(tilts, 3, MPI_DOUBLE, 0, world);

  lammps_reset_box(lmp,boxlo,boxhi,tilts[0],tilts[1],tilts[2]);
}

/* ---------------------------------------------------------------------- */

void MDIEngine2::create_atoms()
{
  int ierr;
  int natoms;
  ierr = MDI_Recv(&natoms, 1, MDI_INT, mdicomm);
  // ierr checks everywhere?
  MPI_Bcast(&natoms, 1, MPI_INT, 0, world);

  tagint *id = nullptr;
  int *type = nullptr;
  double *x = nullptr;
  double *v = nullptr;
  imageint *image = nullptr;

  while (1) {
    char label;
    ierr = MDI_Recv(&label, 1, MDI_CHAR, mdicomm);
    MPI_Bcast(&label, 1, MPI_CHAR, 0, world);

    if (label == '0') break;

    if (label == 'i') {
      id = new tagint[natoms];
      ierr = MDI_Recv(id, natoms, MDI_INT, mdicomm);
      MPI_Bcast(id, natoms, MPI_INT, 0, world);
    } else if (label == 't') {
      type = new int[natoms];
      ierr = MDI_Recv(type, natoms, MDI_INT, mdicomm);
      MPI_Bcast(type, natoms, MPI_INT, 0, world);
    } else if (label == 'x') {
      x = new double[3*natoms];
      ierr = MDI_Recv(x, 3*natoms, MDI_DOUBLE, mdicomm);
      MPI_Bcast(x, 3*natoms, MPI_DOUBLE, 0, world);
    } else if (label == 'v') {
      v = new double[3*natoms];
      ierr = MDI_Recv(v, 3*natoms, MDI_DOUBLE, mdicomm);
      MPI_Bcast(v, 3*natoms, MPI_DOUBLE, 0, world);
    } else if (label == 'i') { 
      image = new imageint[natoms];
      ierr = MDI_Recv(image, natoms, MDI_INT, mdicomm);
      MPI_Bcast(image, natoms, MPI_INT, 0, world);
    }
  }

  if (!x || !type) 
    error->all(FLERR,"MDI create_atoms: did not receive atom coords or types");

  int ncreate = lammps_create_atoms(lmp,natoms,id,type,x,v,image,1);

  if (ncreate != natoms) 
    error->all(FLERR, "MDI create_atoms: created atoms != sent atoms");
}

/* ---------------------------------------------------------------------- */

void MDIEngine2::send_pressure()
{
  double pressure = press->compute_scalar();
  pressure *= lmp2mdi_pressure;

  int ierr = MDI_Send(&pressure, 1, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: Unable to send pressure to driver");
}

/* ---------------------------------------------------------------------- */

void MDIEngine2::send_ptensor()
{
  double pvector[6];
  press->compute_vector();
  for (int i = 0; i < 6; i++)
    pvector[i] = press->vector[i] * lmp2mdi_pressure;

  int ierr = MDI_Send(pvector,6,MDI_DOUBLE,mdicomm);
  if (ierr) error->all(FLERR, "MDI: Unable to send ptensor to driver");
}

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// utility methods
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

void MDIEngine2::reallocate()
{
  if (maxbuf > atom->natoms) return;

  if (3*atom->natoms > MAXSMALLINT)
    error->all(FLERR,"Natoms is too large to use with mdi engine");

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


void MDIEngine2::unit_conversions()
{
  double angstrom_to_bohr,kelvin_to_hartree,ev_to_hartree;

  MDI_Conversion_factor("angstrom","bohr",&angstrom_to_bohr);
  MDI_Conversion_factor("kelvin_energy","hartree",&kelvin_to_hartree);
  MDI_Conversion_factor("electron_volt","hartree",&ev_to_hartree);

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

  // force units

  mdi2lmp_force = 1.0;
  lmp2mdi_force = 1.0;

  if (lmpunits == REAL) {
    lmp2mdi_force = (kelvin_to_hartree / force->boltz) / angstrom_to_bohr;
    mdi2lmp_force = 1.0 / lmp2mdi_force;
  } else if (lmpunits == METAL) {
    lmp2mdi_force = ev_to_hartree / angstrom_to_bohr;
    mdi2lmp_force = angstrom_to_bohr / ev_to_hartree;
  }

  // pressure units

  mdi2lmp_pressure = 1.0;
  lmp2mdi_pressure = 1.0;

  // velocity units

  mdi2lmp_velocity = 1.0;
  lmp2mdi_velocity = 1.0;
}
