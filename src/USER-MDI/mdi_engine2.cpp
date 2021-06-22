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
#include "library_mdi.h"

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
#include "mdi.h"
#include "memory.h"
#include "min.h"
#include "minimize.h"
#include "modify.h"
#include "output.h"
#include "timer.h"
#include "update.h"
#include "verlet.h"

using namespace LAMMPS_NS;

enum {NONE, REAL, METAL};    // LAMMPS units which MDI supports
enum {DEFAULT, MD, OPT};    // MDI engine mode

/* ----------------------------------------------------------------------
   trigger LAMMPS to start acting as an MDI engine
   either as a standalone code or as a plugin
     MDI_Init() for standalone code is in main.cpp
     MDI_Init() for plugin is in library_mdi.cpp::MDI_Plugin_init_lammps()
   endlessly loop over receiving commands from driver and responding
   when EXIT command is received, mdi/engine command exits
---------------------------------------------------------------------- */

void MDIEngine2::command(int narg, char **arg)
{
  // args to choose what nodes LAMMPS will act with
  // also whether a fix mdi/engine is needed
  // NOTE: add args

  if (narg > 0) error->all(FLERR, "Illegal mdi/engine command");

  // check that LAMMPS defines what an MDI engine needs

  if (atom->tag_enable == 0) 
    error->all(FLERR, "Cannot use mdi/engine without atom IDs");

  if (atom->natoms && atom->tag_consecutive() == 0) 
    error->all(FLERR, "mdi/engine requires consecutive atom IDs");

  // must use real or metal units with MDI (atomic scale)
  // real: coords = Ang, eng = Kcal/mole, force = Kcal/mole/Ang
  // metal: coords = Ang, eng = eV, force = eV/Ang

  lmpunits = NONE;
  if (strcmp(update->unit_style, "real") == 0) lmpunits = REAL;
  if (strcmp(update->unit_style, "metal") == 0) lmpunits = METAL;
  if (lmpunits == NONE) error->all(FLERR, "MDI requires real or metal units");

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
  exit_flag = false;
  local_exit_flag = false;

  cmd = new char[MDI_COMMAND_LENGTH];
  node_driver = new char[MDI_COMMAND_LENGTH];
  node_engine = new char[MDI_COMMAND_LENGTH];
  strncpy(node_driver, "\0", MDI_COMMAND_LENGTH);
  strncpy(node_engine, "@DEFAULT", MDI_COMMAND_LENGTH);

  // create computes for KE. PE, pressure
 
  id_ke = utils::strdup(std::string("MDI_ENGINE") + "_ke");
  modify->add_compute(fmt::format("{} all ke", id_ke));

  id_pe = utils::strdup(std::string("MDI_ENGINE") + "_pe");
  modify->add_compute(fmt::format("{} all pe", id_pe));

  id_press = utils::strdup(std::string("MDI_ENGINE") + "_press");
  modify->add_compute(fmt::format("{} all pressure thermo_temp", id_press));

  // store pointers to the new computes

  int icompute_ke = modify->find_compute(id_ke);
  int icompute_pe = modify->find_compute(id_pe);
  int icompute_press = modify->find_compute(id_press);

  ke = modify->compute[icompute_ke];
  pe = modify->compute[icompute_pe];
  press = modify->compute[icompute_press];

  // irregular class and data structs used by MDI

  irregular = new Irregular(lmp);
  add_force = nullptr;

  // define MDI commands that LAMMPS engine recognizes

  mdi_commands();

  // one-time operation to establish a connection with the driver
  
  MDI_Accept_communicator(&mdicomm);
  if (mdicomm <= 0) error->all(FLERR, "Unable to connect to MDI driver");

  // endless engine loop, responding to driver commands

  while (1) {

    // mdi/engine command only recognizes three nodes
    // DEFAULT, INIT_MD, INIT_OPTG

    engine_node("@DEFAULT");

    // MDI commands for dynamics or minimization

    if (strcmp(cmd, "@INIT_MD") == 0) {
      mdi_md();
      if (strcmp(cmd, "EXIT")) break;

    } else if (strcmp(cmd, "@INIT_OPTG") == 0) {
      mdi_optg();
      if (strcmp(cmd, "EXIT")) break;

    } else if (strcmp(cmd, "EXIT") == 0) {
      break;

    } else
      error->all(FLERR,fmt::format("MDI node exited with invalid command: {}", cmd));
  }

  // clean up
  
  delete [] cmd;
  delete [] node_driver;
  delete [] node_engine;

  //modify->delete_compute(id_pe);
  //modify->delete_compute(id_ke);
  delete irregular;
  memory->destroy(add_force);

  // remove mdi/engine fix that mdi/engine instantiated
  // NOTE: make this optional

  //if (added_mdi_engine_fix) modify->delete_fix("MDI_ENGINE_INTERNAL");
}

/* ---------------------------------------------------------------------- */

void MDIEngine2::engine_node(const char *node)
{
  int ierr;

  // do not process commands if engine and driver are not at same node

  strncpy(node_engine, node, MDI_COMMAND_LENGTH);
  if (strcmp(node_driver, "\0") != 0 && strcmp(node_driver, node_engine) != 0)
    local_exit_flag = true;

  // respond to commands from the driver

  while (not exit_flag and not local_exit_flag) {

    // read the next command from the driver
    // all procs call this, but only proc 0 receives the command

    ierr = MDI_Recv_command(cmd, mdicomm);
    if (ierr) error->all(FLERR, "MDI: Unable to receive command from driver");

    // broadcast command to the other MPI tasks

    MPI_Bcast(cmd, MDI_COMMAND_LENGTH, MPI_CHAR, 0, world);

    // execute the command

    execute_command(cmd,mdicomm);

    // check if driver node is now something other than engine node

    if (strcmp(node_driver, "\0") != 0 && strcmp(node_driver, node_engine) != 0)
      local_exit_flag = true;
  }

  // local exit occured so turn off local exit flag

  local_exit_flag = false;
}


/* ----------------------------------------------------------------------
   process a single driver command
---------------------------------------------------------------------- */

int MDIEngine2::execute_command(const char *command, MDI_Comm mdicomm)
{
  int ierr;

  // confirm this command is supported at this node
  // NOTE: logic with ierr and command exists is faulty 

  int command_exists = 1;
  if (root) {
    ierr = MDI_Check_command_exists(node_engine, command, MDI_COMM_NULL, &command_exists);
    if (ierr) 
      error->all(FLERR, "MDI: Unable to check whether current command is supported");
  }
  if (!command_exists)
    error->all(FLERR, "MDI: Received a command unsupported at engine node");

  // respond to each possible driver command

  if (strcmp(command, ">NATOMS") == 0) {
    // NOTE: needs to be 64-bit value or copied into 64-bit value
    ierr = MDI_Recv((char *) &atom->natoms, 1, MDI_INT, mdicomm);
    if (ierr) error->all(FLERR, "MDI: Unable to receive number of atoms from driver");
    MPI_Bcast(&atom->natoms, 1, MPI_INT, 0, world);

  } else if (strcmp(command, "<NATOMS") == 0) {
    // NOTE: does MDI expect 32-bit value?
    int64_t mdi_natoms = atom->natoms;
    ierr = MDI_Send((char *) &mdi_natoms, 1, MDI_INT64_T, mdicomm);
    if (ierr != 0) error->all(FLERR, "MDI: Unable to send number of atoms to driver");

  } else if (strcmp(command, "<NTYPES") == 0) {
    ierr = MDI_Send((char *) &atom->ntypes, 1, MDI_INT, mdicomm);
    if (ierr != 0) error->all(FLERR, "MDI: Unable to send number of atom types to driver");

  } else if (strcmp(command, "<TYPES") == 0) {
    send_types();

  } else if (strcmp(command, "<LABELS") == 0) {
    send_labels();

  } else if (strcmp(command, "<MASSES") == 0) {
    send_masses();

  } else if (strcmp(command, "<CELL") == 0) {
    send_cell();

  } else if (strcmp(command, ">CELL") == 0) {
    receive_cell();

  } else if (strcmp(command, "<CELL_DISPL") == 0) {
    send_celldispl();

  } else if (strcmp(command, ">CELL_DISPL") == 0) {
    receive_celldispl();

  } else if (strcmp(command, ">COORDS") == 0) {
    receive_coordinates();

  } else if (strcmp(command, "<COORDS") == 0) {
    send_coordinates();

  } else if (strcmp(command, "<CHARGES") == 0) {
    send_charges();

  } else if (strcmp(command, "<ENERGY") == 0) {
    send_energy();

  } else if (strcmp(command, "<FORCES") == 0) {
    send_forces();

  } else if (strcmp(command, ">FORCES") == 0) {
    receive_forces(0);

  } else if (strcmp(command, ">+FORCES") == 0) {
    receive_forces(1);

  } else if (strcmp(command, "@INIT_MD") == 0) {
    if (mode != DEFAULT) 
      error->all(FLERR, "MDI: MDI is already performing a simulation");
    mode = MD;
    local_exit_flag = true;

  } else if (strcmp(command, "@INIT_OPTG") == 0) {
    if (mode != DEFAULT) 
      error->all(FLERR, "MDI: MDI is already performing a simulation");
    mode = OPT;
    local_exit_flag = true;

  } else if (strcmp(command, "@") == 0) {
    strncpy(node_driver, "\0", MDI_COMMAND_LENGTH);
    local_exit_flag = true;

  } else if (strcmp(command, "<@") == 0) {
    ierr = MDI_Send(node_engine, MDI_NAME_LENGTH, MDI_CHAR, mdicomm);
    if (ierr) error->all(FLERR, "MDI: Unable to send node to driver");

  } else if (strcmp(command, "<KE") == 0) {
    send_ke();

  } else if (strcmp(command, "<PE") == 0) {
    send_pe();

  } else if (strcmp(command, "@DEFAULT") == 0) {
    mode = DEFAULT;

    strncpy(node_driver, "@DEFAULT", MDI_COMMAND_LENGTH);
    local_exit_flag = true;

    // are we in the middle of a geometry optimization?
    if (mode == OPT) {
      // ensure that the energy and force tolerances are met
      update->etol = std::numeric_limits<double>::max();
      update->ftol = std::numeric_limits<double>::max();

      // set the maximum number of force evaluations to 0
      update->max_eval = 0;
    }

  } else if (strcmp(command, "@COORDS") == 0) {
    strncpy(node_driver, "@COORDS", MDI_COMMAND_LENGTH);
    local_exit_flag = true;

  } else if (strcmp(command, "@FORCES") == 0) {
    strncpy(node_driver, "@FORCES", MDI_COMMAND_LENGTH);
    local_exit_flag = true;

  } else if (strcmp(command, "EXIT") == 0) {
    // exit the driver code
    exit_flag = true;

    // are we in the middle of a geometry optimization?
    if (mode == OPT) {
      // ensure that the energy and force tolerances are met
      update->etol = std::numeric_limits<double>::max();
      update->ftol = std::numeric_limits<double>::max();

      // set the maximum number of force evaluations to 0
      update->max_eval = 0;
    }

  // -------------------------------------------------------
  // LAMMPS specific commands
  // -------------------------------------------------------

  } else if (strcmp(command, "COMMAND") == 0) {
    single_command();
  } else if (strcmp(command, "COMMANDS") == 0) {
    many_commands();
  } else if (strcmp(command, "INFILE") == 0) {
    infile();
  } else if (strcmp(command, "RESET_BOX") == 0) {
    reset_box();
  } else if (strcmp(command, "CREATE_ATOM") == 0) {
    create_atoms();
  } else if (strcmp(command, "<PRESSURE") == 0) {
    send_pressure();
  } else if (strcmp(command, "<VIRIAL") == 0) {
    send_virial();

  // -------------------------------------------------------
  // unknown command
  // -------------------------------------------------------

  } else {
    error->all(FLERR, "MDI: Unknown command from driver");
  }

  return 0;
}

/* ----------------------------------------------------------------------
   define which MDI commands the LAMMPS engine recognizes
   standard MDI commands and custom LAMMPS commands
---------------------------------------------------------------------- */

void MDIEngine2::mdi_commands()
{
  // NOTE: define which commands are allowed at which node

  // ------------------------------------
  // list of nodes and commands that an MDI-compliant MD code supports
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
  MDI_Register_command("@DEFAULT", "<KE");
  MDI_Register_command("@DEFAULT", "<PE");
  MDI_Register_command("@DEFAULT", "<TYPES");
  MDI_Register_command("@DEFAULT", ">CELL");
  MDI_Register_command("@DEFAULT", ">CELL_DISPL");
  MDI_Register_command("@DEFAULT", ">COORDS");
  MDI_Register_command("@DEFAULT", "@INIT_MD");
  MDI_Register_command("@DEFAULT", "@INIT_OPTG");
  MDI_Register_command("@DEFAULT", "EXIT");

  // node for setting up and running a dynamics simulation
  // NOTE: remove some of these in various nodes below

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
  // list of custom MDI nodes and commands which LAMMPS adds support for
  // ------------------------------------

  MDI_Register_command("@DEFAULT", "COMMAND");
  MDI_Register_command("@DEFAULT", "COMMANDS");
  MDI_Register_command("@DEFAULT", "INFILE");
  MDI_Register_command("@DEFAULT", "RESET_BOX");
  MDI_Register_command("@DEFAULT", "CREATE_ATOM");
  MDI_Register_command("@DEFAULT", "<PRESSURE");
  MDI_Register_command("@DEFAULT", "<VIRIAL");
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
  if (strcmp(cmd, "@DEFAULT") == 0 || strcmp(cmd, "EXIT") == 0) return;

  // setup the MD simulation

  update->integrate->setup(1);

  engine_node("@FORCES");
  if (strcmp(cmd, "@DEFAULT") == 0 || strcmp(cmd, "EXIT") == 0) return;

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

    if (strcmp(cmd, "@DEFAULT") == 0 || strcmp(cmd, "EXIT") == 0) return;
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

  if (strcmp(cmd, "@DEFAULT") == 0 || strcmp(cmd, "EXIT") == 0) return;

  // setup the minimization

  update->minimize->setup();

  if (strcmp(cmd, "@DEFAULT") == 0 || strcmp(cmd, "EXIT") == 0) return;

  // Start a minimization, which is configured to run (essentially)
  //       infinite steps.  When the driver sends the EXIT command,
  //       the minimizer's energy and force tolerances are set to
  //       extremely large values, causing the minimization to end.

  update->minimize->iterate(update->nsteps);

  // return if driver sends @DEFAULT or EXIT

  if (strcmp(cmd, "@DEFAULT") == 0 || strcmp(cmd, "EXIT") == 0) return;

  error->all(FLERR,
             fmt::format("MDI reached end of OPTG simulation "
                         "with invalid command: {}",
                         cmd));
}

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// response to individual MDI driver commands
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

void MDIEngine2::receive_coordinates()
{
  // get conversion factor to atomic units
  double posconv;

  // real: coords = Ang, eng = Kcal/mole, force = Kcal/mole/Ang
  // metal: coords = Ang, eng = eV, force = eV/Ang

  if (lmpunits == REAL) {
    double angstrom_to_bohr;
    MDI_Conversion_factor("angstrom", "bohr", &angstrom_to_bohr);
    posconv = force->angstrom / angstrom_to_bohr;
  } else if (lmpunits == METAL) {
    double angstrom_to_bohr;
    MDI_Conversion_factor("angstrom", "bohr", &angstrom_to_bohr);
    posconv = force->angstrom / angstrom_to_bohr;
  }

  // create buffer to hold all coords

  double *buffer;
  buffer = new double[3 * atom->natoms];

  int ierr = MDI_Recv((char *) buffer, 3 * atom->natoms, MDI_DOUBLE, mdicomm);
  if (ierr != 0) error->all(FLERR, "MDI: Unable to receive coordinates from driver");
  MPI_Bcast(buffer, 3 * atom->natoms, MPI_DOUBLE, 0, world);

  // pick local atoms from the buffer

  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    x[i][0] = buffer[3 * (atom->tag[i] - 1) + 0] * posconv;
    x[i][1] = buffer[3 * (atom->tag[i] - 1) + 1] * posconv;
    x[i][2] = buffer[3 * (atom->tag[i] - 1) + 2] * posconv;
  }

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

  delete[] buffer;
}

/* ---------------------------------------------------------------------- */

void MDIEngine2::send_coordinates()
{
  // get conversion factor to atomic units
  double posconv;
  if (lmpunits == REAL) {
    double angstrom_to_bohr;
    MDI_Conversion_factor("angstrom", "bohr", &angstrom_to_bohr);
    posconv = force->angstrom / angstrom_to_bohr;
  } else if (lmpunits == METAL) {
    double angstrom_to_bohr;
    MDI_Conversion_factor("angstrom", "bohr", &angstrom_to_bohr);
    posconv = force->angstrom / angstrom_to_bohr;
  }

  int64_t ncoords = 3 * atom->natoms;
  double *coords;
  double *coords_reduced;
  memory->create(coords, ncoords, "mdi/engine:coords");
  memory->create(coords_reduced, ncoords, "mdi/engine:coords_reduced");

  // zero coords

  for (int64_t icoord = 0; icoord < ncoords; icoord++) coords[icoord] = 0.0;

  // copy local atoms into buffer at correct locations

  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    coords[3 * (atom->tag[i] - 1) + 0] = x[i][0] / posconv;
    coords[3 * (atom->tag[i] - 1) + 1] = x[i][1] / posconv;
    coords[3 * (atom->tag[i] - 1) + 2] = x[i][2] / posconv;
  }

  MPI_Reduce(coords, coords_reduced, 3 * atom->natoms, MPI_DOUBLE, MPI_SUM, 0, world);

  int ierr = MDI_Send((char *) coords_reduced, 3 * atom->natoms, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: Unable to send coordinates to driver");

  memory->destroy(coords);
  memory->destroy(coords_reduced);
}

/* ---------------------------------------------------------------------- */

void MDIEngine2::send_charges()
{
  double *charges;
  double *charges_reduced;

  memory->create(charges, atom->natoms, "mdi/engine:charges");
  memory->create(charges_reduced, atom->natoms, "mdi/engine:charges_reduced");

  // zero the charges array

  for (int icharge = 0; icharge < atom->natoms; icharge++) charges[icharge] = 0.0;

  // pick local atoms from the buffer

  double *charge = atom->q;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) { charges[atom->tag[i] - 1] = charge[i]; }

  MPI_Reduce(charges, charges_reduced, atom->natoms, MPI_DOUBLE, MPI_SUM, 0, world);

  int ierr = MDI_Send((char *) charges_reduced, atom->natoms, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: Unable to send charges to driver");

  memory->destroy(charges);
  memory->destroy(charges_reduced);
}

/* ---------------------------------------------------------------------- */

void MDIEngine2::send_energy()
{
  // get conversion factor to atomic units
  double energy_conv;
  if (lmpunits == REAL) {
    double kelvin_to_hartree;
    MDI_Conversion_factor("kelvin_energy", "hartree", &kelvin_to_hartree);
    energy_conv = kelvin_to_hartree / force->boltz;
  } else if (lmpunits == METAL) {
    double ev_to_hartree;
    MDI_Conversion_factor("electron_volt", "hartree", &ev_to_hartree);
    energy_conv = ev_to_hartree;
  }

  double kelvin_to_hartree;
  MDI_Conversion_factor("kelvin_energy", "hartree", &kelvin_to_hartree);

  double potential_energy = pe->compute_scalar();
  double kinetic_energy = ke->compute_scalar();
  double total_energy;
  double *send_energy = &total_energy;

  // convert the energy to atomic units
  potential_energy *= energy_conv;
  kinetic_energy *= energy_conv;
  total_energy = potential_energy + kinetic_energy;

  int ierr = MDI_Send((char *) send_energy, 1, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: Unable to send potential energy to driver");
}

/* ---------------------------------------------------------------------- */

void MDIEngine2::send_pe()
{
  // get conversion factor to atomic units

  double energy_conv = 1.0;
  /*
  if (lmpunits == REAL) {
    double kelvin_to_hartree;
    MDI_Conversion_factor("kelvin_energy", "hartree", &kelvin_to_hartree);
    energy_conv = kelvin_to_hartree / force->boltz;
  } else if (lmpunits == METAL) {
    double ev_to_hartree;
    MDI_Conversion_factor("electron_volt", "hartree", &ev_to_hartree);
    energy_conv = ev_to_hartree;
  }
  */

  double potential_energy = pe->compute_scalar();

  // convert the energy to atomic units

  potential_energy *= energy_conv;

  int ierr = MDI_Send((char *) &potential_energy, 1, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: Unable to send potential energy to driver");
}

/* ---------------------------------------------------------------------- */

void MDIEngine2::send_ke()
{
  // get conversion factor to atomic units
  double energy_conv;
  if (lmpunits == REAL) {
    double kelvin_to_hartree;
    MDI_Conversion_factor("kelvin_energy", "hartree", &kelvin_to_hartree);
    energy_conv = kelvin_to_hartree / force->boltz;
  } else if (lmpunits == METAL) {
    double ev_to_hartree;
    MDI_Conversion_factor("electron_volt", "hartree", &ev_to_hartree);
    energy_conv = ev_to_hartree;
  }

  double kinetic_energy = ke->compute_scalar();

  // convert the energy to atomic units

  kinetic_energy *= energy_conv;

  int ierr = MDI_Send((char *) &kinetic_energy, 1, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: Unable to send potential energy to driver");
}

/* ---------------------------------------------------------------------- */

void MDIEngine2::send_types()
{
  int *const type = atom->type;

  int ierr = MDI_Send((char *) type, atom->natoms, MDI_INT, mdicomm);
  if (ierr) error->all(FLERR, "MDI: Unable to send atom types to driver");
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
  if (ierr) error->all(FLERR, "MDI: Unable to send atom types to driver");

  delete[] labels;
}

/* ---------------------------------------------------------------------- */

void MDIEngine2::send_masses()
{
  double *const rmass = atom->rmass;
  double *const mass = atom->mass;
  int *const type = atom->type;
  int nlocal = atom->nlocal;

  double *mass_by_atom;
  double *mass_by_atom_reduced;
  memory->create(mass_by_atom, atom->natoms, "mdi/engine:mass_by_atom");
  memory->create(mass_by_atom_reduced, atom->natoms, "mdi/engine:mass_by_atom_reduced");
  for (int iatom = 0; iatom < atom->natoms; iatom++) { mass_by_atom[iatom] = 0.0; }

  // determine the atomic masses

  if (rmass) {
    for (int iatom = 0; iatom < nlocal; iatom++) {
      mass_by_atom[atom->tag[iatom] - 1] = rmass[iatom];
    }
  } else {
    for (int iatom = 0; iatom < nlocal; iatom++) {
      mass_by_atom[atom->tag[iatom] - 1] = mass[type[iatom]];
    }
  }

  MPI_Reduce(mass_by_atom, mass_by_atom_reduced, atom->natoms, MPI_DOUBLE, MPI_SUM, 0, world);

  // send the atomic masses to the driver

  int ierr = MDI_Send((char *) mass_by_atom_reduced, atom->natoms, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: Unable to send atom masses to driver");

  memory->destroy(mass_by_atom);
  memory->destroy(mass_by_atom_reduced);
}

/* ---------------------------------------------------------------------- */

void MDIEngine2::send_forces()
{
  // get conversion factor to atomic units
  double force_conv;
  if (lmpunits == REAL) {
    double kelvin_to_hartree;
    double angstrom_to_bohr;
    MDI_Conversion_factor("kelvin_energy", "hartree", &kelvin_to_hartree);
    MDI_Conversion_factor("angstrom", "bohr", &angstrom_to_bohr);
    force_conv = (kelvin_to_hartree / force->boltz) * (force->angstrom / angstrom_to_bohr);
  } else if (lmpunits == METAL) {
    double ev_to_hartree;
    double angstrom_to_bohr;
    MDI_Conversion_factor("electron_volt", "hartree", &ev_to_hartree);
    MDI_Conversion_factor("angstrom", "bohr", &angstrom_to_bohr);
    force_conv = ev_to_hartree / angstrom_to_bohr;
  }

  double *forces;
  double *forces_reduced;
  double *x_buf;

  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int64_t ncoords = 3 * atom->natoms;

  memory->create(forces, ncoords, "mdi/engine:forces");
  memory->create(forces_reduced, ncoords, "mdi/engine:forces_reduced");
  x_buf = new double[3 * nlocal];

  // zero the forces array

  for (int iforce = 0; iforce < 3 * atom->natoms; iforce++) forces[iforce] = 0.0;

  // if not at a node, calculate the forces

  if (strcmp(node_engine, "@DEFAULT") == 0) {
    // certain fixes, such as shake, move the coordinates
    // to ensure that the coordinates do not change, store a copy
    double **x = atom->x;
    for (int i = 0; i < nlocal; i++) {
      x_buf[3 * i + 0] = x[i][0];
      x_buf[3 * i + 1] = x[i][1];
      x_buf[3 * i + 2] = x[i][2];
    }

    // calculate the forces
    update->whichflag = 1;    // 1 for dynamics
    update->nsteps = 1;
    lmp->init();
    update->integrate->setup_minimal(1);

    if (strcmp(node_engine, "@DEFAULT") == 0) {
      // restore the original set of coordinates
      double **x_new = atom->x;
      for (int i = 0; i < nlocal; i++) {
        x_new[i][0] = x_buf[3 * i + 0];
        x_new[i][1] = x_buf[3 * i + 1];
        x_new[i][2] = x_buf[3 * i + 2];
      }
    }
  }

  // pick local atoms from the buffer
  double **f = atom->f;
  for (int i = 0; i < nlocal; i++) {
    forces[3 * (atom->tag[i] - 1) + 0] = f[i][0] * force_conv;
    forces[3 * (atom->tag[i] - 1) + 1] = f[i][1] * force_conv;
    forces[3 * (atom->tag[i] - 1) + 2] = f[i][2] * force_conv;
  }

  // reduce the forces onto rank 0

  MPI_Reduce(forces, forces_reduced, 3 * atom->natoms, MPI_DOUBLE, MPI_SUM, 0, world);

  // send the forces through MDI

  int ierr = MDI_Send((char *) forces_reduced, 3 * atom->natoms, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: Unable to send atom forces to driver");

  memory->destroy(forces);
  memory->destroy(forces_reduced);
  delete[] x_buf;
}

/* ---------------------------------------------------------------------- */

// Receive forces from the driver
//    mode = 0: replace current forces with forces from driver
//    mode = 1: add forces from driver to current forces

void MDIEngine2::receive_forces(int mode)
{
  // get conversion factor to atomic units
  double force_conv;
  if (lmpunits == REAL) {
    double kelvin_to_hartree;
    double angstrom_to_bohr;
    MDI_Conversion_factor("kelvin_energy", "hartree", &kelvin_to_hartree);
    MDI_Conversion_factor("angstrom", "bohr", &angstrom_to_bohr);
    force_conv = (kelvin_to_hartree / force->boltz) * (force->angstrom / angstrom_to_bohr);
  } else if (lmpunits == METAL) {
    double ev_to_hartree;
    double angstrom_to_bohr;
    MDI_Conversion_factor("electron_volt", "hartree", &ev_to_hartree);
    MDI_Conversion_factor("angstrom", "bohr", &angstrom_to_bohr);
    force_conv = ev_to_hartree / angstrom_to_bohr;
  }

  int64_t ncoords = 3 * atom->natoms;
  double *forces;
  memory->create(forces, ncoords, "mdi/engine:forces");

  int ierr = MDI_Recv((char *) forces, 3 * atom->natoms, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: Unable to receive atom forces to driver");
  MPI_Bcast(forces, 3 * atom->natoms, MPI_DOUBLE, 0, world);

  // pick local atoms from the buffer
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (mode == 0) {    // Replace
    for (int i = 0; i < nlocal; i++) {
      f[i][0] = forces[3 * (atom->tag[i] - 1) + 0] / force_conv;
      f[i][1] = forces[3 * (atom->tag[i] - 1) + 1] / force_conv;
      f[i][2] = forces[3 * (atom->tag[i] - 1) + 2] / force_conv;
    }
  } else {
    for (int i = 0; i < nlocal; i++) {
      f[i][0] += forces[3 * (atom->tag[i] - 1) + 0] / force_conv;
      f[i][1] += forces[3 * (atom->tag[i] - 1) + 1] / force_conv;
      f[i][2] += forces[3 * (atom->tag[i] - 1) + 2] / force_conv;
    }
  }

  memory->destroy(forces);
}

/* ---------------------------------------------------------------------- */

void MDIEngine2::send_cell()
{
  double angstrom_to_bohr;
  MDI_Conversion_factor("angstrom", "bohr", &angstrom_to_bohr);

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

  // convert the units to bohr

  double unit_conv = force->angstrom * angstrom_to_bohr;
  for (int icell = 0; icell < 9; icell++) { celldata[icell] *= unit_conv; }

  int ierr = MDI_Send((char *) celldata, 9, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: Unable to send cell dimensions to driver");
}

/* ---------------------------------------------------------------------- */

void MDIEngine2::receive_cell()
{
  double celldata[9];

  // receive the new cell vector from the driver

  int ierr = MDI_Recv((char *) celldata, 9, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: Unable to send cell dimensions to driver");
  MPI_Bcast(&celldata[0], 9, MPI_DOUBLE, 0, world);

  double angstrom_to_bohr;
  MDI_Conversion_factor("angstrom", "bohr", &angstrom_to_bohr);
  double unit_conv = force->angstrom * angstrom_to_bohr;
  for (int icell = 0; icell < 9; icell++) { celldata[icell] /= unit_conv; }

  // ensure that the new cell vector is orthogonal

  double small = std::numeric_limits<double>::min();
  if (abs(celldata[1]) > small or abs(celldata[2]) > small or abs(celldata[3]) > small or
      abs(celldata[5]) > small or abs(celldata[6]) > small or abs(celldata[7]) > small) {
    error->all(FLERR,
               "MDI: LAMMPS currently only supports the >CELL command for orthogonal cell vectors");
  }

  // set the new LAMMPS cell dimensions
  //    This only works for orthogonal cell vectors.
  //    Supporting the more general case would be possible,
  //    but considerably more complex.

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
  double angstrom_to_bohr;
  MDI_Conversion_factor("angstrom", "bohr", &angstrom_to_bohr);

  double celldata[3];

  celldata[0] = domain->boxlo[0];
  celldata[1] = domain->boxlo[1];
  celldata[2] = domain->boxlo[2];

  // convert the units to bohr

  double unit_conv = force->angstrom * angstrom_to_bohr;
  for (int icell = 0; icell < 3; icell++) { celldata[icell] *= unit_conv; }

  int ierr = MDI_Send((char *) celldata, 3, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: Unable to send cell displacement to driver");
}

/* ---------------------------------------------------------------------- */

void MDIEngine2::receive_celldispl()
{
  // receive the cell displacement from the driver

  double celldata[3];
  int ierr = MDI_Recv((char *) celldata, 3, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: Unable to receive cell displacement from driver");
  MPI_Bcast(&celldata[0], 3, MPI_DOUBLE, 0, world);

  double angstrom_to_bohr;
  MDI_Conversion_factor("angstrom", "bohr", &angstrom_to_bohr);
  double unit_conv = force->angstrom * angstrom_to_bohr;

  double old_boxlo[3];
  old_boxlo[0] = domain->boxlo[0];
  old_boxlo[1] = domain->boxlo[1];
  old_boxlo[2] = domain->boxlo[2];

  // adjust the values of boxlo and boxhi for the new cell displacement vector
  domain->boxlo[0] = celldata[0] / unit_conv;
  domain->boxlo[1] = celldata[1] / unit_conv;
  domain->boxlo[2] = celldata[2] / unit_conv;
  domain->boxhi[0] += domain->boxlo[0] - old_boxlo[0];
  domain->boxhi[1] += domain->boxlo[1] - old_boxlo[1];
  domain->boxhi[2] += domain->boxlo[2] - old_boxlo[2];
}

/* ---------------------------------------------------------------------- */

void MDIEngine2::exchange_forces()
{
  // one-time allocation of add_force array

  if (!add_force) {
    int64_t ncoords = 3 * atom->natoms;
    memory->create(add_force, ncoords, "mdi/engine:add_force");
    for (int64_t i = 0; i < ncoords; i++) add_force[i] = 0.0;
  }

  double **f = atom->f;
  const int *const mask = atom->mask;
  const int nlocal = atom->nlocal;

  // add forces from the driver
  // NOTE: how is group for fix invoked now
  int groupbit;

  for (int i = 0; i < nlocal; ++i) {
    if (mask[i] & groupbit) {
      f[i][0] += add_force[3 * (atom->tag[i] - 1) + 0];
      f[i][1] += add_force[3 * (atom->tag[i] - 1) + 1];
      f[i][2] += add_force[3 * (atom->tag[i] - 1) + 2];
    }
  }
}

/* ---------------------------------------------------------------------- */

void MDIEngine2::single_command()
{
  int length;
  int ierr = MDI_Recv(&length, 1, MDI_INT, mdicomm);
  if (ierr) 
    error->all(FLERR, "MDI single_command: no length");
  MPI_Bcast(&length, 1, MPI_INT, 0, world);

  char *cmd = new char[length+1];
  ierr = MDI_Recv(cmd, length, MDI_CHAR, mdicomm);
  if (ierr) error->all(FLERR, "MDI single_command: did not receive command");
  MPI_Bcast(cmd, length, MPI_CHAR, 0, world);
  cmd[length] = '\0';

  lammps_command(lmp,cmd);

  delete [] cmd;
}

/* ---------------------------------------------------------------------- */

void MDIEngine2::many_commands()
{
  int length;
  int ierr = MDI_Recv(&length, 1, MDI_INT, mdicomm);
  if (ierr) 
    error->all(FLERR, "MDI many_commands: no length");
  MPI_Bcast(&length, 1, MPI_INT, 0, world);

  char *cmds = new char[length+1];
  ierr = MDI_Recv(cmds, length, MDI_CHAR, mdicomm);
  if (ierr) error->all(FLERR, "MDI many_commands: did not receive commands");
  MPI_Bcast(cmds, length, MPI_CHAR, 0, world);
  cmds[length] = '\0';

  char *ptr;
  char *cmd = cmds;

  while (*cmd) {
    ptr = strchr(cmd,'\n');
    if (ptr) *ptr = '\0';
    lammps_command(lmp,cmd);
    if (!ptr) break;
    cmd = ptr+1;
  }

  delete [] cmds;
}

/* ---------------------------------------------------------------------- */

void MDIEngine2::infile()
{
  int length;
  int ierr = MDI_Recv(&length, 1, MDI_INT, mdicomm);
  if (ierr) error->all(FLERR, "MDI infile: no length");
  MPI_Bcast(&length, 1, MPI_INT, 0, world);

  char *infile = new char[length+1];
  ierr = MDI_Recv(infile, length, MDI_CHAR, mdicomm);
  if (ierr) error->all(FLERR, "MDI infile: did not receive filename");
  MPI_Bcast(infile, length, MPI_CHAR, 0, world);
  cmd[length] = '\0';

  lammps_file(lmp,infile);

  delete [] infile;
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
  // get conversion factor to atomic units

  double pressure_conv = 1.0;

  /*
  if (lmpunits == REAL) {
    double kelvin_to_hartree;
    MDI_Conversion_factor("kelvin_energy", "hartree", &kelvin_to_hartree);
    pressure_conv = kelvin_to_hartree / force->boltz;
  } else if (lmpunits == METAL) {
    double ev_to_hartree;
    MDI_Conversion_factor("electron_volt", "hartree", &ev_to_hartree);
    pressure_conv = ev_to_hartree;
  }
  */

  double pressure = press->compute_scalar();

  // convert the pressure to atomic units

  pressure *= pressure_conv;

  int ierr = MDI_Send(&pressure, 1, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: Unable to send pressure to driver");
}

/* ---------------------------------------------------------------------- */

void MDIEngine2::send_virial()
{
  // get conversion factor to atomic units

  double pressure_conv = 1.0;

  /*
  if (lmpunits == REAL) {
    double kelvin_to_hartree;
    MDI_Conversion_factor("kelvin_energy", "hartree", &kelvin_to_hartree);
    pressure_conv = kelvin_to_hartree / force->boltz;
  } else if (lmpunits == METAL) {
    double ev_to_hartree;
    MDI_Conversion_factor("electron_volt", "hartree", &ev_to_hartree);
    pressure_conv = ev_to_hartree;
  }
  */

  press->compute_vector();

  // convert the pressure to atomic units

  //pressure *= pressure_conv;

  int ierr = MDI_Send(press->vector, 6, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: Unable to send pressure to driver");
}
