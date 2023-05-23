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

#include "mdi_engine.h"

#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "fix_mdi_engine.h"
#include "force.h"
#include "group.h"
#include "input.h"
#include "integrate.h"
#include "irregular.h"
#include "library.h"
#include "library_mdi.h"
#include "memory.h"
#include "min.h"
#include "modify.h"
#include "neighbor.h"
#include "output.h"
#include "thermo.h"
#include "timer.h"
#include "update.h"

#include <cstring>
#include <limits>

#include <mdi.h>

using namespace LAMMPS_NS;

enum { NATIVE, REAL, METAL };    // LAMMPS units which MDI supports
enum { DEFAULT, MD, OPT };       // top-level MDI engine modes

// per-atom data which engine commands access

enum { TYPE, CHARGE, MASS, COORD, VELOCITY, FORCE, ADDFORCE };

#define MAXELEMENT 118

/* ----------------------------------------------------------------------
   trigger LAMMPS to start acting as an MDI engine
   either in standalone mode or plugin mode
     MDI_Init() for standalone mode is in main.cpp
     MDI_Init() for plugin mode is in library_mdi.cpp::MDI_Plugin_init_lammps()
   endlessly loop over receiving commands from driver and responding
   when EXIT command is received, mdi engine command exits
---------------------------------------------------------------------- */

MDIEngine::MDIEngine(LAMMPS *_lmp, int narg, char **arg) : Pointers(_lmp)
{
  // check requirements for LAMMPS to work with MDI as an engine

  if (atom->tag_enable == 0) error->all(FLERR, "MDI engine requires atom IDs");

  if (atom->natoms && atom->tag_consecutive() == 0)
    error->all(FLERR, "MDI engine requires consecutive atom IDs");

  // optional args

  elements = nullptr;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "elements") == 0) {
      const char *symbols[] = {
        "H" , "He", "Li", "Be", "B" , "C" , "N" , "O" , "F" , "Ne",
        "Na", "Mg", "Al", "Si", "P" , "S" , "Cl", "Ar", "K" , "Ca",
        "Sc", "Ti", "V" , "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
        "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y" , "Zr",
        "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
        "Sb", "Te", "I" , "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
        "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
        "Lu", "Hf", "Ta", "W" , "Re", "Os", "Ir", "Pt", "Au", "Hg",
        "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
        "Pa", "U" , "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
        "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
        "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
      };

      int ntypes = atom->ntypes;
      delete[] elements;
      elements = new int[ntypes + 1];
      if (iarg + ntypes + 1 > narg) error->all(FLERR, "Illegal mdi engine command");
      for (int i = 1; i <= ntypes; i++) {
        int anum;
        for (anum = 0; anum < MAXELEMENT; anum++)
          if (strcmp(arg[iarg + i],symbols[anum]) == 0) break;
        if (anum == MAXELEMENT)
          error->all(FLERR,"Invalid chemical element in mdi engine command");
        elements[i] = anum + 1;
      }
      iarg += ntypes + 1;
    } else
      error->all(FLERR, "Illegal mdi engine command");
  }

  // error check an MDI element does not map to multiple atom types

  if (elements) {
    int ntypes = atom->ntypes;
    for (int i = 1; i < ntypes; i++)
      for (int j = i + 1; j <= ntypes; j++) {
        if (elements[i] == 0 || elements[j] == 0) continue;
        if (elements[i] == elements[j])
          error->all(FLERR, "MDI engine element cannot map to multiple types");
      }
  }

  // confirm LAMMPS is being run as an engine

  int role;
  MDI_Get_role(&role);
  if (role != MDI_ENGINE)
    error->all(FLERR, "Must invoke LAMMPS as an MDI engine to use mdi engine");

  // root = 1 for proc 0, otherwise 0

  root = (comm->me == 0) ? 1 : 0;

  // MDI setup

  mdicmd = new char[MDI_COMMAND_LENGTH];
  node_engine = new char[MDI_COMMAND_LENGTH];
  strncpy(node_engine, "@DEFAULT", MDI_COMMAND_LENGTH);
  node_driver = new char[MDI_COMMAND_LENGTH];
  strncpy(node_driver, "\0", MDI_COMMAND_LENGTH);

  // create computes for KE. PE, pressure
  // pressure compute only calculates virial, no kinetic term

  id_ke = utils::strdup(std::string("MDI_ENGINE") + "_ke");
  ke = modify->add_compute(fmt::format("{} all ke", id_ke));

  id_pe = utils::strdup(std::string("MDI_ENGINE") + "_pe");
  pe = modify->add_compute(fmt::format("{} all pe", id_pe));

  id_press = utils::strdup(std::string("MDI_ENGINE") + "_press");
  press = modify->add_compute(fmt::format("{} all pressure NULL virial", id_press));

  // irregular class used if >COORDS change dramatically

  irregular = new Irregular(lmp);

  // set unit conversion factors

  if (strcmp(update->unit_style, "real") == 0)
    lmpunits = REAL;
  else if (strcmp(update->unit_style, "metal") == 0)
    lmpunits = METAL;
  else
    lmpunits = NATIVE;

  unit_conversions();

  // internal state of engine

  flag_natoms = 0;
  flag_types = flag_charges = flag_coords = flag_velocities = 0;
  flag_cell = flag_cell_displ = 0;

  sys_types = nullptr;
  sys_coords = sys_velocities = nullptr;
  sys_charges = nullptr;

  buf1 = buf1all = nullptr;
  buf3 = buf3all = nullptr;
  ibuf1 = ibuf1all = nullptr;

  maxatom = 0;
  sys_natoms = static_cast<int>(atom->natoms);
  reallocate();

  nsteps = 0;

  etol = ftol = 1.0e-6;
  niterate = -1;
  max_eval = std::numeric_limits<int>::max();

  nbytes = -1;

  actionflag = 0;

  // define MDI commands that LAMMPS engine recognizes

  mdi_commands();

  // register a callback function with MDI used when engine runs in plugin mode
  // execute_command_plugin_wrapper() must be a static method

  MDI_Set_execute_command_func(execute_command_plugin_wrapper, this);

  // one-time operation to establish a connection with the driver

  MDI_Accept_communicator(&mdicomm);
  if (mdicomm <= 0) error->all(FLERR, "Unable to connect to MDI driver");

  // endless engine loop, responding to driver commands

  mode = DEFAULT;
  node_match = true;
  exit_command = false;

  while (true) {

    // top-level mdi engine only recognizes three nodes
    // DEFAULT, INIT_MD, INIT_OPTG

    engine_node("@DEFAULT");

    // MDI commands for dynamics or minimization

    if (strcmp(mdicmd, "@INIT_MD") == 0) {
      mdi_md();
      if (exit_command) break;

    } else if (strcmp(mdicmd, "@INIT_OPTG") == 0) {
      mdi_optg();
      if (exit_command) break;

    } else if (exit_command) {
      break;

    } else
      error->all(FLERR, fmt::format("MDI engine exited with invalid command: {}", mdicmd));
  }

  // clean up

  delete[] elements;

  delete[] mdicmd;
  delete[] node_engine;
  delete[] node_driver;

  modify->delete_compute(id_ke);
  modify->delete_compute(id_pe);
  modify->delete_compute(id_press);

  delete[] id_ke;
  delete[] id_pe;
  delete[] id_press;

  delete irregular;

  // delete buffers

  deallocate();
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

  strncpy(node_engine, node, MDI_COMMAND_LENGTH);

  if (strcmp(node_driver, "\0") != 0 && strcmp(node_driver, node_engine) != 0) node_match = false;

  // respond to commands from the driver

  while (!exit_command && node_match) {

    // read the next command from the driver
    // all procs call this, but only proc 0 receives the command

    ierr = MDI_Recv_command(mdicmd, mdicomm);
    if (ierr) error->all(FLERR, "MDI: Unable to receive command from driver");

    // broadcast command to the other MPI tasks

    MPI_Bcast(mdicmd, MDI_COMMAND_LENGTH, MPI_CHAR, 0, world);

    // execute the command

    execute_command(mdicmd, mdicomm);

    // check if driver request is now different than engine node

    if (strcmp(node_driver, "\0") != 0 && strcmp(node_driver, node_engine) != 0) node_match = false;
  }

  // node exit was triggered so reset node_match

  node_match = true;
}

/* ----------------------------------------------------------------------
   wrapper function on execute_command()
   invoked as callback by MDI when engine operates in plugin mode
   this is a static method in mdi_engine.h
---------------------------------------------------------------------- */

int MDIEngine::execute_command_plugin_wrapper(const char *command, MDI_Comm comm, void *class_obj)
{
  auto mdi_engine = (MDIEngine *) class_obj;
  return mdi_engine->execute_command(command, comm);
}

/* ----------------------------------------------------------------------
   process a single driver command
   called by engine_node() in loop when engine runs as stand-alone code
   called by execute_command_plugin_wrapper() when engine runs as plugin lib
---------------------------------------------------------------------- */

int MDIEngine::execute_command(const char *command, MDI_Comm mdicomm)
{
  int ierr;

  // confirm this command is supported at this node
  // otherwise is error

  int command_exists;
  if (root) {
    ierr = MDI_Check_command_exists(node_engine, command, MDI_COMM_NULL, &command_exists);
    if (ierr) error->one(FLERR, "MDI: Cannot confirm that command '{}' is supported", command);
  }

  MPI_Bcast(&command_exists, 1, MPI_INT, 0, world);
  if (!command_exists)
    error->all(FLERR, "MDI: Received command '{}' unsupported by engine node {}", command,
               node_engine);

  // ---------------------------------------
  // respond to MDI standard commands
  // receives first, sends second, node commands third
  // ---------------------------------------

  if (strcmp(command, ">CELL") == 0) {
    receive_cell();

  } else if (strcmp(command, ">CELL_DISPL") == 0) {
    receive_cell_displ();

  } else if (strcmp(command, ">CHARGES") == 0) {
    receive_charges();

  } else if (strcmp(command, ">COORDS") == 0) {
    receive_coords();

  } else if (strcmp(command, ">ELEMENTS") == 0) {
    if (!elements) error->all(FLERR, "MDI engine command did not define element list");
    receive_elements();

  } else if (strcmp(command, ">FORCES") == 0) {
    receive_double3(FORCE);

  } else if (strcmp(command, ">+FORCES") == 0) {
    receive_double3(ADDFORCE);

  } else if (strcmp(command, ">NATOMS") == 0) {
    receive_natoms();

  } else if (strcmp(command, ">NSTEPS") == 0) {
    receive_nsteps();

  } else if (strcmp(command, ">TOLERANCE") == 0) {
    receive_tolerance();

  } else if (strcmp(command, ">TYPES") == 0) {
    receive_types();

  } else if (strcmp(command, ">VELOCITIES") == 0) {
    if (strcmp(node_engine, "@DEFAULT") == 0)
      receive_velocities();
    else
      receive_double3(VELOCITY);

    // -----------------------------------------------

  } else if (strcmp(command, "<@") == 0) {
    ierr = MDI_Send(node_engine, MDI_NAME_LENGTH, MDI_CHAR, mdicomm);
    if (ierr) error->all(FLERR, "MDI: <@ data");

  } else if (strcmp(command, "<CELL") == 0) {
    send_cell();

  } else if (strcmp(command, "<CELL_DISPL") == 0) {
    send_cell_displ();

  } else if (strcmp(command, "<CHARGES") == 0) {
    send_double1(CHARGE);

  } else if (strcmp(command, "<COORDS") == 0) {
    send_double3(COORD);

  } else if (strcmp(command, "<ENERGY") == 0) {
    if (!actionflag && strcmp(node_engine, "@DEFAULT") == 0) evaluate();
    send_total_energy();

  } else if (strcmp(command, "<FORCES") == 0) {
    if (!actionflag && strcmp(node_engine, "@DEFAULT") == 0) evaluate();
    send_double3(FORCE);

  } else if (strcmp(command, "<LABELS") == 0) {
    send_labels();

  } else if (strcmp(command, "<MASSES") == 0) {
    send_double1(MASS);

  } else if (strcmp(command, "<NATOMS") == 0) {
    send_natoms();

  } else if (strcmp(command, "<PE") == 0) {
    if (!actionflag && strcmp(node_engine, "@DEFAULT") == 0) evaluate();
    send_pe();

  } else if (strcmp(command, "<STRESS") == 0) {
    if (!actionflag && strcmp(node_engine, "@DEFAULT") == 0) evaluate();
    send_stress();

  } else if (strcmp(command, "<TYPES") == 0) {
    send_int1(TYPE);

  } else if (strcmp(command, "<VELOCITIES") == 0) {
    send_double3(VELOCITY);

    // -----------------------------------------------

    // MDI action commands at @DEFAULT node

  } else if (strcmp(command, "MD") == 0) {
    md();

  } else if (strcmp(command, "OPTG") == 0) {
    optg();

    // -----------------------------------------------

    // MDI node commands

  } else if (strcmp(command, "@INIT_MD") == 0) {
    if (mode != DEFAULT) error->all(FLERR, "MDI: MDI engine is already performing a simulation");
    mode = MD;
    strncpy(node_driver, command, MDI_COMMAND_LENGTH);
    node_match = false;

  } else if (strcmp(command, "@INIT_OPTG") == 0) {
    if (mode != DEFAULT) error->all(FLERR, "MDI: MDI engine is already performing a simulation");
    mode = OPT;
    strncpy(node_driver, command, MDI_COMMAND_LENGTH);
    node_match = false;

  } else if (strcmp(command, "@") == 0) {
    strncpy(node_driver, "\0", MDI_COMMAND_LENGTH);
    node_match = false;

  } else if (strcmp(command, "@DEFAULT") == 0) {
    mode = DEFAULT;
    strncpy(node_driver, command, MDI_COMMAND_LENGTH);
    node_match = false;

  } else if (strcmp(command, "@COORDS") == 0) {
    strncpy(node_driver, command, MDI_COMMAND_LENGTH);
    node_match = false;

  } else if (strcmp(command, "@FORCES") == 0) {
    strncpy(node_driver, command, MDI_COMMAND_LENGTH);
    node_match = false;

  } else if (strcmp(command, "@ENDSTEP") == 0) {
    strncpy(node_driver, command, MDI_COMMAND_LENGTH);
    node_match = false;

    // exit command

  } else if (strcmp(command, "EXIT") == 0) {
    exit_command = true;

    // -------------------------------------------------------
    // custom LAMMPS commands
    // -------------------------------------------------------

  } else if (strcmp(command, "NBYTES") == 0) {
    nbytes_command();
  } else if (strcmp(command, "COMMAND") == 0) {
    single_command();
  } else if (strcmp(command, "COMMANDS") == 0) {
    many_commands();
  } else if (strcmp(command, "INFILE") == 0) {
    infile();
  } else if (strcmp(command, "<KE") == 0) {
    send_ke();

    // -------------------------------------------------------
    // unknown command
    // -------------------------------------------------------

  } else {
    error->all(FLERR, "MDI: Unknown command {} received from driver", command);
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
  MDI_Register_command("@DEFAULT", ">ELEMENTS");
  MDI_Register_command("@DEFAULT", ">NATOMS");
  MDI_Register_command("@DEFAULT", ">NSTEPS");
  MDI_Register_command("@DEFAULT", ">TOLERANCE");
  MDI_Register_command("@DEFAULT", ">TYPES");
  MDI_Register_command("@DEFAULT", ">VELOCITIES");
  MDI_Register_command("@DEFAULT", "MD");
  MDI_Register_command("@DEFAULT", "OPTG");
  MDI_Register_command("@DEFAULT", "@INIT_MD");
  MDI_Register_command("@DEFAULT", "@INIT_OPTG");
  MDI_Register_command("@DEFAULT", "EXIT");

  // default node, custom commands added by LAMMPS

  MDI_Register_command("@DEFAULT", "NBYTES");
  MDI_Register_command("@DEFAULT", "COMMAND");
  MDI_Register_command("@DEFAULT", "COMMANDS");
  MDI_Register_command("@DEFAULT", "INFILE");
  MDI_Register_command("@DEFAULT", "<KE");

  // node for setting up and running a dynamics simulation

  MDI_Register_node("@INIT_MD");
  MDI_Register_command("@INIT_MD", "<@");
  MDI_Register_command("@INIT_MD", "@");
  MDI_Register_command("@INIT_MD", "@DEFAULT");
  MDI_Register_command("@INIT_MD", "@COORDS");
  MDI_Register_command("@INIT_MD", "@FORCES");
  MDI_Register_command("@INIT_MD", "@ENDSTEP");
  MDI_Register_command("@INIT_MD", "EXIT");

  // node for setting up and running a minimization

  MDI_Register_node("@INIT_OPTG");
  MDI_Register_command("@INIT_OPTG", "<@");
  MDI_Register_command("@INIT_OPTG", "@");
  MDI_Register_command("@INIT_OPTG", "@DEFAULT");
  MDI_Register_command("@INIT_OPTG", "@COORDS");
  MDI_Register_command("@INIT_OPTG", "@FORCES");
  MDI_Register_command("@INIT_OPTG", "EXIT");

  // node at POST_INTEGRATE location in timestep
  // only used if fix MDI/ENGINE is instantiated

  MDI_Register_node("@COORDS");
  MDI_Register_command("@COORDS", "<@");
  MDI_Register_command("@COORDS", "<COORDS");
  MDI_Register_command("@COORDS", "<VELOCITIES");
  MDI_Register_command("@COORDS", ">COORDS");
  MDI_Register_command("@COORDS", ">VELOCITIES");
  MDI_Register_command("@COORDS", "@");
  MDI_Register_command("@COORDS", "@DEFAULT");
  MDI_Register_command("@COORDS", "@COORDS");
  MDI_Register_command("@COORDS", "@FORCES");
  MDI_Register_command("@COORDS", "@ENDSTEP");
  MDI_Register_command("@COORDS", "EXIT");

  // node at POST_FORCE location in timestep
  // only used if fix MDI/ENGINE is instantiated
  // two register callbacks allow LAMMPS to interact more easily
  //   with drivers which don't know LAMMPS control flow

  MDI_Register_node("@FORCES");
  MDI_Register_callback("@FORCES", ">FORCES");
  MDI_Register_callback("@FORCES", ">+FORCES");
  MDI_Register_command("@FORCES", "<@");
  MDI_Register_command("@FORCES", "<COORDS");
  MDI_Register_command("@FORCES", "<ENERGY");
  MDI_Register_command("@FORCES", "<FORCES");
  MDI_Register_command("@FORCES", "<KE");
  MDI_Register_command("@FORCES", "<PE");
  MDI_Register_command("@FORCES", "<STRESS");
  MDI_Register_command("@FORCES", "<VELOCITIES");
  MDI_Register_command("@FORCES", ">FORCES");
  MDI_Register_command("@FORCES", ">+FORCES");
  MDI_Register_command("@FORCES", ">VELOCITIES");
  MDI_Register_command("@FORCES", "@");
  MDI_Register_command("@FORCES", "@DEFAULT");
  MDI_Register_command("@FORCES", "@COORDS");
  MDI_Register_command("@FORCES", "@FORCES");
  MDI_Register_command("@FORCES", "@ENDSTEP");
  MDI_Register_command("@FORCES", "EXIT");

  // node at END_OF_STEP location in timestep
  // only used if fix MDI/ENGINE is instantiated

  MDI_Register_node("@ENDSTEP");
  MDI_Register_command("@ENDSTEP", "<@");
  MDI_Register_command("@ENDSTEP", "<ENERGY");
  MDI_Register_command("@ENDSTEP", "<FORCES");
  MDI_Register_command("@ENDSTEP", "<KE");
  MDI_Register_command("@ENDSTEP", "<PE");
  MDI_Register_command("@ENDSTEP", "<STRESS");
  MDI_Register_command("@ENDSTEP", "@");
  MDI_Register_command("@ENDSTEP", "@DEFAULT");
  MDI_Register_command("@ENDSTEP", "@COORDS");
  MDI_Register_command("@ENDSTEP", "@FORCES");
  MDI_Register_command("@ENDSTEP", "@ENDSTEP");
  MDI_Register_command("@ENDSTEP", "EXIT");
}

/* ----------------------------------------------------------------------
   run MD simulation one step at a time, controlled by driver
---------------------------------------------------------------------- */

void MDIEngine::mdi_md()
{
  // create or update system if requested prior to @INIT_MD

  int flag_create = flag_natoms | flag_types;
  if (flag_create)
    create_system();
  else {
    if (flag_cell || flag_cell_displ) adjust_box();
    if (flag_charges) adjust_charges();
    if (flag_coords) adjust_coords();
    if (flag_velocities) adjust_velocities();
  }

  // add an instance of fix MDI/ENGINE
  // delete the instance before this method returns

  modify->add_fix("MDI_ENGINE_INTERNAL all MDI/ENGINE");
  FixMDIEngine *mdi_fix =
      dynamic_cast<FixMDIEngine *>(modify->get_fix_by_id("MDI_ENGINE_INTERNAL"));
  mdi_fix->mdi_engine = this;

  // initialize LAMMPS and setup() the simulation

  update->whichflag = 1;
  timer->init_timeout();

  update->nsteps = 1;
  update->beginstep = update->firststep = update->ntimestep;
  update->endstep = update->laststep = update->ntimestep + update->nsteps;

  lmp->init();

  // engine is now at @INIT_MD node
  // any @ command from driver will start the simulation

  engine_node("@INIT_MD");
  if (strcmp(mdicmd, "EXIT") == 0) return;

  // run one step at a time forever
  // driver triggers exit with @ command other than @COORDS,@FORCES,@ENDSTEP

  update->integrate->setup(1);

  timer->init();
  timer->barrier_start();

  while (true) {
    update->nsteps += 1;
    update->laststep += 1;
    update->endstep = update->laststep;
    output->next = update->ntimestep + 1;

    update->integrate->run(1);

    if (strcmp(mdicmd, "@COORDS") != 0 && strcmp(mdicmd, "@FORCES") != 0 &&
        strcmp(mdicmd, "@ENDSTEP") != 0)
      break;
  }

  timer->barrier_stop();
  update->integrate->cleanup();

  modify->delete_fix("MDI_ENGINE_INTERNAL");

  // clear flags

  flag_natoms = flag_types = 0;
  flag_cell = flag_cell_displ = flag_charges = flag_coords = flag_velocities = 0;
}

/* ----------------------------------------------------------------------
   run MD simulation for >NSTEPS
---------------------------------------------------------------------- */

void MDIEngine::md()
{
  // create or update system if requested prior to MD command

  int flag_create = flag_natoms | flag_types;
  if (flag_create)
    create_system();
  else {
    if (flag_cell || flag_cell_displ) adjust_box();
    if (flag_charges) adjust_charges();
    if (flag_coords) adjust_coords();
    if (flag_velocities) adjust_velocities();
  }

  // initialize LAMMPS and setup() the simulation
  // run the simulation for nsteps
  // clean up

  update->whichflag = 1;
  timer->init_timeout();

  update->nsteps = nsteps;
  update->beginstep = update->firststep = update->ntimestep;
  update->endstep = update->laststep = update->ntimestep + update->nsteps;

  lmp->init();
  update->integrate->setup(1);

  timer->init();
  timer->barrier_start();

  update->integrate->run(nsteps);

  timer->barrier_stop();
  update->integrate->cleanup();

  // clear flags

  flag_natoms = flag_types = 0;
  flag_cell = flag_cell_displ = flag_charges = flag_coords = flag_velocities = 0;

  actionflag = 1;
}

/* ----------------------------------------------------------------------
   perform minimization one iteration at a time, controlled by driver
---------------------------------------------------------------------- */

void MDIEngine::mdi_optg()
{
  // create or update system if requested prior to @INIT_OPTG

  int flag_create = flag_natoms | flag_types;
  if (flag_create)
    create_system();
  else {
    if (flag_cell || flag_cell_displ) adjust_box();
    if (flag_charges) adjust_charges();
    if (flag_coords) adjust_coords();
    if (flag_velocities) adjust_velocities();
  }

  // add an instance of fix MDI/ENGINE
  // delete the instance before this method returns

  modify->add_fix("MDI_ENGINE_INTERNAL all MDI/ENGINE");
  FixMDIEngine *mdi_fix =
      dynamic_cast<FixMDIEngine *>(modify->get_fix_by_id("MDI_ENGINE_INTERNAL"));
  mdi_fix->mdi_engine = this;

  // initialize LAMMPS and setup() the simulation

  update->whichflag = 2;
  timer->init_timeout();

  update->nsteps = (niterate >= 0) ? niterate : max_eval;
  update->beginstep = update->firststep = update->ntimestep;
  update->endstep = update->laststep = update->firststep + update->nsteps;

  update->etol = 0.0;
  update->ftol = 0.0;
  update->max_eval = std::numeric_limits<int>::max();

  lmp->init();

  // engine is now at @INIT_OPTG node
  // any @ command from driver will start the minimization

  engine_node("@INIT_OPTG");
  if (strcmp(mdicmd, "EXIT") == 0) return;

  // run one iteration at a time forever
  // driver triggers exit with @ command other than @COORDS,@FORCES
  // two issues with running in this mode:
  //   @COORDS and @FORCES are not just triggered per min iteration
  //      but also for line search eng/force evals
  //   if driver triggers exit on step that is not multiple of thermo output
  //      then energy/virial not computed, and <PE, <STRESS will fail

  update->minimize->setup();

  timer->init();
  timer->barrier_start();

  while (true) {
    update->minimize->run(1);

    if (strcmp(mdicmd, "@COORDS") != 0 && strcmp(mdicmd, "@FORCES") != 0) break;
  }

  timer->barrier_stop();
  update->minimize->cleanup();

  modify->delete_fix("MDI_ENGINE_INTERNAL");

  // clear flags

  flag_natoms = flag_types = 0;
  flag_cell = flag_cell_displ = flag_charges = flag_coords = flag_velocities = 0;
}

/* ----------------------------------------------------------------------
   perform minimization to convergence using >TOLERANCE settings
---------------------------------------------------------------------- */

void MDIEngine::optg()
{
  // create or update system if requested prior to OPTG command

  int flag_create = flag_natoms | flag_types;
  if (flag_create)
    create_system();
  else {
    if (flag_cell || flag_cell_displ) adjust_box();
    if (flag_charges) adjust_charges();
    if (flag_coords) adjust_coords();
    if (flag_velocities) adjust_velocities();
  }

  // initialize LAMMPS and setup() the simulation
  // run the minimization using 4 >TOLERANCE parameters
  // clean up

  update->whichflag = 2;
  timer->init_timeout();

  update->nsteps = (niterate >= 0) ? niterate : max_eval;
  update->beginstep = update->firststep = update->ntimestep;
  update->endstep = update->laststep = update->firststep + update->nsteps;

  update->etol = etol;
  update->ftol = ftol;
  update->max_eval = max_eval;

  lmp->init();
  update->minimize->setup();

  timer->init();
  timer->barrier_start();

  update->minimize->run(niterate);

  timer->barrier_stop();
  update->minimize->cleanup();

  // clear flags

  flag_natoms = flag_types = 0;
  flag_cell = flag_cell_displ = flag_charges = flag_coords = flag_velocities = 0;

  actionflag = 1;
}

/* ----------------------------------------------------------------------
   evaluate() invoked by <ENERGY, <FORCES, <PE, <STRESS
   if flag_natoms or flag_types set, create a new system
   if any receive flags set, evaulate eng/forces/stress
---------------------------------------------------------------------- */

void MDIEngine::evaluate()
{
  int flag_create = flag_natoms | flag_types;
  int flag_other = flag_cell | flag_cell_displ | flag_charges | flag_coords | flag_velocities;

  // create or update system if requested

  if (flag_create)
    create_system();
  else if (flag_other) {
    if (flag_cell || flag_cell_displ) adjust_box();
    if (flag_charges) adjust_charges();
    if (flag_coords) adjust_coords();
    if (flag_velocities) adjust_velocities();
  }

  // if new system created or first-time eval:
  //   init LAMMPS and eval eng/force/virial via setup(1)
  // else:
  //   atom coords may or may not be updated incrementally
  //     incremental: timstepping an MD simulation
  //     non-incremental: e.g. processing snapshots from a dump file
  //   advance system by single step
  //   ensure potential energy and virial are tallied on new step
  //   check if reneighboing needed
  //   if no, just invoke setup_minimal(0)
  //   if yes, do an irregular->migrate_check() and migrate_atoms() if needed
  //     this can only be done if comm->style is not tiled
  //     also requires atoms be in box and lamda coords (for triclinic)
  //   finally invoke setup_minimal(1) to trigger exchange() & reneigh()
  // NOTE: what this logic still lacks for comm->style tiled,
  //       is when to invoke migrate_atoms() if necessary, not easy to detect

  if (flag_create || neighbor->ago < 0) {
    update->whichflag = 1;
    lmp->init();
    update->integrate->setup(1);
    update->whichflag = 0;

  } else if (flag_other) {
    update->ntimestep++;
    pe->addstep(update->ntimestep);
    press->addstep(update->ntimestep);

    int nflag = neighbor->decide();

    if (nflag == 0) {
      comm->forward_comm();
      update->integrate->setup_minimal(0);
      modify->clearstep_compute();
      output->thermo->compute(1);

    } else {
      if (comm->style == Comm::BRICK) {
        if (domain->triclinic) domain->x2lamda(atom->nlocal);
        domain->pbc();
        domain->reset_box();
        if (irregular->migrate_check()) irregular->migrate_atoms();
        if (domain->triclinic) domain->lamda2x(atom->nlocal);
      }

      update->integrate->setup_minimal(1);
      modify->clearstep_compute();
      output->thermo->compute(1);
    }

    modify->addstep_compute(update->ntimestep + 1);
  }

  // clear flags that trigger next eval

  flag_natoms = flag_types = 0;
  flag_cell = flag_cell_displ = flag_charges = flag_coords = flag_velocities = 0;
}

/* ----------------------------------------------------------------------
   create a new system
   >CELL, >NATOMS, >TYPES or >ELEMENTS, >COORDS commands are required
   >CELL_DISPL, >CHARGES, >VELOCITIES commands are optional
---------------------------------------------------------------------- */

void MDIEngine::create_system()
{
  // check requirements

  if (flag_cell == 0 || flag_natoms == 0 || flag_types == 0 || flag_coords == 0)
    error->all(FLERR,
               "MDI create_system requires >CELL, >NATOMS, "
               ">TYPES or >ELEMENTS, >COORDS MDI commands");

  // remove all existing atoms via delete_atoms command

  lmp->input->one("delete_atoms group all");

  // invoke lib->reset_box()

  double boxlo[3], boxhi[3];
  double xy, yz, xz;

  if (flag_cell_displ) {
    boxlo[0] = sys_cell_displ[0];
    boxlo[1] = sys_cell_displ[1];
    boxlo[2] = sys_cell_displ[2];
  } else {
    boxlo[0] = boxlo[1] = boxlo[2] = 0.0;
  }

  boxhi[0] = boxlo[0] + sys_cell[0];
  boxhi[1] = boxlo[1] + sys_cell[4];
  boxhi[2] = boxlo[2] + sys_cell[8];

  xy = sys_cell[3];
  yz = sys_cell[7];
  xz = sys_cell[6];

  lammps_reset_box(lmp, boxlo, boxhi, xy, yz, xz);

  // invoke lib->create_atoms()
  // create list of 1 to sys_natoms IDs
  // optionally set charges if specified by ">CHARGES"

  tagint *sys_ids;
  memory->create(sys_ids, sys_natoms, "mdi:sys_ids");
  for (int i = 0; i < sys_natoms; i++) sys_ids[i] = i + 1;

  if (flag_velocities)
    lammps_create_atoms(lmp, sys_natoms, sys_ids, sys_types, sys_coords, sys_velocities, nullptr,
                        1);
  else
    lammps_create_atoms(lmp, sys_natoms, sys_ids, sys_types, sys_coords, nullptr, nullptr, 1);

  if (flag_charges) lammps_scatter_atoms(lmp, (char *) "q", 1, 1, sys_charges);

  memory->destroy(sys_ids);

  // new system

  update->ntimestep = 0;
}

/* ----------------------------------------------------------------------
   adjust simulation box
---------------------------------------------------------------------- */

void MDIEngine::adjust_box()
{
  // convert atoms to lamda coords before changing box

  domain->x2lamda(atom->nlocal);

  // if >CELL command received,
  // convert celldata to new boxlo, boxhi, and tilt factors

  if (flag_cell) {
    domain->boxhi[0] = sys_cell[0] + domain->boxlo[0];
    domain->boxhi[1] = sys_cell[4] + domain->boxlo[1];
    domain->boxhi[2] = sys_cell[8] + domain->boxlo[2];
    domain->xy = sys_cell[3];
    domain->xz = sys_cell[6];
    domain->yz = sys_cell[7];
  }

  // if >CELL_DISPL command received,
  // convert cell_displ to new boxlo and boxhi

  if (flag_cell_displ) {
    double old_boxlo[3];
    old_boxlo[0] = domain->boxlo[0];
    old_boxlo[1] = domain->boxlo[1];
    old_boxlo[2] = domain->boxlo[2];

    domain->boxlo[0] = sys_cell_displ[0];
    domain->boxlo[1] = sys_cell_displ[1];
    domain->boxlo[2] = sys_cell_displ[2];

    domain->boxhi[0] += domain->boxlo[0] - old_boxlo[0];
    domain->boxhi[1] += domain->boxlo[1] - old_boxlo[1];
    domain->boxhi[2] += domain->boxlo[2] - old_boxlo[2];
  }

  // reset all Domain variables that depend on box size/shape
  // convert atoms from lamda coords back to new box coords

  domain->set_global_box();
  domain->set_local_box();
  domain->lamda2x(atom->nlocal);
}

/* ----------------------------------------------------------------------
   overwrite charges
---------------------------------------------------------------------- */

void MDIEngine::adjust_charges()
{
  double *q = atom->q;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  int ilocal;

  for (int i = 0; i < nlocal; i++) {
    ilocal = static_cast<int>(tag[i]) - 1;
    q[i] = sys_charges[ilocal];
  }
}

/* ----------------------------------------------------------------------
   overwrite coords
---------------------------------------------------------------------- */

void MDIEngine::adjust_coords()
{
  double **x = atom->x;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  int ilocal;

  for (int i = 0; i < nlocal; i++) {
    ilocal = static_cast<int>(tag[i]) - 1;
    x[i][0] = sys_coords[3 * ilocal + 0];
    x[i][1] = sys_coords[3 * ilocal + 1];
    x[i][2] = sys_coords[3 * ilocal + 2];
  }
}

/* ----------------------------------------------------------------------
   overwrite velocities
---------------------------------------------------------------------- */

void MDIEngine::adjust_velocities()
{
  double **v = atom->v;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  int ilocal;

  for (int i = 0; i < nlocal; i++) {
    ilocal = static_cast<int>(tag[i]) - 1;
    v[i][0] = sys_velocities[3 * ilocal + 0];
    v[i][1] = sys_velocities[3 * ilocal + 1];
    v[i][2] = sys_velocities[3 * ilocal + 2];
  }
}

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------/
// MDI ">" driver commands that send data
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   >CELL command
   reset simulation box edge vectors
   in conjunction with >CELL_DISPL this can change box arbitrarily
   can be done to create a new box
   can be done incrementally during MD or OPTG
---------------------------------------------------------------------- */

void MDIEngine::receive_cell()
{
  actionflag = 0;
  flag_cell = 1;
  int ierr = MDI_Recv(sys_cell, 9, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: >CELL data");
  MPI_Bcast(sys_cell, 9, MPI_DOUBLE, 0, world);

  for (int icell = 0; icell < 9; icell++) sys_cell[icell] *= mdi2lmp_length;

  // error check that edge vectors match LAMMPS triclinic requirement
  // 3,7,6 = xy, yz, xz tilt factors

  if (sys_cell[1] != 0.0 || sys_cell[2] != 0.0 || sys_cell[5] != 0.0)
    error->all(FLERR, "MDI: Received cell edges are not an upper triangular matrix");

  if (sys_cell[3] != 0.0 || sys_cell[7] != 0.0 || sys_cell[6] != 0.0)
    if (!domain->triclinic)
      error->all(FLERR,
                 "MDI: Received cell edges are for a triclinic box, "
                 "but LAMMPS is using an orthogonal box");
}

/* ----------------------------------------------------------------------
   >CELL_DISPL command
   reset simulation box lower left corner
   in conjunction with >CELL this can change box arbitrarily
   can be done to create a new box
   can be done incrementally during MD or OPTG
---------------------------------------------------------------------- */

void MDIEngine::receive_cell_displ()
{
  actionflag = 0;
  flag_cell_displ = 1;
  int ierr = MDI_Recv(sys_cell_displ, 3, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: >CELL_DISPLS data");
  MPI_Bcast(sys_cell_displ, 3, MPI_DOUBLE, 0, world);

  for (int icell = 0; icell < 3; icell++) sys_cell_displ[icell] *= mdi2lmp_length;
}

/* ----------------------------------------------------------------------
   >CHARGES command
---------------------------------------------------------------------- */

void MDIEngine::receive_charges()
{
  actionflag = 0;
  flag_charges = 1;
  int ierr = MDI_Recv(sys_charges, sys_natoms, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: >CHARGES data");
  MPI_Bcast(sys_charges, sys_natoms, MPI_DOUBLE, 0, world);
}

/* ----------------------------------------------------------------------
   >COORDS command
---------------------------------------------------------------------- */

void MDIEngine::receive_coords()
{
  actionflag = 0;
  flag_coords = 1;
  int n = 3 * sys_natoms;
  int ierr = MDI_Recv(sys_coords, n, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: >COORDS data");
  MPI_Bcast(sys_coords, n, MPI_DOUBLE, 0, world);
  for (int i = 0; i < n; i++) sys_coords[i] *= mdi2lmp_length;
}

/* ----------------------------------------------------------------------
   >ELEMENTS command
   receive elements for each atom = atomic numbers
   convert to LAMMPS atom types and store in sys_types
---------------------------------------------------------------------- */

void MDIEngine::receive_elements()
{
  actionflag = 0;
  flag_types = 1;
  int ierr = MDI_Recv(sys_types, sys_natoms, MDI_INT, mdicomm);
  if (ierr) error->all(FLERR, "MDI: >ELEMENTS data");
  MPI_Bcast(sys_types, sys_natoms, MPI_INT, 0, world);

  // convert from element atomic numbers to LAMMPS atom types
  // use mapping provided by mdi engine command

  int ntypes = atom->ntypes;
  int itype;

  for (int i = 0; i < sys_natoms; i++) {
    for (itype = 1; itype <= ntypes; itype++) {
      if (sys_types[i] == elements[itype]) {
        sys_types[i] = itype;
        break;
      }
    }
    if (itype > ntypes) error->all(FLERR, "MDI element not found in element list");
  }
}

/* ----------------------------------------------------------------------
   >NATOMS command
   natoms cannot exceed 32-bit int for use with MDI
---------------------------------------------------------------------- */

void MDIEngine::receive_natoms()
{
  actionflag = 0;
  flag_natoms = 1;
  int ierr = MDI_Recv(&sys_natoms, 1, MDI_INT, mdicomm);
  if (ierr) error->all(FLERR, "MDI: >NATOMS data");
  MPI_Bcast(&sys_natoms, 1, MPI_INT, 0, world);
  if (sys_natoms < 0) error->all(FLERR, "MDI received natoms < 0");
  reallocate();
}

/* ----------------------------------------------------------------------
   >NSTEPS command
   receive nsteps for timestepping
---------------------------------------------------------------------- */

void MDIEngine::receive_nsteps()
{
  int ierr = MDI_Recv(&nsteps, 1, MDI_INT, mdicomm);
  if (ierr) error->all(FLERR, "MDI: >NSTEPS data");
  MPI_Bcast(&nsteps, 1, MPI_INT, 0, world);
  if (nsteps < 0) error->all(FLERR, "MDI received nsteps < 0");
}

/* ----------------------------------------------------------------------
   >TOLERANCE command
   receive 4 minimization tolerance params
---------------------------------------------------------------------- */

void MDIEngine::receive_tolerance()
{
  double params[4];
  int ierr = MDI_Recv(params, 4, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: >TOLERANCE data");
  MPI_Bcast(params, 4, MPI_INT, 0, world);

  etol = params[0];
  ftol = params[1];
  niterate = static_cast<int>(params[2]);
  max_eval = static_cast<int>(params[3]);

  if (etol < 0.0 || ftol < 0.0 || niterate < 0 || max_eval < 0)
    error->all(FLERR, "MDI received invalid toleranace parameters");
}

/* ----------------------------------------------------------------------
   >TYPES command
---------------------------------------------------------------------- */

void MDIEngine::receive_types()
{
  actionflag = 0;
  flag_types = 1;
  int ierr = MDI_Recv(sys_types, sys_natoms, MDI_INT, mdicomm);
  if (ierr) error->all(FLERR, "MDI: >TYPES data");
  MPI_Bcast(sys_types, sys_natoms, MPI_INT, 0, world);
}

/* ----------------------------------------------------------------------
   >VELOCITIES command
---------------------------------------------------------------------- */

void MDIEngine::receive_velocities()
{
  actionflag = 0;
  flag_velocities = 1;
  int n = 3 * sys_natoms;
  int ierr = MDI_Recv(sys_velocities, n, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: >VELOCITIES data");
  MPI_Bcast(sys_velocities, n, MPI_DOUBLE, 0, world);
  for (int i = 0; i < n; i++) sys_coords[i] *= mdi2lmp_velocity;
}

/* ----------------------------------------------------------------------
   receive vector of 3 doubles for all atoms
   atoms are ordered by atomID, 1 to Natoms
   used by >FORCES command
---------------------------------------------------------------------- */

void MDIEngine::receive_double3(int which)
{
  int n = 3 * sys_natoms;
  int ierr = MDI_Recv(buf3, n, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: <double3 data");
  MPI_Bcast(buf3, n, MPI_DOUBLE, 0, world);

  // use atomID to index into ordered buf3

  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  int ilocal;

  if (which == FORCE) {
    double **f = atom->f;
    for (int i = 0; i < nlocal; i++) {
      ilocal = static_cast<int>(tag[i]) - 1;
      f[i][0] = buf3[3 * ilocal + 0] * mdi2lmp_force;
      f[i][1] = buf3[3 * ilocal + 1] * mdi2lmp_force;
      f[i][2] = buf3[3 * ilocal + 2] * mdi2lmp_force;
    }
  } else if (which == ADDFORCE) {
    double **f = atom->f;
    for (int i = 0; i < nlocal; i++) {
      ilocal = static_cast<int>(tag[i]) - 1;
      f[i][0] += buf3[3 * ilocal + 0] * mdi2lmp_force;
      f[i][1] += buf3[3 * ilocal + 1] * mdi2lmp_force;
      f[i][2] += buf3[3 * ilocal + 2] * mdi2lmp_force;
    }
  } else if (which == VELOCITY) {
    double **v = atom->v;
    for (int i = 0; i < nlocal; i++) {
      ilocal = static_cast<int>(tag[i]) - 1;
      v[i][0] = buf3[3 * ilocal + 0] * mdi2lmp_velocity;
      v[i][1] = buf3[3 * ilocal + 1] * mdi2lmp_velocity;
      v[i][2] = buf3[3 * ilocal + 2] * mdi2lmp_velocity;
    }
  }
}

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// MDI "<" driver commands that request data
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

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

  for (int icell = 0; icell < 9; icell++) celldata[icell] *= lmp2mdi_length;

  int ierr = MDI_Send(celldata, 9, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: <CELL data");
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

  for (int icell = 0; icell < 3; icell++) celldata[icell] *= lmp2mdi_length;

  int ierr = MDI_Send(celldata, 3, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: <CELL_DISPLS data");
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

  int ierr = MDI_Send(&total_energy, 1, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: <ENERGY data");
}

/* ----------------------------------------------------------------------
   <LABELS command
   convert numeric atom type to string for each atom
   atoms are ordered by atomID, 1 to Natoms
---------------------------------------------------------------------- */

void MDIEngine::send_labels()
{
  auto labels = new char[sys_natoms * MDI_LABEL_LENGTH];
  memset(labels, ' ', sys_natoms * MDI_LABEL_LENGTH);

  memset(ibuf1, 0, sys_natoms * sizeof(int));

  // use atomID to index into ordered ibuf1

  tagint *tag = atom->tag;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  int ilocal;

  for (int i = 0; i < nlocal; i++) {
    ilocal = static_cast<int>(tag[i]) - 1;
    ibuf1[ilocal] = type[i];
  }

  MPI_Reduce(ibuf1, ibuf1all, sys_natoms, MPI_INT, MPI_SUM, 0, world);

  if (comm->me == 0) {
    for (int iatom = 0; iatom < sys_natoms; iatom++) {
      std::string label = std::to_string(ibuf1all[iatom]);
      int label_len = std::min(int(label.length()), MDI_LABEL_LENGTH);
      strncpy(&labels[iatom * MDI_LABEL_LENGTH], label.c_str(), label_len);
    }
  }

  int ierr = MDI_Send(labels, sys_natoms * MDI_LABEL_LENGTH, MDI_CHAR, mdicomm);
  if (ierr) error->all(FLERR, "MDI: <LABELS data");

  delete[] labels;
}

/* ----------------------------------------------------------------------
   <NATOMS command
   natoms cannot exceed 32-bit int for use with MDI
---------------------------------------------------------------------- */

void MDIEngine::send_natoms()
{
  int ierr = MDI_Send(&sys_natoms, 1, MDI_INT, mdicomm);
  if (ierr != 0) error->all(FLERR, "MDI: <NATOMS data");
}

/* ----------------------------------------------------------------------
   <PE command
   send potential energy
---------------------------------------------------------------------- */

void MDIEngine::send_pe()
{
  double potential_energy = pe->compute_scalar();
  potential_energy *= lmp2mdi_energy;

  int ierr = MDI_Send(&potential_energy, 1, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: <PE data");
}

/* ----------------------------------------------------------------------
   <STRESS command
   send 9-component stress tensor (no kinetic energy term)
   should be intensive quantity (divided by volume in pressure compute)
   MDI stress tensor units are energy/volume,
     so conversion factor includes nktv2p to convert pressure back to virial
---------------------------------------------------------------------- */

void MDIEngine::send_stress()
{
  double vtensor_full[9];
  press->compute_vector();
  vtensor_full[0] = press->vector[0] * lmp2mdi_pressure;
  vtensor_full[4] = press->vector[1] * lmp2mdi_pressure;
  vtensor_full[8] = press->vector[2] * lmp2mdi_pressure;
  vtensor_full[1] = vtensor_full[3] = press->vector[3] * lmp2mdi_pressure;
  vtensor_full[2] = vtensor_full[6] = press->vector[4] * lmp2mdi_pressure;
  vtensor_full[5] = vtensor_full[7] = press->vector[5] * lmp2mdi_pressure;

  int ierr = MDI_Send(vtensor_full, 9, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: <STRESS data");
}

/* ----------------------------------------------------------------------
   send vector of 1 double for all atoms
   atoms are ordered by atomID, 1 to Natoms
   used by <CHARGE, <MASSES commands
---------------------------------------------------------------------- */

void MDIEngine::send_double1(int which)
{
  memset(buf1, 0, sys_natoms * sizeof(double));

  // use atomID to index into ordered buf1

  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  int ilocal;

  if (which == CHARGE) {
    double *q = atom->q;
    for (int i = 0; i < nlocal; i++) {
      ilocal = static_cast<int>(tag[i]) - 1;
      buf1[ilocal] = q[i];
    }
  } else if (which == MASS) {
    double *mass = atom->mass;
    double *rmass = atom->rmass;
    if (rmass) {
      for (int i = 0; i < nlocal; i++) {
        ilocal = static_cast<int>(tag[i]) - 1;
        buf1[ilocal] = rmass[i];
      }
    } else {
      int *type = atom->type;
      for (int i = 0; i < nlocal; i++) {
        ilocal = static_cast<int>(tag[i]) - 1;
        buf1[ilocal] = mass[type[i]];
      }
    }
  }

  MPI_Reduce(buf1, buf1all, sys_natoms, MPI_DOUBLE, MPI_SUM, 0, world);

  int ierr = MDI_Send(buf1all, sys_natoms, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: <double1 data");
}

/* ----------------------------------------------------------------------
   send vector of 1 int for all atoms
   atoms are ordered by atomID, 1 to Natoms
   use by <TYPES command
---------------------------------------------------------------------- */

void MDIEngine::send_int1(int which)
{
  memset(ibuf1, 0, sys_natoms * sizeof(int));

  // use atomID to index into ordered ibuf1

  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  int ilocal;

  if (which == TYPE) {
    int *type = atom->type;
    for (int i = 0; i < nlocal; i++) {
      ilocal = static_cast<int>(tag[i]) - 1;
      ibuf1[ilocal] = type[i];
    }
  }

  MPI_Reduce(ibuf1, ibuf1all, sys_natoms, MPI_INT, MPI_SUM, 0, world);

  int ierr = MDI_Send(ibuf1all, sys_natoms, MDI_INT, mdicomm);
  if (ierr) error->all(FLERR, "MDI: <int1 data");
}

/* ----------------------------------------------------------------------
   <COORDS, <FORCES, <VELOCITIES commands
   send vector of 3 doubles for all atoms
   atoms are ordered by atomID, 1 to Natoms
---------------------------------------------------------------------- */

void MDIEngine::send_double3(int which)
{
  memset(buf3, 0, 3 * sys_natoms * sizeof(double));

  // use atomID to index into ordered buf3

  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  int ilocal;

  if (which == COORD) {
    double **x = atom->x;
    for (int i = 0; i < nlocal; i++) {
      ilocal = static_cast<int>(tag[i]) - 1;
      buf3[3 * ilocal + 0] = x[i][0] * lmp2mdi_length;
      buf3[3 * ilocal + 1] = x[i][1] * lmp2mdi_length;
      buf3[3 * ilocal + 2] = x[i][2] * lmp2mdi_length;
    }
  } else if (which == FORCE) {
    double **f = atom->f;
    for (int i = 0; i < nlocal; i++) {
      ilocal = static_cast<int>(tag[i]) - 1;
      buf3[3 * ilocal + 0] = f[i][0] * lmp2mdi_force;
      buf3[3 * ilocal + 1] = f[i][1] * lmp2mdi_force;
      buf3[3 * ilocal + 2] = f[i][2] * lmp2mdi_force;
    }
  } else if (which == VELOCITY) {
    double **v = atom->v;
    for (int i = 0; i < nlocal; i++) {
      ilocal = static_cast<int>(tag[i]) - 1;
      buf3[3 * ilocal + 0] = v[i][0] * lmp2mdi_velocity;
      buf3[3 * ilocal + 1] = v[i][1] * lmp2mdi_velocity;
      buf3[3 * ilocal + 2] = v[i][2] * lmp2mdi_velocity;
    }
  }

  MPI_Reduce(buf3, buf3all, 3 * sys_natoms, MPI_DOUBLE, MPI_SUM, 0, world);

  int ierr = MDI_Send(buf3all, 3 * sys_natoms, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: <double3 data");
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
  int ierr = MDI_Recv(&nbytes, 1, MDI_INT, mdicomm);
  if (ierr) error->all(FLERR, "MDI: NBYTES data");
  MPI_Bcast(&nbytes, 1, MPI_INT, 0, world);
}

/* ----------------------------------------------------------------------
   COMMAND command
   store received value as string of length nbytes
   invoke as a LAMMPS command
---------------------------------------------------------------------- */

void MDIEngine::single_command()
{
  if (nbytes < 0) error->all(FLERR, "MDI: COMMAND nbytes has not been set");

  auto cmd = new char[nbytes + 1];
  int ierr = MDI_Recv(cmd, nbytes + 1, MDI_CHAR, mdicomm);
  if (ierr) error->all(FLERR, "MDI: COMMAND data");
  MPI_Bcast(cmd, nbytes + 1, MPI_CHAR, 0, world);
  cmd[nbytes] = '\0';

  lammps_command(lmp, cmd);

  delete[] cmd;
}

/* ----------------------------------------------------------------------
   COMMANDS command
   store received value as multi-line string of length nbytes
   invoke as multiple LAMMPS commands
---------------------------------------------------------------------- */

void MDIEngine::many_commands()
{
  if (nbytes < 0) error->all(FLERR, "MDI: COMMANDS nbytes has not been set");

  auto cmds = new char[nbytes + 1];
  int ierr = MDI_Recv(cmds, nbytes + 1, MDI_CHAR, mdicomm);
  if (ierr) error->all(FLERR, "MDI: COMMANDS data");
  MPI_Bcast(cmds, nbytes + 1, MPI_CHAR, 0, world);
  cmds[nbytes] = '\0';

  lammps_commands_string(lmp, cmds);

  delete[] cmds;
}

/* ----------------------------------------------------------------------
   INFILE command
   store received value as infile of length length_param
   invoke as a LAMMPS input script
---------------------------------------------------------------------- */

void MDIEngine::infile()
{
  if (nbytes < 0) error->all(FLERR, "MDI: INFILE nbytes has not been set");

  auto infile = new char[nbytes + 1];
  int ierr = MDI_Recv(infile, nbytes + 1, MDI_CHAR, mdicomm);
  if (ierr) error->all(FLERR, "MDI: INFILE data for {}", infile);
  MPI_Bcast(infile, nbytes + 1, MPI_CHAR, 0, world);
  infile[nbytes] = '\0';

  lammps_file(lmp, infile);

  delete[] infile;
}

/* ----------------------------------------------------------------------
   <KE command
   send kinetic energy
---------------------------------------------------------------------- */

void MDIEngine::send_ke()
{
  double kinetic_energy = ke->compute_scalar();
  kinetic_energy *= lmp2mdi_energy;

  int ierr = MDI_Send(&kinetic_energy, 1, MDI_DOUBLE, mdicomm);
  if (ierr) error->all(FLERR, "MDI: <KE data");
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
  if (sys_natoms <= maxatom) return;

  bigint nsize = (bigint) sys_natoms * 3;
  if (nsize > MAXSMALLINT) error->all(FLERR, "Natoms too large to use with mdi engine");

  maxatom = sys_natoms;

  deallocate();
  allocate();
}

void MDIEngine::deallocate()
{
  memory->destroy(sys_types);
  memory->destroy(sys_coords);
  memory->destroy(sys_velocities);
  memory->destroy(sys_charges);

  memory->destroy(ibuf1);
  memory->destroy(buf1);
  memory->destroy(buf3);

  memory->destroy(ibuf1all);
  memory->destroy(buf1all);
  memory->destroy(buf3all);
}

void MDIEngine::allocate()
{
  memory->create(sys_types, maxatom, "mdi:sys_types");
  memory->create(sys_coords, 3 * maxatom, "mdi:sys_coords");
  memory->create(sys_velocities, 3 * maxatom, "mdi:sys_velocities");
  memory->create(sys_charges, maxatom, "mdi:sys_charges");

  memory->create(ibuf1, maxatom, "mdi:ibuf1");
  memory->create(buf1, maxatom, "mdi:buf1");
  memory->create(buf3, 3 * maxatom, "mdi:buf3");

  memory->create(ibuf1all, maxatom, "mdi:ibuf1all");
  memory->create(buf1all, maxatom, "mdi:buf1all");
  memory->create(buf3all, 3 * maxatom, "mdi:buf3all");
}

/* ----------------------------------------------------------------------
   MDI to/from LAMMPS conversion factors
------------------------------------------------------------------------- */

void MDIEngine::unit_conversions()
{
  double angstrom_to_bohr, kelvin_to_hartree, ev_to_hartree, second_to_aut;

  MDI_Conversion_factor("angstrom", "bohr", &angstrom_to_bohr);
  MDI_Conversion_factor("kelvin_energy", "hartree", &kelvin_to_hartree);
  MDI_Conversion_factor("electron_volt", "hartree", &ev_to_hartree);
  MDI_Conversion_Factor("second", "atomic_unit_of_time", &second_to_aut);

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
  // MDI energy/volume = Hartree/Bohr^3,
  //   so need to remove LAMMPS nktv2p from pressure

  mdi2lmp_pressure = 1.0;
  lmp2mdi_pressure = 1.0;

  if (lmpunits == REAL) {
    lmp2mdi_pressure = (kelvin_to_hartree / force->boltz) /
        (angstrom_to_bohr * angstrom_to_bohr * angstrom_to_bohr) / force->nktv2p;
    mdi2lmp_pressure = 1.0 / lmp2mdi_pressure;
  } else if (lmpunits == METAL) {
    lmp2mdi_pressure =
        ev_to_hartree / (angstrom_to_bohr * angstrom_to_bohr * angstrom_to_bohr) / force->nktv2p;
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
