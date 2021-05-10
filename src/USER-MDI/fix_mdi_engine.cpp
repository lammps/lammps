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

#include "fix_mdi_engine.h"

#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "irregular.h"
#include "library.h"
#include "memory.h"
#include "modify.h"
#include "update.h"
#include "verlet.h"

#include <limits>
#include <string.h>

enum{NONE,REAL,METAL};    // LAMMPS units which MDI supports

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixMDIEngine::FixMDIEngine(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  id_pe(NULL), pe(NULL),
  id_ke(NULL), ke(NULL)
{
  if (narg != 3) error->all(FLERR,"Illegal fix mdi command");

  // NOTE: real & metal the 2 atomic-scale units LAMMPS has
  //       I suggest LAMMPS for MDI support both
  // real: coords = Ang, eng = Kcal/mole, force = Kcal/mole/Ang
  // metal: coords = Ang, eng = eV, force = eV/Ang

  lmpunits = NONE;
  if (strcmp(update->unit_style,"real") == 0) lmpunits = REAL;
  if (strcmp(update->unit_style,"metal") == 0) lmpunits = METAL;
  if (lmpunits == NONE) error->all(FLERR,"MDI requires real or metal units");

  // MDI setup

  most_recent_init = 0;
  exit_flag = false;
  local_exit_flag = false;
  target_command = new char[MDI_COMMAND_LENGTH+1];
  command = new char[MDI_COMMAND_LENGTH+1];
  current_node = new char[MDI_COMMAND_LENGTH];
  target_node = new char[MDI_COMMAND_LENGTH];
  strncpy(target_node, "\0", MDI_COMMAND_LENGTH);
  strncpy(current_node, "@DEFAULT", MDI_COMMAND_LENGTH);

  // accept a communicator to the driver
  // master = 1 for proc 0, otherwise 0

  master = (comm->me==0) ? 1 : 0;

  if (master) {
    MDI_Accept_Communicator(&driver_socket);
    if (driver_socket <= 0) error->all(FLERR,"Unable to connect to driver");
  } else driver_socket = 0;

  // create computes for KE and PE

  id_pe = utils::strdup(std::string(id) + "_pe");
  modify->add_compute(fmt::format("{} all pe",id_pe));

  id_ke = utils::strdup(std::string(id) + "_ke");
  modify->add_compute(fmt::format("{} all ke",id_ke));

  // irregular class and data structs used by MDI

  irregular = new Irregular(lmp);
  add_force = NULL;
}

/* ---------------------------------------------------------------------- */

FixMDIEngine::~FixMDIEngine()
{
  delete [] target_command;
  delete [] command;
  delete [] current_node;
  delete [] target_node;

  modify->delete_compute(id_pe);
  modify->delete_compute(id_ke);
  delete irregular;
  memory->destroy(add_force);
}

/* ---------------------------------------------------------------------- */

int FixMDIEngine::setmask()
{
  int mask = 0;

  mask |= POST_INTEGRATE;
  mask |= PRE_REVERSE;
  mask |= POST_FORCE;
  mask |= MIN_PRE_FORCE;   // NOTE: whack this?
  mask |= MIN_POST_FORCE;

  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::exchange_forces()
{
  double **f = atom->f;
  const int * const mask  = atom->mask;
  const int nlocal = atom->nlocal;

  // add forces from the driver

  for (int i = 0; i < nlocal; ++i) {
    if (mask[i] & groupbit) {
      f[i][0] += add_force[3*(atom->tag[i]-1)+0];
      f[i][1] += add_force[3*(atom->tag[i]-1)+1];
      f[i][2] += add_force[3*(atom->tag[i]-1)+2];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::init()
{
  // confirm that two required computes are still available

  int icompute_pe = modify->find_compute(id_pe);
  if (icompute_pe < 0)
    error->all(FLERR,"Potential energy ID for fix mdi/engine does not exist");
  int icompute_ke = modify->find_compute(id_ke);
  if (icompute_pe < 0)
    error->all(FLERR,"Kinetic energy ID for fix mdi/engine does not exist");

  pe = modify->compute[icompute_pe];
  ke = modify->compute[icompute_ke];

  // one-time allocation of add_force array
  // NOTE: moved this here b/c natoms may not be defined when Fix is constructed
  // NOTE: check that 3*natoms does not overflow a 32-bit int

  if (!add_force) {
    memory->create(add_force,3*atom->natoms,"mdi/engine:add_force");
    for (int i = 0; i < 3*atom->natoms; i++) add_force[i] = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::setup(int vflag)
{
  // NOTE: this seems an odd place to compute these
  //       I think it would be better to compute these on-demand
  //         in response to a driver request?

  potential_energy = pe->compute_scalar();
  kinetic_energy = ke->compute_scalar();

  // trigger potential energy computation on next timestep
  // NOTE: there is no triggering needed for KE

  pe->addstep(update->ntimestep+1);

  if (most_recent_init == 1) engine_mode("@PRE-FORCES");
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::min_setup(int vflag)
{
  // NOTE: this seems an odd place to compute these

  potential_energy = pe->compute_scalar();
  kinetic_energy = ke->compute_scalar();

  engine_mode("@FORCES");

  // trigger potential energy computation on next timestep

  pe->addstep(update->ntimestep+1);
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::post_integrate()
{
  engine_mode("@COORDS");
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::pre_reverse(int eflag, int vflag)
{
  // NOTE: this seems an odd place to compute these

  potential_energy = pe->compute_scalar();
  kinetic_energy = ke->compute_scalar();

  engine_mode("@PRE-FORCES");

  // trigger potential energy computation on next timestep

  pe->addstep(update->ntimestep+1);
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::min_pre_force(int vflag)
{
  engine_mode("@COORDS");
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::min_post_force(int vflag)
{
  // NOTE: this seems an odd place to compute these

  potential_energy = pe->compute_scalar();
  kinetic_energy = ke->compute_scalar();

  engine_mode("@FORCES");

  // trigger potential energy computation on next timestep

  pe->addstep(update->ntimestep+1);
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::post_force(int vflag)
{
  if (most_recent_init == 1) engine_mode("@FORCES");
  else if (most_recent_init == 2) engine_mode("@FORCES");

  // NOTE: should this also be done in this method?
  // trigger potential energy computation on next timestep
  // NOTE: in general, forcing the pair styles to compute PE every step
  //       is inefficient, would be better to think of another way to do this,
  //       e.g. at very end of a step, 
  //       possibly by mdi_engine when knows it is needed

  pe->addstep(update->ntimestep+1);
}

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// rest of file processes and responds to MDI driver commands
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   process a single command from driver
---------------------------------------------------------------------- */

int FixMDIEngine::execute_command(const char *command, MDI_Comm mdicomm)
{
  /*
  if (screen)
    fprintf(screen,"MDI command: %s\n",command);
  if (logfile)
    fprintf(logfile,"MDI command: %s:\n",command);
  */

  // confirm this command is supported at this node

  int command_exists = 0;
  ierr = MDI_Check_Command_Exists(current_node, command, 
                                  MDI_COMM_NULL, &command_exists);
  if (ierr != 0)
    error->all(FLERR,"MDI: Unable to check whether current command is supported");
  if (command_exists != 1)
    error->all(FLERR,"MDI: Received an unsupported at current node");

  // respond to any driver command

  // send calculation status to the driver
  // NOTE: why does STATUS have extra spaces?

  if (strcmp(command,"STATUS      ") == 0 ) {
    if (master) {
      ierr = MDI_Send_Command("READY", mdicomm);
      if (ierr != 0)
	error->all(FLERR,"MDI: Unable to return status to driver");
    }

  } else if (strcmp(command,">NATOMS") == 0 ) {
    if (master) {
      ierr = MDI_Recv((char*) &atom->natoms, 1, MDI_INT, mdicomm);
      if (ierr != 0)
	error->all(FLERR,"MDI: Unable to receive number of atoms from driver");
      }
    MPI_Bcast(&atom->natoms,1,MPI_INT,0,world);

  } else if (strcmp(command,"<NATOMS") == 0 ) {
    if (master) {
      int64_t mdi_natoms = atom->natoms;
      ierr = MDI_Send((char*) &mdi_natoms, 1, MDI_INT64_T, mdicomm);
      if (ierr != 0)
	error->all(FLERR,"MDI: Unable to send number of atoms to driver");
    }

  } else if (strcmp(command,"<NTYPES") == 0 ) {
    if (master) {
      ierr = MDI_Send((char*) &atom->ntypes, 1, MDI_INT, mdicomm);
      if (ierr != 0)
	error->all(FLERR,"MDI: Unable to send number of atom types to driver");
    }
  
  } else if (strcmp(command,"<TYPES") == 0 ) {
    send_types(error);

  } else if (strcmp(command,"<LABELS") == 0 ) {
    send_labels(error);

  } else if (strcmp(command,"<MASSES") == 0 ) {
    send_masses(error);

  // NOTE: ParSplice is going to definitely also need a >CELL command

  } else if (strcmp(command,"<CELL") == 0 ) {
    send_cell(error);

  } else if (strcmp(command,">COORDS") == 0 ) {
    receive_coordinates(error);

  } else if (strcmp(command,"<COORDS") == 0 ) {
    send_coordinates(error);

  } else if (strcmp(command,"<CHARGES") == 0 ) {
    send_charges(error);

  } else if (strcmp(command,"<ENERGY") == 0 ) {
    send_energy(error);

  } else if (strcmp(command,"<FORCES") == 0 ) {
    send_forces(error);

  } else if (strcmp(command,">FORCES") == 0 ) {
    receive_forces(error);

  // receive additional forces from the driver
  // these can be added prior to SHAKE or other post-processing
  // NOTE: maybe this is now not necessary?

  } else if (strcmp(command,"+FORCES") == 0 ) {
    add_forces(error);

  // initialize new MD simulation or minimization
  // return control to return to mdi_engine

  } else if (strcmp(command,"@INIT_MD") == 0 ) {
    if (most_recent_init != 0)
      error->all(FLERR,"MDI: MDI is already performing a simulation");
    most_recent_init = 1;
    local_exit_flag = true;

  // initialize new energy minimization
  // return control to return to mdi_engine

  } else if (strcmp(command,"@INIT_OPTG") == 0 ) {
    if ( most_recent_init != 0 )
      error->all(FLERR,"MDI: MDI is already performing a simulation");
    most_recent_init = 2;
    local_exit_flag = true;

  } else if (strcmp(command,"@") == 0 ) {
    strncpy(target_node, "\0", MDI_COMMAND_LENGTH);
    local_exit_flag = true;

  } else if (strcmp(command,"<@") == 0 ) {
    if (master) {
      ierr = MDI_Send(current_node, MDI_NAME_LENGTH, MDI_CHAR, mdicomm);
      if (ierr != 0) error->all(FLERR,"MDI: Unable to send node to driver");
    }

  } else if (strcmp(command,"<KE") == 0 ) {
    send_ke(error);

  } else if (strcmp(command,"<PE") == 0 ) {
    send_pe(error);

  // NOTE: should there be a test here for @DEFAULT ??

  } else if (strcmp(command,"@COORDS") == 0 ) {
    strncpy(target_node, "@COORDS", MDI_COMMAND_LENGTH);
    local_exit_flag = true;

  } else if (strcmp(command,"@PRE-FORCES") == 0 ) {
    strncpy(target_node, "@PRE-FORCES", MDI_COMMAND_LENGTH);
    local_exit_flag = true;

  } else if (strcmp(command,"@FORCES") == 0 ) {
    strncpy(target_node, "@FORCES", MDI_COMMAND_LENGTH);
    local_exit_flag = true;

  } else if (strcmp(command,"EXIT_SIM") == 0 ) {
    most_recent_init = 0;

    local_exit_flag = true;

    // are we in the middle of a geometry optimization?
    if ( most_recent_init == 2 ) {
      // ensure that the energy and force tolerances are met
      update->etol = std::numeric_limits<double>::max();
      update->ftol = std::numeric_limits<double>::max();

      // set the maximum number of force evaluations to 0
      update->max_eval = 0;
    }

  } else if (strcmp(command,"EXIT") == 0 ) {
    // exit the driver code
    exit_flag = true;

    // are we in the middle of a geometry optimization?
    if ( most_recent_init == 2 ) {
      // ensure that the energy and force tolerances are met
      update->etol = std::numeric_limits<double>::max();
      update->ftol = std::numeric_limits<double>::max();

      // set the maximum number of force evaluations to 0
      update->max_eval = 0;
    }

  } else {
    error->all(FLERR,"MDI: Unknown command from driver");
  }

  return 0;
}

/* ---------------------------------------------------------------------- */

char *FixMDIEngine::engine_mode(const char *node)
{
  /*
  if (screen)
    fprintf(screen,"MDI ENGINE MODE: %i\n",node);
  if (logfile)
    fprintf(logfile,"MDI ENGINE MODE: %i\n",node);
  */

  // do not process commands if engine and driver are not at same node
  // target_node = node that driver has set via a @ command
  // current_node = node that engine (LAMMPS) has set

  strncpy(current_node,node,MDI_COMMAND_LENGTH);
  if (strcmp(target_node,"\0") != 0 && strcmp(target_node,current_node) != 0)
    local_exit_flag = true;

  // register the execute_command function with MDI
  // NOTE: does this need to be done multiple times ??

  MDI_Set_execute_command_func(lammps_execute_mdi_command, this);

  // respond to commands from the driver

  while (not exit_flag and not local_exit_flag) {

    // read the next command from the driver
    // NOTE: all procs call this, but only proc 0 receives command?

    ierr = MDI_Recv_Command(command, driver_socket);
    if (ierr != 0)
      error->all(FLERR,"MDI: Unable to receive command from driver");

    // broadcast command to the other MPI tasks

    MPI_Bcast(command,MDI_COMMAND_LENGTH,MPI_CHAR,0,world);

    // execute the command

    this->execute_command(command, driver_socket);

    // check if the target node is something other than the current node

    if (strcmp(target_node,"\0") != 0 && strcmp(target_node,current_node) != 0 )
      local_exit_flag = true;
  }

  // local exit occured so turn off local exit flag

  local_exit_flag = false;

  return command;
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::receive_coordinates(Error* error)
{
  double posconv;

  // NOTE: logic like this everywhere else

  if (lmpunits == REAL) {
    double angstrom_to_bohr;
    MDI_Conversion_Factor("angstrom", "bohr", &angstrom_to_bohr);
    posconv=force->angstrom/angstrom_to_bohr;
  } else if (lmpunits == METAL) {
    // ??
  }

  // create buffer to hold all coords

  double *buffer;
  buffer = new double[3*atom->natoms];

  if (master) {
    ierr = MDI_Recv((char*) buffer, 3*atom->natoms, MDI_DOUBLE, driver_socket);
    if (ierr != 0)
      error->all(FLERR,"MDI: Unable to receive coordinates from driver");
  }
  MPI_Bcast(buffer,3*atom->natoms,MPI_DOUBLE,0,world);

  // pick local atoms from the buffer

  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    x[i][0]=buffer[3*(atom->tag[i]-1)+0]*posconv;
    x[i][1]=buffer[3*(atom->tag[i]-1)+1]*posconv;
    x[i][2]=buffer[3*(atom->tag[i]-1)+2]*posconv;
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

  delete [] buffer;
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::send_coordinates(Error* error)
{
  double posconv;
  double angstrom_to_bohr;
  MDI_Conversion_Factor("angstrom", "bohr", &angstrom_to_bohr);
  posconv=force->angstrom/angstrom_to_bohr;

  double *coords;
  double *coords_reduced;

  // NOTE: I suggest
  // double *coords;
  // memory->create(coords,3*atom->natoms,"mdi:coords");

  coords = new double[3*atom->natoms];
  coords_reduced = new double[3*atom->natoms];

  // zero coords

  for (int icoord = 0; icoord < 3*atom->natoms; icoord++) coords[icoord] = 0.0;

  // copy local atoms into buffer at correct locations

  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    coords[3*(atom->tag[i]-1)+0] = x[i][0]/posconv;
    coords[3*(atom->tag[i]-1)+1] = x[i][1]/posconv;
    coords[3*(atom->tag[i]-1)+2] = x[i][2]/posconv;
  }

  MPI_Reduce(coords,coords_reduced,3*atom->natoms,MPI_DOUBLE,MPI_SUM,0,world);

  if (master) {
    ierr = MDI_Send((char*) coords_reduced,3*atom->natoms,MDI_DOUBLE,
                    driver_socket);
    if (ierr != 0)
      error->all(FLERR,"MDI: Unable to send coordinates to driver");
  }

  delete [] coords;
  delete [] coords_reduced;
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::send_charges(Error* error)
{
  double *charges;
  double *charges_reduced;

  charges = new double[atom->natoms];
  charges_reduced = new double[atom->natoms];

  // zero the charges array

  for (int icharge = 0; icharge < atom->natoms; icharge++) charges[icharge] = 0.0;

  // pick local atoms from the buffer

  double *charge = atom->q;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    charges[atom->tag[i]-1] = charge[i];
  }

  MPI_Reduce(charges, charges_reduced, atom->natoms, MPI_DOUBLE, MPI_SUM, 0, world);

  if (master) { 
    ierr = MDI_Send((char*) charges_reduced, atom->natoms, MDI_DOUBLE, driver_socket);
    if (ierr != 0)
      error->all(FLERR,"MDI: Unable to send charges to driver");
  }

  delete [] charges;
  delete [] charges_reduced;
}


/* ---------------------------------------------------------------------- */

void FixMDIEngine::send_energy(Error* error)
{
  double kelvin_to_hartree;
  MDI_Conversion_Factor("kelvin_energy", "hartree", &kelvin_to_hartree);

  // NOTE: I suggest you invoke the 2 computes here

  double pe = potential_energy;
  double ke = kinetic_energy;
  double total_energy;
  double *send_energy = &total_energy;

  // convert the energy to atomic units
  pe *= kelvin_to_hartree/force->boltz;
  ke *= kelvin_to_hartree/force->boltz;
  total_energy = pe + ke;

  if (master) {
    ierr = MDI_Send((char*) send_energy, 1, MDI_DOUBLE, driver_socket);
    if (ierr != 0)
      error->all(FLERR,"MDI: Unable to send potential energy to driver");
  }
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::send_pe(Error* error)
{
  double kelvin_to_hartree;
  MDI_Conversion_Factor("kelvin_energy", "hartree", &kelvin_to_hartree);

  // NOTE: I suggest you invoke the PE compute here

  double pe = potential_energy;
  double *send_energy = &pe;

  // convert the energy to atomic units
  pe *= kelvin_to_hartree/force->boltz;

  if (master) {
    ierr = MDI_Send((char*) send_energy, 1, MDI_DOUBLE, driver_socket);
    if (ierr != 0)
      error->all(FLERR,"MDI: Unable to send potential energy to driver");
  }
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::send_ke(Error* error)
{
  double kelvin_to_hartree;
  MDI_Conversion_Factor("kelvin_energy", "hartree", &kelvin_to_hartree);

  // NOTE: I suggest you invoke the KE compute here

  double ke = kinetic_energy;
  double *send_energy = &ke;

  // convert the energy to atomic units
  ke *= kelvin_to_hartree/force->boltz;

  if (master) {
    ierr = MDI_Send((char*) send_energy, 1, MDI_DOUBLE, driver_socket);
    if (ierr != 0)
      error->all(FLERR,"MDI: Unable to send potential energy to driver");
  }
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::send_types(Error* error)
{
  int * const type = atom->type;

  // NOTE: why is this not supported?
  //       maybe MDI labels = LAMMPS types?

  if (master) { 
    ierr = MDI_Send((char*) type, atom->natoms, MDI_INT, driver_socket);
    if (ierr != 0)
      error->all(FLERR,"MDI: Unable to send atom types to driver");
  }
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::send_labels(Error* error)
{
  char *labels = new char[atom->natoms * MDI_LABEL_LENGTH];
  memset(labels, ' ', atom->natoms * MDI_LABEL_LENGTH);

  for (int iatom=0; iatom < atom->natoms; iatom++) {
    std::string label = std::to_string( atom->type[iatom] );
    int label_len = std::min( int(label.length()), MDI_LABEL_LENGTH );
    strncpy(&labels[iatom * MDI_LABEL_LENGTH], label.c_str(), label_len);
  }

  if (master) { 
    ierr = MDI_Send( labels, atom->natoms * MDI_LABEL_LENGTH, MDI_CHAR, driver_socket);
    if (ierr != 0)
      error->all(FLERR,"MDI: Unable to send atom types to driver");
  }

  delete [] labels;
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::send_masses(Error* error)
{
  double * const rmass = atom->rmass;
  double * const mass = atom->mass;
  int * const type = atom->type;
  int nlocal = atom->nlocal;

  double *mass_by_atom = new double[atom->natoms];
  double *mass_by_atom_reduced = new double[atom->natoms];
  for (int iatom=0; iatom < atom->natoms; iatom++) {
    mass_by_atom[iatom] = 0.0;
  }

  // determine the atomic masses

  if (rmass) {
    for (int iatom=0; iatom < nlocal; iatom++) {
      mass_by_atom[ atom->tag[iatom] - 1 ] = rmass[iatom];
    }
  }
  else {
    for (int iatom=0; iatom < nlocal; iatom++) {
      mass_by_atom[ atom->tag[iatom] - 1 ] = mass[ type[iatom] ];
    }
  }

  MPI_Reduce(mass_by_atom, mass_by_atom_reduced, atom->natoms, MPI_DOUBLE, MPI_SUM, 0, world);

  // send the atomic masses to the driver

  if (master) {
    ierr = MDI_Send((char*) mass_by_atom_reduced, atom->natoms, MDI_DOUBLE, driver_socket);
    if (ierr != 0)
      error->all(FLERR,"MDI: Unable to send atom masses to driver");
  }
  
  delete [] mass_by_atom;
  delete [] mass_by_atom_reduced;
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::send_forces(Error* error)
{
  double angstrom_to_bohr;
  MDI_Conversion_Factor("angstrom", "bohr", &angstrom_to_bohr);
  double kelvin_to_hartree;
  MDI_Conversion_Factor("kelvin_energy", "hartree", &kelvin_to_hartree);

  double potconv, posconv, forceconv;
  potconv=kelvin_to_hartree/force->boltz;
  posconv=force->angstrom/angstrom_to_bohr;
  forceconv=potconv*posconv;

  double *forces;
  double *forces_reduced;
  double *x_buf;

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  forces = new double[3*atom->natoms];
  forces_reduced = new double[3*atom->natoms];
  x_buf = new double[3*nlocal];

  // zero the forces array

  for (int iforce = 0; iforce < 3*atom->natoms; iforce++) forces[iforce] = 0.0;

  // if not at a node, calculate the forces

  if ( strcmp(current_node, "@DEFAULT") == 0 ) {
    // certain fixes, such as shake, move the coordinates
    // to ensure that the coordinates do not change, store a copy
    double **x = atom->x;
    for (int i = 0; i < nlocal; i++) {
      x_buf[3*i+0] = x[i][0];
      x_buf[3*i+1] = x[i][1];
      x_buf[3*i+2] = x[i][2];
    }

    // calculate the forces
    update->whichflag = 1; // 1 for dynamics
    update->nsteps = 1;
    lmp->init();
    update->integrate->setup_minimal(1);

    // NOTE: can this be done here instead of below?

    if ( strcmp(current_node, "@DEFAULT") == 0 ) {
      // restore the original set of coordinates
      double **x_new = atom->x;
      for (int i = 0; i < nlocal; i++) {
        x_new[i][0] = x_buf[3*i+0];
        x_new[i][1] = x_buf[3*i+1];
        x_new[i][2] = x_buf[3*i+2];
      }
    }
  }

  // pick local atoms from the buffer
  double **f = atom->f;
  for (int i = 0; i < nlocal; i++) {
    forces[3*(atom->tag[i]-1)+0] = f[i][0]*forceconv;
    forces[3*(atom->tag[i]-1)+1] = f[i][1]*forceconv;
    forces[3*(atom->tag[i]-1)+2] = f[i][2]*forceconv;
  }

  // reduce the forces onto rank 0
  MPI_Reduce(forces, forces_reduced, 3*atom->natoms, MPI_DOUBLE, MPI_SUM, 0, world);

  // send the forces through MDI
  if (master) {
    ierr = MDI_Send((char*) forces_reduced, 3*atom->natoms, MDI_DOUBLE, driver_socket);
    if (ierr != 0)
      error->all(FLERR,"MDI: Unable to send atom forces to driver");
  }

  delete [] forces;
  delete [] forces_reduced;
  delete [] x_buf;

}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::receive_forces(Error* error)
{
  double angstrom_to_bohr;
  MDI_Conversion_Factor("angstrom", "bohr", &angstrom_to_bohr);
  double kelvin_to_hartree;
  MDI_Conversion_Factor("kelvin_energy", "hartree", &kelvin_to_hartree);

  double potconv, posconv, forceconv;
  potconv=kelvin_to_hartree/force->boltz;
  posconv=force->angstrom/angstrom_to_bohr;
  forceconv=potconv*posconv;

  double *forces;
  forces = new double[3*atom->natoms];

  if (master) {
    ierr = MDI_Recv((char*) forces, 3*atom->natoms, MDI_DOUBLE, driver_socket);
    if (ierr != 0)
      error->all(FLERR,"MDI: Unable to receive atom forces to driver");
  }
  MPI_Bcast(forces,3*atom->natoms,MPI_DOUBLE,0,world);

  // pick local atoms from the buffer
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    f[i][0] = forces[3*(atom->tag[i]-1)+0]/forceconv;
    f[i][1] = forces[3*(atom->tag[i]-1)+1]/forceconv;
    f[i][2] = forces[3*(atom->tag[i]-1)+2]/forceconv;
  }

  delete [] forces;
}

/* ---------------------------------------------------------------------- */

// NOTE: if keeping add_forces (see NOTE above)
// then could make one replace_add_forces method with an extra "mode" arg
//   for replace or add
// since these 2 methods are nearly identical

void FixMDIEngine::add_forces(Error* error)
{
  double angstrom_to_bohr;
  MDI_Conversion_Factor("angstrom", "bohr", &angstrom_to_bohr);
  double kelvin_to_hartree;
  MDI_Conversion_Factor("kelvin_energy", "hartree", &kelvin_to_hartree);

  double potconv, posconv, forceconv;
  potconv=kelvin_to_hartree/force->boltz;
  posconv=force->angstrom * angstrom_to_bohr;
  forceconv=potconv*posconv;

  double *forces;
  forces = new double[3*atom->natoms];

  if (master) {
    ierr = MDI_Recv((char*) forces, 3*atom->natoms, MDI_DOUBLE, driver_socket);
    if (ierr != 0)
      error->all(FLERR,"MDI: Unable to receive atom forces to driver");
  }
  MPI_Bcast(forces,3*atom->natoms,MPI_DOUBLE,0,world);

  // pick local atoms from the buffer
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    f[i][0] += forces[3*(atom->tag[i]-1)+0]/forceconv;
    f[i][1] += forces[3*(atom->tag[i]-1)+1]/forceconv;
    f[i][2] += forces[3*(atom->tag[i]-1)+2]/forceconv;
  }

  delete [] forces;
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::send_cell(Error* error)
{
  double angstrom_to_bohr;
  MDI_Conversion_Factor("angstrom", "bohr", &angstrom_to_bohr);

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
  /*
  celldata[9 ] = domain->boxlo[0];
  celldata[10] = domain->boxlo[1];
  celldata[11] = domain->boxlo[2];
  */

  // convert the units to bohr

  double unit_conv = force->angstrom * angstrom_to_bohr;
  for (int icell=0; icell < 9; icell++) {
    celldata[icell] *= unit_conv;
  }

  if (master) { 
    ierr = MDI_Send((char*) celldata, 9, MDI_DOUBLE, driver_socket);
    if (ierr != 0)
      error->all(FLERR,"MDI: Unable to send cell dimensions to driver");
  }
}
