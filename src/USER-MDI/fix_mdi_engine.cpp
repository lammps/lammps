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
#include "mdi.h"
#include "memory.h"
#include "modify.h"
#include "update.h"
#include "verlet.h"

#include <limits>
#include <string.h>

using namespace LAMMPS_NS;
using namespace FixConst;

/****************************************************************************/


/***************************************************************
 * create class and parse arguments in LAMMPS script. Syntax:
 * fix ID group-ID mdi_engine [couple <group-ID>]
 ***************************************************************/
FixMDIEngine::FixMDIEngine(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  id_pe(NULL), pe(NULL),
  id_ke(NULL), ke(NULL)
{

  if (narg > 3)
    error->all(FLERR,"Illegal fix mdi command");

  // allocate arrays
  memory->create(add_force,3*atom->natoms,"mdi:add_force");
  for (int i=0; i< 3*atom->natoms; i++) {
    add_force[i] = 0.0;
  }

  // create a new compute pe style
  // id = fix-ID + pe, compute group = all

  master = (comm->me==0) ? 1 : 0;

  // create instance of the Irregular class
  irregular = new Irregular(lmp);

  most_recent_init = 0;
  exit_flag = false;
  local_exit_flag = false;
  target_command = new char[MDI_COMMAND_LENGTH+1];
  command = new char[MDI_COMMAND_LENGTH+1];
  current_node = new char[MDI_COMMAND_LENGTH];
  target_node = new char[MDI_COMMAND_LENGTH];
  strncpy(target_node, "\0", MDI_COMMAND_LENGTH);
  strncpy(current_node, "@DEFAULT", MDI_COMMAND_LENGTH);

  // create and add a compute for PE calculation
  int n_pe = strlen(id) + 4;
  id_pe = new char[n_pe];
  strcpy(id_pe,id);
  strcat(id_pe,"_pe");
  char **newarg_pe = new char*[3];
  newarg_pe[0] = id_pe;
  newarg_pe[1] = (char *) "all";
  newarg_pe[2] = (char *) "pe";
  modify->add_compute(3,newarg_pe);
  delete [] newarg_pe;

  // create and add a compute for KE calculation
  int n_ke = strlen(id) + 4;
  id_ke = new char[n_ke];
  strcpy(id_ke,id);
  strcat(id_ke,"_ke");
  char **newarg_ke = new char*[3];
  newarg_ke[0] = id_ke;
  newarg_ke[1] = (char *) "all";
  newarg_ke[2] = (char *) "ke";
  modify->add_compute(3,newarg_ke);
  delete [] newarg_ke;

  // accept a communicator to the driver
  int ierr;
  if (master) {
    MDI_Accept_Communicator(&driver_socket);
    if (driver_socket <= 0)
      error->all(FLERR,"Unable to connect to driver");
  } else driver_socket=0;

}

/*********************************
 * Clean up on deleting the fix. *
 *********************************/
FixMDIEngine::~FixMDIEngine()
{
  modify->delete_compute(id_pe);
  modify->delete_compute(id_ke);
  delete irregular;
  delete [] id_pe;
  delete [] id_ke;
  delete [] target_command;
  delete [] command;
}

/* ---------------------------------------------------------------------- */
int FixMDIEngine::setmask()
{
  int mask = 0;

  // MD masks
  mask |= POST_INTEGRATE;
  mask |= PRE_REVERSE;
  mask |= POST_FORCE;

  // Minimizer masks
  mask |= MIN_PRE_FORCE;
  mask |= MIN_POST_FORCE;

  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::exchange_forces()
{
  double **f = atom->f;
  const int * const mask  = atom->mask;
  const int nlocal = atom->nlocal;

  // add the forces from the driver
  for (int i=0; i < nlocal; ++i) {
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
  // Confirm that the required computes are available
  int icompute_pe = modify->find_compute(id_pe);
  if (icompute_pe < 0)
    error->all(FLERR,"Potential energy ID for fix mdi does not exist");
  int icompute_ke = modify->find_compute(id_ke);
  if (icompute_pe < 0)
    error->all(FLERR,"Kinetic energy ID for fix mdi does not exist");

  pe = modify->compute[icompute_pe];
  ke = modify->compute[icompute_ke];

  return;

}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::setup(int vflag)
{
  //compute the potential energy
  potential_energy = pe->compute_scalar();
  kinetic_energy = ke->compute_scalar();

  // trigger potential energy computation on next timestep
  pe->addstep(update->ntimestep+1);
  ke->addstep(update->ntimestep+1);

  if ( most_recent_init == 1 ) { // md
    // @PRE-FORCES
    engine_mode("@PRE-FORCES");
  }
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::min_setup(int vflag)
{
  potential_energy = pe->compute_scalar();
  kinetic_energy = ke->compute_scalar();

  // trigger potential energy computation on next timestep
  pe->addstep(update->ntimestep+1);
  ke->addstep(update->ntimestep+1);

  // @FORCES
  engine_mode("@FORCES");
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::post_integrate()
{
  engine_mode("@COORDS");
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::pre_reverse(int eflag, int vflag)
{
  // calculate the energy
  potential_energy = pe->compute_scalar();
  kinetic_energy = ke->compute_scalar();

  // @PRE-FORCES
  engine_mode("@PRE-FORCES");

  // trigger potential energy computation on next timestep
  pe->addstep(update->ntimestep+1);
  ke->addstep(update->ntimestep+1);
}


/* ---------------------------------------------------------------------- */

void FixMDIEngine::min_pre_force(int vflag)
{
  // @COORDS
  engine_mode("@COORDS");
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::min_post_force(int vflag)
{
  // calculate the energy
  potential_energy = pe->compute_scalar();
  kinetic_energy = ke->compute_scalar();

  // @FORCES
  engine_mode("@FORCES");

  // trigger potential energy computation on next timestep
  pe->addstep(update->ntimestep+1);
  ke->addstep(update->ntimestep+1);
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::post_force(int vflag)
{
  if ( most_recent_init == 1 ) { // md
    // @FORCES
    engine_mode("@FORCES");
  }
  else if ( most_recent_init == 2 ) { // optg
    // @FORCES
    engine_mode("@FORCES");
  }
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

  // flag to indicate whether the engine should continue listening for commands at this node
  strncpy(current_node, node, MDI_COMMAND_LENGTH);
  if ( strcmp(target_node,"\0") != 0 and strcmp(target_node, current_node) != 0 ) {
    local_exit_flag = true;
  }

  /* ----------------------------------------------------------------- */
  // Answer commands from the driver
  /* ----------------------------------------------------------------- */

  while (not exit_flag and not local_exit_flag) {

    if (master) { 
      // read the next command from the driver
      ierr = MDI_Recv_Command(command, driver_socket);
      if (ierr != 0)
        error->all(FLERR,"Unable to receive command from driver");
      command[MDI_COMMAND_LENGTH]=0;
    }
    // broadcast the command to the other tasks
    MPI_Bcast(command,MDI_COMMAND_LENGTH,MPI_CHAR,0,world);

    /*
    if (screen)
      fprintf(screen,"MDI command: %s\n",command);
    if (logfile)
      fprintf(logfile,"MDI command: %s:\n",command);
    */

    // confirm that this command is supported at this node
    int command_exists = 0;
    ierr = MDI_Check_Command_Exists(current_node, command, MDI_COMM_NULL, &command_exists);
    if (ierr != 0)
        error->all(FLERR,"Unable to check whether the current command is supported");
    if ( command_exists != 1 )
        error->all(FLERR,"Received a command that is not supported at the current node");

    if (strcmp(command,"STATUS      ") == 0 ) {
      // send the calculation status to the driver
      if (master) {
	ierr = MDI_Send_Command("READY", driver_socket);
        if (ierr != 0)
          error->all(FLERR,"Unable to return status to driver");
      }
    }
    else if (strcmp(command,">NATOMS") == 0 ) {
      // receive the number of atoms from the driver
      if (master) {
        ierr = MDI_Recv((char*) &atom->natoms, 1, MDI_INT, driver_socket);
        if (ierr != 0)
          error->all(FLERR,"Unable to receive number of atoms from driver");
      }
      MPI_Bcast(&atom->natoms,1,MPI_INT,0,world);
    }
    else if (strcmp(command,"<NATOMS") == 0 ) {
      // send the number of atoms to the driver
      if (master) {
        int64_t mdi_natoms = atom->natoms;
        ierr = MDI_Send((char*) &mdi_natoms, 1, MDI_INT64_T, driver_socket);
        if (ierr != 0)
          error->all(FLERR,"Unable to send number of atoms to driver");
      }
    }
    else if (strcmp(command,"<NTYPES") == 0 ) {
      // send the number of atom types to the driver
      if (master) {
        ierr = MDI_Send((char*) &atom->ntypes, 1, MDI_INT, driver_socket);
        if (ierr != 0)
          error->all(FLERR,"Unable to send number of atom types to driver");
      }
    }
    else if (strcmp(command,"<TYPES") == 0 ) {
      // send the atom types
      send_types(error);
    }
    else if (strcmp(command,"<LABELS") == 0 ) {
      // send the atom labels
      send_labels(error);
    }
    else if (strcmp(command,"<MASSES") == 0 ) {
      // send the atom types
      send_masses(error);
    }
    else if (strcmp(command,"<CELL") == 0 ) {
      // send the cell dimensions to the driver
      send_cell(error);
    }
    else if (strcmp(command,">COORDS") == 0 ) {
      // receive the coordinate information
      receive_coordinates(error);
    }
    else if (strcmp(command,"<COORDS") == 0 ) {
      // send the coordinate information
      send_coordinates(error);
    }
    else if (strcmp(command,"<CHARGES") == 0 ) {
      // send the charges
      send_charges(error);
    }
    else if (strcmp(command,"<ENERGY") == 0 ) {
      // send the total energy to the driver
      send_energy(error);
    }
    else if (strcmp(command,"<FORCES") == 0 ) {
      // send the forces to the driver
      send_forces(error);
    }
    else if (strcmp(command,">FORCES") == 0 ) {
      // receive the forces from the driver
      receive_forces(error);
    }
    else if (strcmp(command,"+FORCES") == 0 ) {
      // receive additional forces from the driver
      // these are added prior to SHAKE or other post-processing
      add_forces(error);
    }
    else if (strcmp(command,"@INIT_MD") == 0 ) {
      if ( most_recent_init != 0 ) {
	error->all(FLERR,"MDI is already performing a simulation");
      }

      // initialize a new MD simulation
      most_recent_init = 1;
      local_exit_flag = true;
    }
    else if (strcmp(command,"@INIT_OPTG") == 0 ) {
      if ( most_recent_init != 0 ) {
	error->all(FLERR,"MDI is already performing a simulation");
      }

      // initialize a new geometry optimization
      most_recent_init = 2;
      local_exit_flag = true;
      //optg_init(error);
    }
    else if (strcmp(command,"@") == 0 ) {
      strncpy(target_node, "\0", MDI_COMMAND_LENGTH);
      local_exit_flag = true;
    }
    else if (strcmp(command,"<@") == 0 ) {
      if (master) {
	ierr = MDI_Send(current_node, MDI_NAME_LENGTH, MDI_CHAR, driver_socket);
	if (ierr != 0)
	  error->all(FLERR,"Unable to send node to driver");
      }
    }
    else if (strcmp(command,"<KE") == 0 ) {
      // send the kinetic energy to the driver
      send_ke(error);
    }
    else if (strcmp(command,"<PE") == 0 ) {
      // send the potential energy to the driver
      send_pe(error);
    }
    else if (strcmp(command,"@COORDS") == 0 ) {
      strncpy(target_node, "@COORDS", MDI_COMMAND_LENGTH);
      local_exit_flag = true;
    }
    else if (strcmp(command,"@PRE-FORCES") == 0 ) {
      strncpy(target_node, "@PRE-FORCES", MDI_COMMAND_LENGTH);
      local_exit_flag = true;
    }
    else if (strcmp(command,"@FORCES") == 0 ) {
      strncpy(target_node, "@FORCES", MDI_COMMAND_LENGTH);
      local_exit_flag = true;
    }
    else if (strcmp(command,"EXIT_SIM") == 0 ) {
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
    }
    else if (strcmp(command,"EXIT") == 0 ) {
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
    }
    else {
      // the command is not supported
      error->all(FLERR,"Unknown command from driver");
    }

    // check if the target node is something other than the current node
    if ( strcmp(target_node,"\0") != 0 and strcmp(target_node, current_node) != 0 ) {
      local_exit_flag = true;
    }

  }

  // a local exit has completed, so turn off the local exit flag
  local_exit_flag = false;

  return command;

}


void FixMDIEngine::receive_coordinates(Error* error)
{
  // unable to convert to atomic units if LAMMPS is using lj units
  if (strcmp(update->unit_style,"lj") == 0)
    error->all(FLERR,"The >COORDS MDI command does not support lj units");

  double posconv;
  double angstrom_to_bohr;
  MDI_Conversion_Factor("angstrom", "bohr", &angstrom_to_bohr);
  posconv=force->angstrom/angstrom_to_bohr;

  // create a buffer to hold the coordinates
  double *buffer;
  buffer = new double[3*atom->natoms];

  if (master) {
    ierr = MDI_Recv((char*) buffer, 3*atom->natoms, MDI_DOUBLE, driver_socket);
    if (ierr != 0)
      error->all(FLERR,"Unable to receive coordinates from driver");
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

  // move atoms to new processors via irregular()
  // only needed if migrate_check() says an atom moves to far
  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  if (irregular->migrate_check()) irregular->migrate_atoms();
  if (domain->triclinic) domain->lamda2x(atom->nlocal);

  delete [] buffer;
}


void FixMDIEngine::send_coordinates(Error* error)
{
  // unable to convert to atomic units if LAMMPS is using lj units
  if (strcmp(update->unit_style,"lj") == 0)
    error->all(FLERR,"The <COORDS MDI command does not support lj units");

  double posconv;
  double angstrom_to_bohr;
  MDI_Conversion_Factor("angstrom", "bohr", &angstrom_to_bohr);
  posconv=force->angstrom/angstrom_to_bohr;

  double *coords;
  double *coords_reduced;

  coords = new double[3*atom->natoms];
  coords_reduced = new double[3*atom->natoms];

  // zero the coordinates array
  for (int icoord = 0; icoord < 3*atom->natoms; icoord++) {
    coords[icoord] = 0.0;
  }

  // pick local atoms from the buffer
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    coords[3*(atom->tag[i]-1)+0] = x[i][0]/posconv;
    coords[3*(atom->tag[i]-1)+1] = x[i][1]/posconv;
    coords[3*(atom->tag[i]-1)+2] = x[i][2]/posconv;
  }

  MPI_Reduce(coords, coords_reduced, 3*atom->natoms, MPI_DOUBLE, MPI_SUM, 0, world);

  if (master) {
    ierr = MDI_Send((char*) coords_reduced, 3*atom->natoms, MDI_DOUBLE, driver_socket);
    if (ierr != 0)
      error->all(FLERR,"Unable to send coordinates to driver");
  }

  delete [] coords;
  delete [] coords_reduced;
}


void FixMDIEngine::send_charges(Error* error)
{
  double *charges;
  double *charges_reduced;

  charges = new double[atom->natoms];
  charges_reduced = new double[atom->natoms];

  // zero the charges array
  for (int icharge = 0; icharge < atom->natoms; icharge++) {
    charges[icharge] = 0.0;
  }

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
      error->all(FLERR,"Unable to send charges to driver");
  }

  delete [] charges;
  delete [] charges_reduced;
}


void FixMDIEngine::send_energy(Error* error)
{
  // unable to convert to atomic units if LAMMPS is using lj units
  if (strcmp(update->unit_style,"lj") == 0)
    error->all(FLERR,"The <ENERGY MDI command does not support lj units");

  double kelvin_to_hartree;
  MDI_Conversion_Factor("kelvin_energy", "hartree", &kelvin_to_hartree);

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
      error->all(FLERR,"Unable to send potential energy to driver");
  }
}


void FixMDIEngine::send_pe(Error* error)
{
  // unable to convert to atomic units if LAMMPS is using lj units
  if (strcmp(update->unit_style,"lj") == 0)
    error->all(FLERR,"The <PE MDI command does not support lj units");

  double kelvin_to_hartree;
  MDI_Conversion_Factor("kelvin_energy", "hartree", &kelvin_to_hartree);

  double pe = potential_energy;
  double *send_energy = &pe;

  // convert the energy to atomic units
  pe *= kelvin_to_hartree/force->boltz;

  if (master) {
    ierr = MDI_Send((char*) send_energy, 1, MDI_DOUBLE, driver_socket);
    if (ierr != 0)
      error->all(FLERR,"Unable to send potential energy to driver");
  }
}


void FixMDIEngine::send_ke(Error* error)
{
  // unable to convert to atomic units if LAMMPS is using lj units
  if (strcmp(update->unit_style,"lj") == 0)
    error->all(FLERR,"The <KE MDI command does not support lj units");

  double kelvin_to_hartree;
  MDI_Conversion_Factor("kelvin_energy", "hartree", &kelvin_to_hartree);

  double ke = kinetic_energy;
  double *send_energy = &ke;

  // convert the energy to atomic units
  ke *= kelvin_to_hartree/force->boltz;

  if (master) {
    ierr = MDI_Send((char*) send_energy, 1, MDI_DOUBLE, driver_socket);
    if (ierr != 0)
      error->all(FLERR,"Unable to send potential energy to driver");
  }
}


void FixMDIEngine::send_types(Error* error)
{
  int * const type = atom->type;

  if (master) { 
    ierr = MDI_Send((char*) type, atom->natoms, MDI_INT, driver_socket);
    if (ierr != 0)
      error->all(FLERR,"Unable to send atom types to driver");
  }
}


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
      error->all(FLERR,"Unable to send atom types to driver");
  }

  delete [] labels;
}


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
      error->all(FLERR,"Unable to send atom masses to driver");
  }
  
  delete [] mass_by_atom;
  delete [] mass_by_atom_reduced;
}


void FixMDIEngine::send_forces(Error* error)
{
  // unable to convert to atomic units if LAMMPS is using lj units
  if (strcmp(update->unit_style,"lj") == 0)
    error->all(FLERR,"The <FORCES MDI command does not support lj units");

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
  for (int iforce = 0; iforce < 3*atom->natoms; iforce++) {
    forces[iforce] = 0.0;
  }

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
      error->all(FLERR,"Unable to send atom forces to driver");
  }

  if ( strcmp(current_node, "@DEFAULT") == 0 ) {
    // restore the original set of coordinates
    double **x_new = atom->x;
    for (int i = 0; i < nlocal; i++) {
      x_new[i][0] = x_buf[3*i+0];
      x_new[i][1] = x_buf[3*i+1];
      x_new[i][2] = x_buf[3*i+2];
    }
  }

  delete [] forces;
  delete [] forces_reduced;
  delete [] x_buf;

}


void FixMDIEngine::receive_forces(Error* error)
{
  // unable to convert to atomic units if LAMMPS is using lj units
  if (strcmp(update->unit_style,"lj") == 0)
    error->all(FLERR,"The >FORCES MDI command does not support lj units");

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
      error->all(FLERR,"Unable to receive atom forces to driver");
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


void FixMDIEngine::add_forces(Error* error)
{
  // unable to convert to atomic units if LAMMPS is using lj units
  if (strcmp(update->unit_style,"lj") == 0)
    error->all(FLERR,"The >+FORCES MDI command does not support lj units");

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
      error->all(FLERR,"Unable to receive atom forces to driver");
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


void FixMDIEngine::send_cell(Error* error)
{
  // unable to convert to atomic units if LAMMPS is using lj units
  if (strcmp(update->unit_style,"lj") == 0)
    error->all(FLERR,"The <CELL MDI command does not support lj units");

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
      error->all(FLERR,"Unable to send cell dimensions to driver");
  }
}
