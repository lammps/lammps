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

#include "fix_mdi_engine.h"
#include "library_mdi.h"

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

enum { NONE, REAL, METAL };    // LAMMPS units which MDI supports

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixMDIEngine::FixMDIEngine(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), id_pe(nullptr), id_ke(nullptr), pe(nullptr), ke(nullptr)
{
  if (narg != 3) error->all(FLERR, "Illegal fix mdi command");

  // The 2 atomic-scale units LAMMPS has are:
  //    real: coords = Ang, eng = Kcal/mole, force = Kcal/mole/Ang
  //    metal: coords = Ang, eng = eV, force = eV/Ang

  lmpunits = NONE;
  if (strcmp(update->unit_style, "real") == 0) lmpunits = REAL;
  if (strcmp(update->unit_style, "metal") == 0) lmpunits = METAL;
  if (lmpunits == NONE) error->all(FLERR, "MDI requires real or metal units");

  // MDI setup

  most_recent_init = 0;
  exit_flag = false;
  local_exit_flag = false;
  target_command = new char[MDI_COMMAND_LENGTH + 1];
  command = new char[MDI_COMMAND_LENGTH + 1];
  current_node = new char[MDI_COMMAND_LENGTH];
  target_node = new char[MDI_COMMAND_LENGTH];
  strncpy(target_node, "\0", MDI_COMMAND_LENGTH);
  strncpy(current_node, "@DEFAULT", MDI_COMMAND_LENGTH);

  // register the execute_command function with MDI

  MDI_Set_execute_command_func(lammps_execute_mdi_command, this);

  // accept a communicator to the driver
  // master = 1 for proc 0, otherwise 0

  master = (comm->me == 0) ? 1 : 0;

  MDI_Accept_communicator(&driver_socket);
  if (driver_socket <= 0) error->all(FLERR, "Unable to connect to driver");

  // create computes for KE and PE

  id_pe = utils::strdup(std::string(id) + "_pe");
  modify->add_compute(fmt::format("{} all pe", id_pe));

  id_ke = utils::strdup(std::string(id) + "_ke");
  modify->add_compute(fmt::format("{} all ke", id_ke));

  // irregular class and data structs used by MDI

  irregular = new Irregular(lmp);
  add_force = nullptr;
}

/* ---------------------------------------------------------------------- */

FixMDIEngine::~FixMDIEngine()
{
  delete[] target_command;
  delete[] command;
  delete[] current_node;
  delete[] target_node;

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
  mask |= POST_FORCE;
  mask |= MIN_PRE_FORCE;
  mask |= MIN_POST_FORCE;

  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::exchange_forces()
{
  double **f = atom->f;
  const int *const mask = atom->mask;
  const int nlocal = atom->nlocal;

  // add forces from the driver

  for (int i = 0; i < nlocal; ++i) {
    if (mask[i] & groupbit) {
      f[i][0] += add_force[3 * (atom->tag[i] - 1) + 0];
      f[i][1] += add_force[3 * (atom->tag[i] - 1) + 1];
      f[i][2] += add_force[3 * (atom->tag[i] - 1) + 2];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::init()
{
  // confirm that two required computes are still available

  int icompute_pe = modify->find_compute(id_pe);
  if (icompute_pe < 0) error->all(FLERR, "Potential energy ID for fix mdi/engine does not exist");
  int icompute_ke = modify->find_compute(id_ke);
  if (icompute_pe < 0) error->all(FLERR, "Kinetic energy ID for fix mdi/engine does not exist");

  pe = modify->compute[icompute_pe];
  ke = modify->compute[icompute_ke];

  // one-time allocation of add_force array

  if (!add_force) {
    int64_t ncoords = 3 * atom->natoms;
    memory->create(add_force, ncoords, "mdi/engine:add_force");
    for (int64_t i = 0; i < ncoords; i++) add_force[i] = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::min_setup(int /* vflag */)
{
  engine_mode("@FORCES");
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::post_integrate()
{
  engine_mode("@COORDS");
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::min_pre_force(int /* vflag */)
{
  engine_mode("@COORDS");
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::min_post_force(int /* vflag */)
{
  engine_mode("@FORCES");
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::post_force(int /* vflag */)
{
  if (most_recent_init == 1)
    engine_mode("@FORCES");
  else if (most_recent_init == 2)
    engine_mode("@FORCES");
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
  // confirm this command is supported at this node

  int command_exists = 1;
  if (master) {
    ierr = MDI_Check_command_exists(current_node, command, MDI_COMM_NULL, &command_exists);
  }
  if (ierr != 0) error->all(FLERR, "MDI: Unable to check whether current command is supported");
  if (command_exists != 1)
    error->all(FLERR, "MDI: Received a command that is unsupported at current node");

  // respond to any driver command

  if (strcmp(command, ">NATOMS") == 0) {
    int64_t mdi_natoms = 0;
    ierr = MDI_Recv((char *) &mdi_natoms, 1, MDI_INT64_T, mdicomm);
    if (ierr != 0) error->all(FLERR, "MDI: Unable to receive number of atoms from driver");
    error->all(FLERR, "MDI: '>NATOMS' driver command not (yet) supported");
    // FIXME: to import the number of atoms, more steps than below are needed for LAMMPS.
    //        also a check for overflow is needed in case natoms is 32-bit
    atom->natoms = mdi_natoms;
    MPI_Bcast(&atom->natoms, 1, MPI_LMP_BIGINT, 0, world);

  } else if (strcmp(command, "<NATOMS") == 0) {
    int64_t mdi_natoms = atom->natoms;
    ierr = MDI_Send((char *) &mdi_natoms, 1, MDI_INT64_T, mdicomm);
    if (ierr != 0) error->all(FLERR, "MDI: Unable to send number of atoms to driver");

  } else if (strcmp(command, "<NTYPES") == 0) {
    ierr = MDI_Send((char *) &atom->ntypes, 1, MDI_INT, mdicomm);
    if (ierr != 0) error->all(FLERR, "MDI: Unable to send number of atom types to driver");

  } else if (strcmp(command, "<TYPES") == 0) {
    send_types(error);

  } else if (strcmp(command, "<LABELS") == 0) {
    send_labels(error);

  } else if (strcmp(command, "<MASSES") == 0) {
    send_masses(error);

  } else if (strcmp(command, "<CELL") == 0) {
    send_cell(error);

  } else if (strcmp(command, ">CELL") == 0) {
    receive_cell(error);

  } else if (strcmp(command, "<CELL_DISPL") == 0) {
    send_celldispl(error);

  } else if (strcmp(command, ">CELL_DISPL") == 0) {
    receive_celldispl(error);

  } else if (strcmp(command, ">COORDS") == 0) {
    receive_coordinates(error);

  } else if (strcmp(command, "<COORDS") == 0) {
    send_coordinates(error);

  } else if (strcmp(command, "<CHARGES") == 0) {
    send_charges(error);

  } else if (strcmp(command, "<ENERGY") == 0) {
    send_energy(error);

  } else if (strcmp(command, "<FORCES") == 0) {
    send_forces(error);

    // replace current forces with forces received from the driver

  } else if (strcmp(command, ">FORCES") == 0) {
    receive_forces(error, 0);

    // add forces received from the driver to current forces

  } else if (strcmp(command, ">+FORCES") == 0) {
    receive_forces(error, 1);

    // initialize new MD simulation or minimization
    // return control to return to mdi/engine

  } else if (strcmp(command, "@INIT_MD") == 0) {
    if (most_recent_init != 0) error->all(FLERR, "MDI: MDI is already performing a simulation");
    most_recent_init = 1;
    local_exit_flag = true;

    // initialize new energy minimization
    // return control to return to mdi/engine

  } else if (strcmp(command, "@INIT_OPTG") == 0) {
    if (most_recent_init != 0) error->all(FLERR, "MDI: MDI is already performing a simulation");
    most_recent_init = 2;
    local_exit_flag = true;

  } else if (strcmp(command, "@") == 0) {
    strncpy(target_node, "\0", MDI_COMMAND_LENGTH);
    local_exit_flag = true;

  } else if (strcmp(command, "<@") == 0) {
    ierr = MDI_Send(current_node, MDI_NAME_LENGTH, MDI_CHAR, mdicomm);
    if (ierr != 0) error->all(FLERR, "MDI: Unable to send node to driver");

  } else if (strcmp(command, "<KE") == 0) {
    send_ke(error);

  } else if (strcmp(command, "<PE") == 0) {
    send_pe(error);

  } else if (strcmp(command, "@DEFAULT") == 0) {
    most_recent_init = 0;

    strncpy(target_node, "@DEFAULT", MDI_COMMAND_LENGTH);
    local_exit_flag = true;

    // are we in the middle of a geometry optimization?
    if (most_recent_init == 2) {
      // ensure that the energy and force tolerances are met
      update->etol = std::numeric_limits<double>::max();
      update->ftol = std::numeric_limits<double>::max();

      // set the maximum number of force evaluations to 0
      update->max_eval = 0;
    }

  } else if (strcmp(command, "@COORDS") == 0) {
    strncpy(target_node, "@COORDS", MDI_COMMAND_LENGTH);
    local_exit_flag = true;

  } else if (strcmp(command, "@FORCES") == 0) {
    strncpy(target_node, "@FORCES", MDI_COMMAND_LENGTH);
    local_exit_flag = true;

  } else if (strcmp(command, "EXIT") == 0) {
    // exit the driver code
    exit_flag = true;

    // are we in the middle of a geometry optimization?
    if (most_recent_init == 2) {
      // ensure that the energy and force tolerances are met
      update->etol = std::numeric_limits<double>::max();
      update->ftol = std::numeric_limits<double>::max();

      // set the maximum number of force evaluations to 0
      update->max_eval = 0;
    }

  } else {
    error->all(FLERR, "MDI: Unknown command from driver");
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

  strncpy(current_node, node, MDI_COMMAND_LENGTH);
  if (strcmp(target_node, "\0") != 0 && strcmp(target_node, current_node) != 0)
    local_exit_flag = true;

  // respond to commands from the driver

  while (not exit_flag and not local_exit_flag) {

    // read the next command from the driver
    // all procs call this, but only proc 0 receives the command

    ierr = MDI_Recv_command(command, driver_socket);
    if (ierr != 0) error->all(FLERR, "MDI: Unable to receive command from driver");

    // broadcast command to the other MPI tasks

    MPI_Bcast(command, MDI_COMMAND_LENGTH, MPI_CHAR, 0, world);

    // execute the command

    this->execute_command(command, driver_socket);

    // check if the target node is something other than the current node

    if (strcmp(target_node, "\0") != 0 && strcmp(target_node, current_node) != 0)
      local_exit_flag = true;
  }

  // local exit occured so turn off local exit flag

  local_exit_flag = false;

  return command;
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::receive_coordinates(Error *error)
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

  ierr = MDI_Recv((char *) buffer, 3 * atom->natoms, MDI_DOUBLE, driver_socket);
  if (ierr != 0) error->all(FLERR, "MDI: Unable to receive coordinates from driver");
  MPI_Bcast(buffer, 3 * atom->natoms, MPI_DOUBLE, 0, world);

  // pick local atoms from the buffer

  double **x = atom->x;
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

void FixMDIEngine::send_coordinates(Error *error)
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
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    coords[3 * (atom->tag[i] - 1) + 0] = x[i][0] / posconv;
    coords[3 * (atom->tag[i] - 1) + 1] = x[i][1] / posconv;
    coords[3 * (atom->tag[i] - 1) + 2] = x[i][2] / posconv;
  }

  MPI_Reduce(coords, coords_reduced, 3 * atom->natoms, MPI_DOUBLE, MPI_SUM, 0, world);

  ierr = MDI_Send((char *) coords_reduced, 3 * atom->natoms, MDI_DOUBLE, driver_socket);
  if (ierr != 0) error->all(FLERR, "MDI: Unable to send coordinates to driver");

  memory->destroy(coords);
  memory->destroy(coords_reduced);
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::send_charges(Error *error)
{
  double *charges;
  double *charges_reduced;

  memory->create(charges, atom->natoms, "mdi/engine:charges");
  memory->create(charges_reduced, atom->natoms, "mdi/engine:charges_reduced");

  // zero the charges array

  for (int icharge = 0; icharge < atom->natoms; icharge++) charges[icharge] = 0.0;

  // pick local atoms from the buffer

  double *charge = atom->q;
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) { charges[atom->tag[i] - 1] = charge[i]; }

  MPI_Reduce(charges, charges_reduced, atom->natoms, MPI_DOUBLE, MPI_SUM, 0, world);

  ierr = MDI_Send((char *) charges_reduced, atom->natoms, MDI_DOUBLE, driver_socket);
  if (ierr != 0) error->all(FLERR, "MDI: Unable to send charges to driver");

  memory->destroy(charges);
  memory->destroy(charges_reduced);
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::send_energy(Error *error)
{
  // get conversion factor to atomic units
  double energy_conv = 1.0;
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

  ierr = MDI_Send((char *) send_energy, 1, MDI_DOUBLE, driver_socket);
  if (ierr != 0) error->all(FLERR, "MDI: Unable to send potential energy to driver");
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::send_pe(Error *error)
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

  double potential_energy = pe->compute_scalar();
  double *send_energy = &potential_energy;

  // convert the energy to atomic units
  potential_energy *= energy_conv;

  ierr = MDI_Send((char *) send_energy, 1, MDI_DOUBLE, driver_socket);
  if (ierr != 0) error->all(FLERR, "MDI: Unable to send potential energy to driver");
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::send_ke(Error *error)
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
  double *send_energy = &kinetic_energy;

  // convert the energy to atomic units
  kinetic_energy *= energy_conv;

  ierr = MDI_Send((char *) send_energy, 1, MDI_DOUBLE, driver_socket);
  if (ierr != 0) error->all(FLERR, "MDI: Unable to send potential energy to driver");
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::send_types(Error *error)
{
  int *const type = atom->type;

  ierr = MDI_Send((char *) type, atom->natoms, MDI_INT, driver_socket);
  if (ierr != 0) error->all(FLERR, "MDI: Unable to send atom types to driver");
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::send_labels(Error *error)
{
  char *labels = new char[atom->natoms * MDI_LABEL_LENGTH];
  memset(labels, ' ', atom->natoms * MDI_LABEL_LENGTH);

  for (int iatom = 0; iatom < atom->natoms; iatom++) {
    std::string label = std::to_string(atom->type[iatom]);
    int label_len = std::min(int(label.length()), MDI_LABEL_LENGTH);
    strncpy(&labels[iatom * MDI_LABEL_LENGTH], label.c_str(), label_len);
  }

  ierr = MDI_Send(labels, atom->natoms * MDI_LABEL_LENGTH, MDI_CHAR, driver_socket);
  if (ierr != 0) error->all(FLERR, "MDI: Unable to send atom types to driver");

  delete[] labels;
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::send_masses(Error *error)
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

  ierr = MDI_Send((char *) mass_by_atom_reduced, atom->natoms, MDI_DOUBLE, driver_socket);
  if (ierr != 0) error->all(FLERR, "MDI: Unable to send atom masses to driver");

  memory->destroy(mass_by_atom);
  memory->destroy(mass_by_atom_reduced);
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::send_forces(Error *error)
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

  int nlocal = atom->nlocal;
  int64_t ncoords = 3 * atom->natoms;

  memory->create(forces, ncoords, "mdi/engine:forces");
  memory->create(forces_reduced, ncoords, "mdi/engine:forces_reduced");
  x_buf = new double[3 * nlocal];

  // zero the forces array

  for (int iforce = 0; iforce < 3 * atom->natoms; iforce++) forces[iforce] = 0.0;

  // if not at a node, calculate the forces

  if (strcmp(current_node, "@DEFAULT") == 0) {
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

    if (strcmp(current_node, "@DEFAULT") == 0) {
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
  ierr = MDI_Send((char *) forces_reduced, 3 * atom->natoms, MDI_DOUBLE, driver_socket);
  if (ierr != 0) error->all(FLERR, "MDI: Unable to send atom forces to driver");

  memory->destroy(forces);
  memory->destroy(forces_reduced);
  delete[] x_buf;
}

/* ---------------------------------------------------------------------- */

// Receive forces from the driver
//    mode = 0: replace current forces with forces from driver
//    mode = 1: add forces from driver to current forces

void FixMDIEngine::receive_forces(Error *error, int mode)
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

  ierr = MDI_Recv((char *) forces, 3 * atom->natoms, MDI_DOUBLE, driver_socket);
  if (ierr != 0) error->all(FLERR, "MDI: Unable to receive atom forces to driver");
  MPI_Bcast(forces, 3 * atom->natoms, MPI_DOUBLE, 0, world);

  // pick local atoms from the buffer
  double **f = atom->f;
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

void FixMDIEngine::send_cell(Error *error)
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

  ierr = MDI_Send((char *) celldata, 9, MDI_DOUBLE, driver_socket);
  if (ierr != 0) error->all(FLERR, "MDI: Unable to send cell dimensions to driver");
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::receive_cell(Error *error)
{
  double celldata[9];

  // receive the new cell vector from the driver
  ierr = MDI_Recv((char *) celldata, 9, MDI_DOUBLE, driver_socket);
  if (ierr != 0) error->all(FLERR, "MDI: Unable to send cell dimensions to driver");
  MPI_Bcast(&celldata[0], 9, MPI_DOUBLE, 0, world);

  double angstrom_to_bohr;
  MDI_Conversion_factor("angstrom", "bohr", &angstrom_to_bohr);
  double unit_conv = force->angstrom * angstrom_to_bohr;
  for (int icell = 0; icell < 9; icell++) { celldata[icell] /= unit_conv; }

  // ensure that the new cell vector is orthogonal
  double small = std::numeric_limits<double>::min();
  if (fabs(celldata[1]) > small or fabs(celldata[2]) > small or fabs(celldata[3]) > small or
      fabs(celldata[5]) > small or fabs(celldata[6]) > small or fabs(celldata[7]) > small) {
    error->all(FLERR, "MDI: LAMMPS currently only supports the >CELL command for orthogonal cell vectors");
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

void FixMDIEngine::send_celldispl(Error *error)
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

  ierr = MDI_Send((char *) celldata, 3, MDI_DOUBLE, driver_socket);
  if (ierr != 0) error->all(FLERR, "MDI: Unable to send cell displacement to driver");
}

/* ---------------------------------------------------------------------- */

void FixMDIEngine::receive_celldispl(Error *error)
{
  // receive the cell displacement from the driver
  double celldata[3];
  ierr = MDI_Recv((char *) celldata, 3, MDI_DOUBLE, driver_socket);
  if (ierr != 0) error->all(FLERR, "MDI: Unable to receive cell displacement from driver");
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
