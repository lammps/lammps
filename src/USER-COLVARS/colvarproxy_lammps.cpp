// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.


#include <mpi.h>
#include "lammps.h"
#include "atom.h"
#include "error.h"
#include "output.h"
#include "random_park.h"

#include "fix_colvars.h"

#include "colvarmodule.h"
#include "colvar.h"
#include "colvarbias.h"
#include "colvaratoms.h"
#include "colvarproxy.h"
#include "colvarproxy_lammps.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <cerrno>
#include <cstdio>
#include <cstring>

#include <iostream>
#include <sstream>
#include <string>

#define HASH_FAIL  -1

////////////////////////////////////////////////////////////////////////
// local helper functions

// safely move filename to filename.extension
static int my_backup_file(const char *filename, const char *extension)
{
  struct stat sbuf;
  if (stat(filename, &sbuf) == 0) {
    if (!extension) extension = ".BAK";
    char *backup = new char[strlen(filename)+strlen(extension)+1];
    strcpy(backup, filename);
    strcat(backup, extension);
#if defined(_WIN32) && !defined(__CYGWIN__)
    remove(backup);
#endif
    if (rename(filename,backup)) {
      char *sys_err_msg = strerror(errno);
      if (!sys_err_msg)  sys_err_msg = (char *) "(unknown error)";
      fprintf(stderr,"Error renaming file %s to %s: %s\n",
              filename, backup, sys_err_msg);
      delete [] backup;
      return COLVARS_ERROR;
    }
    delete [] backup;
  }
  return COLVARS_OK;
}

////////////////////////////////////////////////////////////////////////

colvarproxy_lammps::colvarproxy_lammps(LAMMPS_NS::LAMMPS *lmp,
                                       const char *inp_name,
                                       const char *out_name,
                                       const int seed,
                                       const double temp,
                                       MPI_Comm root2root)
  : _lmp(lmp), inter_comm(root2root)
{
  if (cvm::debug())
    log("Initializing the colvars proxy object.\n");

  _random = new LAMMPS_NS::RanPark(lmp,seed);

  first_timestep=true;
  total_force_requested=false;
  previous_step=-1;
  t_target=temp;
  do_exit=false;
  restart_every=0;

  // User-scripted forces are not available in LAMMPS
  force_script_defined = false;
  have_scripts = false;

  // set input restart name and strip the extension, if present
  input_prefix_str = std::string(inp_name ? inp_name : "");
  if (input_prefix_str.rfind(".colvars.state") != std::string::npos)
    input_prefix_str.erase(input_prefix_str.rfind(".colvars.state"),
                            std::string(".colvars.state").size());

  // output prefix is always given
  output_prefix_str = std::string(out_name);
  // not so for restarts
  restart_output_prefix_str = std::string("rest");

  // check if it is possible to save output configuration
  if ((!output_prefix_str.size()) && (!restart_output_prefix_str.size())) {
    fatal_error("Error: neither the final output state file or "
                "the output restart file could be defined, exiting.\n");
  }

  // try to extract a restart prefix from a potential restart command.
  LAMMPS_NS::Output *outp = _lmp->output;
  if ((outp->restart_every_single > 0) && (outp->restart1 != 0)) {
      restart_output_prefix_str = std::string(outp->restart1);
  } else if  ((outp->restart_every_double > 0) && (outp->restart2a != 0)) {
    restart_output_prefix_str = std::string(outp->restart2a);
  }
  // trim off unwanted stuff from the restart prefix
  if (restart_output_prefix_str.rfind(".*") != std::string::npos)
    restart_output_prefix_str.erase(restart_output_prefix_str.rfind(".*"),2);

  // initialize multi-replica support, if available
  if (replica_enabled()) {
    MPI_Comm_rank(inter_comm, &inter_me);
    MPI_Comm_size(inter_comm, &inter_num);
  }

  if (cvm::debug())
    log("Done initializing the colvars proxy object.\n");
}


void colvarproxy_lammps::init(const char *conf_file)
{
  version_int = get_version_from_string(COLVARPROXY_VERSION);

  // create the colvarmodule instance
  colvars = new colvarmodule(this);

  cvm::log("Using LAMMPS interface, version "+
           cvm::to_str(COLVARPROXY_VERSION)+".\n");

  my_angstrom  = _lmp->force->angstrom;
  my_boltzmann = _lmp->force->boltz;
  my_timestep  = _lmp->update->dt * _lmp->force->femtosecond;

  // TODO move one or more of these to setup() if needed
  colvars->read_config_file(conf_file);
  colvars->setup_input();
  colvars->setup_output();

  if (_lmp->update->ntimestep != 0) {
    cvm::log("Setting initial step number from LAMMPS: "+
             cvm::to_str(_lmp->update->ntimestep)+"\n");
    colvars->it = colvars->it_restart =
      static_cast<cvm::step_number>(_lmp->update->ntimestep);
  }

  if (cvm::debug()) {
    cvm::log("atoms_ids = "+cvm::to_str(atoms_ids)+"\n");
    cvm::log("atoms_ncopies = "+cvm::to_str(atoms_ncopies)+"\n");
    cvm::log("atoms_positions = "+cvm::to_str(atoms_positions)+"\n");
    cvm::log(cvm::line_marker);
    cvm::log("Info: done initializing the colvars proxy object.\n");
  }
}

void colvarproxy_lammps::add_config_file(const char *conf_file)
{
  colvars->read_config_file(conf_file);
}

void colvarproxy_lammps::add_config_string(const std::string &conf)
{
  colvars->read_config_string(conf);
}

colvarproxy_lammps::~colvarproxy_lammps()
{
  delete _random;
  if (colvars != NULL) {
    delete colvars;
    colvars = NULL;
  }
}

// re-initialize data where needed
int colvarproxy_lammps::setup()
{
  my_timestep  = _lmp->update->dt * _lmp->force->femtosecond;
  return colvars->setup();
}

// trigger colvars computation
double colvarproxy_lammps::compute()
{
  if (cvm::debug()) {
    cvm::log(std::string(cvm::line_marker)+
        "colvarproxy_lammps step no. "+
        cvm::to_str(_lmp->update->ntimestep)+" [first - last = "+
        cvm::to_str(_lmp->update->beginstep)+" - "+
        cvm::to_str(_lmp->update->endstep)+"]\n");
  }

  if (first_timestep) {
    first_timestep = false;
  } else {
    // Use the time step number from LAMMPS Update object
    if ( _lmp->update->ntimestep - previous_step == 1 )
      colvars->it++;
    // Other cases could mean:
    // - run 0
    // - beginning of a new run statement
    // then the internal counter should not be incremented
  }
  previous_step = _lmp->update->ntimestep;

  unit_cell_x.set(_lmp->domain->xprd, 0.0, 0.0);
  unit_cell_y.set(0.0, _lmp->domain->yprd, 0.0);
  unit_cell_z.set(0.0, 0.0, _lmp->domain->zprd);

  if (_lmp->domain->xperiodic == 0 && _lmp->domain->yperiodic == 0 &&
      _lmp->domain->zperiodic == 0) {
    boundaries_type = boundaries_non_periodic;
    reset_pbc_lattice();
  } else if ((_lmp->domain->nonperiodic == 0) &&
             (_lmp->domain->dimension == 3) &&
             (_lmp->domain->triclinic == 0)) {
    // Orthogonal unit cell
    boundaries_type = boundaries_pbc_ortho;
    colvarproxy_system::update_pbc_lattice();
    // It is safer to let LAMMPS deal with high-tilt triclinic boxes
  } else {
    boundaries_type = boundaries_unsupported;
  }

  if (cvm::debug()) {
    cvm::log(std::string(cvm::line_marker)+
             "colvarproxy_lammps, step no. "+cvm::to_str(colvars->it)+"\n"+
             "Updating internal data.\n");
  }

  // zero the forces on the atoms, so that they can be accumulated by the colvars
  for (size_t i = 0; i < atoms_new_colvar_forces.size(); i++) {
    atoms_new_colvar_forces[i].reset();
  }

  bias_energy = 0.0;

  if (cvm::debug()) {
    cvm::log("atoms_ids = "+cvm::to_str(atoms_ids)+"\n");
    cvm::log("atoms_ncopies = "+cvm::to_str(atoms_ncopies)+"\n");
    cvm::log("atoms_positions = "+cvm::to_str(atoms_positions)+"\n");
    cvm::log("atoms_new_colvar_forces = "+cvm::to_str(atoms_new_colvar_forces)+"\n");
  }

  // call the collective variable module
  colvars->calc();

  if (cvm::debug()) {
    cvm::log("atoms_ids = "+cvm::to_str(atoms_ids)+"\n");
    cvm::log("atoms_ncopies = "+cvm::to_str(atoms_ncopies)+"\n");
    cvm::log("atoms_positions = "+cvm::to_str(atoms_positions)+"\n");
    cvm::log("atoms_new_colvar_forces = "+cvm::to_str(atoms_new_colvar_forces)+"\n");
  }

  return bias_energy;
}

void colvarproxy_lammps::serialize_status(std::string &rst)
{
  std::ostringstream os;
  colvars->write_restart(os);
  rst = os.str();
}

void colvarproxy_lammps::write_output_files()
{
  // TODO skip output if undefined
  colvars->write_restart_file(cvm::output_prefix()+".colvars.state");
  colvars->write_output_files();
}

// set status from string
bool colvarproxy_lammps::deserialize_status(std::string &rst)
{
  std::istringstream is;
  is.str(rst);

  if (!colvars->read_restart(is)) {
    return false;
  } else {
    return true;
  }
}


cvm::rvector colvarproxy_lammps::position_distance(cvm::atom_pos const &pos1,
                                                   cvm::atom_pos const &pos2)
  const
{
  double xtmp = pos2.x - pos1.x;
  double ytmp = pos2.y - pos1.y;
  double ztmp = pos2.z - pos1.z;
  _lmp->domain->minimum_image(xtmp,ytmp,ztmp);
  return cvm::rvector(xtmp, ytmp, ztmp);
}


void colvarproxy_lammps::log(std::string const &message)
{
  std::istringstream is(message);
  std::string line;
  while (std::getline(is, line)) {
    if (_lmp->screen)
      fprintf(_lmp->screen,"colvars: %s\n",line.c_str());
    if (_lmp->logfile)
      fprintf(_lmp->logfile,"colvars: %s\n",line.c_str());
  }
}


void colvarproxy_lammps::error(std::string const &message)
{
  // In LAMMPS, all errors are fatal
  fatal_error(message);
}


void colvarproxy_lammps::fatal_error(std::string const &message)
{
  log(message);
  _lmp->error->one(FLERR,
                   "Fatal error in the collective variables module.\n");
}


int colvarproxy_lammps::backup_file(char const *filename)
{
  if (std::string(filename).rfind(std::string(".colvars.state"))
      != std::string::npos) {
    return my_backup_file(filename, ".old");
  } else {
    return my_backup_file(filename, ".BAK");
  }
}


// multi-replica support

void colvarproxy_lammps::replica_comm_barrier()
{
  MPI_Barrier(inter_comm);
}


int colvarproxy_lammps::replica_comm_recv(char* msg_data,
                                          int buf_len, int src_rep)
{
  MPI_Status status;
  int retval;

  retval = MPI_Recv(msg_data,buf_len,MPI_CHAR,src_rep,0,inter_comm,&status);
  if (retval == MPI_SUCCESS) {
    MPI_Get_count(&status, MPI_CHAR, &retval);
  } else retval = 0;
  return retval;
}


int colvarproxy_lammps::replica_comm_send(char* msg_data,
                                          int msg_len, int dest_rep)
{
  int retval;
  retval = MPI_Send(msg_data,msg_len,MPI_CHAR,dest_rep,0,inter_comm);
  if (retval == MPI_SUCCESS) {
    retval = msg_len;
  } else retval = 0;
  return retval;
}



int colvarproxy_lammps::check_atom_id(int atom_number)
{
  int const aid = atom_number;

  if (cvm::debug())
    log("Adding atom "+cvm::to_str(atom_number)+
        " for collective variables calculation.\n");

  // TODO add upper boundary check?
  if ( (aid < 0) ) {
    cvm::error("Error: invalid atom number specified, "+
               cvm::to_str(atom_number)+"\n", INPUT_ERROR);
    return INPUT_ERROR;
  }

  return aid;
}


int colvarproxy_lammps::init_atom(int atom_number)
{
  int aid = atom_number;

  for (size_t i = 0; i < atoms_ids.size(); i++) {
    if (atoms_ids[i] == aid) {
      // this atom id was already recorded
      atoms_ncopies[i] += 1;
      return i;
    }
  }

  aid = check_atom_id(atom_number);
  if (aid < 0) {
    return aid;
  }

  int const index = colvarproxy::add_atom_slot(aid);
  // add entries for the LAMMPS-specific fields
  atoms_types.push_back(0);

  return index;
}

