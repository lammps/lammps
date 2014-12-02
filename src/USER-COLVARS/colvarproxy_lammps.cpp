
#include "mpi.h"
#include "lammps.h"
#include "atom.h"
#include "error.h"
#include "output.h"
#include "random_park.h"

#include "fix_colvars.h"

#include "colvarmodule.h"
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
  system_force_requested=false;
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
}


void colvarproxy_lammps::init(const char *conf_file)
{
  // create the colvarmodule instance
  colvars = new colvarmodule(this);

  cvm::log("Using LAMMPS interface, version "+
            cvm::to_str(COLVARPROXY_VERSION)+".\n");

  // TODO move one or more of these to setup() if needed
  colvars->config_file(conf_file);
  colvars->setup_input();
  colvars->setup_output();

  if (_lmp->update->ntimestep != 0) {
    cvm::log("Initializing step number as firstTimestep.\n");
    colvars->it = colvars->it_restart = _lmp->update->ntimestep;
  }

  if (cvm::debug()) {
    log("colvars_atoms = "+cvm::to_str(colvars_atoms)+"\n");
    log("colvars_atoms_ncopies = "+cvm::to_str(colvars_atoms_ncopies)+"\n");
    log("positions = "+cvm::to_str(positions)+"\n");
    log("applied_forces = "+cvm::to_str(applied_forces)+"\n");
    log(cvm::line_marker);
    log("Info: done initializing the colvars proxy object.\n");
  }
}

colvarproxy_lammps::~colvarproxy_lammps()
{
  delete _random;
  if (colvars != NULL) {
    colvars->write_output_files();
    delete colvars;
    colvars = NULL;
  }
}

// re-initialize data where needed
void colvarproxy_lammps::setup()
{
  colvars->setup();
}

// trigger colvars computation
double colvarproxy_lammps::compute()
{
  if (first_timestep) {
    first_timestep = false;
  } else {
    // Use the time step number inherited from LAMMPS
    if ( _lmp->update->ntimestep - previous_step == 1 )
      colvars->it++;
    // Other cases could mean:
    // - run 0
    // - beginning of a new run statement
    // then the internal counter should not be incremented
  }
  previous_step = _lmp->update->ntimestep;

  if (cvm::debug()) {
    cvm::log(cvm::line_marker+
             "colvarproxy_lammps, step no. "+cvm::to_str(colvars->it)+"\n"+
             "Updating internal data.\n");
  }

  // backup applied forces if necessary to calculate system forces
  if (system_force_requested)
    previous_applied_forces = applied_forces;

  // zero the forces on the atoms, so that they can be accumulated by the colvars
  for (size_t i = 0; i < applied_forces.size(); i++) {
    applied_forces[i].x = applied_forces[i].y = applied_forces[i].z = 0.0;
  }

  bias_energy = 0.0;

  // call the collective variable module
  colvars->calc();

#if 0
  for (int i=0; i < colvars_atoms.size(); ++i) {
    fprintf(stderr,"CV: atom %d/%d/%d pos: %g %g %g  for: %g %g %g\n",
            colvars_atoms[i], colvars_atoms_ncopies[i],
            positions[i].type, positions[i].x, positions[i].y, positions[i].z,
            applied_forces[i].x, applied_forces[i].y, applied_forces[i].z);
  }
#endif

  return bias_energy;
}

void colvarproxy_lammps::serialize_status(std::string &rst)
{
  std::ostringstream os;
  colvars->write_restart(os);
  rst = os.str();
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
{
  double xtmp = pos2.x - pos1.x;
  double ytmp = pos2.y - pos1.y;
  double ztmp = pos2.z - pos1.z;
  _lmp->domain->minimum_image(xtmp,ytmp,ztmp);
  return cvm::rvector(xtmp, ytmp, ztmp);
}

cvm::real colvarproxy_lammps::position_dist2(cvm::atom_pos const &pos1,
                                             cvm::atom_pos const &pos2)
{
  double xtmp = pos2.x - pos1.x;
  double ytmp = pos2.y - pos1.y;
  double ztmp = pos2.z - pos1.z;
  _lmp->domain->minimum_image(xtmp,ytmp,ztmp);
  return cvm::real(xtmp*xtmp + ytmp*ytmp + ztmp*ztmp);
}


inline void colvarproxy_lammps::select_closest_image(cvm::atom_pos &pos,
                                                     cvm::atom_pos const &ref)
{
  double xtmp = pos.x - ref.x;
  double ytmp = pos.y - ref.y;
  double ztmp = pos.z - ref.z;
  _lmp->domain->minimum_image(xtmp,ytmp,ztmp);
  pos.x = ref.x + xtmp;
  pos.y = ref.y + ytmp;
  pos.z = ref.z + ztmp;
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
  if (!cvm::debug())
    log("If this error message is unclear, try recompiling the "
         "colvars library and LAMMPS with -DCOLVARS_DEBUG.\n");

  _lmp->error->one(FLERR,
                   "Fatal error in the collective variables module.\n");
}

void colvarproxy_lammps::exit(std::string const &message)
{
  log(message);
  log("Request to exit the simulation made.\n");
  do_exit=true;
}

enum e_pdb_field {
  e_pdb_none,
  e_pdb_occ,
  e_pdb_beta,
  e_pdb_x,
  e_pdb_y,
  e_pdb_z,
  e_pdb_ntot
};


e_pdb_field pdb_field_str2enum(std::string const &pdb_field_str)
{
  e_pdb_field pdb_field = e_pdb_none;

  if (colvarparse::to_lower_cppstr(pdb_field_str) ==
      colvarparse::to_lower_cppstr("O")) {
    pdb_field = e_pdb_occ;
  }

  if (colvarparse::to_lower_cppstr(pdb_field_str) ==
      colvarparse::to_lower_cppstr("B")) {
    pdb_field = e_pdb_beta;
  }

  if (colvarparse::to_lower_cppstr(pdb_field_str) ==
      colvarparse::to_lower_cppstr("X")) {
    pdb_field = e_pdb_x;
  }

  if (colvarparse::to_lower_cppstr(pdb_field_str) ==
      colvarparse::to_lower_cppstr("Y")) {
    pdb_field = e_pdb_y;
  }

  if (colvarparse::to_lower_cppstr(pdb_field_str) ==
      colvarparse::to_lower_cppstr("Z")) {
    pdb_field = e_pdb_z;
  }

  if (pdb_field == e_pdb_none) {
    cvm::fatal_error("Error: unsupported PDB field, \""+
                      pdb_field_str+"\".\n");
  }

  return pdb_field;
}

int colvarproxy_lammps::load_coords(char const *pdb_filename,
                                    std::vector<cvm::atom_pos> &pos,
                                    const std::vector<int> &indices,
                                    std::string const &pdb_field_str,
                                    double const pdb_field_value)
{
  cvm::fatal_error("Reading collective variable coordinates "
                   "from a PDB file is currently not supported.\n");
  return COLVARS_ERROR;
}

int colvarproxy_lammps::load_atoms(char const *pdb_filename,
                                   std::vector<cvm::atom> &atoms,
                                   std::string const &pdb_field_str,
                                   double const pdb_field_value)
{
  cvm::fatal_error("Selecting collective variable atoms "
                    "from a PDB file is currently not supported.\n");
  return COLVARS_ERROR;
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

int colvarproxy_lammps::init_lammps_atom(const int &aid, cvm::atom *atom)
{
  atom->id = aid;
  atom->mass = 0.0;

  for (size_t i = 0; i < colvars_atoms.size(); i++) {
    if (colvars_atoms[i] == aid) {
      // this atom id was already recorded
      colvars_atoms_ncopies[i] += 1;
      return i;
    }
  }

  // allocate a new slot for this atom
  colvars_atoms_ncopies.push_back(1);
  colvars_atoms.push_back(aid);
  struct commdata c;
  c.tag = aid;
  c.type = 0;
  c.x = c.y = c.z = 0.0;
  positions.push_back(c);
  total_forces.push_back(c);
  applied_forces.push_back(c);

  return colvars_atoms.size()-1;
}

// multi-replica support

void colvarproxy_lammps::replica_comm_barrier() {
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

// atom member functions, LAMMPS specific implementations

cvm::atom::atom(const int &id)
{

  if (cvm::debug())
    cvm::log("Adding atom "+cvm::to_str(id)+
             " for collective variables calculation.\n");

  if (id < 0)
    cvm::fatal_error("Error: invalid atom ID specified, "+
                     cvm::to_str(id)+"\n");

  int idx = ((colvarproxy_lammps *) cvm::proxy)->init_lammps_atom(id,this);
  if (idx < 0)
    cvm::fatal_error("Error: atom ID , "+cvm::to_str(id)+" does not exist.\n");

  this->index = idx;
  if (cvm::debug())
    cvm::log("The index of this atom in the colvarproxy_lammps arrays is "+
             cvm::to_str(this->index)+".\n");

  this->reset_data();
}


/// For AMBER topologies, the segment id is automatically set to
/// "MAIN" (the segment id assigned by NAMD's AMBER topology parser),
/// and is therefore optional when an AMBER topology is used
cvm::atom::atom(cvm::residue_id const &residue,
                std::string const     &atom_name,
                std::string const     &segment_id)
{
  cvm::fatal_error("Creating collective variable atoms "
                   "from a PDB file is currently not supported.\n");
}


// copy constructor
cvm::atom::atom(cvm::atom const &a)
  : index(a.index), id(a.id), mass(a.mass)
{
  // init_lammps_atom() has already been called by a's constructor, no
  // need to call it again

  // need to increment the counter anyway
  colvarproxy_lammps *cp = (colvarproxy_lammps *) cvm::proxy;
  cp->colvars_atoms_ncopies[this->index] += 1;
}


cvm::atom::~atom()
{
  if (this->index >= 0) {
    colvarproxy_lammps *cp = (colvarproxy_lammps *) cvm::proxy;
    if (cp->colvars_atoms_ncopies[this->index] > 0)
      cp->colvars_atoms_ncopies[this->index] -= 1;
  }
}

void cvm::atom::read_position()
{
  colvarproxy_lammps const * const cp = (colvarproxy_lammps *) cvm::proxy;
  this->pos.x = cp->positions[this->index].x;
  this->pos.y = cp->positions[this->index].y;
  this->pos.z = cp->positions[this->index].z;
  this->mass = cp->positions[this->index].m;
}

void cvm::atom::read_velocity()
{
  cvm::fatal_error("Error: read_velocity is not yet implemented.\n");
}

void cvm::atom::read_system_force()
{
  colvarproxy_lammps const * const cp = (colvarproxy_lammps *) cvm::proxy;
  this->system_force.x = cp->total_forces[this->index].x
    - cp->previous_applied_forces[this->index].x;
  this->system_force.y = cp->total_forces[this->index].y
    - cp->previous_applied_forces[this->index].y;
  this->system_force.z = cp->total_forces[this->index].z
    - cp->previous_applied_forces[this->index].z;
}

void cvm::atom::apply_force(cvm::rvector const &new_force)
{
  colvarproxy_lammps *cp = (colvarproxy_lammps *) cvm::proxy;
  cp->applied_forces[this->index].x += new_force.x;
  cp->applied_forces[this->index].y += new_force.y;
  cp->applied_forces[this->index].z += new_force.z;
}
