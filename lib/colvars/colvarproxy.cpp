// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#if !defined(WIN32) || defined(__CYGWIN__)
#include <unistd.h>
#endif
#include <cerrno>

#include <sstream>
#include <cstring>
#include <cstdio>

#if defined(_OPENMP)
#include <omp.h>
#endif

#if defined(NAMD_TCL) || defined(VMDTCL)
#define COLVARS_TCL
#include <tcl.h>
#endif

#include "colvarmodule.h"
#include "colvarproxy.h"
#include "colvarscript.h"
#include "colvaratoms.h"



colvarproxy_system::colvarproxy_system()
{
  reset_pbc_lattice();
}


colvarproxy_system::~colvarproxy_system() {}


void colvarproxy_system::add_energy(cvm::real energy) {}


void colvarproxy_system::request_total_force(bool yesno)
{
  if (yesno == true)
    cvm::error("Error: total forces are currently not implemented.\n",
               COLVARS_NOT_IMPLEMENTED);
}


bool colvarproxy_system::total_forces_enabled() const
{
  return false;
}


bool colvarproxy_system::total_forces_same_step() const
{
  return false;
}


inline int round_to_integer(cvm::real x)
{
  return cvm::floor(x+0.5);
}


void colvarproxy_system::update_pbc_lattice()
{
  // Periodicity is assumed in all directions

  if (boundaries_type == boundaries_unsupported ||
      boundaries_type == boundaries_non_periodic) {
    cvm::error("Error: setting PBC lattice with unsupported boundaries.\n",
               BUG_ERROR);
    return;
  }

  {
    cvm::rvector const v = cvm::rvector::outer(unit_cell_y, unit_cell_z);
    reciprocal_cell_x = v/(v*unit_cell_x);
  }
  {
    cvm::rvector const v = cvm::rvector::outer(unit_cell_z, unit_cell_x);
    reciprocal_cell_y = v/(v*unit_cell_y);
  }
  {
    cvm::rvector const v = cvm::rvector::outer(unit_cell_x, unit_cell_y);
    reciprocal_cell_z = v/(v*unit_cell_z);
  }
}


void colvarproxy_system::reset_pbc_lattice()
{
  unit_cell_x.reset();
  unit_cell_y.reset();
  unit_cell_z.reset();
  reciprocal_cell_x.reset();
  reciprocal_cell_y.reset();
  reciprocal_cell_z.reset();
}


cvm::rvector colvarproxy_system::position_distance(cvm::atom_pos const &pos1,
                                                   cvm::atom_pos const &pos2)
  const
{
  if (boundaries_type == boundaries_unsupported) {
    cvm::error("Error: unsupported boundary conditions.\n", INPUT_ERROR);
  }

  cvm::rvector diff = (pos2 - pos1);

  if (boundaries_type == boundaries_non_periodic) return diff;

  cvm::real const x_shift = round_to_integer(reciprocal_cell_x*diff);
  cvm::real const y_shift = round_to_integer(reciprocal_cell_y*diff);
  cvm::real const z_shift = round_to_integer(reciprocal_cell_z*diff);

  diff.x -= x_shift*unit_cell_x.x + y_shift*unit_cell_y.x +
    z_shift*unit_cell_z.x;
  diff.y -= x_shift*unit_cell_x.y + y_shift*unit_cell_y.y +
    z_shift*unit_cell_z.y;
  diff.z -= x_shift*unit_cell_x.z + y_shift*unit_cell_y.z +
    z_shift*unit_cell_z.z;

  return diff;
}



colvarproxy_atoms::colvarproxy_atoms()
{
  updated_masses_ = updated_charges_ = false;
}


colvarproxy_atoms::~colvarproxy_atoms()
{
  reset();
}


int colvarproxy_atoms::reset()
{
  atoms_ids.clear();
  atoms_ncopies.clear();
  atoms_masses.clear();
  atoms_charges.clear();
  atoms_positions.clear();
  atoms_total_forces.clear();
  atoms_new_colvar_forces.clear();
  return COLVARS_OK;
}


int colvarproxy_atoms::add_atom_slot(int atom_id)
{
  atoms_ids.push_back(atom_id);
  atoms_ncopies.push_back(1);
  atoms_masses.push_back(1.0);
  atoms_charges.push_back(0.0);
  atoms_positions.push_back(cvm::rvector(0.0, 0.0, 0.0));
  atoms_total_forces.push_back(cvm::rvector(0.0, 0.0, 0.0));
  atoms_new_colvar_forces.push_back(cvm::rvector(0.0, 0.0, 0.0));
  return (atoms_ids.size() - 1);
}


int colvarproxy_atoms::init_atom(cvm::residue_id const &residue,
                                 std::string const     &atom_name,
                                 std::string const     &segment_id)
{
  cvm::error("Error: initializing an atom by name and residue number is currently not supported.\n",
             COLVARS_NOT_IMPLEMENTED);
  return COLVARS_NOT_IMPLEMENTED;
}


int colvarproxy_atoms::check_atom_id(cvm::residue_id const &residue,
                                     std::string const     &atom_name,
                                     std::string const     &segment_id)
{
  colvarproxy_atoms::init_atom(residue, atom_name, segment_id);
  return COLVARS_NOT_IMPLEMENTED;
}


void colvarproxy_atoms::clear_atom(int index)
{
  if (((size_t) index) >= atoms_ids.size()) {
    cvm::error("Error: trying to disable an atom that was not previously requested.\n",
               INPUT_ERROR);
  }
  if (atoms_ncopies[index] > 0) {
    atoms_ncopies[index] -= 1;
  }
}


int colvarproxy_atoms::load_atoms(char const *filename,
                                  cvm::atom_group &atoms,
                                  std::string const &pdb_field,
                                  double)
{
  return cvm::error("Error: loading atom identifiers from a file "
                    "is currently not implemented.\n",
                    COLVARS_NOT_IMPLEMENTED);
}


int colvarproxy_atoms::load_coords(char const *filename,
                                   std::vector<cvm::atom_pos> &pos,
                                   std::vector<int> const &sorted_ids,
                                   std::string const &pdb_field,
                                   double)
{
  return cvm::error("Error: loading atomic coordinates from a file "
                    "is currently not implemented.\n",
                    COLVARS_NOT_IMPLEMENTED);
}



colvarproxy_atom_groups::colvarproxy_atom_groups() {}


colvarproxy_atom_groups::~colvarproxy_atom_groups()
{
  reset();
}


int colvarproxy_atom_groups::reset()
{
  atom_groups_ids.clear();
  atom_groups_ncopies.clear();
  atom_groups_masses.clear();
  atom_groups_charges.clear();
  atom_groups_coms.clear();
  atom_groups_total_forces.clear();
  atom_groups_new_colvar_forces.clear();
  return COLVARS_OK;
}


int colvarproxy_atom_groups::add_atom_group_slot(int atom_group_id)
{
  atom_groups_ids.push_back(atom_group_id);
  atom_groups_ncopies.push_back(1);
  atom_groups_masses.push_back(1.0);
  atom_groups_charges.push_back(0.0);
  atom_groups_coms.push_back(cvm::rvector(0.0, 0.0, 0.0));
  atom_groups_total_forces.push_back(cvm::rvector(0.0, 0.0, 0.0));
  atom_groups_new_colvar_forces.push_back(cvm::rvector(0.0, 0.0, 0.0));
  return (atom_groups_ids.size() - 1);
}


int colvarproxy_atom_groups::scalable_group_coms()
{
  return COLVARS_NOT_IMPLEMENTED;
}


int colvarproxy_atom_groups::init_atom_group(std::vector<int> const &atoms_ids)
{
  cvm::error("Error: initializing a group outside of the Colvars module "
             "is currently not supported.\n",
             COLVARS_NOT_IMPLEMENTED);
  return COLVARS_NOT_IMPLEMENTED;
}


void colvarproxy_atom_groups::clear_atom_group(int index)
{
  if (((size_t) index) >= atom_groups_ids.size()) {
    cvm::error("Error: trying to disable an atom group "
               "that was not previously requested.\n",
               INPUT_ERROR);
  }
  if (atom_groups_ncopies[index] > 0) {
    atom_groups_ncopies[index] -= 1;
  }
}



colvarproxy_smp::colvarproxy_smp()
{
  b_smp_active = true; // May be disabled by user option
  omp_lock_state = NULL;
#if defined(_OPENMP)
  if (smp_thread_id() == 0) {
    omp_lock_state = reinterpret_cast<void *>(new omp_lock_t);
    omp_init_lock(reinterpret_cast<omp_lock_t *>(omp_lock_state));
  }
#endif
}


colvarproxy_smp::~colvarproxy_smp()
{
#if defined(_OPENMP)
  if (smp_thread_id() == 0) {
    if (omp_lock_state) {
      delete reinterpret_cast<omp_lock_t *>(omp_lock_state);
    }
  }
#endif
}


int colvarproxy_smp::smp_enabled()
{
#if defined(_OPENMP)
  if (b_smp_active) {
    return COLVARS_OK;
  }
  return COLVARS_ERROR;
#else
  return COLVARS_NOT_IMPLEMENTED;
#endif
}


int colvarproxy_smp::smp_colvars_loop()
{
#if defined(_OPENMP)
  colvarmodule *cv = cvm::main();
  colvarproxy *proxy = cv->proxy;
#pragma omp parallel for
  for (size_t i = 0; i < cv->variables_active_smp()->size(); i++) {
    colvar *x = (*(cv->variables_active_smp()))[i];
    int x_item = (*(cv->variables_active_smp_items()))[i];
    if (cvm::debug()) {
      cvm::log("["+cvm::to_str(proxy->smp_thread_id())+"/"+
               cvm::to_str(proxy->smp_num_threads())+
               "]: calc_colvars_items_smp(), i = "+cvm::to_str(i)+", cv = "+
               x->name+", cvc = "+cvm::to_str(x_item)+"\n");
    }
    x->calc_cvcs(x_item, 1);
  }
  return cvm::get_error();
#else
  return COLVARS_NOT_IMPLEMENTED;
#endif
}


int colvarproxy_smp::smp_biases_loop()
{
#if defined(_OPENMP)
  colvarmodule *cv = cvm::main();
#pragma omp parallel
  {
#pragma omp for
    for (size_t i = 0; i < cv->biases_active()->size(); i++) {
      colvarbias *b = (*(cv->biases_active()))[i];
      if (cvm::debug()) {
        cvm::log("Calculating bias \""+b->name+"\" on thread "+
                 cvm::to_str(smp_thread_id())+"\n");
      }
      b->update();
    }
  }
  return cvm::get_error();
#else
  return COLVARS_NOT_IMPLEMENTED;
#endif
}


int colvarproxy_smp::smp_biases_script_loop()
{
#if defined(_OPENMP)
  colvarmodule *cv = cvm::main();
#pragma omp parallel
  {
#pragma omp single nowait
    {
      cv->calc_scripted_forces();
    }
#pragma omp for
    for (size_t i = 0; i < cv->biases_active()->size(); i++) {
      colvarbias *b = (*(cv->biases_active()))[i];
      if (cvm::debug()) {
        cvm::log("Calculating bias \""+b->name+"\" on thread "+
                 cvm::to_str(smp_thread_id())+"\n");
      }
      b->update();
    }
  }
  return cvm::get_error();
#else
  return COLVARS_NOT_IMPLEMENTED;
#endif
}




int colvarproxy_smp::smp_thread_id()
{
#if defined(_OPENMP)
  return omp_get_thread_num();
#else
  return COLVARS_NOT_IMPLEMENTED;
#endif
}


int colvarproxy_smp::smp_num_threads()
{
#if defined(_OPENMP)
  return omp_get_max_threads();
#else
  return COLVARS_NOT_IMPLEMENTED;
#endif
}


int colvarproxy_smp::smp_lock()
{
#if defined(_OPENMP)
  omp_set_lock(reinterpret_cast<omp_lock_t *>(omp_lock_state));
#endif
  return COLVARS_OK;
}


int colvarproxy_smp::smp_trylock()
{
#if defined(_OPENMP)
  return omp_test_lock(reinterpret_cast<omp_lock_t *>(omp_lock_state)) ?
    COLVARS_OK : COLVARS_ERROR;
#else
  return COLVARS_OK;
#endif
}


int colvarproxy_smp::smp_unlock()
{
#if defined(_OPENMP)
  omp_unset_lock(reinterpret_cast<omp_lock_t *>(omp_lock_state));
#endif
  return COLVARS_OK;
}



colvarproxy_replicas::colvarproxy_replicas() {}


colvarproxy_replicas::~colvarproxy_replicas() {}


bool colvarproxy_replicas::replica_enabled()
{
  return false;
}


int colvarproxy_replicas::replica_index()
{
  return 0;
}


int colvarproxy_replicas::replica_num()
{
  return 1;
}


void colvarproxy_replicas::replica_comm_barrier() {}


int colvarproxy_replicas::replica_comm_recv(char* msg_data,
                                            int buf_len,
                                            int src_rep)
{
  return COLVARS_NOT_IMPLEMENTED;
}


int colvarproxy_replicas::replica_comm_send(char* msg_data,
                                            int msg_len,
                                            int dest_rep)
{
  return COLVARS_NOT_IMPLEMENTED;
}




colvarproxy_script::colvarproxy_script()
{
  script = NULL;
}


colvarproxy_script::~colvarproxy_script() {}


char const *colvarproxy_script::script_obj_to_str(unsigned char *obj)
{
  cvm::error("Error: trying to print a script object without a scripting "
             "language interface.\n", BUG_ERROR);
  return reinterpret_cast<char *>(obj);
}


std::vector<std::string> colvarproxy_script::script_obj_to_str_vector(unsigned char *obj)
{
  cvm::error("Error: trying to print a script object without a scripting "
             "language interface.\n", BUG_ERROR);
  return std::vector<std::string>();
}


int colvarproxy_script::run_force_callback()
{
  return COLVARS_NOT_IMPLEMENTED;
}


int colvarproxy_script::run_colvar_callback(
      std::string const &name,
      std::vector<const colvarvalue *> const &cvcs,
      colvarvalue &value)
{
  return COLVARS_NOT_IMPLEMENTED;
}


int colvarproxy_script::run_colvar_gradient_callback(
      std::string const &name,
      std::vector<const colvarvalue *> const &cvcs,
      std::vector<cvm::matrix2d<cvm::real> > &gradient)
{
  return COLVARS_NOT_IMPLEMENTED;
}



colvarproxy_tcl::colvarproxy_tcl()
{
  tcl_interp_ = NULL;
}


colvarproxy_tcl::~colvarproxy_tcl()
{
}


void colvarproxy_tcl::init_tcl_pointers()
{
  cvm::error("Error: Tcl support is currently unavailable "
             "outside NAMD or VMD.\n", COLVARS_NOT_IMPLEMENTED);
}


char const *colvarproxy_tcl::tcl_obj_to_str(unsigned char *obj)
{
#if defined(COLVARS_TCL)
  return Tcl_GetString(reinterpret_cast<Tcl_Obj *>(obj));
#else
  return NULL;
#endif
}


int colvarproxy_tcl::tcl_run_force_callback()
{
#if defined(COLVARS_TCL)
  Tcl_Interp *const tcl_interp = reinterpret_cast<Tcl_Interp *>(tcl_interp_);
  std::string cmd = std::string("calc_colvar_forces ")
    + cvm::to_str(cvm::step_absolute());
  int err = Tcl_Eval(tcl_interp, cmd.c_str());
  if (err != TCL_OK) {
    cvm::log(std::string("Error while executing calc_colvar_forces:\n"));
    cvm::error(Tcl_GetStringResult(tcl_interp));
    return COLVARS_ERROR;
  }
  return cvm::get_error();
#else
  return COLVARS_NOT_IMPLEMENTED;
#endif
}


int colvarproxy_tcl::tcl_run_colvar_callback(
                         std::string const &name,
                         std::vector<const colvarvalue *> const &cvc_values,
                         colvarvalue &value)
{
#if defined(COLVARS_TCL)

  Tcl_Interp *const tcl_interp = reinterpret_cast<Tcl_Interp *>(tcl_interp_);
  size_t i;
  std::string cmd = std::string("calc_") + name;
  for (i = 0; i < cvc_values.size(); i++) {
    cmd += std::string(" {") + (*(cvc_values[i])).to_simple_string() +
      std::string("}");
  }
  int err = Tcl_Eval(tcl_interp, cmd.c_str());
  const char *result = Tcl_GetStringResult(tcl_interp);
  if (err != TCL_OK) {
    return cvm::error(std::string("Error while executing ")
                      + cmd + std::string(":\n") +
                      std::string(Tcl_GetStringResult(tcl_interp)), COLVARS_ERROR);
  }
  std::istringstream is(result);
  if (value.from_simple_string(is.str()) != COLVARS_OK) {
    cvm::log("Error parsing colvar value from script:");
    cvm::error(result);
    return COLVARS_ERROR;
  }
  return cvm::get_error();

#else

  return COLVARS_NOT_IMPLEMENTED;

#endif
}


int colvarproxy_tcl::tcl_run_colvar_gradient_callback(
                         std::string const &name,
                         std::vector<const colvarvalue *> const &cvc_values,
                         std::vector<cvm::matrix2d<cvm::real> > &gradient)
{
#if defined(COLVARS_TCL)

  Tcl_Interp *const tcl_interp = reinterpret_cast<Tcl_Interp *>(tcl_interp_);
  size_t i;
  std::string cmd = std::string("calc_") + name + "_gradient";
  for (i = 0; i < cvc_values.size(); i++) {
    cmd += std::string(" {") + (*(cvc_values[i])).to_simple_string() +
      std::string("}");
  }
  int err = Tcl_Eval(tcl_interp, cmd.c_str());
  if (err != TCL_OK) {
    return cvm::error(std::string("Error while executing ")
                      + cmd + std::string(":\n") +
                      std::string(Tcl_GetStringResult(tcl_interp)), COLVARS_ERROR);
  }
  Tcl_Obj **list;
  int n;
  Tcl_ListObjGetElements(tcl_interp, Tcl_GetObjResult(tcl_interp),
                         &n, &list);
  if (n != int(gradient.size())) {
    cvm::error("Error parsing list of gradient values from script: found "
               + cvm::to_str(n) + " values instead of " +
               cvm::to_str(gradient.size()));
    return COLVARS_ERROR;
  }
  for (i = 0; i < gradient.size(); i++) {
    std::istringstream is(Tcl_GetString(list[i]));
    if (gradient[i].from_simple_string(is.str()) != COLVARS_OK) {
      cvm::log("Gradient matrix size: " + cvm::to_str(gradient[i].size()));
      cvm::log("Gradient string: " + cvm::to_str(Tcl_GetString(list[i])));
      cvm::error("Error parsing gradient value from script", COLVARS_ERROR);
      return COLVARS_ERROR;
    }
  }

  return cvm::get_error();

#else

  return COLVARS_NOT_IMPLEMENTED;

#endif
}



colvarproxy_io::colvarproxy_io() {}


colvarproxy_io::~colvarproxy_io() {}


int colvarproxy_io::get_frame(long int&)
{
  return COLVARS_NOT_IMPLEMENTED;
}


int colvarproxy_io::set_frame(long int)
{
  return COLVARS_NOT_IMPLEMENTED;
}


std::ostream * colvarproxy_io::output_stream(std::string const &output_name,
                                             std::ios_base::openmode mode)
{
  if (cvm::debug()) {
    cvm::log("Using colvarproxy::output_stream()\n");
  }

  std::ostream *os = get_output_stream(output_name);
  if (os != NULL) return os;

  if (!(mode & (std::ios_base::app | std::ios_base::ate))) {
    backup_file(output_name);
  }
  std::ofstream *osf = new std::ofstream(output_name.c_str(), mode);
  if (!osf->is_open()) {
    cvm::error("Error: cannot write to file/channel \""+output_name+"\".\n",
               FILE_ERROR);
    return NULL;
  }
  output_stream_names.push_back(output_name);
  output_files.push_back(osf);
  return osf;
}


std::ostream *colvarproxy_io::get_output_stream(std::string const &output_name)
{
  std::list<std::ostream *>::iterator osi  = output_files.begin();
  std::list<std::string>::iterator    osni = output_stream_names.begin();
  for ( ; osi != output_files.end(); osi++, osni++) {
    if (*osni == output_name) {
      return *osi;
    }
  }
  return NULL;
}



int colvarproxy_io::flush_output_stream(std::ostream *os)
{
  std::list<std::ostream *>::iterator osi  = output_files.begin();
  std::list<std::string>::iterator    osni = output_stream_names.begin();
  for ( ; osi != output_files.end(); osi++, osni++) {
    if (*osi == os) {
      ((std::ofstream *) (*osi))->flush();
      return COLVARS_OK;
    }
  }
  return cvm::error("Error: trying to flush an output file/channel "
                    "that wasn't open.\n", BUG_ERROR);
}


int colvarproxy_io::close_output_stream(std::string const &output_name)
{
  std::list<std::ostream *>::iterator osi  = output_files.begin();
  std::list<std::string>::iterator    osni = output_stream_names.begin();
  for ( ; osi != output_files.end(); osi++, osni++) {
    if (*osni == output_name) {
      ((std::ofstream *) (*osi))->close();
      delete *osi;
      output_files.erase(osi);
      output_stream_names.erase(osni);
      return COLVARS_OK;
    }
  }
  return cvm::error("Error: trying to close an output file/channel "
                    "that wasn't open.\n", BUG_ERROR);
}


int colvarproxy_io::backup_file(char const *filename)
{
  // TODO implement this using rename_file()
  return COLVARS_NOT_IMPLEMENTED;
}


int colvarproxy_io::remove_file(char const *filename)
{
  if (std::remove(filename)) {
    if (errno != ENOENT) {
      return cvm::error("Error: in removing file \""+std::string(filename)+
                        "\".\n.",
                        FILE_ERROR);
    }
  }
  return COLVARS_OK;
}


int colvarproxy_io::rename_file(char const *filename, char const *newfilename)
{
  int error_code = COLVARS_OK;
#if defined(WIN32) && !defined(__CYGWIN__)
  // On straight Windows, must remove the destination before renaming it
  error_code |= remove_file(newfilename);
#endif
  int rename_exit_code = 0;
  while ((rename_exit_code = std::rename(filename, newfilename)) != 0) {
    if (errno == EINTR) continue;
    // Call log() instead of error to allow the next try
    cvm::log("Error: in renaming file \""+std::string(filename)+"\" to \""+
             std::string(newfilename)+"\".\n.");
    error_code |= FILE_ERROR;
    if (errno == EXDEV) continue;
    break;
  }
  return rename_exit_code ? error_code : COLVARS_OK;
}



colvarproxy::colvarproxy()
{
  colvars = NULL;
  b_simulation_running = true;
  b_delete_requested = false;
}


colvarproxy::~colvarproxy() {}


int colvarproxy::reset()
{
  int error_code = COLVARS_OK;
  error_code |= colvarproxy_atoms::reset();
  error_code |= colvarproxy_atom_groups::reset();
  return error_code;
}


int colvarproxy::request_deletion()
{
  return cvm::error("Error: \"delete\" command is only available in VMD; "
                    "please use \"reset\" instead.\n",
                    COLVARS_NOT_IMPLEMENTED);
}


int colvarproxy::setup()
{
  return COLVARS_OK;
}


int colvarproxy::update_input()
{
  return COLVARS_OK;
}


int colvarproxy::update_output()
{
  return COLVARS_OK;
}


size_t colvarproxy::restart_frequency()
{
  return 0;
}


int colvarproxy::get_version_from_string(char const *version_string)
{
  std::string const v(version_string);
  std::istringstream is(v.substr(0, 4) + v.substr(5, 2) + v.substr(8, 2));
  int newint;
  is >> newint;
  return newint;
}
