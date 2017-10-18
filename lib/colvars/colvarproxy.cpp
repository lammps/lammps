// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include <sstream>
#include <string.h>

#if defined(_OPENMP)
#include <omp.h>
#endif

#include "colvarmodule.h"
#include "colvarproxy.h"
#include "colvarscript.h"
#include "colvaratoms.h"



colvarproxy_system::colvarproxy_system() {}


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


cvm::real colvarproxy_system::position_dist2(cvm::atom_pos const &pos1,
                                             cvm::atom_pos const &pos2)
{
  return (position_distance(pos1, pos2)).norm2();
}



colvarproxy_atoms::colvarproxy_atoms() {}


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
                                  double const)
{
  return cvm::error("Error: loading atom identifiers from a file "
                    "is currently not implemented.\n",
                    COLVARS_NOT_IMPLEMENTED);
}


int colvarproxy_atoms::load_coords(char const *filename,
                                   std::vector<cvm::atom_pos> &pos,
                                   const std::vector<int> &indices,
                                   std::string const &pdb_field,
                                   double const)
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
    omp_init_lock(reinterpret_cast<omp_lock_t *>(omp_lock_state));
  }
#endif
}


colvarproxy_smp::~colvarproxy_smp() {}


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


char *colvarproxy_script::script_obj_to_str(unsigned char *obj)
{
  return reinterpret_cast<char *>(obj);
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
  std::list<std::ostream *>::iterator osi  = output_files.begin();
  std::list<std::string>::iterator    osni = output_stream_names.begin();
  for ( ; osi != output_files.end(); osi++, osni++) {
    if (*osni == output_name) {
      return *osi;
    }
  }
  if (!(mode & (std::ios_base::app | std::ios_base::ate))) {
    backup_file(output_name);
  }
  std::ofstream *os = new std::ofstream(output_name.c_str(), mode);
  if (!os->is_open()) {
    cvm::error("Error: cannot write to file/channel \""+output_name+"\".\n",
               FILE_ERROR);
    return NULL;
  }
  output_stream_names.push_back(output_name);
  output_files.push_back(os);
  return os;
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
  return COLVARS_NOT_IMPLEMENTED;
}



colvarproxy::colvarproxy()
{
  colvars = NULL;
  b_simulation_running = true;
}


colvarproxy::~colvarproxy() {}


int colvarproxy::reset()
{
  int error_code = COLVARS_OK;
  error_code |= colvarproxy_atoms::reset();
  error_code |= colvarproxy_atom_groups::reset();
  return error_code;
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











