// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

// Using access() to check if a file exists (until we can assume C++14/17)
#if !defined(WIN32) || defined(__CYGWIN__)
#include <unistd.h>
#endif
#if defined(WIN32)
#include <io.h>
#endif

#include <cerrno>

#include <sstream>
#include <cstring>
#include <cstdio>

#if defined(_OPENMP)
#include <omp.h>
#endif

#include "colvarmodule.h"
#include "colvarproxy.h"
#include "colvarscript.h"
#include "colvaratoms.h"
#include "colvarmodule_utils.h"



colvarproxy_system::colvarproxy_system()
{
  angstrom_value = 0.0;
  kcal_mol_value = 0.0;
  boundaries_type = boundaries_unsupported;
  total_force_requested = false;
  indirect_lambda_biasing_force = 0.0;
  cached_alch_lambda_changed = false;
  cached_alch_lambda = -1.0;
  reset_pbc_lattice();
}


colvarproxy_system::~colvarproxy_system() {}


int colvarproxy_system::set_unit_system(std::string const & /* units */,
                                        bool /* check_only */)
{
  return COLVARS_NOT_IMPLEMENTED;
}


cvm::real colvarproxy_system::backend_angstrom_value()
{
  return 1.0;
}


cvm::real colvarproxy_system::boltzmann()
{
  return 0.001987191;
}


cvm::real colvarproxy_system::temperature()
{
  // TODO define, document and implement a user method to set the value of this
  return 300.0;
}


cvm::real colvarproxy_system::dt()
{
  // TODO define, document and implement a user method to set the value of this
  return 1.0;
}


cvm::real colvarproxy_system::rand_gaussian()
{
  // TODO define, document and implement a user method to set the value of this
  return 0.0;
}


void colvarproxy_system::add_energy(cvm::real /* energy */) {}


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
  return int(cvm::floor(x+0.5));
}


void colvarproxy_system::update_pbc_lattice()
{
  // Periodicity is assumed in all directions

  if (boundaries_type == boundaries_unsupported ||
      boundaries_type == boundaries_non_periodic) {
    cvm::error("Error: setting PBC lattice with unsupported boundaries.\n",
               COLVARS_BUG_ERROR);
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
    cvm::error("Error: unsupported boundary conditions.\n", COLVARS_INPUT_ERROR);
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


int colvarproxy_system::get_molid(int &)
{
  cvm::error("Error: only VMD allows the use of multiple \"molecules\", "
             "i.e. multiple molecular systems.", COLVARS_NOT_IMPLEMENTED);
  return -1;
}


int colvarproxy_system::get_alch_lambda(cvm::real * /* lambda */)
{
  return cvm::error("Error in get_alch_lambda: alchemical lambda dynamics is not supported by this build.",
    COLVARS_NOT_IMPLEMENTED);
}


void colvarproxy_system::set_alch_lambda(cvm::real lambda)
{
  cached_alch_lambda = lambda;
  cached_alch_lambda_changed = true;
}


int colvarproxy_system::send_alch_lambda()
{
  return cvm::error("Error in set_alch_lambda: alchemical lambda dynamics is not supported by this build.",
    COLVARS_NOT_IMPLEMENTED);
}


int colvarproxy_system::get_dE_dlambda(cvm::real * /* force */)
{
  return cvm::error("Error in get_dE_dlambda: alchemical lambda dynamics is not supported by this build.",
    COLVARS_NOT_IMPLEMENTED);
}


int colvarproxy_system::apply_force_dE_dlambda(cvm::real* /* force */)
{
  return cvm::error("Error in apply_force_dE_dlambda: function is not implemented by this build.",
    COLVARS_NOT_IMPLEMENTED);
}


int colvarproxy_system::get_d2E_dlambda2(cvm::real*)
{
  return cvm::error("Error in get_d2E_dlambda2: function is not implemented by this build.",
    COLVARS_NOT_IMPLEMENTED);
}


colvarproxy_atoms::colvarproxy_atoms()
{
  atoms_rms_applied_force_ = atoms_max_applied_force_ = 0.0;
  atoms_max_applied_force_id_ = -1;
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


int colvarproxy_atoms::init_atom(int /* atom_number */)
{
  return COLVARS_NOT_IMPLEMENTED;
}


int colvarproxy_atoms::check_atom_id(int /* atom_number */)
{
  return COLVARS_NOT_IMPLEMENTED;
}


int colvarproxy_atoms::init_atom(cvm::residue_id const & /* residue */,
                                 std::string const     & /* atom_name */,
                                 std::string const     & /* segment_id */)
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
               COLVARS_INPUT_ERROR);
  }
  if (atoms_ncopies[index] > 0) {
    atoms_ncopies[index] -= 1;
  }
}


int colvarproxy_atoms::load_atoms(char const * /* filename */,
                                  cvm::atom_group & /* atoms */,
                                  std::string const & /* pdb_field */,
                                  double)
{
  return cvm::error("Error: loading atom identifiers from a file "
                    "is currently not implemented.\n",
                    COLVARS_NOT_IMPLEMENTED);
}


int colvarproxy_atoms::load_coords(char const * /* filename */,
                                   std::vector<cvm::atom_pos> & /* pos */,
                                   std::vector<int> const & /* sorted_ids */,
                                   std::string const & /* pdb_field */,
                                   double)
{
  return cvm::error("Error: loading atomic coordinates from a file "
                    "is currently not implemented.\n",
                    COLVARS_NOT_IMPLEMENTED);
}


void colvarproxy_atoms::compute_rms_atoms_applied_force()
{
  atoms_rms_applied_force_ =
    compute_norm2_stats<cvm::rvector, 0, false>(atoms_new_colvar_forces);
}


void colvarproxy_atoms::compute_max_atoms_applied_force()
{
  int minmax_index = -1;
  size_t const n_atoms_ids = atoms_ids.size();
  if ((n_atoms_ids > 0) && (n_atoms_ids == atoms_new_colvar_forces.size())) {
    atoms_max_applied_force_ =
      compute_norm2_stats<cvm::rvector, 1, true>(atoms_new_colvar_forces,
                                                 &minmax_index);
    if (minmax_index >= 0) {
      atoms_max_applied_force_id_ = atoms_ids[minmax_index];
    } else {
      atoms_max_applied_force_id_ = -1;
    }
  } else {
    atoms_max_applied_force_ =
      compute_norm2_stats<cvm::rvector, 1, false>(atoms_new_colvar_forces);
    atoms_max_applied_force_id_ = -1;
  }
}



colvarproxy_atom_groups::colvarproxy_atom_groups()
{
  atom_groups_rms_applied_force_ = atom_groups_max_applied_force_ = 0.0;
}


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


int colvarproxy_atom_groups::init_atom_group(std::vector<int> const & /* atoms_ids */)
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
               COLVARS_INPUT_ERROR);
  }
  if (atom_groups_ncopies[index] > 0) {
    atom_groups_ncopies[index] -= 1;
  }
}


void colvarproxy_atom_groups::compute_rms_atom_groups_applied_force()
{
  atom_groups_rms_applied_force_ =
    compute_norm2_stats<cvm::rvector, 0, false>(atom_groups_new_colvar_forces);
}


void colvarproxy_atom_groups::compute_max_atom_groups_applied_force()
{
  atom_groups_max_applied_force_ =
    compute_norm2_stats<cvm::rvector, 1, false>(atom_groups_new_colvar_forces);
}



colvarproxy_smp::colvarproxy_smp()
{
  b_smp_active = true; // May be disabled by user option
  omp_lock_state = NULL;
#if defined(_OPENMP)
  if (omp_get_thread_num() == 0) {
    omp_lock_state = reinterpret_cast<void *>(new omp_lock_t);
    omp_init_lock(reinterpret_cast<omp_lock_t *>(omp_lock_state));
  }
#endif
}


colvarproxy_smp::~colvarproxy_smp()
{
#if defined(_OPENMP)
  if (omp_get_thread_num() == 0) {
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



colvarproxy_script::colvarproxy_script()
{
  script = NULL;
  have_scripts = false;
}


colvarproxy_script::~colvarproxy_script()
{
  if (script != NULL) {
    delete script;
    script = NULL;
  }
}


int colvarproxy_script::run_force_callback()
{
  return COLVARS_NOT_IMPLEMENTED;
}


int colvarproxy_script::run_colvar_callback(std::string const & /* name */,
                                            std::vector<const colvarvalue *> const & /* cvcs */,
                                            colvarvalue & /* value */)
{
  return COLVARS_NOT_IMPLEMENTED;
}


int colvarproxy_script::run_colvar_gradient_callback(std::string const & /* name */,
                                                     std::vector<const colvarvalue *> const & /* cvcs */,
                                                     std::vector<cvm::matrix2d<cvm::real> > & /* gradient */)
{
  return COLVARS_NOT_IMPLEMENTED;
}



colvarproxy_io::colvarproxy_io()
{
  input_buffer_ = NULL;
  restart_frequency_engine = 0;
}


colvarproxy_io::~colvarproxy_io() {}


int colvarproxy_io::get_frame(long int&)
{
  return COLVARS_NOT_IMPLEMENTED;
}


int colvarproxy_io::set_frame(long int)
{
  return COLVARS_NOT_IMPLEMENTED;
}


int colvarproxy_io::backup_file(char const *filename)
{
  // Simplified version of NAMD_file_exists()
  int exit_code;
  do {
#if defined(WIN32) && !defined(__CYGWIN__)
    // We could use _access_s here, but it is probably too new
    exit_code = _access(filename, 00);
#else
    exit_code = access(filename, F_OK);
#endif
  } while ((exit_code != 0) && (errno == EINTR));
  if (exit_code != 0) {
    if (errno == ENOENT) {
      // File does not exist
      return COLVARS_OK;
    } else {
      return cvm::error("Unknown error while checking if file \""+
                        std::string(filename)+"\" exists.\n", COLVARS_ERROR);
    }
  }

  // The file exists, then rename it
  if (std::string(filename).rfind(std::string(".colvars.state")) !=
      std::string::npos) {
    return rename_file(filename, (std::string(filename)+".old").c_str());
  } else {
    return rename_file(filename, (std::string(filename)+".BAK").c_str());
  }
}


int colvarproxy_io::remove_file(char const *filename)
{
  int error_code = COLVARS_OK;
#if defined(WIN32) && !defined(__CYGWIN__)
  // Because the file may be open by other processes, rename it to filename.old
  std::string const renamed_file(std::string(filename)+".old");
  // It may still be there from an interrupted run, so remove it to be safe
  std::remove(renamed_file.c_str());
  int rename_exit_code = 0;
  while ((rename_exit_code = std::rename(filename,
                                         renamed_file.c_str())) != 0) {
    if (errno == EINTR) continue;
    error_code |= COLVARS_FILE_ERROR;
    break;
  }
  // Ask to remove filename.old, but ignore any errors raised
  std::remove(renamed_file.c_str());
#else
  if (std::remove(filename)) {
    if (errno != ENOENT) {
      error_code |= COLVARS_FILE_ERROR;
    }
  }
#endif
  if (error_code != COLVARS_OK) {
    return cvm::error("Error: in removing file \""+std::string(filename)+
                      "\".\n.",
                      error_code);
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
    error_code |= COLVARS_FILE_ERROR;
    if (errno == EXDEV) continue;
    break;
  }
  return rename_exit_code ? error_code : COLVARS_OK;
}



colvarproxy::colvarproxy()
{
  colvars = NULL;
  b_simulation_running = true;
  b_simulation_continuing = false;
  b_delete_requested = false;
  version_int = -1;
  features_hash = 0;
}


colvarproxy::~colvarproxy()
{
  close_files();
  if (colvars != NULL) {
    delete colvars;
    colvars = NULL;
  }
}


int colvarproxy::close_files()
{
  if (smp_enabled() == COLVARS_OK && smp_thread_id() > 0) {
    // Nothing to do on non-master threads
    return COLVARS_OK;
  }
  std::list<std::string>::iterator    osni = output_stream_names.begin();
  std::list<std::ostream *>::iterator osi  = output_files.begin();
  for ( ; osi != output_files.end(); osi++, osni++) {
    ((std::ofstream *) (*osi))->close();
    delete *osi;
  }
  output_files.clear();
  output_stream_names.clear();
  return COLVARS_OK;
}


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


int colvarproxy::end_of_step()
{
  // Disable flags that Colvars doesn't need any more
  updated_masses_ = updated_charges_ = false;

  // Compute force statistics
  compute_rms_atoms_applied_force();
  compute_max_atoms_applied_force();
  compute_rms_atom_groups_applied_force();
  compute_max_atom_groups_applied_force();
  compute_rms_volmaps_applied_force();
  compute_max_volmaps_applied_force();

  if (cached_alch_lambda_changed) {
    send_alch_lambda();
    cached_alch_lambda_changed = false;
  }
  return COLVARS_OK;
}


int colvarproxy::post_run()
{
  int error_code = COLVARS_OK;
  if (colvars->output_prefix().size()) {
    error_code |= colvars->write_restart_file(cvm::output_prefix()+".colvars.state");
    error_code |= colvars->write_output_files();
  }
  error_code |= flush_output_streams();
  return error_code;
}


void colvarproxy::print_input_atomic_data()
{
  cvm::log(cvm::line_marker);

  cvm::log("Step "+cvm::to_str(cvm::step_absolute())+", "+
           "atoms_ids = "+cvm::to_str(atoms_ids)+"\n");
  cvm::log("Step "+cvm::to_str(cvm::step_absolute())+", "+
           "atoms_ncopies = "+cvm::to_str(atoms_ncopies)+"\n");
  cvm::log("Step "+cvm::to_str(cvm::step_absolute())+", "+
           "atoms_masses = "+cvm::to_str(atoms_masses)+"\n");
  cvm::log("Step "+cvm::to_str(cvm::step_absolute())+", "+
           "atoms_charges = "+cvm::to_str(atoms_charges)+"\n");
  cvm::log("Step "+cvm::to_str(cvm::step_absolute())+", "+
           "atoms_positions = "+cvm::to_str(atoms_positions,
                                            cvm::cv_width,
                                            cvm::cv_prec)+"\n");
  cvm::log("Step "+cvm::to_str(cvm::step_absolute())+", "+
           "atoms_total_forces = "+cvm::to_str(atoms_total_forces,
                                               cvm::cv_width,
                                               cvm::cv_prec)+"\n");

  cvm::log(cvm::line_marker);

  cvm::log("Step "+cvm::to_str(cvm::step_absolute())+", "+
           "atom_groups_ids = "+cvm::to_str(atom_groups_ids)+"\n");
  cvm::log("Step "+cvm::to_str(cvm::step_absolute())+", "+
           "atom_groups_ncopies = "+cvm::to_str(atom_groups_ncopies)+"\n");
  cvm::log("Step "+cvm::to_str(cvm::step_absolute())+", "+
           "atom_groups_masses = "+cvm::to_str(atom_groups_masses)+"\n");
  cvm::log("Step "+cvm::to_str(cvm::step_absolute())+", "+
           "atom_groups_charges = "+cvm::to_str(atom_groups_charges)+"\n");
  cvm::log("Step "+cvm::to_str(cvm::step_absolute())+", "+
           "atom_groups_coms = "+cvm::to_str(atom_groups_coms,
                                             cvm::cv_width,
                                             cvm::cv_prec)+"\n");
  cvm::log("Step "+cvm::to_str(cvm::step_absolute())+", "+
           "atom_groups_total_forces = "+cvm::to_str(atom_groups_total_forces,
                                                     cvm::cv_width,
                                                     cvm::cv_prec)+"\n");

  cvm::log(cvm::line_marker);

  cvm::log("Step "+cvm::to_str(cvm::step_absolute())+", "+
           "volmaps_ids = "+cvm::to_str(volmaps_ids)+"\n");
  cvm::log("Step "+cvm::to_str(cvm::step_absolute())+", "+
           "volmaps_values = "+cvm::to_str(volmaps_values)+"\n");

  cvm::log(cvm::line_marker);
}


void colvarproxy::print_output_atomic_data()
{
  cvm::log(cvm::line_marker);
  cvm::log("Step "+cvm::to_str(cvm::step_absolute())+", "+
           "atoms_new_colvar_forces = "+cvm::to_str(atoms_new_colvar_forces,
                                                    colvarmodule::cv_width,
                                                    colvarmodule::cv_prec)+"\n");
  cvm::log(cvm::line_marker);

  cvm::log("Step "+cvm::to_str(cvm::step_absolute())+", "+
           "atom_groups_new_colvar_forces = "+
           cvm::to_str(atom_groups_new_colvar_forces,
                       colvarmodule::cv_width,
                       colvarmodule::cv_prec)+"\n");

  cvm::log(cvm::line_marker);

  cvm::log("Step "+cvm::to_str(cvm::step_absolute())+", "+
           "volmaps_new_colvar_forces = "+
           cvm::to_str(volmaps_new_colvar_forces)+"\n");

  cvm::log(cvm::line_marker);
}


void colvarproxy::log(std::string const &message)
{
  fprintf(stdout, "colvars: %s", message.c_str());
}


void colvarproxy::error(std::string const &message)
{
  // TODO handle errors?
  colvarproxy::log(message);
}


void colvarproxy::add_error_msg(std::string const &message)
{
  std::istringstream is(message);
  std::string line;
  while (std::getline(is, line)) {
    error_output += line+"\n";
  }
}


void colvarproxy::clear_error_msgs()
{
  error_output.clear();
}


std::string const & colvarproxy::get_error_msgs()
{
  return error_output;
}


int colvarproxy::get_version_from_string(char const *version_string)
{
  std::string const v(version_string);
  std::istringstream is(v.substr(0, 4) + v.substr(5, 2) + v.substr(8, 2));
  int newint;
  is >> newint;
  return newint;
}


void colvarproxy::smp_stream_error()
{
  cvm::error("Error: trying to access an output stream from a "
             "multi-threaded region (bug).  For a quick workaround, use "
             "\"smp off\" in the Colvars config.\n", COLVARS_BUG_ERROR);
}


std::ostream * colvarproxy::output_stream(std::string const &output_name,
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
               COLVARS_FILE_ERROR);
    return NULL;
  }
  output_stream_names.push_back(output_name);
  output_files.push_back(osf);
  return osf;
}


std::ostream *colvarproxy::get_output_stream(std::string const &output_name)
{
  if (smp_enabled() == COLVARS_OK) {
    if (smp_thread_id() > 0) smp_stream_error();
  }
  std::list<std::ostream *>::iterator osi  = output_files.begin();
  std::list<std::string>::iterator    osni = output_stream_names.begin();
  for ( ; osi != output_files.end(); osi++, osni++) {
    if (*osni == output_name) {
      return *osi;
    }
  }
  return NULL;
}



int colvarproxy::flush_output_stream(std::ostream *os)
{
  if (smp_enabled() == COLVARS_OK) {
    if (smp_thread_id() > 0) smp_stream_error();
  }
  std::list<std::ostream *>::iterator osi  = output_files.begin();
  std::list<std::string>::iterator    osni = output_stream_names.begin();
  for ( ; osi != output_files.end(); osi++, osni++) {
    if (*osi == os) {
      ((std::ofstream *) (*osi))->flush();
      return COLVARS_OK;
    }
  }
  return cvm::error("Error: trying to flush an output file/channel "
                    "that wasn't open.\n", COLVARS_BUG_ERROR);
}


int colvarproxy::flush_output_streams()
{
  if (smp_enabled() == COLVARS_OK && smp_thread_id() > 0)
    return COLVARS_OK;

  std::list<std::ostream *>::iterator osi  = output_files.begin();
  for ( ; osi != output_files.end(); osi++) {
    ((std::ofstream *) (*osi))->flush();
  }
  return COLVARS_OK;
}


int colvarproxy::close_output_stream(std::string const &output_name)
{
  if (smp_enabled() == COLVARS_OK) {
    if (smp_thread_id() > 0) smp_stream_error();
  }
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
                    "that wasn't open.\n", COLVARS_BUG_ERROR);
}
