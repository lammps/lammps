// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#ifndef COLVARPROXY_H
#define COLVARPROXY_H

#include <fstream>
#include <list>

#include "colvarmodule.h"
#include "colvartypes.h"
#include "colvarvalue.h"


/// \file colvarproxy.h
/// \brief Colvars proxy classes
///
/// This file declares the class for the object responsible for interfacing
/// Colvars with other codes (MD engines, VMD, Python).  The \link colvarproxy
/// \endlink class is a derivative of multiple classes, each devoted to a
/// specific task (e.g. \link colvarproxy_atoms \endlink to access data for
/// individual atoms).
///
/// To interface to a new MD engine, the simplest solution is to derive a new
/// class from \link colvarproxy \endlink.  Currently implemented are: \link
/// colvarproxy_lammps \endlink, \link colvarproxy_namd \endlink, \link
/// colvarproxy_vmd \endlink.


// forward declarations
class colvarscript;


/// Methods for accessing the simulation system (PBCs, integrator, etc)
class colvarproxy_system {

public:

  /// Constructor
  colvarproxy_system();

  /// Destructor
  virtual ~colvarproxy_system();

  /// \brief Value of the unit for atomic coordinates with respect to
  /// angstroms (used by some variables for hard-coded default values)
  virtual cvm::real unit_angstrom() = 0;

  /// \brief Boltzmann constant
  virtual cvm::real boltzmann() = 0;

  /// \brief Target temperature of the simulation (K units)
  virtual cvm::real temperature() = 0;

  /// \brief Time step of the simulation (fs)
  virtual cvm::real dt() = 0;

  /// \brief Pseudo-random number with Gaussian distribution
  virtual cvm::real rand_gaussian(void) = 0;

  /// Pass restraint energy value for current timestep to MD engine
  virtual void add_energy(cvm::real energy) = 0;

  /// \brief Get the PBC-aware distance vector between two positions
  virtual cvm::rvector position_distance(cvm::atom_pos const &pos1,
                                         cvm::atom_pos const &pos2) const;

  /// Recompute PBC reciprocal lattice (assumes XYZ periodicity)
  void update_pbc_lattice();

  /// Set the lattice vectors to zero
  void reset_pbc_lattice();

  /// \brief Tell the proxy whether total forces are needed (they may not
  /// always be available)
  virtual void request_total_force(bool yesno);

  /// Are total forces being used?
  virtual bool total_forces_enabled() const;

  /// Are total forces from the current step available?
  virtual bool total_forces_same_step() const;

protected:

  /// \brief Type of boundary conditions
  ///
  /// Orthogonal and triclinic cells are made available to objects.
  /// For any other conditions (mixed periodicity, triclinic cells in LAMMPS)
  /// minimum-image distances are computed by the host engine regardless.
  enum Boundaries_type {
    boundaries_non_periodic,
    boundaries_pbc_ortho,
    boundaries_pbc_triclinic,
    boundaries_unsupported
  };

  /// Type of boundary conditions
  Boundaries_type boundaries_type;

  /// Bravais lattice vectors
  cvm::rvector unit_cell_x, unit_cell_y, unit_cell_z;

  /// Reciprocal lattice vectors
  cvm::rvector reciprocal_cell_x, reciprocal_cell_y, reciprocal_cell_z;
};


/// \brief Container of atomic data for processing by Colvars
class colvarproxy_atoms {

public:

  /// Constructor
  colvarproxy_atoms();

  /// Destructor
  virtual ~colvarproxy_atoms();

  /// Prepare this atom for collective variables calculation, selecting it by
  /// numeric index (1-based)
  virtual int init_atom(int atom_number) = 0;

  /// Check that this atom number is valid, but do not initialize the
  /// corresponding atom yet
  virtual int check_atom_id(int atom_number) = 0;

  /// Select this atom for collective variables calculation, using name and
  /// residue number.  Not all programs support this: leave this function as
  /// is in those cases.
  virtual int init_atom(cvm::residue_id const &residue,
                        std::string const     &atom_name,
                        std::string const     &segment_id);

  /// Check that this atom is valid, but do not initialize it yet
  virtual int check_atom_id(cvm::residue_id const &residue,
                            std::string const     &atom_name,
                            std::string const     &segment_id);

  /// \brief Used by the atom class destructor: rather than deleting the array slot
  /// (costly) set the corresponding atoms_ncopies to zero
  virtual void clear_atom(int index);

  /// \brief Select atom IDs from a file (usually PDB) \param filename name of
  /// the file \param atoms array to which atoms read from "filename" will be
  /// appended \param pdb_field (optional) if the file is a PDB and this
  /// string is non-empty, select atoms for which this field is non-zero
  /// \param pdb_field_value (optional) if non-zero, select only atoms whose
  /// pdb_field equals this
  virtual int load_atoms(char const *filename,
                         cvm::atom_group &atoms,
                         std::string const &pdb_field,
                         double pdb_field_value = 0.0);

  /// \brief Load a set of coordinates from a file (usually PDB); if "pos" is
  /// already allocated, the number of its elements must match the number of
  /// entries in "filename" \param filename name of the file \param pos array
  /// of coordinates \param sorted_ids array of sorted internal IDs, used to
  /// loop through the file only once \param pdb_field (optional) if the file
  /// is a PDB and this string is non-empty, select atoms for which this field
  /// is non-zero \param pdb_field_value (optional) if non-zero, select only
  /// atoms whose pdb_field equals this
  virtual int load_coords(char const *filename,
                          std::vector<cvm::atom_pos> &pos,
                          std::vector<int> const &sorted_ids,
                          std::string const &pdb_field,
                          double pdb_field_value = 0.0);

  /// Clear atomic data
  int reset();

  /// Get the numeric ID of the given atom (for the program)
  inline int get_atom_id(int index) const
  {
    return atoms_ids[index];
  }

  /// Get the mass of the given atom
  inline cvm::real get_atom_mass(int index) const
  {
    return atoms_masses[index];
  }

  /// Get the charge of the given atom
  inline cvm::real get_atom_charge(int index) const
  {
    return atoms_charges[index];
  }

  /// Read the current position of the given atom
  inline cvm::rvector get_atom_position(int index) const
  {
    return atoms_positions[index];
  }

  /// Read the current total force of the given atom
  inline cvm::rvector get_atom_total_force(int index) const
  {
    return atoms_total_forces[index];
  }

  /// Request that this force is applied to the given atom
  inline void apply_atom_force(int index, cvm::rvector const &new_force)
  {
    atoms_new_colvar_forces[index] += new_force;
  }

  /// Read the current velocity of the given atom
  inline cvm::rvector get_atom_velocity(int index)
  {
    cvm::error("Error: reading the current velocity of an atom "
               "is not yet implemented.\n",
               COLVARS_NOT_IMPLEMENTED);
    return cvm::rvector(0.0);
  }

  inline std::vector<int> *modify_atom_ids()
  {
    return &atoms_ids;
  }

  inline std::vector<cvm::real> *modify_atom_masses()
  {
    // assume that we are requesting masses to change them
    updated_masses_ = true;
    return &atoms_masses;
  }

  inline std::vector<cvm::real> *modify_atom_charges()
  {
    // assume that we are requesting charges to change them
    updated_charges_ = true;
    return &atoms_charges;
  }

  inline std::vector<cvm::rvector> *modify_atom_positions()
  {
    return &atoms_positions;
  }

  inline std::vector<cvm::rvector> *modify_atom_total_forces()
  {
    return &atoms_total_forces;
  }

  inline std::vector<cvm::rvector> *modify_atom_new_colvar_forces()
  {
    return &atoms_new_colvar_forces;
  }

  /// Record whether masses have been updated
  inline bool updated_masses() const
  {
    return updated_masses_;
  }

  /// Record whether masses have been updated
  inline bool updated_charges() const
  {
    return updated_charges_;
  }

protected:

  /// \brief Array of 0-based integers used to uniquely associate atoms
  /// within the host program
  std::vector<int>          atoms_ids;
  /// \brief Keep track of how many times each atom is used by a separate colvar object
  std::vector<size_t>       atoms_ncopies;
  /// \brief Masses of the atoms (allow redefinition during a run, as done e.g. in LAMMPS)
  std::vector<cvm::real>    atoms_masses;
  /// \brief Charges of the atoms (allow redefinition during a run, as done e.g. in LAMMPS)
  std::vector<cvm::real>    atoms_charges;
  /// \brief Current three-dimensional positions of the atoms
  std::vector<cvm::rvector> atoms_positions;
  /// \brief Most recent total forces on each atom
  std::vector<cvm::rvector> atoms_total_forces;
  /// \brief Forces applied from colvars, to be communicated to the MD integrator
  std::vector<cvm::rvector> atoms_new_colvar_forces;

  /// Whether the masses and charges have been updated from the host code
  bool updated_masses_, updated_charges_;

  /// Used by all init_atom() functions: create a slot for an atom not
  /// requested yet; returns the index in the arrays
  int add_atom_slot(int atom_id);

};


/// \brief Container of atom group data (allow collection of aggregated atomic
/// data)
class colvarproxy_atom_groups {

public:

  /// Contructor
  colvarproxy_atom_groups();

  /// Destructor
  virtual ~colvarproxy_atom_groups();

  /// Clear atom group data
  int reset();

  /// \brief Whether this proxy implementation has capability for scalable groups
  virtual int scalable_group_coms();

  /// Prepare this group for collective variables calculation, selecting atoms by internal ids (0-based)
  virtual int init_atom_group(std::vector<int> const &atoms_ids);

  /// \brief Used by the atom_group class destructor
  virtual void clear_atom_group(int index);

  /// Get the numeric ID of the given atom group (for the MD program)
  inline int get_atom_group_id(int index) const
  {
    return atom_groups_ids[index];
  }

  /// Get the mass of the given atom group
  inline cvm::real get_atom_group_mass(int index) const
  {
    return atom_groups_masses[index];
  }

  /// Get the charge of the given atom group
  inline cvm::real get_atom_group_charge(int index) const
  {
    return atom_groups_charges[index];
  }

  /// Read the current position of the center of mass given atom group
  inline cvm::rvector get_atom_group_com(int index) const
  {
    return atom_groups_coms[index];
  }

  /// Read the current total force of the given atom group
  inline cvm::rvector get_atom_group_total_force(int index) const
  {
    return atom_groups_total_forces[index];
  }

  /// Request that this force is applied to the given atom group
  inline void apply_atom_group_force(int index, cvm::rvector const &new_force)
  {
    atom_groups_new_colvar_forces[index] += new_force;
  }

  /// Read the current velocity of the given atom group
  inline cvm::rvector get_atom_group_velocity(int index)
  {
    cvm::error("Error: reading the current velocity of an atom group is not yet implemented.\n",
               COLVARS_NOT_IMPLEMENTED);
    return cvm::rvector(0.0);
  }

protected:

  /// \brief Array of 0-based integers used to uniquely associate atom groups
  /// within the host program
  std::vector<int>          atom_groups_ids;
  /// \brief Keep track of how many times each group is used by a separate cvc
  std::vector<size_t>       atom_groups_ncopies;
  /// \brief Total masses of the atom groups
  std::vector<cvm::real>    atom_groups_masses;
  /// \brief Total charges of the atom groups (allow redefinition during a run, as done e.g. in LAMMPS)
  std::vector<cvm::real>    atom_groups_charges;
  /// \brief Current centers of mass of the atom groups
  std::vector<cvm::rvector> atom_groups_coms;
  /// \brief Most recently updated total forces on the com of each group
  std::vector<cvm::rvector> atom_groups_total_forces;
  /// \brief Forces applied from colvars, to be communicated to the MD integrator
  std::vector<cvm::rvector> atom_groups_new_colvar_forces;

  /// Used by all init_atom_group() functions: create a slot for an atom group not requested yet
  int add_atom_group_slot(int atom_group_id);
};


/// \brief Methods for SMP parallelization
class colvarproxy_smp {

public:

  /// Constructor
  colvarproxy_smp();

  /// Destructor
  virtual ~colvarproxy_smp();

  /// Whether threaded parallelization should be used (TODO: make this a
  /// cvm::deps feature)
  bool b_smp_active;

  /// Whether threaded parallelization is available (TODO: make this a cvm::deps feature)
  virtual int smp_enabled();

  /// Distribute calculation of colvars (and their components) across threads
  virtual int smp_colvars_loop();

  /// Distribute calculation of biases across threads
  virtual int smp_biases_loop();

  /// Distribute calculation of biases across threads 2nd through last, with all scripted biased on 1st thread
  virtual int smp_biases_script_loop();

  /// Index of this thread
  virtual int smp_thread_id();

  /// Number of threads sharing this address space
  virtual int smp_num_threads();

  /// Lock the proxy's shared data for access by a thread, if threads are implemented; if not implemented, does nothing
  virtual int smp_lock();

  /// Attempt to lock the proxy's shared data
  virtual int smp_trylock();

  /// Release the lock
  virtual int smp_unlock();

protected:

  /// Lock state for OpenMP
  void *omp_lock_state;
};


/// \brief Methods for multiple-replica communication
class colvarproxy_replicas {

public:

  /// Constructor
  colvarproxy_replicas();

  /// Destructor
  virtual ~colvarproxy_replicas();

  /// \brief Indicate if multi-replica support is available and active
  virtual bool replica_enabled();

  /// \brief Index of this replica
  virtual int replica_index();

  /// \brief Total number of replica
  virtual int replica_num();

  /// \brief Synchronize replica
  virtual void replica_comm_barrier();

  /// \brief Receive data from other replica
  virtual int replica_comm_recv(char* msg_data, int buf_len, int src_rep);

  /// \brief Send data to other replica
  virtual int replica_comm_send(char* msg_data, int msg_len, int dest_rep);

};


/// Methods for scripting language interface (Tcl or Python)
class colvarproxy_script {

public:

  /// Constructor
  colvarproxy_script();

  /// Destructor
  virtual ~colvarproxy_script();

  /// Convert a script object (Tcl or Python call argument) to a C string
  virtual char const *script_obj_to_str(unsigned char *obj);

  /// Convert a script object (Tcl or Python call argument) to a vector of strings
  virtual std::vector<std::string> script_obj_to_str_vector(unsigned char *obj);

  /// Pointer to the scripting interface object
  /// (does not need to be allocated in a new interface)
  colvarscript *script;

  /// is a user force script defined?
  bool force_script_defined;

  /// Do we have a scripting interface?
  bool have_scripts;

  /// Run a user-defined colvar forces script
  virtual int run_force_callback();

  virtual int run_colvar_callback(
                std::string const &name,
                std::vector<const colvarvalue *> const &cvcs,
                colvarvalue &value);

  virtual int run_colvar_gradient_callback(
                std::string const &name,
                std::vector<const colvarvalue *> const &cvcs,
                std::vector<cvm::matrix2d<cvm::real> > &gradient);
};


/// Methods for using Tcl within Colvars
class colvarproxy_tcl {

public:

  /// Constructor
  colvarproxy_tcl();

  /// Destructor
  virtual ~colvarproxy_tcl();

  /// Is Tcl available? (trigger initialization if needed)
  int tcl_available();

  /// Tcl implementation of script_obj_to_str()
  char const *tcl_obj_to_str(unsigned char *obj);

  /// Run a user-defined colvar forces script
  int tcl_run_force_callback();

  int tcl_run_colvar_callback(
              std::string const &name,
              std::vector<const colvarvalue *> const &cvcs,
              colvarvalue &value);

  int tcl_run_colvar_gradient_callback(
              std::string const &name,
              std::vector<const colvarvalue *> const &cvcs,
              std::vector<cvm::matrix2d<cvm::real> > &gradient);

protected:

  /// Pointer to Tcl interpreter object
  void *tcl_interp_;

  /// Set Tcl pointers
  virtual void init_tcl_pointers();
};


/// Methods for data input/output
class colvarproxy_io {

public:

  /// Constructor
  colvarproxy_io();

  /// Destructor
  virtual ~colvarproxy_io();

  /// \brief Save the current frame number in the argument given
  // Returns error code
  virtual int get_frame(long int &);

  /// \brief Set the current frame number (as well as colvarmodule::it)
  // Returns error code
  virtual int set_frame(long int);

  /// \brief Returns a reference to the given output channel;
  /// if this is not open already, then open it
  virtual std::ostream *output_stream(std::string const &output_name,
                                      std::ios_base::openmode mode =
                                        std::ios_base::out);

  /// \brief Flushes the given output channel
  virtual int flush_output_stream(std::ostream *os);

  /// \brief Closes the given output channel
  virtual int close_output_stream(std::string const &output_name);

  /// \brief Rename the given file, before overwriting it
  virtual int backup_file(char const *filename);

  /// \brief Rename the given file, before overwriting it
  inline int backup_file(std::string const &filename)
  {
    return backup_file(filename.c_str());
  }

  /// \brief Prefix of the input state file
  inline std::string & input_prefix()
  {
    return input_prefix_str;
  }

  /// \brief Prefix to be used for output restart files
  inline std::string & restart_output_prefix()
  {
    return restart_output_prefix_str;
  }

  /// \brief Prefix to be used for output files (final system
  /// configuration)
  inline std::string & output_prefix()
  {
    return output_prefix_str;
  }

protected:

  /// \brief Prefix to be used for input files (restarts, not
  /// configuration)
  std::string input_prefix_str, output_prefix_str, restart_output_prefix_str;

  /// \brief Currently opened output files: by default, these are ofstream objects.
  /// Allows redefinition to implement different output mechanisms
  std::list<std::ostream *> output_files;
  /// \brief Identifiers for output_stream objects: by default, these are the names of the files
  std::list<std::string>    output_stream_names;

};



/// \brief Interface between the collective variables module and
/// the simulation or analysis program (NAMD, VMD, LAMMPS...).
/// This is the base class: each interfaced program is supported by a derived class.
/// Only pure virtual functions ("= 0") must be reimplemented to ensure baseline functionality.
class colvarproxy
  : public colvarproxy_system,
    public colvarproxy_atoms,
    public colvarproxy_atom_groups,
    public colvarproxy_smp,
    public colvarproxy_replicas,
    public colvarproxy_script,
    public colvarproxy_tcl,
    public colvarproxy_io
{

public:

  /// Pointer to the main object
  colvarmodule *colvars;

  /// Constructor
  colvarproxy();

  /// Destructor
  virtual ~colvarproxy();

  /// Request deallocation of the module (currently only implemented by VMD)
  virtual int request_deletion();

  /// Whether deallocation was requested
  inline bool delete_requested()
  {
    return b_delete_requested;
  }

  /// \brief Reset proxy state, e.g. requested atoms
  virtual int reset();

  /// (Re)initialize required member data after construction
  virtual int setup();

  /// \brief Update data required by the colvars module (e.g. cache atom positions)
  ///
  /// TODO Break up colvarproxy_namd and colvarproxy_lammps function into these
  virtual int update_input();

  /// \brief Update data based from the results of a module update (e.g. send forces)
  virtual int update_output();

  /// Print a message to the main log
  virtual void log(std::string const &message) = 0;

  /// Print a message to the main log and let the rest of the program handle the error
  virtual void error(std::string const &message) = 0;

  /// Print a message to the main log and exit with error code
  virtual void fatal_error(std::string const &message) = 0;

  /// \brief Restarts will be written each time this number of steps has passed
  virtual size_t restart_frequency();

  /// Whether a simulation is running (warn against irrecovarable errors)
  inline bool simulation_running() const
  {
    return b_simulation_running;
  }

  /// Convert a version string "YYYY-MM-DD" into an integer
  int get_version_from_string(char const *version_string);

  /// Get the version number (higher = more recent)
  int version_number() const
  {
    return version_int;
  }

protected:

  /// Whether a simulation is running (warn against irrecovarable errors)
  bool b_simulation_running;

  /// Whether the entire module should be deallocated by the host engine
  bool b_delete_requested;

  /// Integer representing the version string (allows comparisons)
  int version_int;

};


#endif
