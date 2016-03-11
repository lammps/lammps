/// -*- c++ -*-

#ifndef COLVARPROXY_H
#define COLVARPROXY_H

#include <fstream>
#include <list>

#include "colvarmodule.h"
#include "colvarvalue.h"

// forward declarations
class colvarscript;

/// \brief Interface between the collective variables module and
/// the simulation or analysis program (NAMD, VMD, LAMMPS...).
/// This is the base class: each interfaced program is supported by a derived class.
/// Only pure virtual functions ("= 0") must be reimplemented to ensure baseline functionality.

class colvarproxy {

public:

  /// Pointer to the main object
  colvarmodule *colvars;

  /// Default constructor
  inline colvarproxy() : script(NULL) {}

  /// Default destructor
  virtual ~colvarproxy() {}

  /// (Re)initialize required member data after construction
  virtual int setup()
  {
    return COLVARS_OK;
  }

  /// \brief Update data required by the colvars module (e.g. cache atom positions)
  ///
  /// TODO Break up colvarproxy_namd and colvarproxy_lammps function into these
  virtual int update_input()
  {
    return COLVARS_OK;
  }

  /// \brief Update data based from the results of a module update (e.g. send forces)
  virtual int update_output()
  {
    return COLVARS_OK;
  }

  // **************** SIMULATION PARAMETERS ****************

  /// \brief Value of the unit for atomic coordinates with respect to
  /// angstroms (used by some variables for hard-coded default values)
  virtual cvm::real unit_angstrom() = 0;

  /// \brief Boltzmann constant
  virtual cvm::real boltzmann() = 0;

  /// \brief Temperature of the simulation (K)
  virtual cvm::real temperature() = 0;

  /// \brief Time step of the simulation (fs)
  virtual cvm::real dt() = 0;

  /// \brief Pseudo-random number with Gaussian distribution
  virtual cvm::real rand_gaussian(void) = 0;

  /// \brief Get the current frame number
  virtual int frame() { return COLVARS_NOT_IMPLEMENTED; }

  /// \brief Set the current frame number
  // return 0 on success, -1 on failure
  virtual int frame(int) { return COLVARS_NOT_IMPLEMENTED; }

  /// \brief Prefix to be used for input files (restarts, not
  /// configuration)
  std::string input_prefix_str, output_prefix_str, restart_output_prefix_str;

  inline std::string input_prefix()
  {
    return input_prefix_str;
  }

  /// \brief Prefix to be used for output restart files
  inline std::string restart_output_prefix()
  {
    return restart_output_prefix_str;
  }

  /// \brief Prefix to be used for output files (final system
  /// configuration)
  inline std::string output_prefix()
  {
    return output_prefix_str;
  }

  /// \brief Restarts will be written each time this number of steps has passed
  virtual size_t restart_frequency()
  {
    return 0;
  }

protected:

  /// \brief Currently opened output files: by default, these are ofstream objects.
  /// Allows redefinition to implement different output mechanisms
  std::list<std::ostream *> output_files;
  /// \brief Identifiers for output_stream objects: by default, these are the names of the files
  std::list<std::string>    output_stream_names;

public:

  // **************** MULTIPLE REPLICAS COMMUNICATION ****************

  // Replica exchange commands:

  /// \brief Indicate if multi-replica support is available and active
  virtual bool replica_enabled() { return false; }

  /// \brief Index of this replica
  virtual int replica_index() { return 0; }

  /// \brief Total number of replica
  virtual int replica_num() { return 1; }

  /// \brief Synchronize replica
  virtual void replica_comm_barrier() {}

  /// \brief Receive data from other replica
  virtual int replica_comm_recv(char* msg_data, int buf_len, int src_rep) {
    return COLVARS_NOT_IMPLEMENTED;
  }

  /// \brief Send data to other replica
  virtual int replica_comm_send(char* msg_data, int msg_len, int dest_rep) {
    return COLVARS_NOT_IMPLEMENTED;
  }


  // **************** SCRIPTING INTERFACE ****************

  /// Pointer to the scripting interface object
  /// (does not need to be allocated in a new interface)
  colvarscript *script;

  /// is a user force script defined?
  bool force_script_defined;

  /// Do we have a scripting interface?
  bool have_scripts;

  /// Run a user-defined colvar forces script
  virtual int run_force_callback() { return COLVARS_NOT_IMPLEMENTED; }

  virtual int run_colvar_callback(std::string const &name,
                                  std::vector<const colvarvalue *> const &cvcs,
                                  colvarvalue &value)
  { return COLVARS_NOT_IMPLEMENTED; }

  virtual int run_colvar_gradient_callback(std::string const &name,
                                           std::vector<const colvarvalue *> const &cvcs,
                                           std::vector<cvm::matrix2d<cvm::real> > &gradient)
  { return COLVARS_NOT_IMPLEMENTED; }


  // **************** INPUT/OUTPUT ****************

  /// Print a message to the main log
  virtual void log(std::string const &message) = 0;

  /// Print a message to the main log and let the rest of the program handle the error
  virtual void error(std::string const &message) = 0;

  /// Print a message to the main log and exit with error code
  virtual void fatal_error(std::string const &message) = 0;

  /// Print a message to the main log and exit normally
  virtual void exit(std::string const &message)
  {
    cvm::error("Error: exiting without error is not implemented, returning error code.\n",
               COLVARS_NOT_IMPLEMENTED);
  }

  // TODO the following definitions may be moved to a .cpp file

  /// \brief Returns a reference to the given output channel;
  /// if this is not open already, then open it
  virtual std::ostream * output_stream(std::string const &output_name)
  {
    std::list<std::ostream *>::iterator osi  = output_files.begin();
    std::list<std::string>::iterator    osni = output_stream_names.begin();
    for ( ; osi != output_files.end(); osi++, osni++) {
      if (*osni == output_name) {
        return *osi;
      }
    }
    output_stream_names.push_back(output_name);
    std::ofstream * os = new std::ofstream(output_name.c_str());
    if (!os->is_open()) {
      cvm::error("Error: cannot write to file \""+output_name+"\".\n",
                 FILE_ERROR);
    }
    output_files.push_back(os);
    return os;
  }

  /// \brief Closes the given output channel
  virtual int close_output_stream(std::string const &output_name)
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
    cvm::error("Error: trying to close an output file or stream that wasn't open.\n",
               BUG_ERROR);
    return COLVARS_ERROR;
  }

  /// \brief Rename the given file, before overwriting it
  virtual int backup_file(char const *filename)
  {
    return COLVARS_NOT_IMPLEMENTED;
  }



  // **************** ACCESS SYSTEM DATA ****************

  /// Pass restraint energy value for current timestep to MD engine
  virtual void add_energy(cvm::real energy) = 0;

  /// Tell the proxy whether system forces are needed (may not always be available)
  virtual void request_system_force(bool yesno)
  {
    if (yesno == true)
      cvm::error("Error: system forces are currently not implemented.\n",
                 COLVARS_NOT_IMPLEMENTED);
  }

  /// \brief Get the PBC-aware distance vector between two positions
  virtual cvm::rvector position_distance(cvm::atom_pos const &pos1,
                                         cvm::atom_pos const &pos2) = 0;

  /// \brief Get the PBC-aware square distance between two positions;
  /// may need to be reimplemented independently from position_distance() for optimization purposes
  virtual cvm::real position_dist2(cvm::atom_pos const &pos1,
                                   cvm::atom_pos const &pos2)
  {
    return (position_distance(pos1, pos2)).norm2();
  }

  /// \brief Get the closest periodic image to a reference position
  /// \param pos The position to look for the closest periodic image
  /// \param ref_pos The reference position
  virtual void select_closest_image(cvm::atom_pos &pos,
                                    cvm::atom_pos const &ref_pos)
  {
    pos = position_distance(ref_pos, pos) + ref_pos;
  }

  /// \brief Perform select_closest_image() on a set of atomic positions
  ///
  /// After that, distance vectors can then be calculated directly,
  /// without using position_distance()
  void select_closest_images(std::vector<cvm::atom_pos> &pos,
                             cvm::atom_pos const &ref_pos)
  {
    for (std::vector<cvm::atom_pos>::iterator pi = pos.begin();
         pi != pos.end(); ++pi) {
      select_closest_image(*pi, ref_pos);
    }
  }


  // **************** ACCESS ATOMIC DATA ****************
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
  /// \brief Most recent forces applied by external potentials onto each atom
  std::vector<cvm::rvector> atoms_applied_forces;
  /// \brief Forces applied from colvars, to be communicated to the MD integrator
  std::vector<cvm::rvector> atoms_new_colvar_forces;

  /// Used by all init_atom() functions: create a slot for an atom not requested yet
  inline int add_atom_slot(int atom_id)
  {
    atoms_ids.push_back(atom_id);
    atoms_ncopies.push_back(1);
    atoms_masses.push_back(1.0);
    atoms_charges.push_back(0.0);
    atoms_positions.push_back(cvm::rvector(0.0, 0.0, 0.0));
    atoms_total_forces.push_back(cvm::rvector(0.0, 0.0, 0.0));
    atoms_applied_forces.push_back(cvm::rvector(0.0, 0.0, 0.0));
    atoms_new_colvar_forces.push_back(cvm::rvector(0.0, 0.0, 0.0));
    return (atoms_ids.size() - 1);
  }

public:

  /// Prepare this atom for collective variables calculation, selecting it by numeric index (1-based)
  virtual int init_atom(int atom_number) = 0;

  /// Check that this atom number is valid, but do not initialize the corresponding atom yet
  virtual int check_atom_id(int atom_number) = 0;

  /// Select this atom for collective variables calculation, using name and residue number.
  /// Not all programs support this: leave this function as is in those cases.
  virtual int init_atom(cvm::residue_id const &residue,
                        std::string const     &atom_name,
                        std::string const     &segment_id)
  {
    cvm::error("Error: initializing an atom by name and residue number is currently not supported.\n",
               COLVARS_NOT_IMPLEMENTED);
    return COLVARS_NOT_IMPLEMENTED;
  }

  /// Check that this atom is valid, but do not initialize it yet
  virtual int check_atom_id(cvm::residue_id const &residue,
                            std::string const     &atom_name,
                            std::string const     &segment_id)
  {
    cvm::error("Error: initializing an atom by name and residue number is currently not supported.\n",
               COLVARS_NOT_IMPLEMENTED);
    return COLVARS_NOT_IMPLEMENTED;
  }

  /// \brief Used by the atom class destructor: rather than deleting the array slot
  /// (costly) set the corresponding atoms_ncopies to zero
  virtual void clear_atom(int index)
  {
    if (((size_t) index) >= atoms_ids.size()) {
      cvm::error("Error: trying to disable an atom that was not previously requested.\n",
                 INPUT_ERROR);
    }
    if (atoms_ncopies[index] > 0) {
      atoms_ncopies[index] -= 1;
    }
  }

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

  /// Read the current total system force of the given atom
  inline cvm::rvector get_atom_system_force(int index) const
  {
    return atoms_total_forces[index] - atoms_applied_forces[index];
  }

  /// Request that this force is applied to the given atom
  inline void apply_atom_force(int index, cvm::rvector const &new_force)
  {
    atoms_new_colvar_forces[index] += new_force;
  }

  /// Read the current velocity of the given atom
  virtual cvm::rvector get_atom_velocity(int index)
  {
    cvm::error("Error: reading the current velocity of an atom is not yet implemented.\n",
               COLVARS_NOT_IMPLEMENTED);
    return cvm::rvector(0.0);
  }

  // useful functions for data management outside this class
  inline std::vector<int> *modify_atom_ids() { return &atoms_ids; }
  inline std::vector<cvm::real> *modify_atom_masses() { return &atoms_masses; }
  inline std::vector<cvm::real> *modify_atom_charges() { return &atoms_charges; }
  inline std::vector<cvm::rvector> *modify_atom_positions() { return &atoms_positions; }
  inline std::vector<cvm::rvector> *modify_atom_total_forces() { return &atoms_total_forces; }
  inline std::vector<cvm::rvector> *modify_atom_applied_forces() { return &atoms_applied_forces; }
  inline std::vector<cvm::rvector> *modify_atom_new_colvar_forces() { return &atoms_new_colvar_forces; }

  /// \brief Read atom identifiers from a file \param filename name of
  /// the file (usually a PDB) \param atoms array to which atoms read
  /// from "filename" will be appended \param pdb_field (optiona) if
  /// "filename" is a PDB file, use this field to determine which are
  /// the atoms to be set
  virtual int load_atoms(char const *filename,
                         cvm::atom_group &atoms,
                         std::string const &pdb_field,
                         double const pdb_field_value = 0.0)
  {
    cvm::error("Error: loading atom identifiers from a file is currently not implemented.\n",
               COLVARS_NOT_IMPLEMENTED);
    return COLVARS_NOT_IMPLEMENTED;
  }

  /// \brief Load the coordinates for a group of atoms from a file
  /// (usually a PDB); if "pos" is already allocated, the number of its
  /// elements must match the number of atoms in "filename"
  virtual int load_coords(char const *filename,
                          std::vector<cvm::atom_pos> &pos,
                          const std::vector<int> &indices,
                          std::string const &pdb_field,
                          double const pdb_field_value = 0.0)
  {
    cvm::error("Error: loading atomic coordinates from a file is currently not implemented.\n");
    return COLVARS_NOT_IMPLEMENTED;
  }

  // **************** ACCESS GROUP DATA ****************

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
  /// \brief Most recent forces applied by external potentials onto each group
  std::vector<cvm::rvector> atom_groups_applied_forces;
  /// \brief Forces applied from colvars, to be communicated to the MD integrator
  std::vector<cvm::rvector> atom_groups_new_colvar_forces;

  /// TODO Add here containers of handles to cvc objects that are computed in parallel

public:

  /// \brief Whether this proxy implementation has capability for scalable groups
  virtual bool has_scalable_groups() const
  {
    return false;
  }

  /// Used by all init_atom_group() functions: create a slot for an atom group not requested yet
  // TODO Add a handle to cvc objects
  inline int add_atom_group_slot(int atom_group_id)
  {
    atom_groups_ids.push_back(atom_group_id);
    atom_groups_masses.push_back(1.0);
    atom_groups_charges.push_back(0.0);
    atom_groups_coms.push_back(cvm::rvector(0.0, 0.0, 0.0));
    atom_groups_total_forces.push_back(cvm::rvector(0.0, 0.0, 0.0));
    atom_groups_applied_forces.push_back(cvm::rvector(0.0, 0.0, 0.0));
    atom_groups_new_colvar_forces.push_back(cvm::rvector(0.0, 0.0, 0.0));
    return (atom_groups_ids.size() - 1);
  }

  /// Prepare this group for collective variables calculation, selecting atoms by internal ids (0-based)
  virtual int init_atom_group(std::vector<int> const &atoms_ids) // TODO Add a handle to cvc objects
  {
    cvm::error("Error: initializing a group outside of the colvars module is currently not supported.\n",
               COLVARS_NOT_IMPLEMENTED);
    return COLVARS_NOT_IMPLEMENTED;
  }

  /// \brief Used by the atom_group class destructor
  virtual void clear_atom_group(int index)
  {
    if (((size_t) index) >= atom_groups_ids.size()) {
      cvm::error("Error: trying to disable an atom group that was not previously requested.\n",
                 INPUT_ERROR);
    }

    atom_groups_ids.erase(atom_groups_ids.begin()+index);
    atom_groups_masses.erase(atom_groups_masses.begin()+index);
    atom_groups_charges.erase(atom_groups_charges.begin()+index);
    atom_groups_coms.erase(atom_groups_coms.begin()+index);
    atom_groups_total_forces.erase(atom_groups_total_forces.begin()+index);
    atom_groups_applied_forces.erase(atom_groups_applied_forces.begin()+index);
    atom_groups_new_colvar_forces.erase(atom_groups_new_colvar_forces.begin()+index);
  }

  /// Get the numeric ID of the given atom group (for the MD program)
  inline cvm::real get_atom_group_id(int index) const
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

  /// Read the current total system force of the given atom group
  inline cvm::rvector get_atom_group_system_force(int index) const
  {
    return atom_groups_total_forces[index] - atom_groups_applied_forces[index];
  }

  /// Request that this force is applied to the given atom group
  inline void apply_atom_group_force(int index, cvm::rvector const &new_force)
  {
    atom_groups_new_colvar_forces[index] += new_force;
  }

  /// Read the current velocity of the given atom group
  virtual cvm::rvector get_atom_group_velocity(int index)
  {
    cvm::error("Error: reading the current velocity of an atom group is not yet implemented.\n",
               COLVARS_NOT_IMPLEMENTED);
    return cvm::rvector(0.0);
  }

};


#endif
