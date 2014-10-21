/// -*- c++ -*-

#ifndef COLVARPROXY_H
#define COLVARPROXY_H


#ifndef COLVARPROXY_VERSION
#define COLVARPROXY_VERSION "2014-10-21"
#endif


#include "colvarmodule.h"
#include "colvarvalue.h"


// return values for the frame() routine
#define COLVARS_NO_SUCH_FRAME -1
#define COLVARS_NOT_IMPLEMENTED -2

// forward declarations
class colvarscript;

/// \brief Interface between the collective variables module and
/// the simulation or analysis program.
/// This is the base class: each program is supported by a derived class.
/// Only pure virtual functions ("= 0") must be reimplemented in a new interface.

class colvarproxy {

public:

  /// Pointer to the main object
  colvarmodule *colvars;

  /// Default constructor
  inline colvarproxy() : script (NULL) {}

  /// Default destructor
  virtual inline ~colvarproxy() {}

  /// (Re)initialize member data after construction
  virtual void setup() {}


  // **************** SYSTEM-WIDE PHYSICAL QUANTITIES ****************

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
  virtual cvm::real rand_gaussian (void) = 0;

  /// \brief Get the current frame number
  virtual int frame() { return COLVARS_NOT_IMPLEMENTED; }

  /// \brief Set the current frame number
  // return 0 on success, -1 on failure
  virtual int frame (int) { return COLVARS_NOT_IMPLEMENTED; }


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

  // **************** SIMULATION PARAMETERS ****************

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
  virtual size_t restart_frequency() = 0;



  // **************** ACCESS ATOMIC DATA ****************

  /// Pass restraint energy value for current timestep to MD engine
  virtual void add_energy (cvm::real energy) = 0;

  /// Tell the proxy whether system forces are needed (may not always be available)
  virtual void request_system_force (bool yesno) = 0;



  // **************** PERIODIC BOUNDARY CONDITIONS ****************

  /// \brief Get the PBC-aware distance vector between two positions
  virtual cvm::rvector position_distance (cvm::atom_pos const &pos1,
                                          cvm::atom_pos const &pos2) = 0;

  /// \brief Get the PBC-aware square distance between two positions;
  /// may be implemented independently from position_distance() for optimization purposes
  virtual cvm::real position_dist2 (cvm::atom_pos const &pos1,
                                    cvm::atom_pos const &pos2);

  /// \brief Get the closest periodic image to a reference position
  /// \param pos The position to look for the closest periodic image
  /// \param ref_pos The reference position
  virtual void select_closest_image (cvm::atom_pos &pos,
                                     cvm::atom_pos const &ref_pos) = 0;

  /// \brief Perform select_closest_image() on a set of atomic positions
  ///
  /// After that, distance vectors can then be calculated directly,
  /// without using position_distance()
  void select_closest_images (std::vector<cvm::atom_pos> &pos,
                              cvm::atom_pos const &ref_pos);

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
                                         std::vector<colvarvalue> &gradient)
  { return COLVARS_NOT_IMPLEMENTED; }

  // **************** INPUT/OUTPUT ****************

  /// Print a message to the main log
  virtual void log (std::string const &message) = 0;

  /// Print a message to the main log and let the rest of the program handle the error
  virtual void error (std::string const &message) = 0;

  /// Print a message to the main log and exit with error code
  virtual void fatal_error (std::string const &message) = 0;

  /// Print a message to the main log and exit normally
  virtual void exit (std::string const &message) = 0;

  /// \brief Read atom identifiers from a file \param filename name of
  /// the file (usually a PDB) \param atoms array to which atoms read
  /// from "filename" will be appended \param pdb_field (optiona) if
  /// "filename" is a PDB file, use this field to determine which are
  /// the atoms to be set
  virtual int load_atoms (char const *filename,
                           std::vector<cvm::atom> &atoms,
                           std::string const &pdb_field,
                           double const pdb_field_value = 0.0) = 0;

  /// \brief Load the coordinates for a group of atoms from a file
  /// (usually a PDB); if "pos" is already allocated, the number of its
  /// elements must match the number of atoms in "filename"
  virtual int load_coords (char const *filename,
                            std::vector<cvm::atom_pos> &pos,
                            const std::vector<int> &indices,
                            std::string const &pdb_field,
                            double const pdb_field_value = 0.0) = 0;

  /// \brief Rename the given file, before overwriting it
  virtual int backup_file (char const *filename)
  { return COLVARS_NOT_IMPLEMENTED; }
};


inline void colvarproxy::select_closest_images (std::vector<cvm::atom_pos> &pos,
                                                cvm::atom_pos const &ref_pos)
{
  for (std::vector<cvm::atom_pos>::iterator pi = pos.begin();
       pi != pos.end(); ++pi) {
    select_closest_image (*pi, ref_pos);
  }
}

inline cvm::real colvarproxy::position_dist2 (cvm::atom_pos const &pos1,
                                              cvm::atom_pos const &pos2)
{
  return (position_distance (pos1, pos2)).norm2();
}

#endif
