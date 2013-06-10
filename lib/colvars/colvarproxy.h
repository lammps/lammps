#ifndef COLVARPROXY_H
#define COLVARPROXY_H



#include "colvarmodule.h"


/// \brief Interface class between the collective variables module and
/// the simulation program

class colvarproxy {

public:

  /// Pointer to the instance of colvarmodule
  colvarmodule *colvars;

  /// Default destructor
  virtual inline ~colvarproxy() {}


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


  // **************** SIMULATION PARAMETERS ****************

  /// \brief Prefix to be used for input files (restarts, not
  /// configuration)
  virtual std::string input_prefix() = 0;

  /// \brief Prefix to be used for output restart files
  virtual std::string restart_output_prefix() = 0;

  /// \brief Prefix to be used for output files (final system
  /// configuration)
  virtual std::string output_prefix() = 0;

  /// \brief Restarts will be fritten each time this number of steps has passed
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



  // **************** INPUT/OUTPUT ****************

  /// Print a message to the main log
  virtual void log (std::string const &message) = 0;

  /// Print a message to the main log and exit with error code
  virtual void fatal_error (std::string const &message) = 0;

  /// Print a message to the main log and exit normally
  virtual void exit (std::string const &message) = 0;

  /// \brief Read atom identifiers from a file \param filename name of
  /// the file (usually a PDB) \param atoms array to which atoms read
  /// from "filename" will be appended \param pdb_field (optiona) if
  /// "filename" is a PDB file, use this field to determine which are
  /// the atoms to be set
  virtual void load_atoms (char const *filename,
                           std::vector<cvm::atom> &atoms,
                           std::string const pdb_field,
                           double const pdb_field_value = 0.0) {}

  /// \brief Load the coordinates for a group of atoms from a file
  /// (usually a PDB); if "pos" is already allocated, the number of its
  /// elements must match the number of atoms in "filename"
  virtual void load_coords (char const *filename,
                            std::vector<cvm::atom_pos> &pos,
                            const std::vector<int> &indices,
                            std::string const pdb_field,
                            double const pdb_field_value = 0.0) = 0;

  /// \brief Rename the given file, before overwriting it
  virtual void backup_file (char const *filename) {}

};


inline void colvarproxy::select_closest_images (std::vector<cvm::atom_pos> &pos,
                                                cvm::atom_pos const &ref_pos)
{
  for (std::vector<cvm::atom_pos>::iterator pi = pos.begin();
       pi != pos.end(); pi++) {
    select_closest_image (*pi, ref_pos);
  }
}

inline cvm::real colvarproxy::position_dist2 (cvm::atom_pos const &pos1,
                                              cvm::atom_pos const &pos2)
{
  return (position_distance (pos1, pos2)).norm2();
}

#endif


// Emacs
// Local Variables:
// mode: C++
// End:
