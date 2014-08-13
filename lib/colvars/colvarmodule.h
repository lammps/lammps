#ifndef COLVARMODULE_H
#define COLVARMODULE_H

#ifndef COLVARS_VERSION
#define COLVARS_VERSION "2014-08-13"
#endif

#ifndef COLVARS_DEBUG
#define COLVARS_DEBUG false
#endif

/// \file colvarmodule.h
/// \brief Collective variables main module
///
/// This file declares the main class for defining and manipulating
/// collective variables: there should be only one instance of this
/// class, because several variables are made static (i.e. they are
/// shared between all object instances) to be accessed from other
/// objects.


#include <iostream>
#include <iomanip>
#include <string>
#include <cstring>
#include <sstream>
#include <fstream>
#include <cmath>
#include <vector>
#include <list>


class colvarparse;
class colvar;
class colvarbias;
class colvarproxy;


/// \brief Collective variables module (main class)
///
/// Class to control the collective variables calculation.  An object
/// (usually one) of this class is spawned from the MD program,
/// containing all i/o routines and general interface.
///
/// At initialization, the colvarmodule object creates a proxy object
/// to provide a transparent interface between the MD program and the
/// child objects
class colvarmodule {

private:

  /// Impossible to initialize the main object without arguments
  colvarmodule();

public:

  friend class colvarproxy;

  /// Defining an abstract real number allows to switch precision
  typedef  double    real;
  /// Residue identifier
  typedef  int       residue_id;

  class rvector;
  template <class T,
            size_t const length> class vector1d;
  template <class T,
            size_t const outer_length,
            size_t const inner_length> class matrix2d;
  class quaternion;
  class rotation;

  /// \brief Atom position (different type name from rvector, to make
  /// possible future PBC-transparent implementations)
  typedef rvector atom_pos;

  /// \brief 3x3 matrix of real numbers
  class rmatrix;

  // allow these classes to access protected data
  class atom;
  class atom_group;
  friend class atom;
  friend class atom_group;
  typedef std::vector<atom>::iterator       atom_iter;
  typedef std::vector<atom>::const_iterator atom_const_iter;

  /// Current step number
  static size_t it;
  /// Starting step number for this run
  static size_t it_restart;

  /// Return the current step number from the beginning of this run
  static inline size_t step_relative()
  {
    return it - it_restart;
  }

  /// Return the current step number from the beginning of the whole
  /// calculation
  static inline size_t step_absolute()
  {
    return it;
  }

  /// If true, get it_restart from the state file; if set to false,
  /// the MD program is providing it
  bool it_restart_from_state_file;

  /// \brief Finite difference step size (if there is no dynamics, or
  /// if gradients need to be tested independently from the size of
  /// dt)
  static real debug_gradients_step_size;

  /// Prefix for all output files for this run
  static std::string output_prefix;

  /// Prefix for files from a previous run (including restart/output)
  static std::string input_prefix;

  /// input restart file name (determined from input_prefix)
  static std::string restart_in_name;


  /// Array of collective variables
  static std::vector<colvar *>     colvars;

  /* TODO: implement named CVCs
  /// Array of named (reusable) collective variable components
  static std::vector<cvc *>     cvcs;
  /// Named cvcs register themselves at initialization time
  inline void register_cvc (cvc *p) {
    cvcs.push_back(p);
  }
  */

  /// Array of collective variable biases
  static std::vector<colvarbias *> biases;
  /// \brief Number of ABF biases initialized (in normal conditions
  /// should be 1)
  static size_t n_abf_biases;
  /// \brief Number of metadynamics biases initialized (in normal
  /// conditions should be 1)
  static size_t n_meta_biases;
  /// \brief Number of restraint biases initialized (no limit on the
  /// number)
  static size_t n_rest_biases;
  /// \brief Number of histograms initialized (no limit on the
  /// number)
  static size_t n_histo_biases;

  /// \brief Whether debug output should be enabled (compile-time option)
  static inline bool debug()
  {
    return COLVARS_DEBUG;
  }


  /// \brief Constructor \param config_name Configuration file name
  /// \param restart_name (optional) Restart file name
  colvarmodule (char const *config_name,
                colvarproxy *proxy_in);

  /// Destructor
  ~colvarmodule();

  /// Initialize collective variables
  void init_colvars (std::string const &conf);

  /// Initialize collective variable biases
  void init_biases (std::string const &conf);

  /// Re-initialize data at the beginning of a run. For use with
  /// MD codes that can change system parameters like atom masses
  /// between run commands.
  void setup();

  /// Load new configuration for the given bias -
  /// currently works for harmonic (force constant and/or centers)
  void change_configuration (std::string const &bias_name, std::string const &conf);

  /// Read a colvar value
  std::string read_colvar(std::string const &name);

  /// Calculate change in energy from using alt. config. for the given bias -
  /// currently works for harmonic (force constant and/or centers)
  real energy_difference (std::string const &bias_name, std::string const &conf);

  /// Calculate collective variables and biases
  void calc();
  /// Read the input restart file
  std::istream & read_restart (std::istream &is);
  /// Write the output restart file
  std::ostream & write_restart (std::ostream &os);
  /// Write all output files (called by the proxy)
  void write_output_files();
  /// \brief Call colvarproxy::backup_file()
  static void backup_file (char const *filename);

  /// Perform analysis
  void analyze();
  /// \brief Read a collective variable trajectory (post-processing
  /// only, not called at runtime)
  bool read_traj (char const *traj_filename,
                  size_t      traj_read_begin,
                  size_t      traj_read_end);

  /// Get the pointer of a colvar from its name (returns NULL if not found)
  static colvar * colvar_p (std::string const &name);

  /// Quick conversion of an object to a string
  template<typename T> static std::string to_str (T const &x,
                                                  size_t const &width = 0,
                                                  size_t const &prec = 0);
  /// Quick conversion of a vector of objects to a string
  template<typename T> static std::string to_str (std::vector<T> const &x,
                                                  size_t const &width = 0,
                                                  size_t const &prec = 0);

  /// Reduce the number of characters in a string
  static inline std::string wrap_string (std::string const &s,
                                         size_t const &nchars)
  {
    if (!s.size())
      return std::string (nchars, ' ');
    else
      return ( (s.size() <= size_t (nchars)) ?
               (s+std::string (nchars-s.size(), ' ')) :
               (std::string (s, 0, nchars)) );
  }

  /// Number of characters to represent a time step
  static size_t const it_width;
  /// Number of digits to represent a collective variables value(s)
  static size_t const cv_prec;
  /// Number of characters to represent a collective variables value(s)
  static size_t const cv_width;
  /// Number of digits to represent the collective variables energy
  static size_t const en_prec;
  /// Number of characters to represent the collective variables energy
  static size_t const en_width;
  /// Line separator in the log output
  static std::string const line_marker;


  // proxy functions

  /// \brief Value of the unit for atomic coordinates with respect to
  /// angstroms (used by some variables for hard-coded default values)
  static real unit_angstrom();

  /// \brief Boltmann constant
  static real boltzmann();

  /// \brief Temperature of the simulation (K)
  static real temperature();

  /// \brief Time step of MD integrator (fs)
  static real dt();

  /// Request calculation of system force from MD engine
  static void request_system_force();

  /// Print a message to the main log
  static void log (std::string const &message);

  /// Print a message to the main log and exit with error code
  static void fatal_error (std::string const &message);

  /// Print a message to the main log and exit normally
  static void exit (std::string const &message);


  /// \brief Get the distance between two atomic positions with pbcs handled
  /// correctly
  static rvector position_distance (atom_pos const &pos1,
                                    atom_pos const &pos2);


  /// \brief Get the square distance between two positions (with
  /// periodic boundary conditions handled transparently)
  ///
  /// Note: in the case of periodic boundary conditions, this provides
  /// an analytical square distance (while taking the square of
  /// position_distance() would produce leads to a cusp)
  static real position_dist2 (atom_pos const &pos1,
                              atom_pos const &pos2);

  /// \brief Get the closest periodic image to a reference position
  /// \param pos The position to look for the closest periodic image
  /// \param ref_pos (optional) The reference position
  static void select_closest_image (atom_pos &pos,
                                    atom_pos const &ref_pos);

  /// \brief Perform select_closest_image() on a set of atomic positions
  ///
  /// After that, distance vectors can then be calculated directly,
  /// without using position_distance()
  static void select_closest_images (std::vector<atom_pos> &pos,
                                     atom_pos const &ref_pos);


  /// \brief Names of groups from a Gromacs .ndx file to be read at startup
  static std::list<std::string> index_group_names;

  /// \brief Groups from a Gromacs .ndx file read at startup
  static std::list<std::vector<int> > index_groups;

  /// \brief Read a Gromacs .ndx file
  static void read_index_file (char const *filename);


  /// \brief Create atoms from a file \param filename name of the file
  /// (usually a PDB) \param atoms array of the atoms to be allocated
  /// \param pdb_field (optiona) if "filename" is a PDB file, use this
  /// field to determine which are the atoms to be set
  static void load_atoms (char const *filename,
                          std::vector<atom> &atoms,
                          std::string const &pdb_field,
                          double const pdb_field_value = 0.0);

  /// \brief Load the coordinates for a group of atoms from a file
  /// (PDB or XYZ)
  static void load_coords (char const *filename,
                           std::vector<atom_pos> &pos,
                           const std::vector<int> &indices,
                           std::string const &pdb_field,
                           double const pdb_field_value = 0.0);

  /// \brief Load the coordinates for a group of atoms from an
  /// XYZ file
  static void load_coords_xyz (char const *filename,
                              std::vector<atom_pos> &pos,
                              const std::vector<int> &indices);

  /// Frequency for collective variables trajectory output
  static size_t cv_traj_freq;

  /// \brief True if only analysis is performed and not a run
  static bool   b_analysis;

  /// Frequency for saving output restarts
  static size_t restart_out_freq;
  /// Output restart file name
  std::string   restart_out_name;

  /// Pseudo-random number with Gaussian distribution
  static real rand_gaussian (void);
protected:

  /// Configuration file
  std::ifstream config_s;

  /// Configuration file parser object
  colvarparse *parse;

  /// Name of the trajectory file
  std::string   cv_traj_name;

  /// Collective variables output trajectory file
  std::ofstream cv_traj_os;

  /// Appending to the existing trajectory file?
  bool          cv_traj_append;

  /// Output restart file
  std::ofstream restart_out_os;

  /// \brief Pointer to the proxy object, used to retrieve atomic data
  /// from the hosting program; it is static in order to be accessible
  /// from static functions in the colvarmodule class
  static colvarproxy *proxy;

  /// \brief Counter for the current depth in the object hierarchy (useg e.g. in outpu
  static size_t depth;

public:

  /// Increase the depth (number of indentations in the output)
  static void increase_depth();

  /// Decrease the depth (number of indentations in the output)
  static void decrease_depth();
};


/// Shorthand for the frequently used type prefix
typedef colvarmodule cvm;


#include "colvartypes.h"


std::ostream & operator << (std::ostream &os, cvm::rvector const &v);
std::istream & operator >> (std::istream &is, cvm::rvector &v);


template<typename T> std::string cvm::to_str (T const &x,
                                              size_t const &width,
                                              size_t const &prec) {
  std::ostringstream os;
  if (width) os.width (width);
  if (prec) {
    os.setf (std::ios::scientific, std::ios::floatfield);
    os.precision (prec);
  }
  os << x;
  return os.str();
}

template<typename T> std::string cvm::to_str (std::vector<T> const &x,
                                              size_t const &width,
                                              size_t const &prec) {
  if (!x.size()) return std::string ("");
  std::ostringstream os;
  if (prec) {
    os.setf (std::ios::scientific, std::ios::floatfield);
  }
  os << "{ ";
  if (width) os.width (width);
  if (prec) os.precision (prec);
  os << x[0];
  for (size_t i = 1; i < x.size(); i++) {
    os << ", ";
    if (width) os.width (width);
    if (prec) os.precision (prec);
    os << x[i];
  }
  os << " }";
  return os.str();
}


#include "colvarproxy.h"


inline cvm::real cvm::unit_angstrom()
{
  return proxy->unit_angstrom();
}

inline cvm::real cvm::boltzmann()
{
  return proxy->boltzmann();
}

inline cvm::real cvm::temperature()
{
  return proxy->temperature();
}

inline cvm::real cvm::dt()
{
  return proxy->dt();
}

inline void cvm::request_system_force()
{
  proxy->request_system_force (true);
}

inline void cvm::select_closest_image (atom_pos &pos,
                                       atom_pos const &ref_pos)
{
  proxy->select_closest_image (pos, ref_pos);
}

inline void cvm::select_closest_images (std::vector<atom_pos> &pos,
                                        atom_pos const &ref_pos)
{
  proxy->select_closest_images (pos, ref_pos);
}

inline cvm::rvector cvm::position_distance (atom_pos const &pos1,
                                            atom_pos const &pos2)
{
  return proxy->position_distance (pos1, pos2);
}

inline cvm::real cvm::position_dist2 (cvm::atom_pos const &pos1,
                                      cvm::atom_pos const &pos2)
{
  return proxy->position_dist2 (pos1, pos2);
}

inline cvm::real cvm::rand_gaussian (void)
{
  return proxy->rand_gaussian();
}

#endif


// Emacs
// Local Variables:
// mode: C++
// End:
