// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#ifndef COLVARATOMS_H
#define COLVARATOMS_H

#include "colvarmodule.h"
#include "colvarproxy.h"
#include "colvarparse.h"
#include "colvardeps.h"


/// \brief Stores numeric id, mass and all mutable data for an atom,
/// mostly used by a \link colvar::cvc \endlink
///
/// This class may be used to keep atomic data such as id, mass,
/// position and collective variable derivatives) altogether.
/// There may be multiple instances with identical
/// numeric id, all acting independently: forces communicated through
/// these instances will be summed together.

class colvarmodule::atom {

protected:

  /// Index in the colvarproxy arrays (\b NOT in the global topology!)
  int index;

public:

  /// Identifier for the MD program (0-based)
  int             id;

  /// Mass
  cvm::real       mass;

  /// Charge
  cvm::real       charge;

  /// \brief Current position (copied from the program, can be
  /// modified if necessary)
  cvm::atom_pos   pos;

  /// \brief Current velocity (copied from the program, can be
  /// modified if necessary)
  cvm::rvector    vel;

  /// \brief System force at the previous step (copied from the
  /// program, can be modified if necessary)
  cvm::rvector    total_force;

  /// \brief Gradient of a scalar collective variable with respect
  /// to this atom
  ///
  /// This can only handle a scalar collective variable (i.e. when
  /// the \link colvarvalue::real_value \endlink member is used
  /// from the \link colvarvalue \endlink class), which is also the
  /// most frequent case. For more complex types of \link
  /// colvarvalue \endlink objects, atomic gradients should be
  /// defined within the specific \link colvar::cvc \endlink
  /// implementation
  cvm::rvector   grad;

  /// \brief Default constructor (sets index and id both to -1)
  atom();

  /// \brief Initialize an atom for collective variable calculation
  /// and get its internal identifier \param atom_number Atom index in
  /// the system topology (1-based)
  atom(int atom_number);

  /// \brief Initialize an atom for collective variable calculation
  /// and get its internal identifier \param residue Residue number
  /// \param atom_name Name of the atom in the residue \param
  /// segment_id For PSF topologies, the segment identifier; for other
  /// type of topologies, may not be required
  atom(cvm::residue_id const &residue,
       std::string const     &atom_name,
       std::string const     &segment_id);

  /// Copy constructor
  atom(atom const &a);

  /// Destructor
  ~atom();

  /// Set mutable data (everything except id and mass) to zero
  inline void reset_data()
  {
    pos = cvm::atom_pos(0.0);
    vel = grad = total_force = cvm::rvector(0.0);
  }

  /// Get the latest value of the mass
  inline void update_mass()
  {
    colvarproxy *p = cvm::proxy;
    mass = p->get_atom_mass(index);
  }

  /// Get the latest value of the charge
  inline void update_charge()
  {
    colvarproxy *p = cvm::proxy;
    charge = p->get_atom_charge(index);
  }

  /// Get the current position
  inline void read_position()
  {
    pos = (cvm::proxy)->get_atom_position(index);
  }

  /// Get the current velocity
  inline void read_velocity()
  {
    vel = (cvm::proxy)->get_atom_velocity(index);
  }

  /// Get the total force
  inline void read_total_force()
  {
    total_force = (cvm::proxy)->get_atom_total_force(index);
  }

  /// \brief Apply a force to the atom
  ///
  /// Note: the force is not applied instantly, but will be used later
  /// by the MD integrator (the colvars module does not integrate
  /// equations of motion.
  ///
  /// Multiple calls to this function by either the same
  /// \link atom \endlink object or different objects with identical
  /// \link id \endlink will all be added together.
  inline void apply_force(cvm::rvector const &new_force) const
  {
    (cvm::proxy)->apply_atom_force(index, new_force);
  }
};



/// \brief Group of \link atom \endlink objects, mostly used by a
/// \link colvar::cvc \endlink object to gather all atomic data
class colvarmodule::atom_group
  : public colvarparse, public colvardeps
{
public:


  /// \brief Default constructor
  atom_group();

  /// \brief Create a group object, assign a name to it
  atom_group(char const *key);

  /// \brief Initialize the group after a (temporary) vector of atoms
  atom_group(std::vector<cvm::atom> const &atoms_in);

  /// \brief Destructor
  ~atom_group();

  /// \brief Optional name to reuse properties of this in other groups
  std::string name;

  /// \brief Keyword used to define the group
  // TODO Make this field part of the data structures that link a group to a CVC
  std::string key;

  /// \brief Set default values for common flags
  int init();

  /// \brief Initialize dependency tree
  virtual int init_dependencies();

  /// \brief Update data required to calculate cvc's
  int setup();

  /// \brief Initialize the group by looking up its configuration
  /// string in conf and parsing it
  int parse(std::string const &conf);

  int add_atom_numbers(std::string const &numbers_conf);
  int add_atoms_of_group(atom_group const * ag);
  int add_index_group(std::string const &index_group_name);
  int add_atom_numbers_range(std::string const &range_conf);
  int add_atom_name_residue_range(std::string const &psf_segid,
                                  std::string const &range_conf);
  int parse_fitting_options(std::string const &group_conf);

  /// \brief Add an atom object to this group
  int add_atom(cvm::atom const &a);

  /// \brief Add an atom ID to this group (the actual atomicdata will be not be handled by the group)
  int add_atom_id(int aid);

  /// \brief Remove an atom object from this group
  int remove_atom(cvm::atom_iter ai);

  /// Set this group as a dummy group (no actual atoms)
  int set_dummy();

  /// If this group is dummy, set the corresponding position
  int set_dummy_pos(cvm::atom_pos const &pos);

  /// \brief Print the updated the total mass and charge of a group.
  /// This is needed in case the hosting MD code has an option to
  /// change atom masses after their initialization.
  void print_properties(std::string const &colvar_name, int i, int j);

  /// \brief Implementation of the feature list for atom group
  static std::vector<feature *> ag_features;

  /// \brief Implementation of the feature list accessor for atom group
  virtual const std::vector<feature *> &features() const
  {
    return ag_features;
  }
  virtual std::vector<feature *> &modify_features()
  {
    return ag_features;
  }
  static void delete_features() {
    for (size_t i=0; i < ag_features.size(); i++) {
      delete ag_features[i];
    }
    ag_features.clear();
  }

protected:

  /// \brief Array of atom objects
  std::vector<cvm::atom> atoms;

  /// \brief Internal atom IDs for host code
  std::vector<int> atoms_ids;

  /// Sorted list of internal atom IDs (populated on-demand by
  /// create_sorted_ids); used to read coordinate files
  std::vector<int> sorted_atoms_ids;

  /// Map entries of sorted_atoms_ids onto the original positions in the group
  std::vector<int> sorted_atoms_ids_map;

  /// \brief Dummy atom position
  cvm::atom_pos dummy_atom_pos;

  /// \brief Index in the colvarproxy arrays (if the group is scalable)
  int index;

public:

  inline cvm::atom & operator [] (size_t const i)
  {
    return atoms[i];
  }

  inline cvm::atom const & operator [] (size_t const i) const
  {
    return atoms[i];
  }

  inline cvm::atom_iter begin()
  {
    return atoms.begin();
  }

  inline cvm::atom_const_iter begin() const
  {
    return atoms.begin();
  }

  inline cvm::atom_iter end()
  {
    return atoms.end();
  }

  inline cvm::atom_const_iter end() const
  {
    return atoms.end();
  }

  inline size_t size() const
  {
    return atoms.size();
  }

  /// \brief If this option is on, this group merely acts as a wrapper
  /// for a fixed position; any calls to atoms within or to
  /// functions that return disaggregated data will fail
  bool b_dummy;

  /// Internal atom IDs (populated during initialization)
  inline std::vector<int> const &ids() const
  {
    return atoms_ids;
  }

  std::string const print_atom_ids() const;

  /// Allocates and populates sorted_ids and sorted_ids_map
  int create_sorted_ids();

  /// Sorted internal atom IDs (populated on-demand by create_sorted_ids);
  /// used to read coordinate files
  inline std::vector<int> const &sorted_ids() const
  {
    return sorted_atoms_ids;
  }

  /// Map entries of sorted_atoms_ids onto the original positions in the group
  inline std::vector<int> const &sorted_ids_map() const
  {
    return sorted_atoms_ids_map;
  }

  /// Detect whether two groups share atoms
  /// If yes, returns 1-based number of a common atom; else, returns 0
  static int overlap(const atom_group &g1, const atom_group &g2);

  /// The rotation calculated automatically if f_ag_rotate is defined
  cvm::rotation rot;

  /// \brief Indicates that the user has explicitly set centerToReference or
  /// rotateReference, and the corresponding reference:
  /// cvc's (eg rmsd, eigenvector) will not override the user's choice
  bool b_user_defined_fit;

  /// \brief use reference coordinates for f_ag_center or f_ag_rotate
  std::vector<cvm::atom_pos> ref_pos;

  /// \brief Center of geometry of the reference coordinates; regardless
  /// of whether f_ag_center is true, ref_pos is centered to zero at
  /// initialization, and ref_pos_cog serves to center the positions
  cvm::atom_pos              ref_pos_cog;

  /// \brief If f_ag_center or f_ag_rotate is true, use this group to
  /// define the transformation (default: this group itself)
  atom_group                *fitting_group;

  /// Total mass of the atom group
  cvm::real total_mass;

  /// Update the total mass of the atom group
  void update_total_mass();

  /// Total charge of the atom group
  cvm::real total_charge;

  /// Update the total mass of the group
  void update_total_charge();

  /// \brief Don't apply any force on this group (use its coordinates
  /// only to calculate a colvar)
  bool noforce;

  /// \brief Get the current positions
  void read_positions();

  /// \brief (Re)calculate the optimal roto-translation
  void calc_apply_roto_translation();

  /// \brief Save aside the center of geometry of the reference positions,
  /// then subtract it from them
  ///
  /// In this way it will be possible to use ref_pos also for the
  /// rotational fit.
  /// This is called either by atom_group::parse or by CVCs that assign
  /// reference positions (eg. RMSD, eigenvector).
  void center_ref_pos();

  /// \brief Move all positions
  void apply_translation(cvm::rvector const &t);

  /// \brief Get the current velocities; this must be called always
  /// *after* read_positions(); if f_ag_rotate is defined, the same
  /// rotation applied to the coordinates will be used
  void read_velocities();

  /// \brief Get the current total_forces; this must be called always
  /// *after* read_positions(); if f_ag_rotate is defined, the same
  /// rotation applied to the coordinates will be used
  void read_total_forces();

  /// Call reset_data() for each atom
  inline void reset_atoms_data()
  {
    for (cvm::atom_iter ai = atoms.begin(); ai != atoms.end(); ai++)
      ai->reset_data();
    if (fitting_group)
      fitting_group->reset_atoms_data();
  }

  /// \brief Recompute all mutable quantities that are required to compute CVCs
  int calc_required_properties();

  /// \brief Return a copy of the current atom positions
  std::vector<cvm::atom_pos> positions() const;

  /// \brief Calculate the center of geometry of the atomic positions, assuming
  /// that they are already pbc-wrapped
  int calc_center_of_geometry();

private:

  /// \brief Center of geometry
  cvm::atom_pos cog;

  /// \brief Center of geometry before any fitting
  cvm::atom_pos cog_orig;

public:

  /// \brief Return the center of geometry of the atomic positions
  inline cvm::atom_pos center_of_geometry() const
  {
    return cog;
  }

  /// \brief Calculate the center of mass of the atomic positions, assuming that
  /// they are already pbc-wrapped
  int calc_center_of_mass();

private:

  /// \brief Center of mass
  cvm::atom_pos com;

  /// \brief The derivative of a scalar variable with respect to the COM
  // TODO for scalable calculations of more complex variables (e.g. rotation),
  // use a colvarvalue of vectors to hold the entire derivative
  cvm::rvector scalar_com_gradient;

public:

  /// \brief Return the center of mass (COM) of the atomic positions
  inline cvm::atom_pos center_of_mass() const
  {
    return com;
  }

  /// \brief Return previously gradient of scalar variable with respect to the
  /// COM
  inline cvm::rvector center_of_mass_scalar_gradient() const
  {
    return scalar_com_gradient;
  }

  /// \brief Return a copy of the current atom positions, shifted by a constant vector
  std::vector<cvm::atom_pos> positions_shifted(cvm::rvector const &shift) const;

  /// \brief Return a copy of the current atom velocities
  std::vector<cvm::rvector> velocities() const;

  ///\brief Calculate the dipole of the atom group around the specified center
  int calc_dipole(cvm::atom_pos const &dipole_center);

private:

  /// Dipole moment of the atom group
  cvm::rvector dip;

public:

  ///\brief Return the (previously calculated) dipole of the atom group
  inline cvm::rvector dipole() const
  {
    return dip;
  }

  /// \brief Return a copy of the total forces
  std::vector<cvm::rvector> total_forces() const;

  /// \brief Return a copy of the aggregated total force on the group
  cvm::rvector total_force() const;


  /// \brief Shorthand: save the specified gradient on each atom,
  /// weighting with the atom mass (mostly used in combination with
  /// \link center_of_mass() \endlink)
  void set_weighted_gradient(cvm::rvector const &grad);

  /// \brief Calculate the derivatives of the fitting transformation
  void calc_fit_gradients();

  /// \brief Derivatives of the fitting transformation
  std::vector<cvm::atom_pos> fit_gradients;

  /// \brief Used by a (scalar) colvar to apply its force on its \link
  /// atom_group \endlink members
  ///
  /// The (scalar) force is multiplied by the colvar gradient for each
  /// atom; this should be used when a colvar with scalar \link
  /// colvarvalue \endlink type is used (this is the most frequent
  /// case: for colvars with a non-scalar type, the most convenient
  /// solution is to sum together the Cartesian forces from all the
  /// colvar components, and use apply_force() or apply_forces()).  If
  /// the group is being rotated to a reference frame (e.g. to express
  /// the colvar independently from the solute rotation), the
  /// gradients are temporarily rotated to the original frame.
  void apply_colvar_force(cvm::real const &force);

  /// \brief Apply a force "to the center of mass", i.e. the force is
  /// distributed on each atom according to its mass
  ///
  /// If the group is being rotated to a reference frame (e.g. to
  /// express the colvar independently from the solute rotation), the
  /// force is rotated back to the original frame.  Colvar gradients
  /// are not used, either because they were not defined (e.g because
  /// the colvar has not a scalar value) or the biases require to
  /// micromanage the force.
  /// This function will be phased out eventually, in favor of
  /// apply_colvar_force() once that is implemented for non-scalar values
  void apply_force(cvm::rvector const &force);

  /// Implements possible actions to be carried out
  /// when a given feature is enabled
  /// This overloads the base function in colvardeps
  void do_feature_side_effects(int id);
};


#endif
