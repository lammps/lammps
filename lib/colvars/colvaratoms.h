// -*- c++ -*-

#ifndef COLVARATOMS_H
#define COLVARATOMS_H

#include "colvarmodule.h"
#include "colvarparse.h"

/// \brief Stores numeric id, mass and all mutable data for an atom,
/// mostly used by a \link cvc \endlink
/// 
/// This class may be used (although not necessarily) to keep atomic
/// data (id, mass, position and collective variable derivatives)
/// altogether.  There may be multiple instances with identical
/// numeric id, all acting independently: forces communicated through
/// these instances will be summed together.
///
/// Read/write operations depend on the underlying code: hence, some
/// member functions are defined in colvarproxy_xxx.h.
class colvarmodule::atom {

protected:

  /// \brief Index in the list of atoms involved by the colvars (\b
  /// NOT in the global topology!)
  size_t        index;

public:

  /// Internal identifier (zero-based)
  int              id;

  /// Mass
  cvm::real      mass;

  /// \brief Current position (copied from the program, can be
  /// manipulated)
  cvm::atom_pos   pos;

  /// \brief Current velocity (copied from the program, can be
  /// manipulated)
  cvm::rvector    vel;

  /// \brief System force at the previous step (copied from the
  /// program, can be manipulated)
  cvm::rvector    system_force;

  /// \brief Gradient of a scalar collective variable with respect
  /// to this atom
  /// 
  /// This can only handle a scalar collective variable (i.e. when
  /// the \link colvarvalue::real_value \endlink member is used
  /// from the \link colvarvalue \endlink class), which is also the
  /// most frequent case. For more complex types of \link
  /// colvarvalue \endlink objects, atomic gradients should be
  /// defined within the specific \link cvc \endlink
  /// implementation
  cvm::rvector   grad;

  /// \brief Default constructor, setting id to a non-valid value
  inline atom() {}

  /// \brief Initialize an atom for collective variable calculation
  /// and get its internal identifier \param atom_number Atom index in
  /// the system topology (starting from 1)
  atom (int const &atom_number);

  /// \brief Initialize an atom for collective variable calculation
  /// and get its internal identifier \param residue Residue number
  /// \param atom_name Name of the atom in the residue \param
  /// segment_id For PSF topologies, the segment identifier; for other
  /// type of topologies, may not be required
  atom (cvm::residue_id const &residue,
        std::string const     &atom_name,
        std::string const     &segment_id = std::string (""));

  /// Copy constructor
  atom (atom const &a);

  /// Destructor
  ~atom();

  /// Set non-constant data (everything except id and mass) to zero
  inline void reset_data() {
    pos = atom_pos (0.0);
    vel = grad = system_force = rvector (0.0);
  }

  /// Get the current position
  void read_position();

  /// Get the current velocity
  void read_velocity();

  /// Get the system force
  void read_system_force();

  /// \brief Apply a force to the atom
  ///
  /// The force will be used later by the MD integrator, the
  /// collective variables module does not integrate equations of
  /// motion.  Multiple calls to this function by either the same
  /// \link atom \endlink object or different objects with identical
  /// \link id \endlink, will all add to the existing MD force.
  void apply_force (cvm::rvector const &new_force);
};




/// \brief Group of \link atom \endlink objects, mostly used by a
/// \link cvc \endlink
///
/// This class inherits from \link colvarparse \endlink and from
/// std::vector<colvarmodule::atom>, and hence all functions and
/// operators (including the bracket operator, group[i]) can be used
/// on an \link atom_group \endlink object.  It can be initialized as
/// a vector, or by parsing a keyword in the configuration.
class colvarmodule::atom_group
  : public std::vector<cvm::atom>,
    public colvarparse
{
public:
  // Note: all members here are kept public, to make possible to any
  // object accessing and manipulating them


  /// \brief If this option is on, this group merely acts as a wrapper
  /// for a fixed position; any calls to atoms within or to
  /// functions that return disaggregated data will fail
  bool b_dummy;
  /// \brief dummy atom position
  cvm::atom_pos dummy_atom_pos;

  /// Sorted list of zero-based (internal) atom ids
  /// (populated on-demand by create_sorted_ids)
  std::vector<int> sorted_ids;

  /// Allocates and populates the sorted list of atom ids
  void create_sorted_ids (void);
  

  /// \brief Before calculating colvars, move the group to overlap the
  /// center of mass of reference coordinates
  bool b_center;

  /// \brief Right after updating atom coordinates (and after
  /// centering coordinates, if b_center is true), rotate the group to
  /// overlap the reference coordinates.  You should not manipulate
  /// atoms individually if you turn on this flag.
  ///
  /// Note: gradients will be calculated in the rotated frame: when
  /// forces will be applied, they will rotated back to the original
  /// frame
  bool b_rotate;
  /// Rotation between the group and its reference coordinates
  cvm::rotation rot;

  /// \brief In case b_center or b_rotate is true, use these reference
  /// coordinates
  std::vector<cvm::atom_pos> ref_pos;
  /// \brief Center of geometry of the reference coordinates; regardless
  /// of whether b_center is true, ref_pos is centered to zero at
  /// initialization, and ref_pos_cog serves to center the positions
  cvm::atom_pos              ref_pos_cog;
  /// \brief In case b_center or b_rotate is true, fit this group to
  /// the reference positions (default: the parent group itself)
  atom_group                *ref_pos_group;


  /// Total mass of the atom group
  cvm::real total_mass;

  /// \brief Don't apply any force on this group (use its coordinates
  /// only to calculate a colvar)
  bool        noforce;


  /// \brief Initialize the group by looking up its configuration
  /// string in conf and parsing it; this is actually done by parse(),
  /// which is a member function so that a group can be initialized
  /// also after construction
  atom_group (std::string const &conf,
              char const        *key,
              atom_group        *ref_pos_group = NULL);

  /// \brief Initialize the group by looking up its configuration
  /// string in conf and parsing it
  void parse (std::string const &conf,
              char const        *key,
              atom_group        *ref_pos_group = NULL);

  /// \brief Initialize the group after a temporary vector of atoms
  atom_group (std::vector<cvm::atom> const &atoms);

  /// \brief Add an atom to this group  
  void add_atom (cvm::atom const &a);

  /// \brief Default constructor
  atom_group();

  /// \brief Destructor
  ~atom_group();

  /// \brief Get the current positions; if b_center or b_rotate are
  /// true, center and/or rotate the coordinates right after reading
  /// them
  void read_positions();

  /// \brief Move all positions
  void apply_translation (cvm::rvector const &t);

  /// \brief Rotate all positions
  void apply_rotation (cvm::rotation const &q);


  /// \brief Get the current velocities; this must be called always
  /// *after* read_positions(); if b_rotate is defined, the same
  /// rotation applied to the coordinates will be used
  void read_velocities();

  /// \brief Get the current system_forces; this must be called always
  /// *after* read_positions(); if b_rotate is defined, the same
  /// rotation applied to the coordinates will be used
  void read_system_forces();

  /// Call reset_data() for each atom
  inline void reset_atoms_data()
  {
    for (cvm::atom_iter ai = this->begin(); ai != this->end(); ai++)
      ai->reset_data();
  }

  /// \brief Return a copy of the current atom positions
  std::vector<cvm::atom_pos> positions() const;

  /// \brief Return a copy of the current atom positions, shifted by a constant vector
  std::vector<cvm::atom_pos> positions_shifted (cvm::rvector const &shift) const;

  /// \brief Return the center of geometry of the positions \param ref_pos
  /// Use the closest periodic images to this position
  cvm::atom_pos center_of_geometry (cvm::atom_pos const &ref_pos);
  /// \brief Return the center of geometry of the positions, assuming
  /// that coordinates are already pbc-wrapped
  cvm::atom_pos center_of_geometry() const;

  /// \brief Return the center of mass of the positions \param ref_pos
  /// Use the closest periodic images to this position
  cvm::atom_pos center_of_mass (cvm::atom_pos const &ref_pos);
  /// \brief Return the center of mass of the positions, assuming that
  /// coordinates are already pbc-wrapped
  cvm::atom_pos center_of_mass() const;

  /// \brief Store atom positions from the previous step
  std::vector<cvm::atom_pos> old_pos;


  /// \brief Return a copy of the current atom velocities
  std::vector<cvm::rvector> velocities() const;


  /// \brief Return a copy of the system forces
  std::vector<cvm::rvector> system_forces() const;
  /// \brief Return a copy of the aggregated total force on the group
  cvm::rvector system_force() const;


  /// \brief Shorthand: save the specified gradient on each atom,
  /// weighting with the atom mass (mostly used in combination with
  /// \link center_of_mass() \endlink)
  void set_weighted_gradient (cvm::rvector const &grad);

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
  /// gradients are temporarily to the original frame.
  void apply_colvar_force (cvm::real const &force);

  /// \brief Apply a force "to the center of mass", i.e. the force is
  /// distributed on each atom according to its mass
  ///
  /// If the group is being rotated to a reference frame (e.g. to
  /// express the colvar independently from the solute rotation), the
  /// force is rotated back to the original frame.  Colvar gradients
  /// are not used, either because they were not defined (e.g because
  /// the colvar has not a scalar value) or the biases require to
  /// micromanage the force.
  void apply_force (cvm::rvector const &force);

  /// \brief Apply an array of forces directly on the individual
  /// atoms; the length of the specified vector must be the same of
  /// this \link atom_group \endlink.
  ///
  /// If the group is being rotated to a reference frame (e.g. to
  /// express the colvar independently from the solute rotation), the
  /// forces are rotated back to the original frame.  Colvar gradients
  /// are not used, either because they were not defined (e.g because
  /// the colvar has not a scalar value) or the biases require to
  /// micromanage the forces.
  void apply_forces (std::vector<cvm::rvector> const &forces);

};


#endif
