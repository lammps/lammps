// -*- c++ -*-

#ifndef COLVARDEPS_H
#define COLVARDEPS_H

#include "colvarmodule.h"
#include "colvarparse.h"

/// Parent class for a member object of a bias, cv or cvc etc. containing dependencies
/// (features) and handling dependency resolution

// Some features like colvar::f_linear have no dependencies, require() doesn't enable anything but fails if unavailable
// Policy: those features are unavailable at all times
// Other features are under user control
// They are unavailable unless requested by the user, then they may be enabled subject to
// satisfied dependencies

// It seems important to have available default to false (for safety) and enabled to false (for efficiency)

class colvardeps {
public:

  colvardeps() {}
  virtual ~colvardeps();

  // Subclasses should initialize the following members:

  std::string description; // reference to object name (cv, cvc etc.)

  /// This contains the current state of each feature for each object
  struct feature_state {
    feature_state(bool a, bool e)
    : available(a), enabled(e) {}

    /// Available means: supported, subject to dependencies as listed,
    /// MAY BE ENABLED AS A RESULT OF DEPENDENCY SOLVING
    /// Remains false for passive flags that are set based on other object properties,
    /// eg. f_cv_linear
    /// Is set to true upon user request for features that are implemented by the user
    /// or under his/her direct control, e.g. f_cv_scripted or f_cv_extended_Lagrangian
    bool available;
    /// Currently enabled - this flag is subject to change dynamically
    /// TODO consider implications for dependency solving: anyone who disables
    /// it should trigger a refresh of parent objects
    bool enabled;     // see if this should be private depending on implementation
    // bool enabledOnce; // this should trigger an update when object is evaluated
  };

  /// List of the state of all features
  std::vector<feature_state *> feature_states;

  /// Describes a feature and its dependecies
  /// used in a static array within each subclass
  class feature {

  public:
    feature() {}
    ~feature() {}

    std::string description; // Set by derived object initializer

    // features that this feature requires in the same object
    // NOTE: we have no safety mechanism against circular dependencies, however, they would have to be internal to an object (ie. requires_self or requires_alt)
    std::vector<int> requires_self;

    // Features that are incompatible, ie. required disabled
    // if enabled, they will cause a dependency failure (they will not be disabled)
    // To enforce these dependencies regardless of the order in which they
    // are enabled, they are always set in a symmetric way, so whichever is enabled
    // second will cause the dependency to fail
    std::vector<int> requires_exclude;

    // sets of features that are required in an alternate way
    // when parent feature is enabled, if none are enabled, the first one listed that is available will be enabled
    std::vector<std::vector<int> > requires_alt;

    // features that this feature requires in children
    std::vector<int> requires_children;
  };

  // Accessor to array of all features with deps, static in most derived classes
  // Subclasses with dynamic dependency trees may override this
  // with a non-static array
  // Intermediate classes (colvarbias and colvarcomp, which are also base classes)
  // implement this as virtual to allow overriding
  virtual std::vector<feature *>&features() = 0;

  void add_child(colvardeps *child);

  void remove_child(colvardeps *child);

  /// Used before deleting an object, if not handled by that object's destructor
  /// (useful for cvcs because their children are member objects)
  void remove_all_children();



private:
  // pointers to objects this object depends on
  // list should be maintained by any code that modifies the object
  // this could be secured by making lists of colvars / cvcs / atom groups private and modified through accessor functions
  std::vector<colvardeps *> children;

  // pointers to objects that depend on this object
  // the size of this array is in effect a reference counter
  std::vector<colvardeps *> parents;

public:
  // disabling a feature f:
  // if parents depend on f, tell them to refresh / check that they are ok?
  // if children provide features to satisfy f ONLY, disable that

  // When the state of this object has changed, recursively tell parents
  // to enforce their dependencies
//   void refresh_parents() {
//
//   }

  // std::vector<colvardeps *> parents; // Needed to trigger a refresh if capabilities of this object change

  // End of members to be initialized by subclasses

  // Checks whether given feature is enabled
  // Defaults to querying f_*_active
  inline bool is_enabled(int f = f_cv_active) const {
    return feature_states[f]->enabled;
  }

  // Checks whether given feature is available
  // Defaults to querying f_*_active
  inline bool is_available(int f = f_cv_active) const {
    return feature_states[f]->available;
  }

  void provide(int feature_id); // set the feature's flag to available in local object

  /// Parse a keyword and enable a feature accordingly
  bool get_keyval_feature(colvarparse *cvp,
                          std::string const &conf, char const *key,
                          int feature_id, bool const &def_value,
                          colvarparse::Parse_Mode const parse_mode = colvarparse::parse_normal);

  int enable(int f, bool dry_run = false, bool toplevel = true);  // enable a feature and recursively solve its dependencies
  // dry_run is set to true to recursively test if a feature is available, without enabling it
//     int disable(int f);


  /// This function is called whenever feature states are changed outside
  /// of the object's control, that is, by parents
  /// Eventually it may also be used when properties of children change
  virtual int refresh_deps() { return COLVARS_OK; }

  // NOTE that all feature enums should start with f_*_active
  enum features_biases {
    /// \brief Bias is active
    f_cvb_active,
    f_cvb_apply_force, // will apply forces
    f_cvb_get_total_force, // requires total forces
    f_cvb_history_dependent, // depends on simulation history
    f_cvb_ntot
  };

  enum features_colvar {
    /// \brief Calculate colvar
    f_cv_active,
    /// \brief Gradients are calculated and temporarily stored, so
    /// that external forces can be applied
    f_cv_gradient,
    /// \brief Collect atomic gradient data from all cvcs into vector
    /// atomic_gradient
    f_cv_collect_gradient,
    /// \brief Calculate the velocity with finite differences
    f_cv_fdiff_velocity,
    /// \brief The total force is calculated, projecting the atomic
    /// forces on the inverse gradient
    f_cv_total_force,
    /// \brief Calculate total force from atomic forces
    f_cv_total_force_calc,
    /// \brief Subtract the applied force from the total force
    f_cv_subtract_applied_force,
    /// \brief Estimate Jacobian derivative
    f_cv_Jacobian,
    /// \brief Do not report the Jacobian force as part of the total force
    /// instead, apply a correction internally to cancel it
    f_cv_hide_Jacobian,
    /// \brief The variable has a harmonic restraint around a moving
    /// center with fictitious mass; bias forces will be applied to
    /// the center
    f_cv_extended_Lagrangian,
    /// \brief The extended system coordinate undergoes Langevin dynamics
    f_cv_Langevin,
    /// \brief Output the potential and kinetic energies
    /// (for extended Lagrangian colvars only)
    f_cv_output_energy,
    /// \brief Output the value to the trajectory file (on by default)
    f_cv_output_value,
    /// \brief Output the velocity to the trajectory file
    f_cv_output_velocity,
    /// \brief Output the applied force to the trajectory file
    f_cv_output_applied_force,
    /// \brief Output the total force to the trajectory file
    f_cv_output_total_force,
    /// \brief A lower boundary is defined
    f_cv_lower_boundary,
    /// \brief An upper boundary is defined
    f_cv_upper_boundary,
    /// \brief Provide a discretization of the values of the colvar to
    /// be used by the biases or in analysis (needs lower and upper
    /// boundary)
    f_cv_grid,
    /// \brief Apply a restraining potential (|x-xb|^2) to the colvar
    /// when it goes below the lower wall
    f_cv_lower_wall,
    /// \brief Apply a restraining potential (|x-xb|^2) to the colvar
    /// when it goes above the upper wall
    f_cv_upper_wall,
    /// \brief Compute running average
    f_cv_runave,
    /// \brief Compute time correlation function
    f_cv_corrfunc,
    /// \brief Value and gradient computed by user script
    f_cv_scripted,
    /// \brief Colvar is periodic
    f_cv_periodic,
    /// \brief is scalar
    f_cv_scalar,
    f_cv_linear,
    f_cv_homogeneous,
    /// \brief Number of colvar features
    f_cv_ntot
  };

  enum features_cvc {
    f_cvc_active,
    f_cvc_scalar,
    f_cvc_gradient,
    f_cvc_inv_gradient,
    /// \brief If enabled, calc_gradients() will call debug_gradients() for every group needed
    f_cvc_debug_gradient,
    f_cvc_Jacobian,
    f_cvc_one_site_total_force,
    f_cvc_com_based,
    f_cvc_scalable,
    f_cvc_scalable_com,
    f_cvc_ntot
  };

  enum features_atomgroup {
    f_ag_active,
    f_ag_center,
    f_ag_rotate,
    f_ag_fitting_group,
    /// Perform a standard minimum msd fit for given atoms
    /// ie. not using refpositionsgroup
//     f_ag_min_msd_fit,
    f_ag_fit_gradient_group,// TODO check that these are sometimes needed separately
                            // maybe for minimum RMSD?
    f_ag_fit_gradient_ref,
    f_ag_atom_forces,
    f_ag_scalable,
    f_ag_scalable_com,
    f_ag_ntot
  };

  void init_cvb_requires();
  void init_cv_requires();
  void init_cvc_requires();
  void init_ag_requires();

  /// \brief print all enabled features and those of children, for debugging
  void print_state();
};

#endif


