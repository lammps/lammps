// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#ifndef COLVARDEPS_H
#define COLVARDEPS_H

#include "colvarmodule.h"
#include "colvarparse.h"

/// \brief Parent class for a member object of a bias, cv or cvc etc. containing features and
/// their dependencies, and handling dependency resolution
///
/// There are 3 kinds of features:
/// 1. Dynamic features are under the control of the dependency resolution
/// system. They may be enabled or disabled depending on dependencies.
/// 2. User features may be enabled based on user input (they may trigger a failure upon dependency resolution, though)
/// 3. Static features are static properties of the object, determined
///   programmatically at initialization time.
///
/// The following diagram summarizes the dependency tree at the bias, colvar, and colvarcomp levels.
/// Isolated and atom group features are not shown to save space.
/// @image html deps_2019.svg
///
/// In all classes, feature 0 is `active`. When an object is inactivated
/// all its children dependencies are dereferenced (free_children_deps)
/// While the object is inactive, no dependency solving is done on children
/// it is done when the object is activated back (restore_children_deps)
class colvardeps {
public:

  colvardeps();
  virtual ~colvardeps();

  // Subclasses should initialize the following members:

  std::string description; // reference to object name (cv, cvc etc.)

  /// This contains the current state of each feature for each object
  // since the feature class only contains static properties
  struct feature_state {
    feature_state(bool a, bool e)
    : available(a), enabled(e), ref_count(0) {}

    /// Feature may be enabled, subject to possible dependencies
    bool available;
    /// Currently enabled - this flag is subject to change dynamically
    /// TODO consider implications for dependency solving: anyone who disables
    /// it should trigger a refresh of parent objects
    bool enabled;     // see if this should be private depending on implementation

    // bool enabledOnce; // this should trigger an update when object is evaluated

    /// Number of features requiring this one as a dependency
    /// When it falls to zero:
    ///  - a dynamic feature is disabled automatically
    ///  - other features may be disabled statically
    int ref_count;
    /// List of features that were enabled by this one
    /// as part of an alternate requirement (for ref counting purposes)
    /// This is necessary because we don't know which feature in the list
    /// we enabled, otherwise
    std::vector<int> alternate_refs;
  };

protected:
  /// Time step multiplier (for coarse-timestep biases & colvars)
  /// Biases and colvars will only be calculated at those times
  /// (f_cvb_awake and f_cv_awake); a
  /// Biases use this to apply "impulse" biasing forces at the outer timestep
  /// Unused by lower-level objects (cvcs and atom groups)
  int   time_step_factor;

  /// List of the states of all features
  std::vector<feature_state> feature_states;

  /// Enum of possible feature types
  enum feature_type {
    f_type_not_set,
    f_type_dynamic,
    f_type_user,
    f_type_static
  };

public:
  /// \brief returns time_step_factor
  inline int get_time_step_factor() const {return time_step_factor;}

  /// Pair a numerical feature ID with a description and type
  void init_feature(int feature_id, const char *description, feature_type type);

  /// Describes a feature and its dependencies
  /// used in a static array within each subclass
  class feature {

  public:
    feature() : type(f_type_not_set) {}
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

    inline bool is_dynamic() { return type == f_type_dynamic; }
    inline bool is_static() { return type == f_type_static; }
    inline bool is_user() { return type == f_type_user; }
    /// Type of this feature, from the enum feature_type
    feature_type type;
  };

  inline bool is_not_set(int id) { return features()[id]->type == f_type_not_set; }
  inline bool is_dynamic(int id) { return features()[id]->type == f_type_dynamic; }
  inline bool is_static(int id) { return features()[id]->type == f_type_static; }
  inline bool is_user(int id) { return features()[id]->type == f_type_user; }

  // Accessor to array of all features with deps, static in most derived classes
  // Subclasses with dynamic dependency trees may override this
  // with a non-static array
  // Intermediate classes (colvarbias and colvarcomp, which are also base classes)
  // implement this as virtual to allow overriding
  virtual const std::vector<feature *> &features() const = 0;
  virtual std::vector<feature *>&modify_features() = 0;

  void add_child(colvardeps *child);

  void remove_child(colvardeps *child);

  /// Used before deleting an object, if not handled by that object's destructor
  /// (useful for cvcs because their children are member objects)
  void remove_all_children();

private:

  /// pointers to objects this object depends on
  /// list should be maintained by any code that modifies the object
  /// this could be secured by making lists of colvars / cvcs / atom groups private and modified through accessor functions
  std::vector<colvardeps *> children;

  /// pointers to objects that depend on this object
  /// the size of this array is in effect a reference counter
  std::vector<colvardeps *> parents;

public:
  // Checks whether given feature is enabled
  // Defaults to querying f_*_active
  inline bool is_enabled(int f = f_cv_active) const {
    return feature_states[f].enabled;
  }

  // Checks whether given feature is available
  // Defaults to querying f_*_active
  inline bool is_available(int f = f_cv_active) const {
    return feature_states[f].available;
  }

  /// Set the feature's available flag, without checking
  /// To be used for dynamic properties
  /// dependencies will be checked by enable()
  void provide(int feature_id, bool truefalse = true);

  /// Enable or disable, depending on flag value
  void set_enabled(int feature_id, bool truefalse = true);

protected:

  /// Parse a keyword and enable a feature accordingly
  bool get_keyval_feature(colvarparse *cvp,
                          std::string const &conf, char const *key,
                          int feature_id, bool const &def_value,
                          colvarparse::Parse_Mode const parse_mode = colvarparse::parse_normal);

public:

  /// Enable a feature and recursively solve its dependencies.
  /// For accurate reference counting, do not add spurious calls to enable()
  /// \param dry_run Recursively test whether a feature is available, without enabling it
  /// \param toplevel False if this is called as part of a chain of dependency resolution.
  /// This is used to diagnose failed dependencies by displaying the full stack:
  /// only the toplevel dependency will throw a fatal error.
  int enable(int f, bool dry_run = false, bool toplevel = true);

  /// Disable a feature, decrease the reference count of its dependencies
  /// and recursively disable them as applicable
  int disable(int f);

  /// disable all enabled features to free their dependencies
  /// to be done when deleting the object
  /// Cannot be in the base class destructor because it needs the derived class features()
  void free_children_deps();

  /// re-enable children features (to be used when object becomes active)
  void restore_children_deps();

  /// Decrement the reference count of a feature
  /// disabling it if it's dynamic and count reaches zero
  int decr_ref_count(int f);

  /// Implements possible actions to be carried out
  /// when a given feature is enabled
  /// Base function does nothing, can be overloaded
  virtual void do_feature_side_effects(int /* id */) {}

  // NOTE that all feature enums should start with f_*_active
  enum features_biases {
    /// \brief Bias is active
    f_cvb_active,
    /// \brief Bias is awake (active on its own accord) this timestep
    f_cvb_awake,
    /// Accumulates data starting from step 0 of a simulation run
    f_cvb_step_zero_data,
    /// \brief will apply forces
    f_cvb_apply_force,
    /// \brief force this bias to act on actual value for extended-Lagrangian coordinates
    f_cvb_bypass_ext_lagrangian,
    /// \brief requires total forces
    f_cvb_get_total_force,
    /// \brief whether this bias should record the accumulated work
    f_cvb_output_acc_work,
    /// \brief depends on simulation history
    f_cvb_history_dependent,
    /// \brief depends on time
    f_cvb_time_dependent,
    /// \brief requires scalar colvars
    f_cvb_scalar_variables,
    /// \brief whether this bias will compute a PMF
    f_cvb_calc_pmf,
    /// \brief whether this bias will compute TI samples
    f_cvb_calc_ti_samples,
    /// \brief whether this bias will write TI samples
    f_cvb_write_ti_samples,
    /// \brief whether this bias should write the TI PMF
    f_cvb_write_ti_pmf,
    /// \brief whether this bias uses an external grid to scale the biasing forces
    f_cvb_scale_biasing_force,
    f_cvb_ntot
  };

  enum features_colvar {
    /// \brief Calculate colvar
    f_cv_active,
    /// \brief Colvar is awake (active on its own accord) this timestep
    f_cv_awake,
    /// \brief Gradients are calculated and temporarily stored, so
    /// that external forces can be applied
    f_cv_gradient,
    /// \brief Collect atomic gradient data from all cvcs into vector
    /// atomic_gradient
    f_cv_collect_gradient,
    /// \brief Build list of atoms involved in CV calculation
    f_cv_collect_atom_ids,
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
    /// \brief An extended variable that sets an external variable in the
    /// back-end (eg. an alchemical coupling parameter for lambda-dynamics)
    /// Can have a single component
    f_cv_external,
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
    /// \brief The lower boundary is not defined from user's choice
    f_cv_hard_lower_boundary,
    /// \brief The upper boundary is not defined from user's choice
    f_cv_hard_upper_boundary,
    /// \brief Reflecting lower boundary condition
    f_cv_reflecting_lower_boundary,
    /// \brief Reflecting upper boundary condition
    f_cv_reflecting_upper_boundary,
    /// \brief Provide a discretization of the values of the colvar to
    /// be used by the biases or in analysis (needs lower and upper
    /// boundary)
    f_cv_grid,
    /// \brief Compute running average
    f_cv_runave,
    /// \brief Compute time correlation function
    f_cv_corrfunc,
    /// \brief Value and gradient computed by user script
    f_cv_scripted,
    /// \brief Value and gradient computed by user function through Lepton
    f_cv_custom_function,
    /// \brief Colvar is periodic
    f_cv_periodic,
    /// \brief The colvar has only one component
    f_cv_single_cvc,
    /// \brief is scalar
    f_cv_scalar,
    f_cv_linear,
    f_cv_homogeneous,
    /// \brief multiple timestep through time_step_factor
    f_cv_multiple_ts,
    /// \brief Number of colvar features
    f_cv_ntot
  };

  enum features_cvc {
    /// Computation of this CVC is enabled
    f_cvc_active,
    /// This CVC computes a scalar value
    f_cvc_scalar,
    /// Values of this CVC lie in a periodic interval
    f_cvc_periodic,
    /// This CVC provides a default value for the colvar's width
    f_cvc_width,
    /// This CVC provides a default value for the colvar's lower boundary
    f_cvc_lower_boundary,
    /// This CVC provides a default value for the colvar's upper boundary
    f_cvc_upper_boundary,
    /// CVC calculates atom gradients
    f_cvc_gradient,
    /// CVC calculates and stores explicit atom gradients on rank 0
    f_cvc_explicit_gradient,
    /// CVC calculates and stores inverse atom gradients (used for total force)
    f_cvc_inv_gradient,
    /// CVC calculates the Jacobian term of the total-force expression
    f_cvc_Jacobian,
    /// The total force for this CVC will be computed from one group only
    f_cvc_one_site_total_force,
    /// calc_gradients() will call debug_gradients() for every group needed
    f_cvc_debug_gradient,
    /// With PBCs, minimum-image convention will be used for distances
    /// (does not affect the periodicity of CVC values, e.g. angles)
    f_cvc_pbc_minimum_image,
    /// This CVC is a function of centers of mass
    f_cvc_com_based,
    /// This CVC can be computed in parallel
    f_cvc_scalable,
    /// Centers-of-mass used in this CVC can be computed in parallel
    f_cvc_scalable_com,
    /// \brief Build list of atoms involved in CVC calculation
    f_cvc_collect_atom_ids,
    /// Number of CVC features
    f_cvc_ntot
  };

  enum features_atomgroup {
    f_ag_active,
    f_ag_center,
    f_ag_center_origin,
    f_ag_rotate,
    f_ag_fitting_group,
    /// Perform a standard minimum msd fit for given atoms
    /// ie. not using refpositionsgroup
//     f_ag_min_msd_fit,
    /// \brief Does not have explicit atom gradients from parent CVC
    f_ag_explicit_gradient,
    f_ag_fit_gradients,
    f_ag_atom_forces,
    f_ag_scalable,
    f_ag_scalable_com,
    /// \brief Build list of atoms involved in atom group
    f_ag_collect_atom_ids,
    f_ag_ntot
  };

  /// Initialize dependency tree for object of a derived class
  virtual int init_dependencies() = 0;

  /// Make feature f require feature g within the same object
  void require_feature_self(int f, int g);

  /// Make features f and g mutually exclusive within the same object
  void exclude_feature_self(int f, int g);

  /// Make feature f require feature g within children
  void require_feature_children(int f, int g);

  /// Make feature f require either g or h within the same object
  void require_feature_alt(int f, int g, int h);

  /// Make feature f require any of g, h, or i within the same object
  void require_feature_alt(int f, int g, int h, int i);

  /// Make feature f require any of g, h, i, or j within the same object
  void require_feature_alt(int f, int g, int h, int i, int j);

  /// \brief print all enabled features and those of children, for debugging
  void print_state();

  /// \brief Check that a feature is enabled, raising COLVARS_BUG_ERROR if not
  inline void check_enabled(int f, std::string const &reason) const
  {
    if (! is_enabled(f)) {
      cvm::error("Error: "+reason+" requires that the feature \""+
                 features()[f]->description+"\" is active.\n", COLVARS_BUG_ERROR);
    }
  }

};

#endif


