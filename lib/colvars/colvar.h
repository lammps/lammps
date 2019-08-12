// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#ifndef COLVAR_H
#define COLVAR_H

#include <iostream>

#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvarparse.h"
#include "colvardeps.h"

#ifdef LEPTON
#include "Lepton.h" // for runtime custom expressions
#endif

/// \brief A collective variable (main class); to be defined, it needs
/// at least one object of a derived class of colvar::cvc; it
/// calculates and returns a \link colvarvalue \endlink object
///
/// This class parses the configuration, defines the behaviour and
/// stores the value (\link colvar::x \endlink) and all related data
/// of a collective variable.  How the value is calculated is defined
/// in \link colvar::cvc \endlink and its derived classes.  The
/// \link colvar \endlink object contains pointers to multiple \link
/// colvar::cvc \endlink derived objects, which can be combined
/// together into one collective variable.  This makes possible to
/// implement new collective variables at runtime based on the
/// existing ones.  Currently, this possibility is limited to a
/// polynomial, using the coefficients cvc::sup_coeff and the
/// exponents cvc::sup_np.  In case of non-scalar variables,
/// only exponents equal to 1 are accepted.
///
/// Please note that most of its members are \link colvarvalue
/// \endlink objects, i.e. they can handle different data types
/// together, and must all be set to the same type of colvar::value()
/// before using them together in assignments or other operations; this is usually done
/// automatically in the constructor.  If you add a new member of
/// \link colvarvalue \endlink type, you should also add its
/// initialization line in the \link colvar \endlink constructor.

class colvar : public colvarparse, public colvardeps {

public:

  /// Name
  std::string name;

  /// \brief Current value (previously set by calc() or by read_traj())
  colvarvalue const & value() const;

  /// \brief Current actual value (not extended DOF)
  colvarvalue const & actual_value() const;

  /// \brief Current running average (if calculated as set by analysis flag)
  colvarvalue const & run_ave() const;

  /// \brief Force constant of the spring
  cvm::real const & force_constant() const;

  /// \brief Current velocity (previously set by calc() or by read_traj())
  colvarvalue const & velocity() const;

  /// \brief Current total force (previously obtained from calc() or
  /// read_traj()).  Note: this is calculated using the atomic forces
  /// from the last simulation step.
  ///
  /// Total atom forces are read from the MD program, the total force
  /// acting on the collective variable is calculated summing those
  /// from all colvar components, the bias and walls forces are
  /// subtracted.
  colvarvalue const & total_force() const;

  /// \brief Typical fluctuation amplitude for this collective
  /// variable (e.g. local width of a free energy basin)
  ///
  /// In metadynamics calculations, \link colvarbias_meta \endlink,
  /// this value is used to calculate the width of a gaussian.  In ABF
  /// calculations, \link colvarbias_abf \endlink, it is used to
  /// calculate the grid spacing in the direction of this collective
  /// variable.
  cvm::real width;

  /// \brief Implementation of the feature list for colvar
  static std::vector<feature *> cv_features;

  /// \brief Implementation of the feature list accessor for colvar
  virtual const std::vector<feature *> &features() const
  {
    return cv_features;
  }
  virtual std::vector<feature *> &modify_features()
  {
    return cv_features;
  }
  static void delete_features() {
    for (size_t i=0; i < cv_features.size(); i++) {
      delete cv_features[i];
    }
    cv_features.clear();
  }

  /// Implements possible actions to be carried out
  /// when a given feature is enabled
  /// This overloads the base function in colvardeps
  void do_feature_side_effects(int id);

  /// List of biases that depend on this colvar
  std::vector<colvarbias *> biases;
protected:


  /*
    extended:
    H = H_{0} + \sum_{i} 1/2*\lambda*(S_i(x(t))-s_i(t))^2 \\
    + \sum_{i} 1/2*m_i*(ds_i(t)/dt)^2 \\
    + \sum_{t'<t} W * exp(-1/2*\sum_{i} (s_i(t')-s_i(t))^2/(\delta{}s_i)^2) \\
    + \sum_{w} (\sum_{i}c_{w,i}s_i(t) - D_w)^M/(\Sigma_w)^M

    normal:
    H = H_{0} + \sum_{t'<t} W * exp(-1/2*\sum_{i} (S_i(x(t'))-S_i(x(t)))^2/(\delta{}S_i)^2) \\
    + \sum_{w} (\sum_{i}c_{w,i}S_i(t) - D_w)^M/(\Sigma_w)^M

    output:
    H = H_{0}   (only output S(x), no forces)

    Here:
    S(x(t)) = x
    s(t)    = x_ext
    DS = Ds = delta
  */


  /// Value of the colvar
  colvarvalue x;

  // TODO: implement functionality to treat these
  // /// Vector of individual values from CVCs
  // colvarvalue x_cvc;

  // /// Jacobian matrix of individual values from CVCs
  // colvarvalue dx_cvc;

  /// Cached reported value (x may be manipulated)
  colvarvalue x_reported;

  /// Finite-difference velocity
  colvarvalue v_fdiff;

  inline colvarvalue fdiff_velocity(colvarvalue const &xold,
                                     colvarvalue const &xnew)
  {
    // using the gradient of the square distance to calculate the
    // velocity (non-scalar variables automatically taken into
    // account)
    cvm::real dt = cvm::dt();
    return ( ( (dt > 0.0) ? (1.0/dt) : 1.0 ) *
             0.5 * dist2_lgrad(xnew, xold) );
  }

  /// Cached reported velocity
  colvarvalue v_reported;

  // Options for extended_lagrangian
  /// Restraint center
  colvarvalue x_ext;
  /// Previous value of the restraint center;
  colvarvalue prev_x_ext;
  /// Velocity of the restraint center
  colvarvalue v_ext;
  /// Previous velocity of the restraint center
  colvarvalue prev_v_ext;
  /// Mass of the restraint center
  cvm::real ext_mass;
  /// Restraint force constant
  cvm::real ext_force_k;
  /// Friction coefficient for Langevin extended dynamics
  cvm::real ext_gamma;
  /// Amplitude of Gaussian white noise for Langevin extended dynamics
  cvm::real ext_sigma;

  /// \brief Harmonic restraint force
  colvarvalue fr;

  /// \brief Jacobian force, when Jacobian_force is enabled
  colvarvalue fj;

  /// Cached reported total force
  colvarvalue ft_reported;

public:


  /// \brief Bias force; reset_bias_force() should be called before
  /// the biases are updated
  colvarvalue fb;

  /// \brief Bias force to the actual value (only useful with extended Lagrangian)
  colvarvalue fb_actual;

  /// \brief Total \em applied force; fr (if extended_lagrangian
  /// is defined), fb (if biases are applied) and the walls' forces
  /// (if defined) contribute to it
  colvarvalue f;

  /// Applied force at the previous step (to be subtracted from total force if needed)
  colvarvalue f_old;

  /// \brief Total force, as derived from the atomic trajectory;
  /// should equal the system force plus \link f \endlink
  colvarvalue ft;


  /// Period, if this variable is periodic
  cvm::real period;
  cvm::real wrap_center;


  /// \brief Expand the boundaries of multiples of width, to keep the
  /// value always within range
  bool expand_boundaries;

  /// \brief Location of the lower boundary
  colvarvalue lower_boundary;
  /// \brief Location of the lower wall
  colvarvalue lower_wall;
  /// \brief Force constant for the lower boundary potential (|x-xb|^2)
  cvm::real   lower_wall_k;
  /// \brief Whether this colvar has a hard lower boundary
  bool        hard_lower_boundary;

  /// \brief Location of the upper boundary
  colvarvalue upper_boundary;
  /// \brief Location of the upper wall
  colvarvalue upper_wall;
  /// \brief Force constant for the upper boundary potential (|x-xb|^2)
  cvm::real   upper_wall_k;
  /// \brief Whether this colvar has a hard upper boundary
  bool        hard_upper_boundary;

  /// \brief Is the interval defined by the two boundaries periodic?
  bool periodic_boundaries() const;

  /// \brief Is the interval defined by the two boundaries periodic?
  bool periodic_boundaries(colvarvalue const &lb, colvarvalue const &ub) const;


  /// Constructor
  colvar();

  /// Main init function
  int init(std::string const &conf);

  /// Parse the CVC configuration and allocate their data
  int init_components(std::string const &conf);

  /// Parse parameters for custom function with Lepton
  int init_custom_function(std::string const &conf);

  /// Init defaults for grid options
  int init_grid_parameters(std::string const &conf);

  /// Init extended Lagrangian parameters
  int init_extended_Lagrangian(std::string const &conf);

  /// Init output flags
  int init_output_flags(std::string const &conf);

  /// \brief Initialize dependency tree
  virtual int init_dependencies();

private:
  /// Parse the CVC configuration for all components of a certain type
  template<typename def_class_name> int init_components_type(std::string const &conf,
                                                             char const *def_desc,
                                                             char const *def_config_key);

public:

  /// Get ready for a run and re-initialize internal data if needed
  void setup();

  /// Destructor
  ~colvar();


  /// \brief Calculate the colvar's value and related quantities
  int calc();

  /// Carry out operations needed before next step is run
  int end_of_step();

  /// \brief Calculate a subset of the colvar components (CVCs) currently active
  /// (default: all active CVCs)
  /// Note: both arguments refer to the sect of *active* CVCs, not all CVCs
  int calc_cvcs(int first = 0, size_t num_cvcs = 0);

  /// Ensure that the selected range of CVCs is consistent
  int check_cvc_range(int first_cvc, size_t num_cvcs);

  /// \brief Calculate the values of the given subset of CVCs
  int calc_cvc_values(int first, size_t num_cvcs);
  /// \brief Same as \link colvar::calc_cvc_values \endlink but for gradients
  int calc_cvc_gradients(int first, size_t num_cvcs);
  /// \brief Same as \link colvar::calc_cvc_values \endlink but for total forces
  int calc_cvc_total_force(int first, size_t num_cvcs);
  /// \brief Same as \link colvar::calc_cvc_values \endlink but for Jacobian derivatives/forces
  int calc_cvc_Jacobians(int first, size_t num_cvcs);

  /// \brief Collect quantities from CVCs and update aggregated data for the colvar
  int collect_cvc_data();

  /// \brief Collect the values of the CVCs
  int collect_cvc_values();
  /// \brief Same as \link colvar::collect_cvc_values \endlink but for gradients
  int collect_cvc_gradients();
  /// \brief Same as \link colvar::collect_cvc_values \endlink but for total forces
  int collect_cvc_total_forces();
  /// \brief Same as \link colvar::collect_cvc_values \endlink but for Jacobian derivatives/forces
  int collect_cvc_Jacobians();
  /// \brief Calculate the quantities associated to the colvar (but not to the CVCs)
  int calc_colvar_properties();

  /// Get the current applied force
  inline colvarvalue const applied_force() const
  {
    if (is_enabled(f_cv_extended_Lagrangian)) {
      return fr;
    }
    return f;
  }

  /// Set the total biasing force to zero
  void reset_bias_force();

  /// Add to the total force from biases
  void add_bias_force(colvarvalue const &force);

  /// Apply a force to the actual value (only meaningful with extended Lagrangian)
  void add_bias_force_actual_value(colvarvalue const &force);

  /// \brief Collect all forces on this colvar, integrate internal
  /// equations of motion of internal degrees of freedom; see also
  /// colvar::communicate_forces()
  /// return colvar energy if extended Lagrandian active
  cvm::real update_forces_energy();

  /// \brief Communicate forces (previously calculated in
  /// colvar::update()) to the external degrees of freedom
  void communicate_forces();

  /// \brief Enables and disables individual CVCs based on flags
  int set_cvc_flags(std::vector<bool> const &flags);

  /// \brief Updates the flags in the CVC objects, and their number
  int update_cvc_flags();

  /// \brief Modify the configuration of CVCs (currently, only base class data)
  int update_cvc_config(std::vector<std::string> const &confs);

protected:
  /// \brief Number of CVC objects with an active flag
  size_t n_active_cvcs;

  /// Sum of square coefficients for active cvcs
  cvm::real active_cvc_square_norm;

  /// Update the sum of square coefficients for active cvcs
  void update_active_cvc_square_norm();

  /// \brief Absolute timestep number when this colvar was last updated
  cvm::step_number prev_timestep;

public:

  /// \brief Return the number of CVC objects defined
  inline size_t num_cvcs() const { return cvcs.size(); }

  /// \brief Return the number of CVC objects with an active flag (as set by update_cvc_flags)
  inline size_t num_active_cvcs() const { return n_active_cvcs; }

  /// \brief Use the internal metrics (as from \link colvar::cvc
  /// \endlink objects) to calculate square distances and gradients
  ///
  /// Handles correctly symmetries and periodic boundary conditions
  cvm::real dist2(colvarvalue const &x1,
                  colvarvalue const &x2) const;

  /// \brief Use the internal metrics (as from \link colvar::cvc
  /// \endlink objects) to calculate square distances and gradients
  ///
  /// Handles correctly symmetries and periodic boundary conditions
  colvarvalue dist2_lgrad(colvarvalue const &x1,
                          colvarvalue const &x2) const;

  /// \brief Use the internal metrics (as from \link colvar::cvc
  /// \endlink objects) to calculate square distances and gradients
  ///
  /// Handles correctly symmetries and periodic boundary conditions
  colvarvalue dist2_rgrad(colvarvalue const &x1,
                          colvarvalue const &x2) const;

  /// \brief Use the internal metrics (as from \link colvar::cvc
  /// \endlink objects) to wrap a value into a standard interval
  ///
  /// Handles correctly symmetries and periodic boundary conditions
  void wrap(colvarvalue &x_unwrapped) const;


  /// Read the analysis tasks
  int parse_analysis(std::string const &conf);

  /// Perform analysis tasks
  int analyze();


  /// Read the value from a collective variable trajectory file
  std::istream & read_traj(std::istream &is);
  /// Output formatted values to the trajectory file
  std::ostream & write_traj(std::ostream &os);
  /// Write a label to the trajectory file (comment line)
  std::ostream & write_traj_label(std::ostream &os);


  /// Read the collective variable from a restart file
  std::istream & read_restart(std::istream &is);
  /// Write the collective variable to a restart file
  std::ostream & write_restart(std::ostream &os);

  /// Write output files (if defined, e.g. in analysis mode)
  int write_output_files();


protected:
  /// Previous value (to calculate velocities during analysis)
  colvarvalue            x_old;

  /// Value read from the most recent state file (if any)
  colvarvalue            x_restart;

  /// True if a state file was just read
  bool                   after_restart;

  /// Time series of values and velocities used in correlation
  /// functions
  std::list< std::list<colvarvalue> > acf_x_history, acf_v_history;
  /// Time series of values and velocities used in correlation
  /// functions (pointers)x
  std::list< std::list<colvarvalue> >::iterator acf_x_history_p, acf_v_history_p;

  /// Time series of values and velocities used in running averages
  std::list< std::list<colvarvalue> > x_history;
  /// Time series of values and velocities used in correlation
  /// functions (pointers)x
  std::list< std::list<colvarvalue> >::iterator x_history_p;

  /// \brief Collective variable with which the correlation is
  /// calculated (default: itself)
  std::string            acf_colvar_name;
  /// Length of autocorrelation function (ACF)
  size_t                 acf_length;
  /// After how many steps the ACF starts
  size_t                 acf_offset;
  /// How many timesteps separate two ACF values
  size_t                 acf_stride;
  /// Number of frames for each ACF point
  size_t                 acf_nframes;
  /// Normalize the ACF to a maximum value of 1?
  bool                   acf_normalize;
  /// ACF values
  std::vector<cvm::real> acf;
  /// Name of the file to write the ACF
  std::string            acf_outfile;

  /// Type of autocorrelation function (ACF)
  enum acf_type_e {
    /// Unset type
    acf_notset,
    /// Velocity ACF, scalar product between v(0) and v(t)
    acf_vel,
    /// Coordinate ACF, scalar product between x(0) and x(t)
    acf_coor,
    /// \brief Coordinate ACF, second order Legendre polynomial
    /// between x(0) and x(t) (does not work with scalar numbers)
    acf_p2coor
  };

  /// Type of autocorrelation function (ACF)
  acf_type_e             acf_type;

  /// \brief Velocity ACF, scalar product between v(0) and v(t)
  void calc_vel_acf(std::list<colvarvalue> &v_history,
                    colvarvalue const      &v);

  /// \brief Coordinate ACF, scalar product between x(0) and x(t)
  /// (does not work with scalar numbers)
  void calc_coor_acf(std::list<colvarvalue> &x_history,
                     colvarvalue const      &x);

  /// \brief Coordinate ACF, second order Legendre polynomial between
  /// x(0) and x(t) (does not work with scalar numbers)
  void calc_p2coor_acf(std::list<colvarvalue> &x_history,
                       colvarvalue const      &x);

  /// Calculate the auto-correlation function (ACF)
  int calc_acf();
  /// Save the ACF to a file
  int write_acf(std::ostream &os);

  /// Length of running average series
  size_t         runave_length;
  /// Timesteps to skip between two values in the running average series
  size_t         runave_stride;
  /// Name of the file to write the running average
  std::string    runave_outfile;
  /// File to write the running average
  std::ostream  *runave_os;
  /// Current value of the running average
  colvarvalue    runave;
  /// Current value of the square deviation from the running average
  cvm::real      runave_variance;

  /// Calculate the running average and its standard deviation
  int calc_runave();

  /// If extended Lagrangian active: colvar energies (kinetic and harmonic potential)
  cvm::real kinetic_energy;
  cvm::real potential_energy;
public:


  // collective variable component base class
  class cvc;

  // list of available collective variable components

  // scalar colvar components
  class distance;
  class distance_z;
  class distance_xy;
  class polar_theta;
  class polar_phi;
  class distance_inv;
  class distance_pairs;
  class dipole_magnitude;
  class angle;
  class dipole_angle;
  class dihedral;
  class coordnum;
  class selfcoordnum;
  class groupcoordnum;
  class h_bond;
  class rmsd;
  class orientation_angle;
  class orientation_proj;
  class tilt;
  class spin_angle;
  class gyration;
  class inertia;
  class inertia_z;
  class eigenvector;
  class alpha_dihedrals;
  class alpha_angles;
  class dihedPC;
  class componentDisabled;
  class CartesianBasedPath;
  class gspath;
  class gzpath;
  class linearCombination;
  class CVBasedPath;
  class gspathCV;
  class gzpathCV;

  // non-scalar components
  class distance_vec;
  class distance_dir;
  class cartesian;
  class orientation;

protected:

  /// \brief Array of \link colvar::cvc \endlink objects
  std::vector<cvc *> cvcs;

  /// \brief Flags to enable or disable cvcs at next colvar evaluation
  std::vector<bool> cvc_flags;

  /// \brief Initialize the sorted list of atom IDs for atoms involved
  /// in all cvcs (called when enabling f_cv_collect_gradients)
  void build_atom_list(void);

private:
  /// Name of scripted function to be used
  std::string scripted_function;

  /// Current cvc values in the order requested by script
  /// when using scriptedFunction
  std::vector<const colvarvalue *> sorted_cvc_values;

#ifdef LEPTON
  /// Vector of evaluators for custom functions using Lepton
  std::vector<Lepton::CompiledExpression *> value_evaluators;

  /// Vector of evaluators for gradients of custom functions
  std::vector<Lepton::CompiledExpression *> gradient_evaluators;

  /// Vector of references to cvc values to be passed to Lepton evaluators
  std::vector<double *> value_eval_var_refs;
  std::vector<double *> grad_eval_var_refs;

  /// Unused value that is written to when a variable simplifies out of a Lepton expression
  double dev_null;
#endif

public:
  /// \brief Sorted array of (zero-based) IDs for all atoms involved
  std::vector<int> atom_ids;

  /// \brief Array of atomic gradients collected from all cvcs
  /// with appropriate components, rotations etc.
  /// For scalar variables only!
  std::vector<cvm::rvector> atomic_gradients;

  inline size_t n_components() const {
    return cvcs.size();
  }

  /// \brief Get vector of vectors of atom IDs for all atom groups
  virtual std::vector<std::vector<int> > get_atom_lists();
};

inline cvm::real const & colvar::force_constant() const
{
  return ext_force_k;
}

inline colvarvalue const & colvar::value() const
{
  return x_reported;
}

inline colvarvalue const & colvar::actual_value() const
{
  return x;
}

inline colvarvalue const & colvar::run_ave() const
{
  return runave;
}

inline colvarvalue const & colvar::velocity() const
{
  return v_reported;
}


inline colvarvalue const & colvar::total_force() const
{
  return ft_reported;
}


inline void colvar::add_bias_force(colvarvalue const &force)
{
  check_enabled(f_cv_gradient,
                std::string("applying a force to the variable \""+name+"\""));
  if (cvm::debug()) {
    cvm::log("Adding biasing force "+cvm::to_str(force)+" to colvar \""+name+"\".\n");
  }
  fb += force;
}


inline void colvar::add_bias_force_actual_value(colvarvalue const &force)
{
  if (cvm::debug()) {
    cvm::log("Adding biasing force "+cvm::to_str(force)+" to colvar \""+name+"\".\n");
  }
  fb_actual += force;
}


inline void colvar::reset_bias_force() {
  fb.type(value());
  fb.reset();
  fb_actual.type(value());
  fb_actual.reset();
}

#endif

