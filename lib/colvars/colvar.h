/// -*- c++ -*-

#ifndef COLVAR_H
#define COLVAR_H

#include <iostream>
#include <iomanip>
#include <cmath>

#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvarparse.h"


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
/// together, and must all be set to the same type of colvar::x by
/// using the colvarvalue::type() member function before using them
/// together in assignments or other operations; this is usually done
/// automatically in the constructor.  If you add a new member of
/// \link colvarvalue \endlink type, you should also add its
/// initialization line in the \link colvar \endlink constructor.

class colvar : public colvarparse {

public:

  /// Name
  std::string name;

  /// Type of value
  colvarvalue::Type type() const;

  /// \brief Current value (previously obtained from calc() or read_traj())
  colvarvalue const & value() const;

  /// \brief Current actual value (not extended DOF)
  colvarvalue const & actual_value() const;

  /// \brief Current velocity (previously obtained from calc() or read_traj())
  colvarvalue const & velocity() const;

  /// \brief Current system force (previously obtained from calc() or
  /// read_traj()).  Note: this is calculated using the atomic forces
  /// from the last simulation step.
  ///
  /// Total atom forces are read from the MD program, the total force
  /// acting on the collective variable is calculated summing those
  /// from all colvar components, the bias and walls forces are
  /// subtracted.
  colvarvalue const & system_force() const;

  /// \brief

  /// \brief Typical fluctuation amplitude for this collective
  /// variable (e.g. local width of a free energy basin)
  ///
  /// In metadynamics calculations, \link colvarbias_meta \endlink,
  /// this value is used to calculate the width of a gaussian.  In ABF
  /// calculations, \link colvarbias_abf \endlink, it is used to
  /// calculate the grid spacing in the direction of this collective
  /// variable.
  cvm::real width;

  /// \brief True if this \link colvar \endlink is a linear
  /// combination of \link cvc \endlink elements
  bool b_linear;

  /// \brief True if this \link colvar \endlink is equal to
  /// its only constituent cvc
  bool b_single_cvc;

  /// \brief True if all \link cvc \endlink objects are capable
  /// of calculating inverse gradients
  bool b_inverse_gradients;

  /// \brief True if all \link cvc \endlink objects are capable
  /// of calculating Jacobian forces
  bool b_Jacobian_force;

  /// \brief Options controlling the behaviour of the colvar during
  /// the simulation, which are set from outside the cvcs
  enum task {
    /// \brief Gradients are calculated and temporarily stored, so
    /// that external forces can be applied
    task_gradients,
    /// \brief Collect atomic gradient data from all cvcs into vector
    /// atomic_gradients
    task_collect_gradients,
    /// \brief Calculate the velocity with finite differences
    task_fdiff_velocity,
    /// \brief The system force is calculated, projecting the atomic
    /// forces on the inverse gradients
    task_system_force,
    /// \brief The variable has a harmonic restraint around a moving
    /// center with fictitious mass; bias forces will be applied to
    /// the center
    task_extended_lagrangian,
    /// \brief The extended system coordinate undergoes Langevin
    /// dynamics
    task_langevin,
    /// \brief Output the potential and kinetic energies
    /// (for extended Lagrangian colvars only)
    task_output_energy,
    /// \brief Compute analytically the "force" arising from the
    /// geometric entropy component (for example, that from the angular
    /// states orthogonal to a distance vector)
    task_Jacobian_force,
    /// \brief Report the Jacobian force as part of the system force
    /// if disabled, apply a correction internally to cancel it
    task_report_Jacobian_force,
    /// \brief Output the value to the trajectory file (on by default)
    task_output_value,
    /// \brief Output the velocity to the trajectory file
    task_output_velocity,
    /// \brief Output the applied force to the trajectory file
    task_output_applied_force,
    /// \brief Output the system force to the trajectory file
    task_output_system_force,
    /// \brief A lower boundary is defined
    task_lower_boundary,
    /// \brief An upper boundary is defined
    task_upper_boundary,
    /// \brief Provide a discretization of the values of the colvar to
    /// be used by the biases or in analysis (needs lower and upper
    /// boundary)
    task_grid,
    /// \brief Apply a restraining potential (|x-xb|^2) to the colvar
    /// when it goes below the lower wall
    task_lower_wall,
    /// \brief Apply a restraining potential (|x-xb|^2) to the colvar
    /// when it goes above the upper wall
    task_upper_wall,
    /// \brief Compute running average
    task_runave,
    /// \brief Compute time correlation function
    task_corrfunc,
    /// \brief Value and gradients computed by user script
    task_scripted,
    /// \brief Number of possible tasks
    task_ntot
  };

  /// Tasks performed by this colvar
  bool tasks[task_ntot];

  /// List of biases that depend on this colvar
  std::vector<colvarbias *> biases;
protected:


  /*
    extended:
    H = H_{0} + \sum_{i} 1/2*\lambda*(S_i(x(t))-s_i(t))^2 \\
    + \sum_{i} 1/2*m_i*(ds_i(t)/dt)^2 \\
    + \sum_{t'<t} W * exp (-1/2*\sum_{i} (s_i(t')-s_i(t))^2/(\delta{}s_i)^2) \\
    + \sum_{w} (\sum_{i}c_{w,i}s_i(t) - D_w)^M/(\Sigma_w)^M

    normal:
    H = H_{0} + \sum_{t'<t} W * exp (-1/2*\sum_{i} (S_i(x(t'))-S_i(x(t)))^2/(\delta{}S_i)^2) \\
    + \sum_{w} (\sum_{i}c_{w,i}S_i(t) - D_w)^M/(\Sigma_w)^M

    output:
    H = H_{0}   (only output S(x), no forces)

    Here:
    S(x(t)) = x
    s(t)    = xr
    DS = Ds = delta
  */


  /// Value of the colvar
  colvarvalue x;

  /// Cached reported value (x may be manipulated)
  colvarvalue x_reported;

  /// Finite-difference velocity
  colvarvalue v_fdiff;

  inline colvarvalue fdiff_velocity (colvarvalue const &xold,
                                     colvarvalue const &xnew)
  {
    // using the gradient of the square distance to calculate the
    // velocity (non-scalar variables automatically taken into
    // account)
    cvm::real dt = cvm::dt();
    return ( ( (dt > 0.0) ? (1.0/dt) : 1.0 ) *
             0.5 * dist2_lgrad (xnew, xold) );
  }

  /// Cached reported velocity
  colvarvalue v_reported;

  // Options for task_extended_lagrangian
  /// Restraint center
  colvarvalue xr;
  /// Velocity of the restraint center
  colvarvalue vr;
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

  /// \brief Jacobian force, when task_Jacobian_force is enabled
  colvarvalue fj;

  /// Cached reported system force
  colvarvalue ft_reported;

public:


  /// \brief Bias force; reset_bias_force() should be called before
  /// the biases are updated
  colvarvalue fb;

  /// \brief Total \em applied force; fr (if task_extended_lagrangian
  /// is defined), fb (if biases are applied) and the walls' forces
  /// (if defined) contribute to it
  colvarvalue f;

  /// \brief Total force, as derived from the atomic trajectory;
  /// should equal the total system force plus \link f \endlink
  colvarvalue ft;


  /// Period, if it is a constant
  cvm::real period;

  /// \brief Same as above, but also takes into account components
  /// with a variable period, such as distanceZ
  bool b_periodic;


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
  bool periodic_boundaries (colvarvalue const &lb, colvarvalue const &ub) const;


  /// Constructor
  colvar (std::string const &conf);

  /// Enable the specified task
  int enable (colvar::task const &t);

  /// Disable the specified task
  void disable (colvar::task const &t);

  /// Get ready for a run and possibly re-initialize internal data
  void setup();

  /// Destructor
  ~colvar();


  /// \brief Calculate the colvar value and all the other requested
  /// quantities
  void calc();

  /// \brief Calculate just the value, and store it in the argument
  void calc_value (colvarvalue &xn);

  /// Set the total biasing force to zero
  void reset_bias_force();

  /// Add to the total force from biases
  void add_bias_force (colvarvalue const &force);

  /// \brief Collect all forces on this colvar, integrate internal
  /// equations of motion of internal degrees of freedom; see also
  /// colvar::communicate_forces()
  /// return colvar energy if extended Lagrandian active
  cvm::real update();

  /// \brief Communicate forces (previously calculated in
  /// colvar::update()) to the external degrees of freedom
  void communicate_forces();


  /// \brief Use the internal metrics (as from \link cvc
  /// \endlink objects) to calculate square distances and gradients
  ///
  /// Handles correctly symmetries and periodic boundary conditions
  cvm::real dist2 (colvarvalue const &x1,
                   colvarvalue const &x2) const;

  /// \brief Use the internal metrics (as from \link cvc
  /// \endlink objects) to calculate square distances and gradients
  ///
  /// Handles correctly symmetries and periodic boundary conditions
  colvarvalue dist2_lgrad (colvarvalue const &x1,
                           colvarvalue const &x2) const;

  /// \brief Use the internal metrics (as from \link cvc
  /// \endlink objects) to calculate square distances and gradients
  ///
  /// Handles correctly symmetries and periodic boundary conditions
  colvarvalue dist2_rgrad (colvarvalue const &x1,
                           colvarvalue const &x2) const;

  /// \brief Use the internal metrics (as from \link cvc
  /// \endlink objects) to wrap a value into a standard interval
  ///
  /// Handles correctly symmetries and periodic boundary conditions
  void wrap (colvarvalue &x) const;


  /// Read the analysis tasks
  int parse_analysis (std::string const &conf);
  /// Perform analysis tasks
  void analyse();


  /// Read the value from a collective variable trajectory file
  std::istream & read_traj (std::istream &is);
  /// Output formatted values to the trajectory file
  std::ostream & write_traj (std::ostream &os);
  /// Write a label to the trajectory file (comment line)
  std::ostream & write_traj_label (std::ostream &os);


  /// Read the collective variable from a restart file
  std::istream & read_restart (std::istream &is);
  /// Write the collective variable to a restart file
  std::ostream & write_restart (std::ostream &os);

  /// Write output files (if defined, e.g. in analysis mode)
  int write_output_files();


protected:

  /// Previous value (to calculate velocities during analysis)
  colvarvalue            x_old;

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
  int calc_vel_acf (std::list<colvarvalue> &v_history,
                     colvarvalue const      &v);

  /// \brief Coordinate ACF, scalar product between x(0) and x(t)
  /// (does not work with scalar numbers)
  void calc_coor_acf (std::list<colvarvalue> &x_history,
                      colvarvalue const      &x);

  /// \brief Coordinate ACF, second order Legendre polynomial between
  /// x(0) and x(t) (does not work with scalar numbers)
  void calc_p2coor_acf (std::list<colvarvalue> &x_history,
                        colvarvalue const      &x);

  /// Calculate the auto-correlation function (ACF)
  int calc_acf();
  /// Save the ACF to a file
  void write_acf (std::ostream &os);

  /// Length of running average series
  size_t         runave_length;
  /// Timesteps to skip between two values in the running average series
  size_t         runave_stride;
  /// Name of the file to write the running average
  std::ofstream  runave_os;
  /// Current value of the running average
  colvarvalue    runave;
  /// Current value of the square deviation from the running average
  cvm::real      runave_variance;

  /// Calculate the running average and its standard deviation
  void calc_runave();

  /// If extended Lagrangian active: colvar energies (kinetic and harmonic potential)
  cvm::real kinetic_energy;
  cvm::real potential_energy;
public:


  // collective variable component base class
  class cvc;

  // currently available collective variable components

  // scalar colvar components
  class distance;
  class distance_z;
  class distance_xy;
  class distance_inv;
  class angle;
  class dihedral;
  class coordnum;
  class selfcoordnum;
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

  // non-scalar components
  class distance_vec;
  class distance_dir;
  class orientation;

protected:

  /// \brief Array of \link cvc \endlink objects
  std::vector<cvc *> cvcs;

  /// \brief Initialize the sorted list of atom IDs for atoms involved
  /// in all cvcs (called when enabling task_collect_gradients)
  void build_atom_list (void);

private:
  /// Name of scripted function to be used
  std::string scripted_function;

  /// Current cvc values in the order requested by script
  /// when using scriptedFunction
  std::vector<const colvarvalue *> sorted_cvc_values;

public:
  /// \brief Sorted array of (zero-based) IDs for all atoms involved
  std::vector<int> atom_ids;

  /// \brief Array of atomic gradients collected from all cvcs
  /// with appropriate components, rotations etc.
  /// For scalar variables only!
  std::vector<cvm::rvector> atomic_gradients;

  inline size_t n_components () const {
    return cvcs.size();
  }
};

inline colvarvalue::Type colvar::type() const
{
  return x.type();
}


inline colvarvalue const & colvar::value() const
{
  return x_reported;
}


inline colvarvalue const & colvar::actual_value() const
{
  return x;
}


inline colvarvalue const & colvar::velocity() const
{
  return v_reported;
}


inline colvarvalue const & colvar::system_force() const
{
  return ft_reported;
}


inline void colvar::add_bias_force (colvarvalue const &force)
{
  fb += force;
}


inline void colvar::reset_bias_force() {
  fb.reset();
}

#endif

