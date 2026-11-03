// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#ifndef COLVARBIAS_RESTRAINT_H
#define COLVARBIAS_RESTRAINT_H

#include "colvarbias.h"

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4250) // Silence diamond inheritance warning
#endif

/// \brief Most general definition of a colvar restraint:
/// see derived classes for specific types
/// (implementation of \link colvarbias \endlink)
class colvarbias_restraint
  : public virtual colvarbias,
    public virtual colvarbias_ti
{

public:

  /// Retrieve colvar values and calculate their biasing forces
  virtual int update();

  /// Load new configuration - force constant and/or centers only
  virtual int change_configuration(std::string const & /* conf */) { return COLVARS_NOT_IMPLEMENTED; }

  /// Calculate change in energy from using alternate configuration
  virtual cvm::real energy_difference(std::string const & /* conf */) { return 0.0; }

  virtual std::string const get_state_params() const;
  virtual int set_state_params(std::string const &conf);
  virtual std::ostream & write_traj_label(std::ostream &os);
  virtual std::ostream & write_traj(std::ostream &os);

  /// \brief Constructor
  colvarbias_restraint(char const *key);

  virtual int init(std::string const &conf);
  virtual ~colvarbias_restraint();


protected:

  /// \brief Potential function for the i-th colvar
  virtual cvm::real restraint_potential(size_t i) const = 0;

  /// \brief Force function for the i-th colvar
  virtual colvarvalue const restraint_force(size_t i) const = 0;

  /// \brief Derivative of the potential function with respect to the force constant
  virtual cvm::real d_restraint_potential_dk(size_t i) const = 0;
};


/// Definition and parsing of the restraint centers
class colvarbias_restraint_centers
  : public virtual colvarbias_restraint
{
public:

  colvarbias_restraint_centers(char const *key);
  virtual int init(std::string const &conf);
  virtual int change_configuration(std::string const &conf);

protected:

  /// \brief Restraint centers
  std::vector<colvarvalue> colvar_centers;
};


/// Definition and parsing of the force constant
class colvarbias_restraint_k
  : public virtual colvarbias_restraint
{
public:

  colvarbias_restraint_k(char const *key);
  virtual int init(std::string const &conf);
  virtual int change_configuration(std::string const &conf);

protected:

  /// \brief Restraint force constant
  cvm::real force_k;

  /// \brief Whether the force constant should be positive
  bool check_positive_k;
};


/// Options to change the restraint configuration over time (shared between centers and k moving)
class colvarbias_restraint_moving
  : public virtual colvarbias_restraint,
    public virtual colvardeps {
public:

  colvarbias_restraint_moving(char const *key)
    : colvarbias_ti(key),
      colvarbias_restraint(key) {}
  // Note: despite the diamond inheritance, most of this function gets only executed once
  virtual int init(std::string const &conf) override;
  virtual int update() override;
  virtual int change_configuration(std::string const & /* conf */) override { return COLVARS_NOT_IMPLEMENTED; }

  virtual std::string const get_state_params() const override;
  virtual int set_state_params(std::string const &conf) override;

protected:

  /// \brief Moving target?
  bool b_chg_centers = false;

  /// \brief Changing force constant?
  bool b_chg_force_k = false;

  /// \brief Changing wall locations?
  bool b_chg_walls = false;

  /// @brief Update the force constant by interpolating between initial and target
  virtual void update_k(cvm::real /* lambda */) {}
  /// @brief Update the centers by interpolating between initial and target
  virtual void update_centers(cvm::real /* lambda */) {}
  /// @brief Update the walls by interpolating between initial and target
  virtual void update_walls(cvm::real /* lambda */) {}
  /// \brief Compute the derivative of the free energy wrt lambda due to changing k
  virtual cvm::real dU_dlambda_k() const { return 0.0; }
  /// \brief Compute the derivative of the free energy wrt lambda due to changing centers
  virtual cvm::real dU_dlambda_centers() const { return 0.0; }
  /// \brief Compute the derivative of the free energy wrt lambda due to changing walls
  virtual cvm::real dU_dlambda_walls() const { return 0.0; }

  /// \brief Perform decoupling of the restraint? If yes, lambda goes from 1 to 0
  bool b_decoupling;

  /// \brief Number of stages over which to perform the change
  /// If zero, perform a continuous change
  /// First step is ignored, then steps within each run of target_nsteps
  /// count for one stage.
  int target_nstages = 0;

  /// \brief Number of current stage of the perturbation
  /// Starts at 0, goes up to target_nstages during the perturbation
  int stage = 0;

  /// \brief Accumulating restraint FE derivative wrt lambda
  cvm::real dA_dlambda = 0.0;

  /// \brief Update the stage number based on the current step
  /// Note: this is idempotent so multiple calls are safe
  inline void update_stage() {
    stage = (cvmodule->step_absolute() - first_step) / target_nsteps;
    if (stage > target_nstages) {
      stage = target_nstages;
    }
  }

  /// \brief Get lambda value for the current stage
  cvm::real current_lambda() const;

  /// \brief Lambda-schedule for custom varying force constant
  std::vector<cvm::real> lambda_schedule;

  /// \brief Number of steps required to reach the target force constant
  /// or restraint centers
  cvm::step_number target_nsteps = 0L;

  /// \brief Equilibration steps for restraint FE calculation through TI
  cvm::step_number target_equil_steps = 0L;

  /// \brief Timestep at which the restraint starts moving
  cvm::step_number first_step = 0L;

  /// \brief Accumulated work (computed when outputAccumulatedWork == true)
  cvm::real acc_work = 0.0;
};


/// Options to change the restraint centers over time
class colvarbias_restraint_centers_moving
  : public virtual colvarbias_restraint_centers,
    public virtual colvarbias_restraint_moving
{
public:

  colvarbias_restraint_centers_moving(char const *key);
  virtual int init(std::string const &conf) override;
  virtual int change_configuration(std::string const & /* conf */) override { return COLVARS_NOT_IMPLEMENTED; }

  virtual std::string const get_state_params() const override;
  virtual int set_state_params(std::string const &conf) override;
  virtual std::ostream & write_traj_label(std::ostream &os) override;
  virtual std::ostream & write_traj(std::ostream &os) override;

protected:

  /// \brief New restraint centers
  std::vector<colvarvalue> target_centers;

  /// \brief Initial value of the restraint centers
  std::vector<colvarvalue> initial_centers;

  /// \brief Increment of the restraint centers at each step
  std::vector<colvarvalue> centers_incr;

  /// \brief Update the centers by interpolating between initial and target
  void update_centers(cvm::real lambda) override;

  /// \brief Compute the derivative of the free energy wrt lambda due to changing centers
  cvm::real dU_dlambda_centers() const override;

  /// Whether to write the current restraint centers to the trajectory file
  bool b_output_centers;

  /// Update the accumulated work
  int update_acc_work();
};


/// Options to change the restraint force constant over time
class colvarbias_restraint_k_moving
  : public virtual colvarbias_restraint_k,
    public virtual colvarbias_restraint_moving
{
public:
  // Do not initialize colvarbias here to avoid diamond inheritance issues
  // because this is an intermediate class
  colvarbias_restraint_k_moving(char const *key);
  virtual int init(std::string const &conf) override;
  virtual int change_configuration(std::string const & /* conf */) override { return COLVARS_NOT_IMPLEMENTED; }

  virtual std::string const get_state_params() const override;
  virtual int set_state_params(std::string const &conf) override;
  virtual std::ostream & write_traj_label(std::ostream &os) override;
  virtual std::ostream & write_traj(std::ostream &os) override;

protected:

  /// \brief Restraint force constant (target value at lambda = 1)
  cvm::real target_force_k;

  /// \brief Restraint force constant at lambda = 0
  /// (generally starting value, but is final value in the decoupling case)
  cvm::real starting_force_k;

  /// \brief Exponent for varying the force constant
  cvm::real lambda_exp;

  /// \brief Increment of the force constant at each step
  cvm::real force_k_incr;

  /// \brief Update the force constant by interpolating between initial and target
  void update_k(cvm::real lambda) override;

  /// \brief Compute the derivative of the free energy wrt lambda due to changing k
  cvm::real dU_dlambda_k() const override;

  /// Update the accumulated work
  int update_acc_work();
};


/// \brief Harmonic bias restraint
/// (implementation of \link colvarbias_restraint \endlink)
class colvarbias_restraint_harmonic
  : public colvarbias_restraint_centers_moving,
    public colvarbias_restraint_k_moving
{
public:
  colvarbias_restraint_harmonic(colvarmodule *cvmodule_in, char const *key);
  virtual int init(std::string const &conf);
  virtual int update();
  virtual std::string const get_state_params() const;
  virtual int set_state_params(std::string const &conf);
  virtual std::ostream & write_traj_label(std::ostream &os);
  virtual std::ostream & write_traj(std::ostream &os);
  virtual int change_configuration(std::string const &conf);
  virtual cvm::real energy_difference(std::string const &conf);

protected:

  virtual cvm::real restraint_potential(size_t i) const;
  virtual colvarvalue const restraint_force(size_t i) const;
  virtual cvm::real d_restraint_potential_dk(size_t i) const;
};


/// \brief Wall restraint
/// (implementation of \link colvarbias_restraint \endlink)
class colvarbias_restraint_harmonic_walls
  : public colvarbias_restraint_k_moving
{
public:

  colvarbias_restraint_harmonic_walls(colvarmodule *cvmodule_in, char const *key);
  virtual int init(std::string const &conf) override;
  virtual int update() override;
  virtual int update_acc_work();
  virtual std::string const get_state_params() const override;
  virtual int set_state_params(std::string const &conf) override;
  virtual int change_configuration(std::string const &conf) override;
  virtual std::ostream & write_traj_label(std::ostream &os) override;
  virtual std::ostream & write_traj(std::ostream &os) override;

protected:

  /// \brief Location of the lower walls
  std::vector<colvarvalue> lower_walls;

  /// \brief Location of the upper walls
  std::vector<colvarvalue> upper_walls;

  /// \brief If both walls are defined, use this k for the lower
  cvm::real lower_wall_k;

  /// \brief If both walls are defined, use this k for the upper
  cvm::real upper_wall_k;

  /// \brief Location of the target lower wall, for moving the wall
  std::vector<colvarvalue> target_lower_walls;

    /// \brief Initial value of the lower walls
  std::vector<colvarvalue> initial_lower_walls;

  /// \brief Increment of the lower walls at each step
  std::vector<colvarvalue> lower_walls_incr;

  /// \brief Location of the target upper wall, for moving the wall
  std::vector<colvarvalue> target_upper_walls;

    /// \brief Initial value of upper lower walls
  std::vector<colvarvalue> initial_upper_walls;

  /// \brief Increment of the upper walls at each step
  std::vector<colvarvalue> upper_walls_incr;

  /// \brief Signed distance to relevant wall (lower if below, upper if above)
  /// In PBC, only sees the closer of the two walls
  /// If negative, relates to lower wall - if positive, to upper wall
  cvm::real colvar_distance(size_t i) const;

  /// \brief Update the walls by interpolating between initial and target
  void update_walls(cvm::real lambda) override;

  /// \brief Compute the derivative of the free energy wrt lambda due to changing walls
  cvm::real dU_dlambda_walls() const override;

  virtual cvm::real restraint_potential(size_t i) const override;
  virtual colvarvalue const restraint_force(size_t i) const override;
  virtual cvm::real d_restraint_potential_dk(size_t i) const override;
};


/// \brief Linear bias restraint
/// (implementation of \link colvarbias_restraint \endlink)
class colvarbias_restraint_linear
  : public colvarbias_restraint_centers_moving,
    public colvarbias_restraint_k_moving
{

public:
  colvarbias_restraint_linear(colvarmodule *cvmodule_in, char const *key);
  virtual int init(std::string const &conf);
  virtual int update();
  virtual int change_configuration(std::string const &conf);
  virtual cvm::real energy_difference(std::string const &conf);

  virtual std::string const get_state_params() const;
  virtual int set_state_params(std::string const &conf);
  virtual std::ostream & write_traj_label(std::ostream &os);
  virtual std::ostream & write_traj(std::ostream &os);

protected:

  virtual cvm::real restraint_potential(size_t i) const;
  virtual colvarvalue const restraint_force(size_t i) const;
  virtual cvm::real d_restraint_potential_dk(size_t i) const;
};


/// Restrain the 1D histogram of a set of variables (or of a multidimensional one)
// TODO this could be reimplemented more cleanly as a derived class of both restraint and histogram
class colvarbias_restraint_histogram : public colvarbias {

public:

  colvarbias_restraint_histogram(colvarmodule *cvmodule_in, char const *key);
  int init(std::string const &conf);
  ~colvarbias_restraint_histogram();

  virtual int update();

  virtual int write_output_files();
  virtual std::ostream & write_traj_label(std::ostream &os);
  virtual std::ostream & write_traj(std::ostream &os);

protected:

  /// Probability density
  cvm::vector1d<cvm::real> p;

  /// Reference probability density
  cvm::vector1d<cvm::real> ref_p;

  /// Difference between probability density and reference
  cvm::vector1d<cvm::real> p_diff;

  /// Lower boundary of the grid
  cvm::real lower_boundary;

  /// Upper boundary of the grid
  cvm::real upper_boundary;

  /// Resolution of the grid
  cvm::real width;

  /// Width of the Gaussians
  cvm::real gaussian_width;

  /// Restraint force constant
  cvm::real force_k;

  /// Write the histogram to a file
  bool b_write_histogram;
};

#ifdef _MSC_VER
#pragma warning(pop)
#endif

#endif
