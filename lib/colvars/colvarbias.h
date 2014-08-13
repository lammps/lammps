#ifndef COLVARBIAS_H
#define COLVARBIAS_H

#include "colvar.h"
#include "colvarparse.h"


/// \brief Collective variable bias, base class
class colvarbias : public colvarparse {
public:

  /// Numeric id of this bias
  int            id;

  /// Name of this bias
  std::string    name;

  /// Add a new collective variable to this bias
  void add_colvar (std::string const &cv_name);

  /// Retrieve colvar values and calculate their biasing forces
  /// Return bias energy
  virtual cvm::real update() = 0;

  /// Load new configuration - force constant and/or centers only
  virtual void change_configuration(std::string const &conf);

  /// Calculate change in energy from using alternate configuration
  virtual cvm::real energy_difference(std::string const &conf);

  /// Perform analysis tasks
  virtual inline void analyse() {}

  /// Send forces to the collective variables
  void communicate_forces();

  /// \brief Constructor
  ///
  /// The constructor of the colvarbias base class is protected, so
  /// that it can only be called from inherited classes
  colvarbias (std::string const &conf, char const *key);

  /// Default constructor
  colvarbias();

  /// Destructor
  virtual inline ~colvarbias() {}

  /// Read the bias configuration from a restart file
  virtual std::istream & read_restart (std::istream &is) = 0;

  /// Write the bias configuration to a restart file
  virtual std::ostream & write_restart (std::ostream &os) = 0;

  /// Write a label to the trajectory file (comment line)
  virtual std::ostream & write_traj_label (std::ostream &os);

  /// Output quantities such as the bias energy to the trajectory file
  virtual std::ostream & write_traj (std::ostream &os);


protected:

  /// \brief Pointers to collective variables to which the bias is
  /// applied; current values and metric functions will be obtained
  /// through each colvar object
  std::vector<colvar *>    colvars;

  /// \brief Current forces from this bias to the colvars
  std::vector<colvarvalue> colvar_forces;

  /// \brief Current energy of this bias (colvar_forces should be obtained by deriving this)
  cvm::real                bias_energy;

  /// Whether to write the current bias energy from this bias to the trajectory file
  bool                     b_output_energy;

  /// \brief Whether this bias has already accumulated information
  /// (when relevant)
  bool                     has_data;

};


/// \brief Bias restraint, optionally moving towards a target
/// (implementation of \link colvarbias \endlink)
class colvarbias_restraint : public colvarbias {

public:

  /// Retrieve colvar values and calculate their biasing forces
  virtual cvm::real update();

  /// Load new configuration - force constant and/or centers only
  virtual void change_configuration(std::string const &conf);

  /// Calculate change in energy from using alternate configuration
  virtual cvm::real energy_difference(std::string const &conf);

  /// Read the bias configuration from a restart file
  virtual std::istream & read_restart (std::istream &is);

  /// Write the bias configuration to a restart file
  virtual std::ostream & write_restart (std::ostream &os);

  /// Write a label to the trajectory file (comment line)
  virtual std::ostream & write_traj_label (std::ostream &os);

  /// Output quantities such as the bias energy to the trajectory file
  virtual std::ostream & write_traj (std::ostream &os);

  /// \brief Constructor
  colvarbias_restraint (std::string const &conf, char const *key);

  /// Destructor
  virtual ~colvarbias_restraint();


protected:

  /// \brief Potential function
  virtual cvm::real restraint_potential(cvm::real k, colvar* x, const colvarvalue& xcenter) const = 0;

  /// \brief Force function
  virtual colvarvalue restraint_force(cvm::real k, colvar* x, const colvarvalue& xcenter) const = 0;

  ///\brief Unit scaling
  virtual cvm::real restraint_convert_k(cvm::real k, cvm::real dist_measure) const = 0;

  /// \brief Restraint centers
  std::vector<colvarvalue> colvar_centers;

  /// \brief Restraint centers without wrapping or constraints applied
  std::vector<colvarvalue> colvar_centers_raw;

  /// \brief Moving target?
  bool b_chg_centers;

  /// \brief New restraint centers
  std::vector<colvarvalue> target_centers;

  /// \brief Amplitude of the restraint centers' increment at each step
  /// (or stage) towards the new values (calculated from target_nsteps)
  std::vector<colvarvalue> centers_incr;

  /// Whether to write the current restraint centers to the trajectory file
  bool b_output_centers;

  /// Whether to write the current accumulated work to the trajectory file
  bool b_output_acc_work;

  /// \brief Accumulated work
  cvm::real acc_work;

  /// \brief Restraint force constant
  cvm::real force_k;

  /// \brief Changing force constant?
  bool b_chg_force_k;

  /// \brief Restraint force constant (target value)
  cvm::real target_force_k;

  /// \brief Restraint force constant (starting value)
  cvm::real starting_force_k;

  /// \brief Lambda-schedule for custom varying force constant
  std::vector<cvm::real> lambda_schedule;

  /// \brief Exponent for varying the force constant
  cvm::real force_k_exp;

  /// \brief Intermediate quantity to compute the restraint free energy
  /// (in TI, would be the accumulating FE derivative)
  cvm::real restraint_FE;


  /// \brief Equilibration steps for restraint FE calculation through TI
  cvm::real target_equil_steps;

  /// \brief Number of stages over which to perform the change
  /// If zero, perform a continuous change
  int target_nstages;

  /// \brief Number of current stage of the perturbation
  int stage;

  /// \brief Number of steps required to reach the target force constant
  /// or restraint centers
  size_t target_nsteps;
};

/// \brief Harmonic bias restraint
/// (implementation of \link colvarbias_restraint \endlink)
class colvarbias_restraint_harmonic : public colvarbias_restraint {

public:
  colvarbias_restraint_harmonic(std::string const &conf, char const *key);

protected: /// \brief Potential function
  virtual cvm::real restraint_potential(cvm::real k,  colvar*  x, const colvarvalue& xcenter) const;

  /// \brief Force function
  virtual colvarvalue restraint_force(cvm::real k,  colvar* x,  const colvarvalue& xcenter) const;

  ///\brief Unit scaling
  virtual cvm::real restraint_convert_k(cvm::real k, cvm::real dist_measure) const;

};

/// \brief Linear bias restraint
/// (implementation of \link colvarbias_restraint \endlink)
class colvarbias_restraint_linear : public colvarbias_restraint {

public:
  colvarbias_restraint_linear(std::string const &conf, char const *key);

protected: /// \brief Potential function
  virtual cvm::real restraint_potential(cvm::real k,  colvar*  x, const colvarvalue& xcenter) const;

  /// \brief Force function
  virtual colvarvalue restraint_force(cvm::real k,  colvar* x,  const colvarvalue& xcenter) const;

  ///\brief Unit scaling
  virtual cvm::real restraint_convert_k(cvm::real k, cvm::real dist_measure) const;

};


#endif



// Emacs
// Local Variables:
// mode: C++
// End:
