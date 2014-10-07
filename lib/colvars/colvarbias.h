/// -*- c++ -*-

#ifndef COLVARBIAS_H
#define COLVARBIAS_H

#include "colvar.h"
#include "colvarparse.h"


/// \brief Collective variable bias, base class
class colvarbias : public colvarparse {
public:

  /// Name of this bias
  std::string    name;

  /// Add a new collective variable to this bias
  void add_colvar (std::string const &cv_name);

  /// Retrieve colvar values and calculate their biasing forces
  /// Return bias energy
  virtual cvm::real update() = 0;

  // TODO: move update_bias here (share with metadynamics)

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
  virtual ~colvarbias();

  /// Read the bias configuration from a restart file
  virtual std::istream & read_restart (std::istream &is) = 0;

  /// Write the bias configuration to a restart file
  virtual std::ostream & write_restart (std::ostream &os) = 0;

  /// Write a label to the trajectory file (comment line)
  virtual std::ostream & write_traj_label (std::ostream &os);

  /// Output quantities such as the bias energy to the trajectory file
  virtual std::ostream & write_traj (std::ostream &os);

  inline cvm::real get_energy () {
    return bias_energy;
  }
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

#endif
