// -*- c++ -*-

#ifndef COLVARBIAS_H
#define COLVARBIAS_H

#include "colvar.h"
#include "colvarparse.h"
#include "colvardeps.h"


/// \brief Collective variable bias, base class
class colvarbias : public colvarparse, public cvm::deps {
public:

  /// Name of this bias
  std::string name;

  /// Type of this bias
  std::string bias_type;

  /// If there is more than one bias of this type, record its rank
  int rank;

  /// Add a new collective variable to this bias
  int add_colvar(std::string const &cv_name);

  /// Retrieve colvar values and calculate their biasing forces
  /// Return bias energy
  virtual int update();

  // TODO: move update_bias here (share with metadynamics)

  /// Load new configuration - force constant and/or centers only
  virtual void change_configuration(std::string const &conf);

  /// Calculate change in energy from using alternate configuration
  virtual cvm::real energy_difference(std::string const &conf);

  /// Give the total number of bins for a given bias.
  virtual int bin_num();
  /// Calculate the bin index for a given bias.
  virtual int current_bin();
  //// Give the count at a given bin index.
  virtual int bin_count(int bin_index);
  //// Share information between replicas, whatever it may be.
  virtual int replica_share();

  /// Perform analysis tasks
  virtual void analyze() {}

  /// Send forces to the collective variables
  void communicate_forces();

  /// \brief Constructor
  colvarbias(char const *key);

  /// \brief Parse config string and (re)initialize
  virtual int init(std::string const &conf);

  /// \brief Set to zero all mutable data
  virtual int reset();

protected:

  /// Default constructor
  colvarbias();

private:

  /// Copy constructor
  colvarbias(colvarbias &);

public:

  /// \brief Delete everything
  virtual int clear();

  /// Destructor
  virtual ~colvarbias();

  /// Read the bias configuration from a restart file
  virtual std::istream & read_restart(std::istream &is) = 0;

  /// Write the bias configuration to a restart file
  virtual std::ostream & write_restart(std::ostream &os) = 0;

  /// Write a label to the trajectory file (comment line)
  virtual std::ostream & write_traj_label(std::ostream &os);

  /// (Re)initialize the output files (does not write them yet)
  virtual int setup_output() { return COLVARS_OK; }

  /// Output quantities such as the bias energy to the trajectory file
  virtual std::ostream & write_traj(std::ostream &os);

  /// Write output files (if defined, e.g. in analysis mode)
  virtual int write_output_files()
  {
    return COLVARS_OK;
  }

  inline cvm::real get_energy() {
    return bias_energy;
  }

  /// \brief Implementation of the feature list for colvarbias
  static std::vector<feature *> cvb_features;

  /// \brief Implementation of the feature list accessor for colvarbias
  virtual std::vector<feature *> &features() {
    return cvb_features;
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
  /// (for history-dependent biases)
  bool                     has_data;

};

#endif
