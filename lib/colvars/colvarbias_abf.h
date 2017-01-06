// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#ifndef COLVARBIAS_ABF_H
#define COLVARBIAS_ABF_H

#include <vector>
#include <list>
#include <sstream>
#include <iomanip>

#include "colvarbias.h"
#include "colvargrid.h"

typedef cvm::real* gradient_t;


/// ABF bias
class colvarbias_abf : public colvarbias {

public:

  colvarbias_abf(char const *key);
  virtual int init(std::string const &conf);
  virtual ~colvarbias_abf();
  virtual int update();

private:

  /// Filename prefix for human-readable gradient/sample count output
  std::string  output_prefix;

  /// Base filename(s) for reading previous gradient data (replaces data from restart file)
  std::vector<std::string> input_prefix;

  bool		update_bias;
  bool		hide_Jacobian;
  size_t	full_samples;
  size_t	min_samples;
  /// frequency for updating output files
  int		output_freq;
  /// Write combined files with a history of all output data?
  bool      b_history_files;
  /// Write CZAR output file for stratified eABF (.zgrad)
  bool      b_czar_window_file;
  size_t    history_freq;

  /// Cap applied biasing force?
  bool                    cap_force;
  std::vector<cvm::real>  max_force;

  // Internal data and methods

  std::vector<int>  bin, force_bin, z_bin;
  gradient_t    system_force, applied_force;

  /// n-dim grid of free energy gradients
  colvar_grid_gradient  *gradients;
  /// n-dim grid of number of samples
  colvar_grid_count     *samples;
  /// n-dim grid: average force on "real" coordinate for eABF z-based estimator
  colvar_grid_gradient  *z_gradients;
  /// n-dim grid of number of samples on "real" coordinate for eABF z-based estimator
  colvar_grid_count     *z_samples;
  /// n-dim grid contining CZAR estimator of "real" free energy gradients
  colvar_grid_gradient  *czar_gradients;

  inline int update_system_force(size_t i)
  {
    if (colvars[i]->is_enabled(f_cv_subtract_applied_force)) {
      // this colvar is already subtracting the ABF force
      system_force[i] = colvars[i]->total_force().real_value;
    } else {
      system_force[i] = colvars[i]->total_force().real_value
        - colvar_forces[i].real_value;
    }
    if (cvm::debug())
      cvm::log("ABF System force calc: cv " + cvm::to_str(i) +
               " fs " + cvm::to_str(system_force[i]) +
               " = ft " + cvm::to_str(colvars[i]->total_force().real_value) +
               " - fa " + cvm::to_str(colvar_forces[i].real_value));
    return COLVARS_OK;
  }

  // shared ABF
  bool     shared_on;
  size_t   shared_freq;
  int   shared_last_step;
  // Share between replicas -- may be called independently of update
  virtual int replica_share();

  // Store the last set for shared ABF
  colvar_grid_gradient  *last_gradients;
  colvar_grid_count     *last_samples;

  // For Tcl implementation of selection rules.
  /// Give the total number of bins for a given bias.
  virtual int bin_num();
  /// Calculate the bin index for a given bias.
  virtual int current_bin();
  //// Give the count at a given bin index.
  virtual int bin_count(int bin_index);

  /// Write human-readable FE gradients and sample count
  void		  write_gradients_samples(const std::string &prefix, bool append = false);
  void		  write_last_gradients_samples(const std::string &prefix, bool append = false);

  /// Read human-readable FE gradients and sample count (if not using restart)
  void		  read_gradients_samples();

  std::istream& read_state_data(std::istream&);
  std::ostream& write_state_data(std::ostream&);
};

#endif

