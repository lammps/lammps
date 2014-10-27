/************************************************************************
 * Headers for the ABF and histogram biases                             *
 ************************************************************************/

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

  colvarbias_abf(std::string const &conf, char const *key);
  ~colvarbias_abf();

  cvm::real update();

private:

  /// Filename prefix for human-readable gradient/sample count output
  std::string	output_prefix;

  /// Base filename(s) for reading previous gradient data (replaces data from restart file)
  std::vector<std::string> input_prefix;

  bool		apply_bias;
  bool		update_bias;
  bool		hide_Jacobian;
  size_t	full_samples;
  size_t	min_samples;
  /// frequency for updating output files (default: same as restartFreq?)
  int		output_freq;
  /// Write combined files with a history of all output data?
  bool      b_history_files;
  size_t    history_freq;

  /// Cap applied biasing force?
  bool                    cap_force;
  std::vector<cvm::real>  max_force;

  // Internal data and methods

  std::vector<int>  bin, force_bin;
  gradient_t	    force;

  /// n-dim grid of free energy gradients
  colvar_grid_gradient  *gradients;
  /// n-dim grid of number of samples
  colvar_grid_count     *samples;

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

  std::istream& read_restart(std::istream&);
  std::ostream& write_restart(std::ostream&);
};


/// Histogram "bias" (does as the name says)
class colvarbias_histogram : public colvarbias {

public:

  colvarbias_histogram(std::string const &conf, char const *key);
  ~colvarbias_histogram();

  cvm::real update();

private:

  /// n-dim histogram
  colvar_grid_count    *grid;
  std::vector<int>  bin;
  std::string	  out_name;

  int		  output_freq;
  void		  write_grid();
  std::ofstream	  grid_os;  /// Stream for writing grid to disk

  std::istream& read_restart(std::istream&);
  std::ostream& write_restart(std::ostream&);
};

#endif
