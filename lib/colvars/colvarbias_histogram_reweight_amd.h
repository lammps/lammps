// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#ifndef COLVARBIAS_HISTOGRAM_REWEIGHT_AMD
#define COLVARBIAS_HISTOGRAM_REWEIGHT_AMD

#include "colvarbias_histogram.h"

/// Reweighted histogram for accelerated molecular dynamics (aMD) or
/// Gaussian aMD (GaMD)
class colvarbias_reweightaMD : public colvarbias_histogram {
public:
  colvarbias_reweightaMD(char const *key);
  virtual ~colvarbias_reweightaMD();
#if (__cplusplus >= 201103L)
  virtual int init(std::string const &conf) override;
  virtual int update() override;
  virtual int write_output_files() override;
#else
  virtual int init(std::string const &conf);
  virtual int update();
  virtual int write_output_files();
#endif

  /// @brief convert histogram to PMF by taking logarithm and multiplying
  ///        it with -1/beta
  /// @param[in,out] hist the origin histogram and also the output PMF
  /// @param[in] hist_count the sampling or biased histogram
  void hist_to_pmf(
    colvar_grid_scalar* hist,
    const colvar_grid_scalar* hist_count) const;

  /// @brief calculate the cumulant expansion to second order
  /// @param[in] hist_dV the histogram of the boosting potential, ΔV(ξ)
  /// @param[in] hist_dV_square the histogram of the square of boosting
  ///                           potential
  /// @param[in] hist_count the sampling or biased histogram
  /// @param[out] cumulant_expansion_factor the factor of the cumulant
  ///                                       expansion to second order
  void compute_cumulant_expansion_factor(
    const colvar_grid_scalar* hist_dV,
    const colvar_grid_scalar* hist_dV_square,
    const colvar_grid_scalar* hist_count,
    colvar_grid_scalar* cumulant_expansion_factor) const;

  /// @brief output the PMF by the exponential average estimator
  /// @param[in] p_output_prefix the prefix of the output file
  /// @param[in] keep_open Allow writing the history of the PMF
  virtual int write_exponential_reweighted_pmf(
    const std::string& p_output_prefix, bool keep_open = false);

  /// @brief output the PMF by the cumulant expansion estimator
  /// @param[in] p_output_prefix the prefix of the output file
  /// @param[in] keep_open Allow writing the history of the expansion
  virtual int write_cumulant_expansion_pmf(
    const std::string& p_output_prefix, bool keep_open = false);

  /// @brief output the biased sampling
  /// @param[in] p_output_prefix the prefix of the output file
  /// @param[in] keep_open Allow writing the history of the samples
  virtual int write_count(
    const std::string& p_output_prefix, bool keep_open = false);
protected:
  /// Current accelMD factor is the from previous frame
  std::vector<int> previous_bin;
  /// Start collecting samples after N steps
  colvarmodule::step_number start_after_steps;

  /// Use cumulant expansion to second order?
  bool b_use_cumulant_expansion;
  colvar_grid_scalar* grid_count;
  colvar_grid_scalar* grid_dV;
  colvar_grid_scalar* grid_dV_square;

  /// Number of timesteps between recording data in history files (if non-zero)
  size_t history_freq;
  bool b_history_files;

  /// Write gradients of the PMF?
  bool b_write_gradients;

  /// save and restore
#if (__cplusplus >= 201103L)
  virtual std::istream & read_state_data(std::istream &is) override;
  virtual std::ostream & write_state_data(std::ostream &os) override;
#else
  virtual std::istream & read_state_data(std::istream &is);
  virtual std::ostream & write_state_data(std::ostream &os);
#endif
private:
  /// temporary grids for evaluating PMFs
  colvar_grid_scalar  *pmf_grid_exp_avg;
  colvar_grid_scalar  *pmf_grid_cumulant;
  colvar_grid_gradient *grad_grid_exp_avg;
  colvar_grid_gradient *grad_grid_cumulant;
};

#endif // COLVARBIAS_HISTOGRAM_REWEIGHT_AMD
