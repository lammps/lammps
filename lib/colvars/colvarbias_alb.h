#ifndef COLVARBIAS_ALB_H
#define COLVARBIAS_ALB_H

#include "colvar.h"
#include "colvarbias.h"

class colvarbias_alb : public colvarbias {

public:
  colvarbias_alb(std::string const &conf, char const *key);
  virtual ~colvarbias_alb();


  virtual cvm::real update();

  /// Read the bias configuration from a restart file
  virtual std::istream & read_restart (std::istream &is);

  /// Write the bias configuration to a restart file
  virtual std::ostream & write_restart (std::ostream &os);

  /// Write a label to the trajectory file (comment line)
  virtual std::ostream & write_traj_label (std::ostream &os);

  /// Output quantities such as the bias energy to the trajectory file
  virtual std::ostream & write_traj (std::ostream &os);

protected:

  /// \brief Restraint centers
  std::vector<colvarvalue> colvar_centers;

  /// \brief colvar parameters, used for calculating the gradient/variance
  std::vector<cvm::real> means;
  std::vector<cvm::real> ssd; // SSD = sum of squares of differences from mean
  int update_calls;

  ///\brief how often to update coupling constant
  int update_freq;

  ///\brief Estimated range of coupling constant values in kT
  std::vector<cvm::real> max_coupling_range;

  //\brief Estimated max for how quickly the rate can change in kT / time
  std::vector<cvm::real> max_coupling_rate;

  /// \brief accumated couping force; used in stochastic online gradient descent algorithm
  std::vector<cvm::real> coupling_accum;

  /// \brief coupling constant
  std::vector<cvm::real> set_coupling;

  /// \brief current coupling constant, which is ramped up during equilibration to coupling
  std::vector<cvm::real> current_coupling;

  /// \brief how quickly to change the coupling constant
  std::vector<cvm::real> coupling_rate;

  // \brief if we're equilibrating our estimates or collecting data
  bool b_equilibration;

  // \brief If the coupling range should be increased
  bool b_hard_coupling_range;


  /// \brief flag for outputting colvar centers
  bool b_output_centers;

  /// \brief flag for outputting current gradient
  bool b_output_grad;

  /// \brief flag for outputting coupling constant
  bool b_output_coupling;

  cvm::real restraint_potential(cvm::real k,  const colvar*  x, const colvarvalue& xcenter) const;

  /// \brief Force function
  colvarvalue restraint_force(cvm::real k,  const colvar* x,  const colvarvalue& xcenter) const;

  ///\brief Unit scaling
  cvm::real restraint_convert_k(cvm::real k, cvm::real dist_measure) const;

};

#endif
