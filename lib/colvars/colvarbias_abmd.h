// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#ifndef COLVARBIAS_ABMD_H
#define COLVARBIAS_ABMD_H

#include "colvarbias_restraint.h"


/// \brief Adiabatic Bias MD
class colvarbias_abmd
  : public colvarbias_ti
{
public:

  colvarbias_abmd(char const *key);
  virtual int init(std::string const &conf);
  virtual int update();
  virtual std::string const get_state_params() const;
  virtual int set_state_params(std::string const &conf);
  virtual std::ostream & write_traj_label(std::ostream &os);
  virtual std::ostream & write_traj(std::ostream &os);

protected:

  /// \brief Location of the moving wall
  cvm::real ref_val = 0.;

  /// \brief Has ref_val already been set?
  bool ref_initialized = false;

  /// \brief Value of the reference where it stops moving
  cvm::real stopping_val = 0.;

  /// \brief Is the target moving down?
  bool decreasing = false;

  /// \brief Restraint force constant
  cvm::real k = 0.;
};


#endif
