// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#ifndef COLVARPROXY_TCL_H
#define COLVARPROXY_TCL_H

#if defined(NAMD_TCL) || defined(VMDTCL)
#define COLVARS_TCL
#endif

#ifdef COLVARS_TCL
#include <tcl.h>
#else
// Allow for placeholders Tcl_Interp* variables
typedef void Tcl_Interp;
#endif

#include <vector>


/// Methods for using Tcl within Colvars
class colvarproxy_tcl {

public:

  /// Constructor
  colvarproxy_tcl();

  /// Destructor
  virtual ~colvarproxy_tcl();

  /// Is Tcl available? (trigger initialization if needed)
  inline bool tcl_available() {
#if defined(COLVARS_TCL)
    return true;
#else
    return false;
#endif
  }

  /// Get a string representation of the Tcl object pointed to by obj
  char const *tcl_get_str(void *obj);

  int tcl_run_script(std::string const &script);

  int tcl_run_file(std::string const &fileName);

  /// Tcl implementation of run_force_callback()
  int tcl_run_force_callback();

  /// Tcl implementation of run_colvar_callback()
  int tcl_run_colvar_callback(
              std::string const &name,
              std::vector<const colvarvalue *> const &cvcs,
              colvarvalue &value);

  /// Tcl implementation of run_colvar_gradient_callback()
  int tcl_run_colvar_gradient_callback(
              std::string const &name,
              std::vector<const colvarvalue *> const &cvcs,
              std::vector<cvm::matrix2d<cvm::real> > &gradient);

  /// Get a pointer to the Tcl interpreter
  inline Tcl_Interp *get_tcl_interp()
  {
    return tcl_interp_;
  }

  /// Set the pointer to the Tcl interpreter
  inline void set_tcl_interp(Tcl_Interp *interp)
  {
    tcl_interp_ = interp;
  }

  /// Set Tcl pointers
  virtual void init_tcl_pointers();

protected:
  /// Pointer to Tcl interpreter object
  Tcl_Interp *tcl_interp_;
};

#endif
