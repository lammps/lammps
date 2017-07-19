// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#ifndef COLVARSCRIPT_H
#define COLVARSCRIPT_H

#include <string>
#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvarbias.h"
#include "colvarproxy.h"

// Only these error values are part of the scripting interface
#define COLVARSCRIPT_ERROR -1
#define COLVARSCRIPT_OK 0

class colvarscript  {

private:
  colvarproxy *proxy;
  colvarmodule *colvars;

  inline colvarscript() {} // no-argument construction forbidden

public:

  friend class colvarproxy;

  colvarscript(colvarproxy * p);
  inline ~colvarscript() {}

  /// If an error is caught by the proxy through fatal_error(), this is set to COLVARSCRIPT_ERROR
  int proxy_error;

  /// If an error is returned by one of the methods, it should set this to the error message
  std::string result;

  /// Run script command with given positional arguments (objects)
  int run(int objc, unsigned char *const objv[]);

private:
  /// Run subcommands on colvar
  int proc_colvar(colvar *cv, int argc, unsigned char *const argv[]);

  /// Run subcommands on bias
  int proc_bias(colvarbias *b, int argc, unsigned char *const argv[]);

  /// Run subcommands on base colvardeps object (colvar, bias, ...)
  int proc_features(colvardeps *obj,
                    int argc, unsigned char *const argv[]);

  /// Build and return a short help
  std::string help_string(void);

public:

  inline char const *obj_to_str(unsigned char *const obj)
  {
    return cvm::proxy->script_obj_to_str(obj);
  }

};


#endif
