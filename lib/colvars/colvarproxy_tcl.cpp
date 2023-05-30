// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include <sstream>
#include <iostream>

#include "colvarmodule.h"
#include "colvarproxy.h"
#include "colvarproxy_tcl.h"
#include "colvaratoms.h"

#ifdef COLVARS_TCL
#include <tcl.h>
#endif



colvarproxy_tcl::colvarproxy_tcl()
{
  tcl_interp_ = NULL;
}


colvarproxy_tcl::~colvarproxy_tcl()
{
}


void colvarproxy_tcl::init_tcl_pointers()
{
  // This is overloaded by NAMD and VMD proxies to use the local interpreters
#if defined(COLVARS_TCL)
  if (tcl_interp_ == NULL) {
    // Allocate a dedicated Tcl interpreter for Colvars
    std::cout << "colvars: Allocating Tcl interpreter." << std::endl;
    set_tcl_interp(Tcl_CreateInterp());
  } else {
    std::cerr << "Error: init_tcl_pointers called with non-NULL tcl_interp_" << std::endl;
  }
#else
  std::cerr << "Error: Tcl support is not available in this build." << std::endl;
#endif
}


char const *colvarproxy_tcl::tcl_get_str(void *obj)
{
#if defined(COLVARS_TCL)
  return Tcl_GetString(reinterpret_cast<Tcl_Obj *>(obj));
#else
  (void) obj;
  return NULL;
#endif
}


int colvarproxy_tcl::tcl_run_script(std::string const &script)
{
#if defined(COLVARS_TCL)
  Tcl_Interp *const interp = get_tcl_interp();
  int err = Tcl_Eval(interp, script.c_str());
  if (err != TCL_OK) {
    cvm::log("Error while executing Tcl script:\n");
    cvm::error(Tcl_GetStringResult(interp));
    return COLVARS_ERROR;
  }
  return cvm::get_error();
#else
  return COLVARS_NOT_IMPLEMENTED;
#endif
}


int colvarproxy_tcl::tcl_run_file(std::string const &fileName)
{
#if defined(COLVARS_TCL)
  Tcl_Interp *const interp = get_tcl_interp();
  int err = Tcl_EvalFile(interp, fileName.c_str());
  if (err != TCL_OK) {
    cvm::log("Error while executing Tcl script file" + fileName + ":\n");
    cvm::error(Tcl_GetStringResult(interp));
    return COLVARS_ERROR;
  }
  return cvm::get_error();
#else
  return COLVARS_NOT_IMPLEMENTED;
#endif
}


int colvarproxy_tcl::tcl_run_force_callback()
{
#if defined(COLVARS_TCL)
  Tcl_Interp *const interp = get_tcl_interp();
  if (Tcl_FindCommand(interp, "calc_colvar_forces", NULL, 0) == NULL) {
    cvm::error("Error: Colvars force procedure calc_colvar_forces is not defined.\n");
    return COLVARS_ERROR;
  }

  std::string cmd = std::string("calc_colvar_forces ")
    + cvm::to_str(cvm::step_absolute());
  int err = Tcl_Eval(interp, cmd.c_str());
  if (err != TCL_OK) {
    cvm::log("Error while executing calc_colvar_forces:\n");
    cvm::error(Tcl_GetStringResult(interp));
    return COLVARS_ERROR;
  }
  return cvm::get_error();
#else
  return COLVARS_NOT_IMPLEMENTED;
#endif
}


int colvarproxy_tcl::tcl_run_colvar_callback(
                         std::string const &name,
                         std::vector<const colvarvalue *> const &cvc_values,
                         colvarvalue &value)
{
#if defined(COLVARS_TCL)

  Tcl_Interp *const interp = get_tcl_interp();
  size_t i;

  std::string cmd = std::string("calc_") + name;
  if (Tcl_FindCommand(interp, cmd.c_str(), NULL, 0) == NULL) {
    cvm::error("Error: scripted colvar procedure \"" + cmd + "\" is not defined.\n");
    return COLVARS_ERROR;
  }

  for (i = 0; i < cvc_values.size(); i++) {
    cmd += std::string(" {") + (*(cvc_values[i])).to_simple_string() +
      std::string("}");
  }
  int err = Tcl_Eval(interp, cmd.c_str());
  const char *result = Tcl_GetStringResult(interp);
  if (err != TCL_OK) {
    return cvm::error(std::string("Error while executing ")
                      + cmd + std::string(":\n") +
                      std::string(Tcl_GetStringResult(interp)),
                      COLVARS_ERROR);
  }
  std::istringstream is(result);
  if (value.from_simple_string(is.str()) != COLVARS_OK) {
    cvm::log("Error parsing colvar value from script:");
    cvm::error(result);
    return COLVARS_ERROR;
  }
  return cvm::get_error();

#else

  (void) name;
  (void) cvc_values;
  (void) value;
  return COLVARS_NOT_IMPLEMENTED;

#endif
}


int colvarproxy_tcl::tcl_run_colvar_gradient_callback(
                         std::string const &name,
                         std::vector<const colvarvalue *> const &cvc_values,
                         std::vector<cvm::matrix2d<cvm::real> > &gradient)
{
#if defined(COLVARS_TCL)

  Tcl_Interp *const interp = get_tcl_interp();
  size_t i;

  std::string cmd = std::string("calc_") + name + "_gradient";
  if (Tcl_FindCommand(interp, cmd.c_str(), NULL, 0) == NULL) {
    cvm::error("Error: scripted colvar gradient procedure \"" + cmd + "\" is not defined.\n");
    return COLVARS_ERROR;
  }

  for (i = 0; i < cvc_values.size(); i++) {
    cmd += std::string(" {") + (*(cvc_values[i])).to_simple_string() +
      std::string("}");
  }
  int err = Tcl_Eval(interp, cmd.c_str());
  if (err != TCL_OK) {
    return cvm::error(std::string("Error while executing ")
                      + cmd + std::string(":\n") +
                      std::string(Tcl_GetStringResult(interp)),
                      COLVARS_ERROR);
  }
  Tcl_Obj **list;
  int n;
  Tcl_ListObjGetElements(interp, Tcl_GetObjResult(interp),
                         &n, &list);
  if (n != int(gradient.size())) {
    cvm::error("Error parsing list of gradient values from script: found "
               + cvm::to_str(n) + " values instead of " +
               cvm::to_str(gradient.size()));
    return COLVARS_ERROR;
  }
  for (i = 0; i < gradient.size(); i++) {
    std::istringstream is(Tcl_GetString(list[i]));
    if (gradient[i].from_simple_string(is.str()) != COLVARS_OK) {
      cvm::log("Gradient matrix size: " + cvm::to_str(gradient[i].size()));
      cvm::log("Gradient string: " + cvm::to_str(Tcl_GetString(list[i])));
      cvm::error("Error parsing gradient value from script", COLVARS_ERROR);
      return COLVARS_ERROR;
    }
  }

  return cvm::get_error();

#else

  (void) name;
  (void) cvc_values;
  (void) gradient;
  return COLVARS_NOT_IMPLEMENTED;

#endif
}
