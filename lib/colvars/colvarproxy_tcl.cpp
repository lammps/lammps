// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include <sstream>

#if defined(NAMD_TCL) || defined(VMDTCL)
#define COLVARS_TCL
#endif
#ifdef COLVARS_TCL
#include <tcl.h>
#endif

#include "colvarmodule.h"
#include "colvarproxy.h"
#include "colvarproxy_tcl.h"
#include "colvaratoms.h"



colvarproxy_tcl::colvarproxy_tcl()
{
  tcl_interp_ = NULL;
}


colvarproxy_tcl::~colvarproxy_tcl()
{
}


void colvarproxy_tcl::init_tcl_pointers()
{
  cvm::error("Error: Tcl support is not available in this build.\n",
             COLVARS_NOT_IMPLEMENTED);
}


char const *colvarproxy_tcl::tcl_get_str(void *obj)
{
#if defined(COLVARS_TCL)
  return Tcl_GetString(reinterpret_cast<Tcl_Obj *>(obj));
#else
  return NULL;
#endif
}


int colvarproxy_tcl::tcl_run_force_callback()
{
#if defined(COLVARS_TCL)
  Tcl_Interp *const tcl_interp =
    reinterpret_cast<Tcl_Interp *>(get_tcl_interp());
  std::string cmd = std::string("calc_colvar_forces ")
    + cvm::to_str(cvm::step_absolute());
  int err = Tcl_Eval(tcl_interp, cmd.c_str());
  if (err != TCL_OK) {
    cvm::log(std::string("Error while executing calc_colvar_forces:\n"));
    cvm::error(Tcl_GetStringResult(tcl_interp));
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

  Tcl_Interp *const tcl_interp =
    reinterpret_cast<Tcl_Interp *>(get_tcl_interp());
  size_t i;
  std::string cmd = std::string("calc_") + name;
  for (i = 0; i < cvc_values.size(); i++) {
    cmd += std::string(" {") + (*(cvc_values[i])).to_simple_string() +
      std::string("}");
  }
  int err = Tcl_Eval(tcl_interp, cmd.c_str());
  const char *result = Tcl_GetStringResult(tcl_interp);
  if (err != TCL_OK) {
    return cvm::error(std::string("Error while executing ")
                      + cmd + std::string(":\n") +
                      std::string(Tcl_GetStringResult(tcl_interp)),
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

  return COLVARS_NOT_IMPLEMENTED;

#endif
}


int colvarproxy_tcl::tcl_run_colvar_gradient_callback(
                         std::string const &name,
                         std::vector<const colvarvalue *> const &cvc_values,
                         std::vector<cvm::matrix2d<cvm::real> > &gradient)
{
#if defined(COLVARS_TCL)

  Tcl_Interp *const tcl_interp =
    reinterpret_cast<Tcl_Interp *>(get_tcl_interp());
  size_t i;
  std::string cmd = std::string("calc_") + name + "_gradient";
  for (i = 0; i < cvc_values.size(); i++) {
    cmd += std::string(" {") + (*(cvc_values[i])).to_simple_string() +
      std::string("}");
  }
  int err = Tcl_Eval(tcl_interp, cmd.c_str());
  if (err != TCL_OK) {
    return cvm::error(std::string("Error while executing ")
                      + cmd + std::string(":\n") +
                      std::string(Tcl_GetStringResult(tcl_interp)),
                      COLVARS_ERROR);
  }
  Tcl_Obj **list;
  int n;
  Tcl_ListObjGetElements(tcl_interp, Tcl_GetObjResult(tcl_interp),
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

  return COLVARS_NOT_IMPLEMENTED;

#endif
}
