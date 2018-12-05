// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#ifndef COLVARSCRIPT_H
//#define COLVARSCRIPT_H // Delay definition until later

#include <string>
#include <vector>
#include <map>

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

  /// If an error is caught by the proxy through fatal_error(), this is set to
  /// COLVARSCRIPT_ERROR
  int proxy_error;

  /// If an error is returned by one of the methods, it should set this to the
  /// error message
  std::string result;

  /// Run script command with given positional arguments (objects)
  int run(int objc, unsigned char *const objv[]);

  /// Set the return value of the script command to the given string
  inline void set_str_result(std::string const &s)
  {
    result = s;
  }

  /// Build and return a short help
  std::string help_string(void) const;

  /// Use scripting language to get the string representation of an object
  inline char const *obj_to_str(unsigned char *const obj)
  {
    return cvm::proxy->script_obj_to_str(obj);
  }

  enum command {
    cv_help,
    cv_version,
    cv_config,
    cv_getconfig,
    cv_configfile,
    cv_reset,
    cv_resetindexgroups,
    cv_delete,
    cv_list,
    cv_list_biases,
    cv_load,
    cv_save,
    cv_update,
    cv_addenergy,
    cv_getenergy,
    cv_printframe,
    cv_printframelabels,
    cv_frame,
    cv_colvar,
    cv_colvar_value,
    cv_colvar_update,
    cv_colvar_type,
    cv_colvar_delete,
    cv_colvar_addforce,
    cv_colvar_getappliedforce,
    cv_colvar_gettotalforce,
    cv_colvar_cvcflags,
    cv_colvar_getconfig,
    cv_colvar_get,
    cv_colvar_set,
    cv_bias,
    cv_bias_energy,
    cv_bias_update,
    cv_bias_delete,
    cv_bias_getconfig,
    cv_bias_get,
    cv_bias_set,
    cv_n_commands
  };

  /// Execute a script command
  inline int exec_command(command c,
                          void *pobj,
                          int objc, unsigned char * const *objv)
  {
    return (*(comm_fns[c]))(pobj, objc, objv);
  }

  /// Get help for a command (TODO reformat for each language?)
  inline std::string command_help(colvarscript::command c) const
  {
    return comm_help[c];
  }

  /// Clear all object results
  inline void clear_results()
  {
    result.clear();
  }

private:

  /// Run subcommands on colvar
  int proc_colvar(colvar *cv, int argc, unsigned char *const argv[]);

  /// Run subcommands on bias
  int proc_bias(colvarbias *b, int argc, unsigned char *const argv[]);

  /// Run subcommands on base colvardeps object (colvar, bias, ...)
  int proc_features(colvardeps *obj,
                    int argc, unsigned char *const argv[]);

  /// Internal identifiers of command strings
  std::map<std::string, command> comm_str_map;

  /// Help strings for each command
  std::vector<std::string> comm_help;

  /// Number of arguments for each command
  std::vector<size_t> comm_n_args;

  /// Arguments for each command
  std::vector< std::vector<std::string> > comm_args;

  /// Implementations of each command
  std::vector<int (*)(void *, int, unsigned char * const *)> comm_fns;

};


/// Get a pointer to the main colvarscript object
inline static colvarscript *colvarscript_obj()
{
  return cvm::main()->proxy->script;
}

/// Get a pointer to the colvar object pointed to by pobj
inline static colvar *colvar_obj(void *pobj)
{
  return reinterpret_cast<colvar *>(pobj);
}

/// Get a pointer to the colvarbias object pointed to by pobj
inline static colvarbias *colvarbias_obj(void *pobj)
{
  return reinterpret_cast<colvarbias *>(pobj);
}


#define CVSCRIPT_COMM_FNAME(COMM) cvscript_ ## COMM

#define CVSCRIPT_COMM_PROTO(COMM)                                       \
  int CVSCRIPT_COMM_FNAME(COMM)(void *, int, unsigned char *const *);

#define CVSCRIPT(COMM,HELP,N_ARGS_MIN,N_ARGS_MAX,ARGS,FN_BODY)  \
  CVSCRIPT_COMM_PROTO(COMM)

#undef COLVARSCRIPT_H
#endif // #ifndef COLVARSCRIPT_H


#ifdef COLVARSCRIPT_CPP
#define CVSCRIPT_COMM_FN(COMM,N_ARGS_MIN,N_ARGS_MAX,ARGS,FN_BODY)       \
  int CVSCRIPT_COMM_FNAME(COMM)(void *pobj,                             \
                                int objc, unsigned char *const objv[])  \
  {                                                                     \
    colvarscript *script = colvarscript_obj();                          \
    script->clear_results();                                            \
    if (objc < 2+N_ARGS_MIN) /* "cv" and "COMM" are 1st and 2nd */ {    \
      script->set_str_result("Missing arguments\n" +                    \
                             script->command_help(colvarscript::COMM)); \
      return COLVARSCRIPT_ERROR;                                        \
    }                                                                   \
    if (objc > 2+N_ARGS_MAX) {                                          \
      script->set_str_result("Too many arguments\n" +                   \
                             script->command_help(colvarscript::COMM)); \
      return COLVARSCRIPT_ERROR;                                        \
    }                                                                   \
    FN_BODY;                                                            \
  }
#undef CVSCRIPT
#define CVSCRIPT(COMM,HELP,N_ARGS_MIN,N_ARGS_MAX,ARGS,FN_BODY) \
  CVSCRIPT_COMM_FN(COMM,N_ARGS_MIN,N_ARGS_MAX,ARGS,FN_BODY)
#endif // #ifdef COLVARSCRIPT_CPP


#ifdef COLVARSCRIPT_INIT_FN
#define CVSCRIPT_COMM_INIT(COMM,HELP,ARGS) {                    \
    comm_str_map[#COMM] = COMM;                                 \
    comm_help[COMM] = HELP;                                     \
    comm_fns[COMM] = &(CVSCRIPT_COMM_FNAME(COMM));              \
  }
#undef CVSCRIPT
#define CVSCRIPT(COMM,HELP,N_ARGS_MIN,N_ARGS_MAX,ARGS,FN_BODY)  \
  CVSCRIPT_COMM_INIT(COMM,HELP,ARGS)
#endif


#if !defined(COLVARSCRIPT_H) || defined(COLVARSCRIPT_INIT_FN)
#define COLVARSCRIPT_H

#ifndef COLVARSCRIPT_INIT_FN
#ifdef __cplusplus
extern "C" {
#endif
#endif

  // Add optional arguments for command-specific help?
  CVSCRIPT(cv_help,
           "Print the help message",
           0, 0,
           {},
           script->set_str_result(script->help_string());
           return COLVARS_OK;
           )

  CVSCRIPT(cv_config,
           "Read configuration from the given string",
           1, 1,
           { "conf (str) - Configuration string" },
           std::string const conf(script->obj_to_str(objv[2]));
           if (cvm::main()->read_config_string(conf) == COLVARS_OK) {
             return COLVARS_OK;
           }
           script->set_str_result("Error parsing configuration string");
           return COLVARSCRIPT_ERROR;
           )

  CVSCRIPT(cv_getconfig,
           "Get the module's configuration string read so far",
           0, 0,
           { },
           script->set_str_result(cvm::main()->get_config());
           return COLVARS_OK;
           )

  CVSCRIPT(cv_resetindexgroups,
           "Clear the index groups loaded so far, allowing to replace them",
           0, 0,
           { },
           cvm::main()->index_group_names.clear();
           cvm::main()->index_groups.clear();
           return COLVARS_OK;
           )

  CVSCRIPT(cv_addenergy,
           "Add an energy to the MD engine",
           1, 1,
           { "E (float) - Amount of energy to add" },
           cvm::main()->total_bias_energy +=
             strtod(script->obj_to_str(objv[2]), NULL);
           return COLVARS_OK;
           )

  CVSCRIPT(cv_getenergy,
           "Get the current Colvars energy",
           1, 1,
           { "E (float) - Store the energy in this variable" },
           double *energy = reinterpret_cast<double *>(objv[2]);
           *energy = cvm::main()->total_bias_energy;
           return COLVARS_OK;
           )

#ifndef COLVARSCRIPT_INIT_FN
#ifdef __cplusplus
} // extern "C"
#endif
#endif

#undef CVSCRIPT

#endif // #ifndef COLVARSCRIPT_H
