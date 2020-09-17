// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include <vector>
#include <cstdlib>
#include <string.h>

#include "colvarproxy.h"
#include "colvardeps.h"
#include "colvarscript.h"
#include "colvarscript_commands.h"



extern "C"
int cvscript_n_commands()
{
  return static_cast<int>(colvarscript::cv_n_commands);
}


extern "C"
char const **cvscript_command_names()
{
  colvarscript *script = colvarscript_obj();
  return script->get_command_names();
}


// Instantiate the body of all script commands

#define CVSCRIPT_COMM_FN(COMM,N_ARGS_MIN,N_ARGS_MAX,ARGS,FN_BODY)       \
  int CVSCRIPT_COMM_FNAME(COMM)(void *pobj,                             \
                                int objc, unsigned char *const objv[])  \
  {                                                                     \
    if (cvm::debug()) {                                                 \
      cvm::log("Executing script function \""+std::string(#COMM)+"\""); \
    }                                                                   \
    colvarscript *script = colvarscript_obj();                          \
    script->clear_str_result();                                         \
    if (script->check_module_cmd_nargs(#COMM,                           \
                                       objc, N_ARGS_MIN, N_ARGS_MAX) != \
        COLVARSCRIPT_OK) {                                              \
      return COLVARSCRIPT_ERROR;                                        \
    }                                                                   \
    FN_BODY;                                                            \
  }
#undef CVSCRIPT
#define CVSCRIPT(COMM,HELP,N_ARGS_MIN,N_ARGS_MAX,ARGS,FN_BODY) \
  CVSCRIPT_COMM_FN(COMM,N_ARGS_MIN,N_ARGS_MAX,ARGS,FN_BODY)

// Skips the colvar- and bias- specific commands
#define COLVARSCRIPT_COMMANDS_GLOBAL

#undef COLVARSCRIPT_COMMANDS_H
#include "colvarscript_commands.h"

#undef CVSCRIPT_COMM_FN
#undef CVSCRIPT
