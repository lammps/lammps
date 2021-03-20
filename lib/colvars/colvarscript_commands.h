// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#ifndef COLVARSCRIPT_COMMANDS_H
#define COLVARSCRIPT_COMMANDS_H

// The following is a complete definition of the scripting API.

// The CVSCRIPT macro is used in four distinct contexts:
// 1) Expand to the functions' prototypes (when included generically)
// 2) List colvarscript::command entries (when included in colvarscript.h)
// 3) Implement colvarscript::init() (when included in colvarscript.cpp)
// 4) Define the functions' bodies (when included in colvarscript_commands.cpp)


// Each command is created by an instance of the CVSCRIPT macro

// The arguments of the CVSCRIPT macro are:

// COMM = the id of the command (must be a member of colvarscript::command)

// HELP = a one-line description (C string literal) for the command

// N_ARGS_MIN = the lowest number of arguments allowed

// N_ARGS_MAX = the highest number of arguments allowed

// ARGS = multi-line string literal describing each parameter; each line
//        follows the format "name : type - description"

// FN_BODY = the implementation of the function; this should be a thin wrapper
//           over existing functions; the "script" pointer to the colvarscript
//           object is already set by the CVSCRIPT_COMM_FN macro; see also the
//           functions in colvarscript_commands.h.

#ifndef CVSCRIPT_COMM_FNAME
#define CVSCRIPT_COMM_FNAME(COMM) cvscript_ ## COMM
#endif

// If CVSCRIPT is not defined, this file yields the function prototypes
#ifndef CVSCRIPT

#ifdef __cplusplus
#define CVSCRIPT_COMM_PROTO(COMM)                                       \
  extern "C" int CVSCRIPT_COMM_FNAME(COMM)(void *,                      \
                                           int, unsigned char *const *);
#else
#define CVSCRIPT_COMM_PROTO(COMM)                                       \
  int CVSCRIPT_COMM_FNAME(COMM)(void *, int, unsigned char *const *);
#endif

#define CVSCRIPT(COMM,HELP,N_ARGS_MIN,N_ARGS_MAX,ARGS,FN_BODY)  \
  CVSCRIPT_COMM_PROTO(COMM)


// Utility functions used to query the command database
extern "C" {

  /// Get the number of colvarscript commands
  int cvscript_n_commands();

  /// Get the names of all commands (array of strings)
  char const ** cvscript_command_names();

}

#endif


CVSCRIPT(cv_addenergy,
         "Add an energy to the MD engine (no effect in VMD)",
         1, 1,
         "E : float - Amount of energy to add",
         char const *Earg =
           script->obj_to_str(script->get_module_cmd_arg(0, objc, objv));
         cvm::main()->total_bias_energy += strtod(Earg, NULL);
         return cvm::get_error(); // TODO Make this multi-language
         )

CVSCRIPT(cv_bias,
         "Prefix for bias-specific commands",
         0, 0,
         "",
         // This cannot be executed from a command line
         return COLVARS_OK;
         )

CVSCRIPT(cv_colvar,
         "Prefix for colvar-specific commands",
         0, 0,
         "",
         // This cannot be executed from a command line
         return COLVARS_OK;
         )

CVSCRIPT(cv_config,
         "Read configuration from the given string",
         1, 1,
         "conf : string - Configuration string",
         char const *conf_str =
           script->obj_to_str(script->get_module_cmd_arg(0, objc, objv));
         std::string const conf(conf_str);
         if (cvm::main()->read_config_string(conf) == COLVARS_OK) {
           return COLVARS_OK;
         }
         script->add_error_msg("Error parsing configuration string");
         return COLVARSCRIPT_ERROR;
         )

CVSCRIPT(cv_configfile,
         "Read configuration from a file",
         1, 1,
         "conf_file : string - Path to configuration file",
         char const *conf_file_name =
           script->obj_to_str(script->get_module_cmd_arg(0, objc, objv));
         if (script->module()->read_config_file(conf_file_name) == COLVARS_OK) {
           return COLVARS_OK;
         } else {
           script->add_error_msg("Error parsing configuration file");
           return COLVARSCRIPT_ERROR;
         }
         )

CVSCRIPT(cv_delete,
         "Delete this Colvars module instance (VMD only)",
         0, 0,
         "",
         return script->proxy()->request_deletion();
         )

CVSCRIPT(cv_frame,
         "Get or set current frame number (VMD only)",
         0, 1,
         "frame : integer - Frame number",
         char const *arg =
           script->obj_to_str(script->get_module_cmd_arg(0, objc, objv));
         if (arg == NULL) {
           long int f = -1;
           if (script->proxy()->get_frame(f) == COLVARS_OK) {
             script->set_result_str(cvm::to_str(f));
             return COLVARS_OK;
           } else {
             script->add_error_msg("Frame number is not available");
             return COLVARSCRIPT_ERROR;
           }
         } else {
           int const f = strtol(const_cast<char *>(arg), NULL, 10);
           int error_code = script->proxy()->set_frame(f);
           if (error_code == COLVARS_NO_SUCH_FRAME) {
             script->add_error_msg("Invalid frame number: \""+std::string(arg)+
                                   "\"\n");
           }
           return error_code;
         }
         return COLVARS_OK;
         )

CVSCRIPT(cv_getconfig,
         "Get the module's configuration string read so far",
         0, 0,
         "",
         script->set_result_str(cvm::main()->get_config());
         return COLVARS_OK;
         )

CVSCRIPT(cv_getenergy,
         "Get the current Colvars energy",
         0, 0,
         "",
         script->set_result_str(cvm::to_str(cvm::main()->total_bias_energy));
         return COLVARS_OK;
         )

CVSCRIPT(cv_help,
         "Get the help string of the Colvars scripting interface",
         0, 1,
         "command : string - Get the help string of this specific command",
         unsigned char *const cmdobj =
           script->get_module_cmd_arg(0, objc, objv);
         if (cmdobj) {
           std::string const cmdstr(script->obj_to_str(cmdobj));
           if (cmdstr.size()) {
             if (cmdstr == std::string("colvar")) {
               script->set_result_str(script->get_cmdline_help_summary(colvarscript::use_colvar));
             } else if (cmdstr == std::string("bias")) {
               script->set_result_str(script->get_cmdline_help_summary(colvarscript::use_bias));
             } else {
               script->set_result_str(script->get_command_cmdline_help(colvarscript::use_module,
                                                                       cmdstr));
             }
             return cvm::get_error();
           } else {
             return COLVARSCRIPT_ERROR;
           }
         } else {
           script->set_result_str(script->get_cmdline_help_summary(colvarscript::use_module));
           return COLVARS_OK;
         }
         )

CVSCRIPT(cv_list,
         "Return a list of all variables or biases",
         0, 1,
         "param : string - \"colvars\" or \"biases\"; default is \"colvars\"",
         std::string res;
         unsigned char *const kwarg = script->get_module_cmd_arg(0, objc, objv);
         std::string const kwstr = kwarg ? script->obj_to_str(kwarg) :
           std::string("colvars");
         if (kwstr == "colvars") {
           for (std::vector<colvar *>::iterator cvi = script->module()->variables()->begin();
                cvi != script->module()->variables()->end();
                ++cvi) {
             res += (cvi == script->module()->variables()->begin() ? "" : " ") + (*cvi)->name;
           }
           script->set_result_str(res);
           return COLVARS_OK;
         } else if (kwstr == "biases") {
           for (std::vector<colvarbias *>::iterator bi = script->module()->biases.begin();
                bi != script->module()->biases.end();
                ++bi) {
             res += (bi == script->module()->biases.begin() ? "" : " ") + (*bi)->name;
           }
           script->set_result_str(res);
           return COLVARS_OK;
         } else {
           script->add_error_msg("Wrong arguments to command \"list\"\n");
           return COLVARSCRIPT_ERROR;
         }
         )

CVSCRIPT(cv_listcommands,
         "Get the list of script functions, prefixed with \"cv_\", \"colvar_\" or \"bias_\"",
         0, 0,
         "",
         int const n_commands = cvscript_n_commands();
         char const **command_names = cvscript_command_names();
         std::string result;
         for (int i = 0; i < n_commands; i++) {
           if (i > 0) result.append(1, ' ');
           result.append(std::string(command_names[i]));
         }
         script->set_result_str(result);
         return COLVARS_OK;
         )

CVSCRIPT(cv_load,
         "Load data from a state file into all matching colvars and biases",
         1, 1,
         "prefix : string - Path to existing state file or input prefix",
         char const *arg =
           script->obj_to_str(script->get_module_cmd_arg(0, objc, objv));
         script->proxy()->input_prefix() = cvm::state_file_prefix(arg);
         if (script->module()->setup_input() == COLVARS_OK) {
           return COLVARS_OK;
         } else {
           script->add_error_msg("Error loading state file");
           return COLVARSCRIPT_ERROR;
         }
         )

CVSCRIPT(cv_loadfromstring,
         "Load state data from a string into all matching colvars and biases",
         1, 1,
         "buffer : string - String buffer containing the state information",
         char const *arg =
           script->obj_to_str(script->get_module_cmd_arg(0, objc, objv));
         script->proxy()->input_buffer() = arg;
         if (script->module()->setup_input() == COLVARS_OK) {
           return COLVARS_OK;
         } else {
           script->add_error_msg("Error loading state string");
           return COLVARSCRIPT_ERROR;
         }
         )

CVSCRIPT(cv_molid,
         "Get or set the molecule ID on which Colvars is defined (VMD only)",
         0, 1,
         "molid : integer - Molecule ID; -1 means undefined",
         char const *arg =
           script->obj_to_str(script->get_module_cmd_arg(0, objc, objv));
         if (arg == NULL) {
           int molid = -1;
           script->proxy()->get_molid(molid);
           script->set_result_str(cvm::to_str(molid));
           return COLVARS_OK;
         } else {
           script->add_error_msg("Error: To change the molecule ID in VMD, use cv delete first.");
           return COLVARS_NOT_IMPLEMENTED;
         }
         )

CVSCRIPT(cv_printframe,
         "Return the values that would be written to colvars.traj",
         0, 0,
         "",
         std::ostringstream os;
         script->module()->write_traj(os);
         script->set_result_str(os.str());
         return COLVARS_OK;
         )

CVSCRIPT(cv_printframelabels,
         "Return the labels that would be written to colvars.traj",
         0, 0,
         "",
         std::ostringstream os;
         script->module()->write_traj_label(os);
         script->set_result_str(os.str());
         return COLVARS_OK;
         )

CVSCRIPT(cv_reset,
         "Delete all internal configuration",
         0, 0,
         "",
         return script->module()->reset();
         )

CVSCRIPT(cv_resetindexgroups,
         "Clear the index groups loaded so far, allowing to replace them",
         0, 0,
         "",
         cvm::main()->index_group_names.clear();
         cvm::main()->index_groups.clear();
         return COLVARS_OK;
         )

CVSCRIPT(cv_save,
         "Change the prefix of all output files and save them",
         1, 1,
         "prefix : string - Output prefix with trailing \".colvars.state\" gets removed)",
         std::string const prefix =
           cvm::state_file_prefix(script->obj_to_str(script->get_module_cmd_arg(0, objc, objv)));
         script->proxy()->output_prefix() = prefix;
         int error_code = COLVARS_OK;
         error_code |= script->module()->setup_output();
         error_code |= script->module()->write_restart_file(prefix+
                                                            ".colvars.state");
         error_code |= script->module()->write_output_files();
         return error_code;
         )

CVSCRIPT(cv_savetostring,
         "Write the Colvars state to a string and return it",
         0, 0,
         "",
         return script->module()->write_restart_string(script->modify_str_result());
         )

CVSCRIPT(cv_units,
         "Get or set the current Colvars unit system",
         0, 1,
         "units : string - The new unit system",
         char const *argstr =
           script->obj_to_str(script->get_module_cmd_arg(0, objc, objv));
         if (argstr) {
           return cvm::proxy->set_unit_system(argstr, false);
         } else {
           script->set_result_str(cvm::proxy->units);
           return COLVARS_OK;
         }
         )

CVSCRIPT(cv_update,
         "Recalculate colvars and biases",
         0, 0,
         "",
         int error_code = script->proxy()->update_input();
         if (error_code) {
           script->add_error_msg("Error updating the Colvars module (input)");
           return error_code;
         }
         error_code |= script->module()->calc();
         if (error_code) {
           script->add_error_msg("Error updating the Colvars module (calc)");
           return error_code;
         }
         error_code |= script->proxy()->update_output();
         if (error_code) {
           script->add_error_msg("Error updating the Colvars module (output)");
         }
         return error_code;
         )

CVSCRIPT(cv_version,
         "Get the Colvars Module version number",
         0, 0,
         "",
         script->set_result_str(COLVARS_VERSION);
         return COLVARS_OK;
         )

// This guard allows compiling colvar and bias function bodies in their
// respecitve files instead of colvarscript_commands.o
#ifndef COLVARSCRIPT_COMMANDS_GLOBAL
#include "colvarscript_commands_colvar.h"
#include "colvarscript_commands_bias.h"
#endif

#endif // #ifndef COLVARSCRIPT_COMMANDS_H
