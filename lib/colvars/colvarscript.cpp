// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include <cstdlib>
#include <cstring>
#include <sstream>

#include "colvarproxy.h"
#include "colvardeps.h"
#include "colvarscript.h"
#include "colvarscript_commands.h"



#ifdef COLVARS_TCL
/// Run the script API via Tcl command-line interface
/// \param clientData Not used
/// \param my_interp Pointer to Tcl_Interp object (read from Colvars if NULL)
/// \param objc Number of Tcl command parameters
/// \param objv Array of command parameters
/// \return Result of the script command
extern "C" int tcl_run_colvarscript_command(ClientData clientData,
                                            Tcl_Interp *interp_in,
                                            int objc, Tcl_Obj *const objv[]);
#endif


colvarscript::colvarscript(colvarproxy *p, colvarmodule *m)
 : proxy_(p),
   colvars(m)
{
  cmd_names = NULL;
  init_commands();
#ifdef COLVARS_TCL
  // must be called after constructing derived proxy class to allow for overloading
  proxy()->init_tcl_pointers();
  // TODO put this in backend functions so we don't have to delete
  Tcl_Interp *const interp = proxy()->get_tcl_interp();
  if (interp == NULL) {
    cvm::error("Error: trying to construct colvarscript without a Tcl interpreter.\n");
    return;
  }
  Tcl_DeleteCommand(interp, "cv");
  Tcl_CreateObjCommand(interp, "cv", tcl_run_colvarscript_command,
                       (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  cvm::log("Redefining the Tcl \"cv\" command to the new script interface.\n");
#endif
}


colvarscript::~colvarscript()
{
  if (cmd_names) {
    delete [] cmd_names;
    cmd_names = NULL;
  }
}


int colvarscript::init_commands()
{
  if (cvm::debug()) {
    cvm::log("Called colvarcript::init_commands()\n");
  }

  cmd_help.resize(colvarscript::cv_n_commands);
  cmd_rethelp.resize(colvarscript::cv_n_commands);
  cmd_n_args_min.resize(colvarscript::cv_n_commands);
  cmd_n_args_max.resize(colvarscript::cv_n_commands);
  cmd_arghelp.resize(colvarscript::cv_n_commands);
  cmd_full_help.resize(colvarscript::cv_n_commands);
  cmd_fns.resize(colvarscript::cv_n_commands);

  if (cmd_names) {
    delete [] cmd_names;
    cmd_names = NULL;
  }
  cmd_names = new char const * [colvarscript::cv_n_commands];

#undef COLVARSCRIPT_COMMANDS_H // disable include guard
#if defined(CVSCRIPT)
#undef CVSCRIPT // disable default macro
#endif
#define CVSCRIPT_COMM_INIT(COMM,HELP,N_ARGS_MIN,N_ARGS_MAX,ARGHELP) {   \
    init_command(COMM,#COMM,HELP,N_ARGS_MIN,N_ARGS_MAX,ARGHELP,&(CVSCRIPT_COMM_FNAME(COMM))); \
  }
#define CVSCRIPT(COMM,HELP,N_ARGS_MIN,N_ARGS_MAX,ARGS,FN_BODY)  \
  CVSCRIPT_COMM_INIT(COMM,HELP,N_ARGS_MIN,N_ARGS_MAX,ARGS)

#include "colvarscript_commands.h"

#undef CVSCRIPT_COMM_INIT
#undef CVSCRIPT

  return COLVARS_OK;
}


int colvarscript::init_command(colvarscript::command const &comm,
                               char const *name, char const *help,
                               int n_args_min, int n_args_max,
                               char const *arghelp,
                               int (*fn)(void *, int, unsigned char * const *))
{
  cmd_str_map[std::string(name)] = comm;
  cmd_names[comm] = name;

  // Initialize short help string and return-value help string (if present)
  {
    std::string const help_str(help);
    std::istringstream is(help_str);
    std::string line;
    std::getline(is, line);
    cmd_help[comm] = line;
    cmd_rethelp[comm] = "";
    while (std::getline(is, line)) {
      cmd_rethelp[comm] += line + "\n";
    }
  }

  // Initialize arguments' help strings
  cmd_n_args_min[comm] = n_args_min;
  cmd_n_args_max[comm] = n_args_max;
  {
    std::string const arghelp_str(arghelp);
    std::istringstream is(arghelp_str);
    std::string line;
    for (int iarg = 0; iarg < n_args_max; iarg++) {
      if (! std::getline(is, line)) {
        return cvm::error("Error: could not initialize help string for scripting "
                          "command \""+std::string(name)+"\".\n", COLVARS_BUG_ERROR);
      }
      cmd_arghelp[comm].push_back(line);
    }
  }

  cmd_full_help[comm] = cmd_help[comm]+"\n";
  if (cmd_n_args_min[comm] > 0) {
    cmd_full_help[comm] += "\nParameters\n";
    cmd_full_help[comm] += "----------\n\n";
    size_t i;
    for (i = 0; i < cmd_n_args_min[comm]; i++) {
      cmd_full_help[comm] += cmd_arghelp[comm][i]+"\n";
    }
    for (i = cmd_n_args_min[comm]; i < cmd_n_args_max[comm]; i++) {
      cmd_full_help[comm] += cmd_arghelp[comm][i]+" (optional)\n";
    }
  }
  if (cmd_rethelp[comm].size() > 0) {
    cmd_full_help[comm] += "\nReturns\n";
    cmd_full_help[comm] += "-------\n\n";
    cmd_full_help[comm] += cmd_rethelp[comm]+"\n";
  }

  cmd_fns[comm] = fn;
  if (cvm::debug()) {
    cvm::log("Defined command \""+std::string(name)+"\", with help string:\n");
    cvm::log(get_command_full_help(name));
  }

  return COLVARS_OK;
}


std::string colvarscript::get_cmd_prefix(colvarscript::Object_type t)
{
  switch (t) {
  case use_module:
    return std::string("cv_"); break;
  case use_colvar:
    return std::string("colvar_"); break;
  case use_bias:
    return std::string("bias_"); break;
  default:
    cvm::error("Error: undefined colvarscript object type.", COLVARS_BUG_ERROR);
    return std::string("");
  }
}



char const *colvarscript::get_command_help(char const *cmd)
{
  if (cmd_str_map.count(cmd) > 0) {
    colvarscript::command const c = cmd_str_map[std::string(cmd)];
    return cmd_help[c].c_str();
  }
  cvm::error("Error: command "+std::string(cmd)+
             " is not implemented.\n", COLVARS_INPUT_ERROR);
  return NULL;
}


char const *colvarscript::get_command_rethelp(char const *cmd)
{
  if (cmd_str_map.count(cmd) > 0) {
    colvarscript::command const c = cmd_str_map[std::string(cmd)];
    return cmd_rethelp[c].c_str();
  }
  cvm::error("Error: command "+std::string(cmd)+
             " is not implemented.\n", COLVARS_INPUT_ERROR);
  return NULL;
}


char const *colvarscript::get_command_arghelp(char const *cmd, int i)
{
  if (cmd_str_map.count(cmd) > 0) {
    colvarscript::command const c = cmd_str_map[std::string(cmd)];
    return cmd_arghelp[c][i].c_str();
  }
  cvm::error("Error: command "+std::string(cmd)+
             " is not implemented.\n", COLVARS_INPUT_ERROR);
  return NULL;
}


int colvarscript::get_command_n_args_min(char const *cmd)
{
  if (cmd_str_map.count(cmd) > 0) {
    colvarscript::command const c = cmd_str_map[std::string(cmd)];
    return cmd_n_args_min[c];
  }
  cvm::error("Error: command "+std::string(cmd)+
             " is not implemented.\n", COLVARS_INPUT_ERROR);
  return -1;
}


int colvarscript::get_command_n_args_max(char const *cmd)
{
  if (cmd_str_map.count(cmd) > 0) {
    colvarscript::command const c = cmd_str_map[std::string(cmd)];
    return cmd_n_args_max[c];
  }
  cvm::error("Error: command "+std::string(cmd)+
             " is not implemented.\n", COLVARS_INPUT_ERROR);
  return -1;
}


char const *colvarscript::get_command_full_help(char const *cmd)
{
  if (cmd_str_map.count(cmd) > 0) {
    colvarscript::command const c = cmd_str_map[std::string(cmd)];
    return cmd_full_help[c].c_str();
  }
  cvm::error("Error: command "+std::string(cmd)+
             " is not implemented.\n", COLVARS_INPUT_ERROR);
  return NULL;
}


std::string colvarscript::get_command_cmdline_syntax(colvarscript::Object_type t,
                                                     colvarscript::command cmd)
{
  std::string const prefix = get_cmd_prefix(t);
  std::string const cmdstr(cmd_names[cmd]);

  // Get the sub-command as used in the command line
  std::string const cmdline_cmd(cmdstr, prefix.size());
  std::string cmdline_args;

  size_t i;
  for (i = 0; i < cmd_n_args_min[cmd]; i++) {
    std::string const &arghelp = cmd_arghelp[cmd][i];
    size_t space = arghelp.find(" : ");
    cmdline_args += " <"+cmd_arghelp[cmd][i].substr(0, space)+">";
  }
  for (i = cmd_n_args_min[cmd]; i < cmd_n_args_max[cmd]; i++) {
    std::string const &arghelp = cmd_arghelp[cmd][i];
    size_t space = arghelp.find(" : ");
    cmdline_args += " ["+cmd_arghelp[cmd][i].substr(0, space)+"]";
  }

  switch (t) {
  case use_module:
    return std::string("cv "+cmdline_cmd+cmdline_args); break;
  case use_colvar:
    return std::string("cv colvar name "+cmdline_cmd+cmdline_args); break;
  case use_bias:
    return std::string("cv bias name "+cmdline_cmd+cmdline_args); break;
  default:
    // Already handled, but silence the warning
    return std::string("");
  }

  return std::string("");
}


std::string colvarscript::get_cmdline_help_summary(colvarscript::Object_type t)
{
  std::string output;
  output += "List of commands:\n\n";

  for (size_t i = 0; i < cmd_help.size(); i++) {
    std::string const prefix = get_cmd_prefix(t);
    command const c = cmd_str_map[std::string(cmd_names[i])];
    if (std::string(cmd_names[i], prefix.size()) == prefix) {
      output += get_command_cmdline_syntax(t, c)+std::string("\n");
    }
  }
  if (t == use_module) {
    output += "\nFor detailed help on each command use:\n"
      "    cv help <command>\n";
    output += "\nTo list all commands acting on collective variables use:\n"
      "    cv help colvar\n";
    output += "\nTo list all commands acting on biases use:\n"
      "    cv help bias\n";
  }
  if (t == use_colvar) {
    output += "\nFor detailed help on each command use:\n"
      "    cv colvar name help <command> (\"name\" does not need to exist)\n";
  }
  if (t == use_bias) {
    output += "\nFor detailed help on each command use:\n"
      "    cv bias name help <command> (\"name\" does not need to exist)\n";
  }
  return output;
}


std::string colvarscript::get_command_cmdline_help(colvarscript::Object_type t,
                                                   std::string const &cmd)
{
  std::string const cmdkey(get_cmd_prefix(t)+cmd);
  if (cmd_str_map.count(cmdkey) > 0) {
    command const c = cmd_str_map[cmdkey];
    return get_command_cmdline_syntax(t, c)+"\n\n"+
      get_command_full_help(cmd_names[c]);
  }
  cvm::set_error_bits(COLVARS_INPUT_ERROR);
  return std::string("Could not find scripting command \""+cmd+"\".");
}


int colvarscript::run(int objc, unsigned char *const objv[])
{
  clear_str_result();

  if (cvm::debug()) {
    cvm::log("Called script run with " + cvm::to_str(objc) + " args:");
    for (int i = 0; i < objc; i++) {
      cvm::log(obj_to_str(objv[i]));
    }
  }

  if (objc < 2) {
    set_result_str("No commands given: use \"cv help\" "
                   "for a list of commands.");
    return COLVARSCRIPT_ERROR;
  }

  // Main command; usually "cv"
  std::string const main_cmd(std::string(obj_to_str(objv[0])));

  // Name of the (sub)command
  std::string const cmd(obj_to_str(objv[1]));

  // Build a safe-to-print command line to print in case of error
  std::string cmdline(main_cmd+std::string(" ")+cmd);

  // Pointer to the function implementing it
  int (*cmd_fn)(void *, int, unsigned char * const *) = NULL;

  // Pointer to object handling the command (the functions are C)
  void *obj_for_cmd = NULL;

  if (cmd == "colvar") {

    if (objc < 4) {
      add_error_msg("Missing parameters: use \""+main_cmd+
                    " help colvar\" for a summary");
      return COLVARSCRIPT_ERROR;
    }
    std::string const name(obj_to_str(objv[2]));
    std::string const subcmd(obj_to_str(objv[3]));
    obj_for_cmd = reinterpret_cast<void *>(cvm::colvar_by_name(name));
    if (obj_for_cmd == NULL) {
      if (subcmd != std::string("help")) {
        // Unless asking for help, a valid colvar name must be given
        add_error_msg("Colvar not found: " + name);
        return COLVARSCRIPT_ERROR;
      }
    }
    cmd_fn = get_cmd_fn(get_cmd_prefix(use_colvar)+subcmd);
    cmdline += std::string(" name ")+subcmd;
    if (objc > 4) {
      cmdline += " ...";
    }

  } else if (cmd == "bias") {

    if (objc < 4) {
      add_error_msg("Missing parameters: use \""+main_cmd+
                    " help bias\" for a summary");
      return COLVARSCRIPT_ERROR;
    }
    std::string const name(obj_to_str(objv[2]));
    std::string const subcmd(obj_to_str(objv[3]));
    obj_for_cmd = reinterpret_cast<void *>(cvm::bias_by_name(name));
    if (obj_for_cmd == NULL) {
      if ((subcmd == "") || (subcmd != std::string("help"))) {
        // Unless asking for help, a valid bias name must be given
        add_error_msg("Bias not found: " + name);
        return COLVARSCRIPT_ERROR;
      }
    }
    cmd_fn = get_cmd_fn(get_cmd_prefix(use_bias)+subcmd);
    cmdline += std::string(" name ")+subcmd;
    if (objc > 4) {
      cmdline += " ...";
    }

  } else {

    cmd_fn = get_cmd_fn(get_cmd_prefix(use_module)+cmd);
    obj_for_cmd = reinterpret_cast<void *>(this);

    if (objc > 2) {
      cmdline += " ...";
    }
  }

  int error_code = COLVARS_OK;

  // If command was found in map, execute it
  if (cmd_fn) {
    error_code = (*cmd_fn)(obj_for_cmd, objc, objv);
  } else {
    add_error_msg("Syntax error: "+cmdline+"\n"
                  "  Run \"cv help\" or \"cv help <command>\" "
                  "to get the correct syntax.\n");
    error_code = COLVARSCRIPT_ERROR;
  }

  return error_code;
}


char *colvarscript::obj_to_str(unsigned char *obj)
{
  char *strobj = reinterpret_cast<char *>(obj);
  if (cvm::debug()) {
    cvm::log("Using simple-cast script::obj_to_str(): result = \"" +
             (strobj ? std::string(strobj) : std::string("(null)")) + "\"");
  }
  return strobj;
}


std::vector<std::string> colvarscript::obj_to_str_vector(unsigned char *obj)
{
  if (cvm::debug()) {
    cvm::log("Using simple-cast colvarscript::obj_to_str_vector().\n");
  }

  std::vector<std::string> new_result;
  std::string const str(reinterpret_cast<char *>(obj));

  // TODO get rid of this once colvarscript can handle both fix_modify and Tcl?
  // LAMMPS has a nicer function in the utils class

  for (size_t i = 0; i < str.length(); i++) {
    char const c = str[i];
    if (c == '\"') {
      i++;
      if (i >= str.length()) {
        cvm::error("Error: could not split the following string:\n"+
                   str+"\n", COLVARS_INPUT_ERROR);
        break;
      }
      new_result.push_back(std::string(""));
      while (str[i] != '\"') {
        new_result.back().append(1, str[i]);
        if (i >= str.length()) {
          cvm::error("Error: could not split the following string:\n"+
                     str+"\n", COLVARS_INPUT_ERROR);
          break;
        } else {
          i++;
        }
      }
    }
  }

  if (cvm::debug()) {
    cvm::log("result = "+cvm::to_str(new_result)+".\n");
  }

  return new_result;
}


int colvarscript::proc_features(colvardeps *obj,
                                int objc, unsigned char *const objv[]) {

  // size was already checked before calling
  std::string const subcmd(obj_to_str(objv[3]));

  if (cvm::debug()) {
    cvm::log("Called proc_features() with " + cvm::to_str(objc) + " args:");
    for (int i = 0; i < objc; i++) {
      cvm::log(obj_to_str(objv[i]));
    }
  }

  if ((subcmd == "get") || (subcmd == "set")) {
    std::vector<colvardeps::feature *> const &features = obj->features();
    std::string const req_feature(obj_to_str(objv[4]));
    colvardeps::feature *f = NULL;
    int fid = 0;
    for (fid = 0; fid < int(features.size()); fid++) {
      if (features[fid]->description ==
          colvarparse::to_lower_cppstr(req_feature)) {
        f = features[fid];
        break;
      }
    }

    if (f == NULL) {

      add_error_msg("Error: feature \""+req_feature+"\" does not exist.\n");
      return COLVARSCRIPT_ERROR;

    } else {

      if (! obj->is_available(fid)) {
        add_error_msg("Error: feature \""+req_feature+"\" is unavailable.\n");
        return COLVARSCRIPT_ERROR;
      }

      if (subcmd == "get") {
        set_result_str(cvm::to_str(obj->is_enabled(fid) ? 1 : 0));
        return COLVARS_OK;
      }

      if (subcmd == "set") {
        if (objc == 6) {
          std::string const yesno =
            colvarparse::to_lower_cppstr(std::string(obj_to_str(objv[5])));
          if ((yesno == std::string("yes")) ||
              (yesno == std::string("on")) ||
              (yesno == std::string("1"))) {
            obj->enable(fid);
            return COLVARS_OK;
          } else if ((yesno == std::string("no")) ||
              (yesno == std::string("off")) ||
              (yesno == std::string("0"))) {
            obj->disable(fid);
            return COLVARS_OK;
          }
        }
        add_error_msg("Missing value when setting feature \""+req_feature+
                      "\".\n");
        return COLVARSCRIPT_ERROR;
      }
    }
  }

  // This shouldn't be reached any more
  return COLVARSCRIPT_ERROR;
}


int colvarscript::unsupported_op()
{
  return cvm::error("Error: unsupported script operation.\n",
                    COLVARS_NOT_IMPLEMENTED);
}


int colvarscript::set_result_str(std::string const &s)
{
  if (cvm::get_error() != COLVARS_OK) {
    // Avoid overwriting the error message
    modify_str_result() += s;
  } else {
    modify_str_result() = s;
  }
  return COLVARS_OK;
}


void colvarscript::add_error_msg(std::string const &s)
{
  modify_str_result() += s;
  // Ensure terminating newlines
  if (s[s.size()-1] != '\n') {
    modify_str_result() += "\n";
  }
}


int colvarscript::clear_str_result()
{
  modify_str_result().clear();
  return COLVARS_OK;
}


extern "C"
int run_colvarscript_command(int objc, unsigned char *const objv[])
{
  colvarmodule *cv = cvm::main();
  colvarscript *script = cv ? cv->proxy->script : NULL;
  if (!script) {
    cvm::error("Called run_colvarscript_command without a script object.\n",
               COLVARS_BUG_ERROR);
    return -1;
  }
  int retval = script->run(objc, objv);
  return retval;
}


extern "C"
const char * get_colvarscript_result()
{
  colvarscript *script = colvarscript_obj();
  if (!script) {
    cvm::error("Called get_colvarscript_result without a script object.\n");
    return NULL;
  }
  return script->str_result().c_str();
}


#if defined(COLVARS_TCL)

#if defined(VMDTCL)
// Function used by VMD to set up the module
int tcl_colvars_vmd_init(Tcl_Interp *interp, int molid);
#endif

#if !defined(VMDTCL) && !defined(NAMD_TCL)
// Initialize Colvars when loaded as a shared library into Tcl interpreter
extern "C" {
  int Colvars_Init(Tcl_Interp *interp) {
    colvarproxy *proxy = new colvarproxy();
    colvarmodule *colvars = new colvarmodule(proxy);
    proxy->set_tcl_interp(interp);
    proxy->colvars = colvars;
    Tcl_CreateObjCommand(interp, "cv", tcl_run_colvarscript_command,
                         (ClientData *) NULL, (Tcl_CmdDeleteProc *) NULL);
    Tcl_EvalEx(interp, "package provide colvars", -1, 0);
    return TCL_OK;
  }
}
#endif


extern "C" int tcl_run_colvarscript_command(ClientData /* clientData */,
                                            Tcl_Interp *my_interp,
                                            int objc, Tcl_Obj *const objv[])
{
  colvarmodule *colvars = cvm::main();

  if (!colvars) {
#if defined(VMDTCL)

    if (objc == 2) {
      if (!strcmp(Tcl_GetString(objv[1]), "molid")) {
        // return invalid molid
        Tcl_SetResult(my_interp, (char *) "-1", TCL_STATIC);
      }
      if (!strcmp(Tcl_GetString(objv[1]), "delete") ||
          !strcmp(Tcl_GetString(objv[1]), "reset")) {
        // nothing to delete or reset
        Tcl_SetResult(my_interp, NULL, TCL_STATIC);
      }
      if (!strcmp(Tcl_GetString(objv[1]), "help")) {
        // print message
        Tcl_SetResult(my_interp,
                      (char *) "First, setup the Colvars module with: "
                      "cv molid <id>|top", TCL_STATIC);
      }
      return TCL_OK;
    }

    if (objc >= 3) {
      // require a molid to create the module
      if (!strcmp(Tcl_GetString(objv[1]), "molid")) {
        int molid = -(1<<16); // This value is used to indicate "top"
        if (strcmp(Tcl_GetString(objv[2]), "top")) {
          // If this is not "top", get the integer value
          Tcl_GetIntFromObj(my_interp, objv[2], &molid);
        }
        return tcl_colvars_vmd_init(my_interp, molid);
      } else {
        Tcl_SetResult(my_interp, (char *) "Syntax error.  First, setup the Colvars module with cv molid <id>|top", TCL_STATIC);
        return TCL_ERROR;
      }
    }

    Tcl_SetResult(my_interp, (char *) "First, setup the Colvars module with: "
                  "cv molid <id>|top", TCL_STATIC);

#else
    Tcl_SetResult(my_interp,
                  const_cast<char *>("Error: Colvars module not yet initialized"),
                  TCL_STATIC);
#endif
    return TCL_ERROR;
  }

  colvarproxy *proxy = colvars->proxy;
  Tcl_Interp *interp = my_interp ? my_interp : proxy->get_tcl_interp();
  colvarscript *script = colvarscript_obj();
  if (!script) {
    char const *errstr = "Called tcl_run_colvarscript_command "
      "without a Colvars script interface set up.\n";
    Tcl_SetResult(interp, const_cast<char *>(errstr), TCL_VOLATILE);
    return TCL_ERROR;
  }

  cvm::clear_error();

  unsigned char * arg_pointers_[100];
  if (objc > 100) {
    std::string const errstr = "Too many positional arguments ("+
      cvm::to_str(objc)+") passed to the \"cv\" command.\n";
    Tcl_SetResult(interp, const_cast<char *>(errstr.c_str()), TCL_VOLATILE);
    return TCL_ERROR;
  }
  for (int i = 0; i < objc; i++) {
    arg_pointers_[i] = reinterpret_cast<unsigned char *>(const_cast<char *>(proxy->tcl_get_str(objv[i])));
  }
  int retval = script->run(objc, arg_pointers_);

  std::string result = proxy->get_error_msgs() + script->str_result();

  Tcl_SetResult(interp, const_cast<char *>(result.c_str()),
                TCL_VOLATILE);

  if (proxy->delete_requested()) {
    if (!proxy->simulation_running()) {
      // Running in VMD
      Tcl_SetResult(interp,
                    const_cast<char *>("Deleting Colvars module"
                                       ": to recreate, use cv molid <molecule ID>"),
                    TCL_STATIC);
    }
    delete proxy;
    proxy = NULL;
  }

  return (retval == COLVARS_OK) ? TCL_OK : TCL_ERROR;
}

#endif // #if defined(COLVARS_TCL)




int colvarscript::set_result_text_from_str(std::string const &x_str,
                                           unsigned char *obj) {
  if (obj) {
    strcpy(reinterpret_cast<char *>(obj), x_str.c_str());
  } else {
    set_result_str(x_str);
  }
  return COLVARS_OK;
}

// Template to convert everything to string and use the above

template <typename T>
int colvarscript::set_result_text(T const &x, unsigned char *obj) {
  std::string const x_str = x.to_simple_string();
  return set_result_text_from_str(x_str, obj);
}


template <typename T>
int colvarscript::pack_vector_elements_text(std::vector<T> const &x,
                                            std::string &x_str) {
  x_str.clear();
  for (size_t i = 0; i < x.size(); ++i) {
    if (i > 0) x_str.append(1, ' ');
    x_str += cvm::to_str(x[i]);
  }
  return COLVARS_OK;
}


// Specializations for plain old data types that don't have a stringifier member

template <>
int colvarscript::set_result_text(int const &x, unsigned char *obj) {
  std::string const x_str = cvm::to_str(x);
  return set_result_text_from_str(x_str, obj);
}

template <>
int colvarscript::set_result_text(std::vector<int> const &x,
                                  unsigned char *obj) {
  std::string x_str("");
  pack_vector_elements_text<int>(x, x_str);
  return set_result_text_from_str(x_str, obj);
}


template <>
int colvarscript::set_result_text(long int const &x, unsigned char *obj) {
  std::string const x_str = cvm::to_str(x);
  return set_result_text_from_str(x_str, obj);
}

template <>
int colvarscript::set_result_text(std::vector<long int> const &x,
                                  unsigned char *obj) {
  std::string x_str("");
  pack_vector_elements_text<long int>(x, x_str);
  return set_result_text_from_str(x_str, obj);
}


template <>
int colvarscript::set_result_text(cvm::real const &x, unsigned char *obj) {
  std::string const x_str = cvm::to_str(x);
  return set_result_text_from_str(x_str, obj);
}

template <>
int colvarscript::set_result_text(std::vector<cvm::real> const &x,
                                  unsigned char *obj) {
  std::string x_str("");
  pack_vector_elements_text<cvm::real>(x, x_str);
  return set_result_text_from_str(x_str, obj);
}


// TODO these can be removed after the Tcl backend is ready (otherwise, the
// default template syntax may break scripts or the Dashboard)

template <>
int colvarscript::set_result_text(std::vector<cvm::rvector> const &x,
                                  unsigned char *obj) {
  std::string x_str("");
  for (size_t i = 0; i < x.size(); i++) {
    if (i > 0) x_str.append(1, ' ');
    x_str += "{ "+x[i].to_simple_string()+" }";
  }
  return set_result_text_from_str(x_str, obj);
}

template <>
int colvarscript::set_result_text(std::vector<colvarvalue> const &x,
                                  unsigned char *obj) {
  std::string x_str("");
  for (size_t i = 0; i < x.size(); i++) {
    if (i > 0) x_str.append(1, ' ');
    x_str += "{ "+x[i].to_simple_string()+" }";
  }
  return set_result_text_from_str(x_str, obj);
}


// Member functions to set script results for each typexc

int colvarscript::set_result_int(int const &x, unsigned char *obj) {
  return set_result_text<int>(x, obj);
}

int colvarscript::set_result_int_vec(std::vector<int> const &x,
                                     unsigned char *obj) {
  return set_result_text< std::vector<int> >(x, obj);
}


int colvarscript::set_result_long_int(long int const &x, unsigned char *obj) {
  return set_result_text<long int>(x, obj);
}

int colvarscript::set_result_long_int_vec(std::vector<long int> const &x,
                                          unsigned char *obj) {
  return set_result_text< std::vector<long int> >(x, obj);
}


int colvarscript::set_result_real(cvm::real const &x, unsigned char *obj) {
  return set_result_text<cvm::real>(x, obj);
}

int colvarscript::set_result_real_vec(std::vector<cvm::real> const &x,
                                      unsigned char *obj) {
  return set_result_text< std::vector<cvm::real> >(x, obj);
}


int colvarscript::set_result_rvector(cvm::rvector const &x, unsigned char *obj) {
  return set_result_text<cvm::rvector>(x, obj);
}

int colvarscript::set_result_rvector_vec(std::vector<cvm::rvector> const &x,
                                         unsigned char *obj) {
  return set_result_text< std::vector<cvm::rvector> >(x, obj);
}


int colvarscript::set_result_colvarvalue(colvarvalue const &x,
                                         unsigned char *obj) {
  return set_result_text<colvarvalue>(x, obj);
}

int colvarscript::set_result_colvarvalue_vec(std::vector<colvarvalue> const &x,
                                             unsigned char *obj) {
  return set_result_text< std::vector<colvarvalue> >(x, obj);
}
