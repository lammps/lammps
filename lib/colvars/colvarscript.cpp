// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include <cstdlib>
#include <cstring>

#define COLVARSCRIPT_CPP
#include "colvarscript.h"
#undef COLVARSCRIPT_CPP

#include "colvarproxy.h"
#include "colvardeps.h"


colvarscript::colvarscript(colvarproxy *p)
 : proxy(p),
   colvars(p->colvars),
   proxy_error(0)
{
  comm_help.resize(colvarscript::cv_n_commands);
  comm_fns.resize(colvarscript::cv_n_commands);
#define COLVARSCRIPT_INIT_FN
#include "colvarscript.h"
#undef COLVARSCRIPT_INIT_FN
}


extern "C" {

  // Generic hooks; NAMD and VMD have Tcl-specific versions in the respective proxies

  int run_colvarscript_command(int objc, unsigned char *const objv[])
  {
    colvarproxy *cvp = cvm::proxy;
    if (!cvp) {
      return -1;
    }
    if (!cvp->script) {
      cvm::error("Called run_colvarscript_command without a script object initialized.\n");
      return -1;
    }
    return cvp->script->run(objc, objv);
  }

  const char * get_colvarscript_result()
  {
    colvarproxy *cvp = cvm::proxy;
    if (!cvp->script) {
      cvm::error("Called run_colvarscript_command without a script object initialized.\n");
      return "";
    }
    return cvp->script->result.c_str();
  }
}


/// Run method based on given arguments
int colvarscript::run(int objc, unsigned char *const objv[])
{
  result.clear();

  if (cvm::debug()) {
    cvm::log("Called script run with " + cvm::to_str(objc) + " args:");
    for (int i = 0; i < objc; i++) {
      cvm::log(obj_to_str(objv[i]));
    }
  }

  if (objc < 2) {
    set_str_result("No commands given: use \"cv help\" "
                   "for a list of commands.");
    return COLVARSCRIPT_ERROR;
  }

  std::string const cmd(obj_to_str(objv[1]));

  int error_code = COLVARS_OK;

  // If command is found in map, execute it
  std::string const cmd_key("cv_"+cmd);
  if (comm_str_map.count(cmd_key) > 0) {
    error_code |= (*(comm_fns[comm_str_map[cmd_key]]))(
                      reinterpret_cast<void *>(this), objc, objv);
    return error_code;
  }

  if (cmd == "colvar") {
    if (objc < 3) {
      result = "Missing parameters\n" + help_string();
      return COLVARSCRIPT_ERROR;
    }
    std::string const name(obj_to_str(objv[2]));
    colvar *cv = cvm::colvar_by_name(name);
    if (cv == NULL) {
      result = "Colvar not found: " + name;
      return COLVARSCRIPT_ERROR;
    }
    return proc_colvar(cv, objc-1, &(objv[1]));
  }

  if (cmd == "bias") {
    if (objc < 3) {
      result = "Missing parameters\n" + help_string();
      return COLVARSCRIPT_ERROR;
    }
    std::string const name(obj_to_str(objv[2]));
    colvarbias *b = cvm::bias_by_name(name);
    if (b == NULL) {
      result = "Bias not found: " + name;
      return COLVARSCRIPT_ERROR;
    }
    return proc_bias(b, objc-1, &(objv[1]));
  }

  if (cmd == "version") {
    result = COLVARS_VERSION;
    return COLVARS_OK;
  }

  if (cmd == "reset") {
    /// Delete every child object
    colvars->reset();
    return COLVARS_OK;
  }

  if (cmd == "delete") {
    // Note: the delete bit may be ignored by some backends
    // it is mostly useful in VMD
    return proxy->request_deletion();
  }

  if (cmd == "update") {
    error_code |= proxy->update_input();
    if (error_code) {
      result += "Error updating the Colvars module.\n";
      return error_code;
    }
    error_code |= colvars->calc();
    error_code |= proxy->update_output();
    if (error_code) {
      result += "Error updating the Colvars module.\n";
    }
    return error_code;
  }

  if (cmd == "list") {
    if (objc == 2) {
      for (std::vector<colvar *>::iterator cvi = colvars->colvars.begin();
           cvi != colvars->colvars.end();
           ++cvi) {
        result += (cvi == colvars->colvars.begin() ? "" : " ") + (*cvi)->name;
      }
      return COLVARS_OK;
    } else if (objc == 3 && !strcmp(obj_to_str(objv[2]), "biases")) {
      for (std::vector<colvarbias *>::iterator bi = colvars->biases.begin();
           bi != colvars->biases.end();
           ++bi) {
        result += (bi == colvars->biases.begin() ? "" : " ") + (*bi)->name;
      }
      return COLVARS_OK;
    } else {
      result = "Wrong arguments to command \"list\"\n" + help_string();
      return COLVARSCRIPT_ERROR;
    }
  }

  /// Parse config from file
  if (cmd == "configfile") {
    if (objc < 3) {
      result = "Missing arguments\n" + help_string();
      return COLVARSCRIPT_ERROR;
    }
    if (colvars->read_config_file(obj_to_str(objv[2])) == COLVARS_OK) {
      return COLVARS_OK;
    } else {
      result = "Error parsing configuration file";
      return COLVARSCRIPT_ERROR;
    }
  }

  /// Parse config from string
  if (cmd == "config") {
    return exec_command(cv_config, NULL, objc, objv);
  }

  /// Load an input state file
  if (cmd == "load") {
    if (objc < 3) {
      result = "Missing arguments\n" + help_string();
      return COLVARSCRIPT_ERROR;
    }
    proxy->input_prefix() = obj_to_str(objv[2]);
    if (colvars->setup_input() == COLVARS_OK) {
      return COLVARS_OK;
    } else {
      result = "Error loading state file";
      return COLVARSCRIPT_ERROR;
    }
  }

  /// Save to an output state file
  if (cmd == "save") {
    if (objc < 3) {
      result = "Missing arguments";
      return COLVARSCRIPT_ERROR;
    }
    proxy->output_prefix() = obj_to_str(objv[2]);
    int error = 0;
    error |= colvars->setup_output();
    error |= colvars->write_restart_file(colvars->output_prefix()+
                                         ".colvars.state");
    error |= colvars->write_output_files();
    return error ? COLVARSCRIPT_ERROR : COLVARS_OK;
  }

  /// Print the values that would go on colvars.traj
  if (cmd == "printframelabels") {
    std::ostringstream os;
    colvars->write_traj_label(os);
    result = os.str();
    return COLVARS_OK;
  }
  if (cmd == "printframe") {
    std::ostringstream os;
    colvars->write_traj(os);
    result = os.str();
    return COLVARS_OK;
  }

  if (cmd == "frame") {
    if (objc == 2) {
      long int f;
      int error = proxy->get_frame(f);
      if (error == COLVARS_OK) {
        result = cvm::to_str(f);
        return COLVARS_OK;
      } else {
        result = "Frame number is not available";
        return COLVARSCRIPT_ERROR;
      }
    } else if (objc == 3) {
      // Failure of this function does not trigger an error, but
      // returns nonzero, to let scripts detect available frames
      int error = proxy->set_frame(strtol(obj_to_str(objv[2]), NULL, 10));
      result = cvm::to_str(error == COLVARS_OK ? 0 : -1);
      return COLVARS_OK;
    } else {
      result = "Wrong arguments to command \"frame\"\n" + help_string();
      return COLVARSCRIPT_ERROR;
    }
  }

  if (cmd == "addenergy") {
    if (objc == 3) {
      colvars->total_bias_energy += strtod(obj_to_str(objv[2]), NULL);
      return COLVARS_OK;
    } else {
      result = "Wrong arguments to command \"addenergy\"\n" + help_string();
      return COLVARSCRIPT_ERROR;
    }
  }

  if (cmd == "help") {
    return exec_command(cv_help, NULL, objc, objv);
  }

  result = "Syntax error\n" + help_string();
  return COLVARSCRIPT_ERROR;
}


int colvarscript::proc_colvar(colvar *cv, int objc, unsigned char *const objv[]) {

  if (objc < 3) {
    result = "Missing arguments";
    return COLVARSCRIPT_ERROR;
  }
  std::string const subcmd(obj_to_str(objv[2]));

  if (subcmd == "value") {
    result = (cv->value()).to_simple_string();
    return COLVARS_OK;
  }

  if (subcmd == "run_ave") {
    result = (cv->run_ave()).to_simple_string();
    return COLVARS_OK;
  }

  if (subcmd == "width") {
    result = cvm::to_str(cv->width, 0, cvm::cv_prec);
    return COLVARS_OK;
  }

  if (subcmd == "type") {
    result = cv->value().type_desc(cv->value().value_type);
    return COLVARS_OK;
  }

  if (subcmd == "update") {
    cv->calc();
    cv->update_forces_energy();
    result = (cv->value()).to_simple_string();
    return COLVARS_OK;
  }

  if (subcmd == "delete") {
    while (cv->biases.size() > 0) {
      size_t i = cv->biases.size()-1;
      cvm::log("Warning: before deleting colvar " + cv->name
        + ", deleting parent bias " + cv->biases[i]->name);
      delete cv->biases[i];
    }
    // colvar destructor is tasked with the cleanup
    delete cv;
    // TODO this could be done by the destructors
    if (colvars->cv_traj_os != NULL) {
      colvars->write_traj_label(*(colvars->cv_traj_os));
    }
    return COLVARS_OK;
  }

  if (subcmd == "getconfig") {
    result = cv->get_config();
    return COLVARS_OK;
  }

  if (subcmd == "getatomgroups") {
    std::vector<std::vector<int> > lists = cv->get_atom_lists();
    std::vector<std::vector<int> >::iterator li = lists.begin();

    for ( ; li != lists.end(); ++li) {
      result += "{";
      std::vector<int>::iterator lj = (*li).begin();
      for ( ; lj != (*li).end(); ++lj) {
        result += cvm::to_str(*lj);
        result += " ";
      }
      result += "} ";
    }
    return COLVARS_OK;
  }

  if (subcmd == "getatomids") {
    std::vector<int>::iterator li = cv->atom_ids.begin();

    for ( ; li != cv->atom_ids.end(); ++li) {
      result += cvm::to_str(*li);
      result += " ";
    }
    return COLVARS_OK;
  }

  if (subcmd == "getgradients") {
    std::vector<cvm::rvector>::iterator li = cv->atomic_gradients.begin();

    for ( ; li != cv->atomic_gradients.end(); ++li) {
      result += "{";
      int j;
      for (j = 0; j < 3; ++j) {
        result += cvm::to_str((*li)[j]);
        result += " ";
      }
      result += "} ";
    }
    return COLVARS_OK;
  }

  if (subcmd == "getappliedforce") {
    result = (cv->applied_force()).to_simple_string();
    return COLVARS_OK;
  }

  if (subcmd == "getsystemforce") {
    // TODO warning here
    result = (cv->total_force()).to_simple_string();
    return COLVARS_OK;
  }

  if (subcmd == "gettotalforce") {
    result = (cv->total_force()).to_simple_string();
    return COLVARS_OK;
  }

  if (subcmd == "addforce") {
    if (objc < 4) {
      result = "addforce: missing parameter: force value\n" + help_string();
      return COLVARSCRIPT_ERROR;
    }
    std::string const f_str(obj_to_str(objv[3]));
    std::istringstream is(f_str);
    is.width(cvm::cv_width);
    is.precision(cvm::cv_prec);
    colvarvalue force(cv->value());
    force.is_derivative();
    if (force.from_simple_string(is.str()) != COLVARS_OK) {
      result = "addforce : error parsing force value";
      return COLVARSCRIPT_ERROR;
    }
    cv->add_bias_force(force);
    result = force.to_simple_string();
    return COLVARS_OK;
  }

  if (subcmd == "cvcflags") {
    if (objc < 4) {
      result = "cvcflags: missing parameter: vector of flags";
      return COLVARSCRIPT_ERROR;
    }
    std::string const flags_str(obj_to_str(objv[3]));
    std::istringstream is(flags_str);
    std::vector<bool> flags;

    int flag;
    while (is >> flag) {
      flags.push_back(flag != 0);
    }

    int res = cv->set_cvc_flags(flags);
    if (res != COLVARS_OK) {
      result = "Error setting CVC flags";
      return COLVARSCRIPT_ERROR;
    }
    result = "0";
    return COLVARS_OK;
  }

  if (subcmd == "modifycvcs") {
    if (objc < 4) {
      result = "cvcflags: missing parameter: vector of strings";
      return COLVARSCRIPT_ERROR;
    }
    std::vector<std::string> const confs(proxy->script_obj_to_str_vector(objv[3]));
    cvm::increase_depth();
    int res = cv->update_cvc_config(confs);
    cvm::decrease_depth();
    if (res != COLVARS_OK) {
      result = "Error setting CVC flags";
      return COLVARSCRIPT_ERROR;
    }
    result = "0";
    return COLVARS_OK;
  }

  if ((subcmd == "get") || (subcmd == "set") || (subcmd == "state")) {
    return proc_features(cv, objc, objv);
  }

  result = "Syntax error\n" + help_string();
  return COLVARSCRIPT_ERROR;
}


int colvarscript::proc_bias(colvarbias *b, int objc, unsigned char *const objv[]) {

  if (objc < 3) {
    result = "Missing arguments";
    return COLVARSCRIPT_ERROR;
  }
  std::string const subcmd(obj_to_str(objv[2]));

  if (subcmd == "energy") {
    result = cvm::to_str(b->get_energy());
    return COLVARS_OK;
  }

  if (subcmd == "update") {
    b->update();
    result = cvm::to_str(b->get_energy());
    return COLVARS_OK;
  }

  if (subcmd == "getconfig") {
    result = b->get_config();
    return COLVARS_OK;
  }

  // Subcommands for MW ABF
  if (subcmd == "bin") {
    int r = b->current_bin();
    result = cvm::to_str(r);
    return COLVARS_OK;
  }

  if (subcmd == "binnum") {
    int r = b->bin_num();
    if (r < 0) {
      result = "Error: calling bin_num() for bias " + b->name;
      return COLVARSCRIPT_ERROR;
    }
    result = cvm::to_str(r);
    return COLVARS_OK;
  }

  if (subcmd == "share") {
    int r = b->replica_share();
    if (r < 0) {
      result = "Error: calling replica_share() for bias " + b->name;
      return COLVARSCRIPT_ERROR;
    }
    result = cvm::to_str(r);
    return COLVARS_OK;
  }
  // End commands for MW ABF

  if (subcmd == "delete") {
    // the bias destructor takes care of the cleanup at cvm level
    delete b;
    // TODO this could be done by the destructors
    if (colvars->cv_traj_os != NULL) {
      colvars->write_traj_label(*(colvars->cv_traj_os));
    }
    return COLVARS_OK;
  }

  if ((subcmd == "get") || (subcmd == "set") || (subcmd == "state")) {
    return proc_features(b, objc, objv);
  }

  if (objc >= 4) {
    std::string const param(obj_to_str(objv[3]));
    if (subcmd == "count") {
      int index;
      if (!(std::istringstream(param) >> index)) {
        result = "bin_count: error parsing bin index";
        return COLVARSCRIPT_ERROR;
      }
      result = cvm::to_str(b->bin_count(index));
      return COLVARS_OK;
    }

    result = "Syntax error\n" + help_string();
    return COLVARSCRIPT_ERROR;
  }

  result = "Syntax error\n" + help_string();
  return COLVARSCRIPT_ERROR;
}


int colvarscript::proc_features(colvardeps *obj,
                                int objc, unsigned char *const objv[]) {
  // size was already checked before calling
  std::string const subcmd(obj_to_str(objv[2]));

  if (objc == 3) {
    if (subcmd == "state") {
      // TODO make this returned as result?
      obj->print_state();
      return COLVARS_OK;
    }

    // get and set commands require more arguments
    result = "Syntax error\n" + help_string();
    return COLVARSCRIPT_ERROR;
  }

  if ((subcmd == "get") || (subcmd == "set")) {
    std::vector<colvardeps::feature *> const &features = obj->features();
    std::string const req_feature(obj_to_str(objv[3]));
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

      result = "Error: feature \""+req_feature+"\" does not exist.\n";
      return COLVARSCRIPT_ERROR;

    } else {

      if (! obj->is_available(fid)) {
        result = "Error: feature \""+req_feature+"\" is unavailable.\n";
        return COLVARSCRIPT_ERROR;
      }

      if (subcmd == "get") {
        result = cvm::to_str(obj->is_enabled(fid) ? 1 : 0);
        return COLVARS_OK;
      }

      if (subcmd == "set") {
        if (objc == 5) {
          std::string const yesno =
            colvarparse::to_lower_cppstr(std::string(obj_to_str(objv[4])));
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
        result = "Syntax error\n" + help_string();
        return COLVARSCRIPT_ERROR;
      }
    }
  }

  result = "Syntax error\n" + help_string();
  return COLVARSCRIPT_ERROR;
}


std::string colvarscript::help_string() const
{
  std::string buf;
  buf = "Usage: cv <subcommand> [args...]\n\
\n\
Managing the Colvars module:\n\
  configfile <file name>      -- read configuration from a file\n\
  config <string>             -- read configuration from the given string\n\
  getconfig                   -- get the module's configuration string\n\
  resetindexgroups            -- clear the index groups loaded so far\n\
  reset                       -- delete all internal configuration\n\
  delete                      -- delete this Colvars module instance\n\
  version                     -- return version of Colvars code\n\
  \n\
Input and output:\n\
  list                        -- return a list of all variables\n\
  list biases                 -- return a list of all biases\n\
  load <file name>            -- load a state file (requires configuration)\n\
  save <file name>            -- save a state file (requires configuration)\n\
  update                      -- recalculate colvars and biases\n\
  addenergy <E>               -- add <E> to the total bias energy\n\
  printframe                  -- return a summary of the current frame\n\
  printframelabels            -- return labels to annotate printframe's output\n";

  long int tmp;
  if (proxy->get_frame(tmp) != COLVARS_NOT_IMPLEMENTED) {
      buf += "\
  frame                       -- return current frame number\n\
  frame <new_frame>           -- set frame number\n";
  }

  buf += "\n\
Accessing collective variables:\n\
  colvar <name> value         -- return the current value of colvar <name>\n\
  colvar <name> update        -- recalculate colvar <name>\n\
  colvar <name> type          -- return the type of colvar <name>\n\
  colvar <name> delete        -- delete colvar <name>\n\
  colvar <name> addforce <F>  -- apply given force on colvar <name>\n\
  colvar <name> getappliedforce -- return applied force of colvar <name>\n\
  colvar <name> gettotalforce -- return total force of colvar <name>\n\
  colvar <name> getconfig     -- return config string of colvar <name>\n\
  colvar <name> cvcflags <fl> -- enable or disable cvcs according to 0/1 flags\n\
  colvar <name> modifycvcs <str> -- pass new config strings to each CVC\n\
  colvar <name> get <f>       -- get the value of the colvar feature <f>\n\
  colvar <name> set <f> <val> -- set the value of the colvar feature <f>\n\
\n\
Accessing biases:\n\
  bias <name> energy          -- return the current energy of bias <name>\n\
  bias <name> update          -- recalculate bias <name>\n\
  bias <name> delete          -- delete bias <name>\n\
  bias <name> getconfig       -- return config string of bias <name>\n\
  bias <name> get <f>         -- get the value of the bias feature <f>\n\
  bias <name> set <f> <val>   -- set the value of the bias feature <f>\n\
";

  return buf;
}
