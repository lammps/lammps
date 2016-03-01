// -*- c++ -*-

#include <cstdlib>
#include <stdlib.h>
#include <string.h>

#include "colvarscript.h"


colvarscript::colvarscript(colvarproxy *p)
 : proxy(p),
   colvars(p->colvars),
   proxy_error(0)
{
}

/// Run method based on given arguments
int colvarscript::run(int argc, char const *argv[]) {

  result = "";

  if (cvm::debug()) {
    cvm::log("Called script run with " + cvm::to_str(argc) + " args");
    for (int i = 0; i < argc; i++) { cvm::log(argv[i]); }
  }

  if (argc < 2) {
    result = help_string();
    return COLVARS_OK;
  }

  std::string cmd = argv[1];

  int error_code = COLVARS_OK;

  if (cmd == "colvar") {
    return proc_colvar(argc-1, &(argv[1]));
  }

  if (cmd == "bias") {
    return proc_bias(argc-1, &(argv[1]));
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
    colvars->reset();
    // Note: the delete bit may be ignored by some backends
    // it is mostly useful in VMD
    colvars->set_error_bits(DELETE_COLVARS);
    return COLVARS_OK;
  }

  if (cmd == "update") {
    error_code |= proxy->update_input();
    error_code |= colvars->calc();
    error_code |= proxy->update_output();
    if (error_code) {
      result += "Error updating the colvars module.\n";
    }
    return error_code;
  }

  if (cmd == "list") {
    if (argc == 2) {
      for (std::vector<colvar *>::iterator cvi = colvars->colvars.begin();
           cvi != colvars->colvars.end();
           ++cvi) {
        result += (cvi == colvars->colvars.begin() ? "" : " ") + (*cvi)->name;
      }
      return COLVARS_OK;
    } else if (argc == 3 && !strcmp(argv[2], "biases")) {
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
    if (argc < 3) {
      result = "Missing arguments\n" + help_string();
      return COLVARSCRIPT_ERROR;
    }
    if (colvars->read_config_file(argv[2]) == COLVARS_OK) {
      return COLVARS_OK;
    } else {
      result = "Error parsing configuration file";
      return COLVARSCRIPT_ERROR;
    }
  }

  /// Parse config from string
  if (cmd == "config") {
    if (argc < 3) {
      result = "Missing arguments\n" + help_string();
      return COLVARSCRIPT_ERROR;
    }
    std::string conf = argv[2];
    if (colvars->read_config_string(conf) == COLVARS_OK) {
      return COLVARS_OK;
    } else {
      result = "Error parsing configuration string";
      return COLVARSCRIPT_ERROR;
    }
  }

  /// Load an input state file
  if (cmd == "load") {
    if (argc < 3) {
      result = "Missing arguments\n" + help_string();
      return COLVARSCRIPT_ERROR;
    }
    proxy->input_prefix_str = argv[2];
    if (colvars->setup_input() == COLVARS_OK) {
      return COLVARS_OK;
    } else {
      result = "Error loading state file";
      return COLVARSCRIPT_ERROR;
    }
  }

  /// Save to an output state file
  if (cmd == "save") {
    if (argc < 3) {
      result = "Missing arguments";
      return COLVARSCRIPT_ERROR;
    }
    proxy->output_prefix_str = argv[2];
    int error = 0;
    error |= colvars->setup_output();
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
    if (argc == 2) {
      int f = proxy->frame();
      if (f >= 0) {
        result = cvm::to_str(f);
        return COLVARS_OK;
      } else {
        result = "Frame number is not available";
        return COLVARSCRIPT_ERROR;
      }
    } else if (argc == 3) {
      // Failure of this function does not trigger an error, but
      // returns the plain result to let scripts detect available frames
      long int f = proxy->frame(strtol(argv[2], NULL, 10));
      colvars->it = proxy->frame();
      result = cvm::to_str(f);
      return COLVARS_OK;
    } else {
      result = "Wrong arguments to command \"frame\"\n" + help_string();
      return COLVARSCRIPT_ERROR;
    }
  }

  result = "Syntax error\n" + help_string();
  return COLVARSCRIPT_ERROR;
}


int colvarscript::proc_colvar(int argc, char const *argv[]) {
  if (argc < 3) {
    result = "Missing parameters\n" + help_string();
    return COLVARSCRIPT_ERROR;
  }

  std::string name = argv[1];
  colvar *cv = cvm::colvar_by_name(name);
  if (cv == NULL) {
    result = "Colvar not found: " + name;
    return COLVARSCRIPT_ERROR;
  }
  std::string subcmd = argv[2];

  if (subcmd == "value") {
    result = (cv->value()).to_simple_string();
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
    cv->update();
    result = (cv->value()).to_simple_string();
    return COLVARS_OK;
  }

  if (subcmd == "delete") {
    if (cv->biases.size() > 0) {
      result = "Cannot delete a colvar currently used by biases, delete those biases first";
      return COLVARSCRIPT_ERROR;
    }
    // colvar destructor is tasked with the cleanup
    delete cv;
    // TODO this could be done by the destructors
    colvars->write_traj_label(colvars->cv_traj_os);
    return COLVARS_OK;
  }

  if (subcmd == "getconfig") {
    result = cv->get_config();
    return COLVARS_OK;
  }

  if (subcmd == "addforce") {
    if (argc < 4) {
      result = "addforce: missing parameter: force value\n" + help_string();
      return COLVARSCRIPT_ERROR;
    }
    std::string f_str = argv[3];
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
    if (argc < 4) {
      result = "cvcflags: missing parameter: vector of flags";
      return COLVARSCRIPT_ERROR;
    }
    std::string flags_str = argv[3];
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

  result = "Syntax error\n" + help_string();
  return COLVARSCRIPT_ERROR;
}


int colvarscript::proc_bias(int argc, char const *argv[]) {
  if (argc < 3) {
    result = "Missing parameters\n" + help_string();
    return COLVARSCRIPT_ERROR;
  }

  std::string name = argv[1];
  colvarbias *b = cvm::bias_by_name(name);
  if (b == NULL) {
    result = "Bias not found: " + name;
    return COLVARSCRIPT_ERROR;
  }

  std::string subcmd = argv[2];

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
    colvars->write_traj_label(colvars->cv_traj_os);
    return COLVARS_OK;
  }

  if (argc >= 4) {
    std::string param = argv[3];
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


std::string colvarscript::help_string()
{
  std::string buf;
  buf = "Usage: cv <subcommand> [args...]\n\
\n\
Managing the colvars module:\n\
  configfile <file name>      -- read configuration from a file\n\
  config <string>             -- read configuration from the given string\n\
  reset                       -- delete all internal configuration\n\
  delete                      -- delete this colvars module instance\n\
  version                     -- return version of colvars code\n\
  \n\
Input and output:\n\
  list                        -- return a list of all variables\n\
  list biases                 -- return a list of all biases\n\
  load <file name>            -- load a state file (requires configuration)\n\
  save <file name>            -- save a state file (requires configuration)\n\
  update                      -- recalculate colvars and biases based\n\
  printframe                  -- return a summary of the current frame\n\
  printframelabels            -- return labels to annotate printframe's output\n";

  if (proxy->frame() != COLVARS_NOT_IMPLEMENTED) {
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
  colvar <name> getconfig     -- return config string of colvar <name>\n\
  colvar <name> cvcflags <fl> -- enable or disable cvcs according to 0/1 flags\n\
\n\
Accessing biases:\n\
  bias <name> energy          -- return the current energy of bias <name>\n\
  bias <name> update          -- recalculate bias <name>\n\
  bias <name> delete          -- delete bias <name>\n\
  bias <name> getconfig       -- return config string of bias <name>\n";

  return buf;
}
