// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.


CVSCRIPT(bias_bin,
         "Get the current grid bin index (1D ABF only for now)\n"
         "bin : integer - Bin index",
         0, 0,
         "",
         script->set_result_int(this_bias->current_bin());
         return COLVARS_OK;
         )

CVSCRIPT(bias_bincount,
         "Get the number of samples at the given grid bin (1D ABF only for now)\n"
         "samples : integer - Number of samples",
         0, 1,
         "index : integer - Grid index; defaults to current bin",
         int index = this_bias->current_bin();
         char const *indexarg =
           script->obj_to_str(script->get_bias_cmd_arg(0, objc, objv));
         if (indexarg) {
           std::string const param(indexarg);
           if (!(std::istringstream(param) >> index)) {
             script->add_error_msg("bincount: error parsing bin index");
             return COLVARSCRIPT_ERROR;
           }
         }
         script->set_result_int(this_bias->bin_count(index));
         return COLVARS_OK;
         )

CVSCRIPT(bias_binnum,
         "Get the total number of grid points of this bias (1D ABF only for now)\n"
         "Bins : integer - Number of grid points",
         0, 0,
         "",
         int r = this_bias->bin_num();
         if (r < 0) {
           script->add_error_msg("Error: calling bin_num() for bias " +
                                 this_bias->name);
           return COLVARSCRIPT_ERROR;
         }
         script->set_result_int(r);
         return COLVARS_OK;
         )

CVSCRIPT(bias_delete,
         "Delete this bias",
         0, 0,
         "",
         delete this_bias;
         return COLVARS_OK;
         )

CVSCRIPT(bias_energy,
         "Get the current energy of this bias\n"
         "E : float - Energy value",
         0, 0,
         "",
         script->set_result_real(this_bias->get_energy());
         return COLVARS_OK;
         )

CVSCRIPT(bias_get,
         "Get the value of the given feature for this bias\n"
         "state : 1/0 - State of the given feature",
         1, 1,
         "feature : string - Name of the feature",
         return script->proc_features(this_bias, objc, objv);
         )

CVSCRIPT(bias_getconfig,
         "Return the configuration string of this bias\n"
         "conf : string - Current configuration string",
         0, 0,
         "",
         script->set_result_str(this_bias->get_config());
         return COLVARS_OK;
         )

CVSCRIPT(bias_help,
         "Get a help summary or the help string of one bias subcommand\n"
         "help : string - Help string",
         0, 1,
         "command : string - Get the help string of this specific command",
         unsigned char *const cmdobj =
           script->get_colvar_cmd_arg(0, objc, objv);
         if (this_bias) {
         }
         if (cmdobj) {
           std::string const cmdstr(script->obj_to_str(cmdobj));
           if (cmdstr.size()) {
             script->set_result_str(script->get_command_cmdline_help(colvarscript::use_bias,
                                                                     cmdstr));
             return COLVARS_OK;
           } else {
             return COLVARSCRIPT_ERROR;
           }
         } else {
           script->set_result_str(script->get_cmdline_help_summary(colvarscript::use_bias));
           return COLVARS_OK;
         }
         )

CVSCRIPT(bias_load,
         "Load data into this bias",
         1, 1,
         "prefix : string - Read from a file with this name or prefix",
         char const *arg =
           script->obj_to_str(script->get_bias_cmd_arg(0, objc, objv));
         return this_bias->read_state_prefix(std::string(arg));
         )

CVSCRIPT(bias_loadfromstring,
         "Load state data into this bias from a string",
         1, 1,
         "buffer : string - String buffer containing the state information",
         char const *buffer = script->obj_to_str(script->get_bias_cmd_arg(0, objc, objv));
         return this_bias->read_state_string(buffer);
         )

CVSCRIPT(bias_save,
         "Save data from this bias into a file with the given prefix",
         1, 1,
         "prefix : string - Prefix for the state file of this bias",
         std::string const prefix =
           cvm::state_file_prefix(script->obj_to_str(script->get_bias_cmd_arg(0, objc, objv)));
         return this_bias->write_state_prefix(prefix);
         )

CVSCRIPT(bias_savetostring,
         "Save data from this bias into a string and return it\n"
         "state : string - The bias state",
         0, 0,
         "",
         return this_bias->write_state_string(script->modify_str_result());
         )

CVSCRIPT(bias_set,
         "Set the given feature of this bias to a new value",
         2, 2,
         "feature : string - Name of the feature\n"
         "value : string - String representation of the new feature value",
         return script->proc_features(this_bias, objc, objv);
         )

CVSCRIPT(bias_share,
         "Share bias information with other replicas (multiple-walker scheme)",
         0, 0,
         "",
         if (this_bias->replica_share() != COLVARS_OK) {
           script->add_error_msg("Error: calling replica_share() for bias " +
                                 this_bias->name);
           return COLVARSCRIPT_ERROR;
         }
         return COLVARS_OK;
         )

CVSCRIPT(bias_state,
         "Print a string representation of the feature state of this bias\n"
         "state : string - String representation of the bias features",
         0, 0,
         "",
         this_bias->print_state();
         return COLVARS_OK;
         )

CVSCRIPT(bias_update,
         "Recompute this bias and return its up-to-date energy\n"
         "E : float - Energy value",
         0, 0,
         "",
         this_bias->update();
         script->set_result_colvarvalue(this_bias->get_energy());
         return COLVARS_OK;
         )
