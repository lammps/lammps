// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.


CVSCRIPT(colvar_addforce,
         "Apply the given force onto this colvar and return the same\n"
         "force : float or array - Applied force; matches colvar dimensionality",
         1, 1,
         "force : float or array - Applied force; must match colvar dimensionality",
         std::string const f_str(script->obj_to_str(script->get_colvar_cmd_arg(0, objc, objv)));
         std::istringstream is(f_str);
         is.width(cvm::cv_width);
         is.precision(cvm::cv_prec);
         colvarvalue force(this_colvar->value());
         force.is_derivative();
         if (force.from_simple_string(is.str()) != COLVARS_OK) {
           script->add_error_msg("addforce : error parsing force value");
           return COLVARSCRIPT_ERROR;
         }
         this_colvar->add_bias_force(force);
         script->set_result_colvarvalue(force);
         return COLVARS_OK;
         )

CVSCRIPT(colvar_cvcflags,
         "Enable or disable individual components by setting their active flags",
         1, 1,
         "flags : integer array - Zero/nonzero value disables/enables the CVC",
         std::string const flags_str(script->obj_to_str(script->get_colvar_cmd_arg(0, objc, objv)));
         std::istringstream is(flags_str);
         std::vector<bool> flags;
         int flag;
         while (is >> flag) {
           flags.push_back(flag != 0);
         }
         int res = this_colvar->set_cvc_flags(flags);
         if (res != COLVARS_OK) {
           script->add_error_msg("Error setting CVC flags");
           return COLVARSCRIPT_ERROR;
         }
         script->set_result_str("0");
         return COLVARS_OK;
         )

CVSCRIPT(colvar_delete,
         "Delete this colvar, along with all biases that depend on it",
         0, 0,
         "",
         delete this_colvar;
         return COLVARS_OK;
         )

CVSCRIPT(colvar_get,
         "Get the value of the given feature for this colvar\n"
         "state : 1/0 - State of the given feature",
         1, 1,
         "feature : string - Name of the feature",
         return script->proc_features(this_colvar, objc, objv);
         )

CVSCRIPT(colvar_getappliedforce,
         "Return the total of the forces applied to this colvar\n"
         "force : float - Applied force; matches the colvar dimensionality",
         0, 0,
         "",
         script->set_result_colvarvalue(this_colvar->applied_force());
         return COLVARS_OK;
         )

CVSCRIPT(colvar_getatomgroups,
         "Return the atom indices used by this colvar as a list of lists\n"
         "groups : array of arrays of ints - Atom indices",
         0, 0,
         "",
         std::string result;
         std::vector<std::vector<int> > lists = this_colvar->get_atom_lists();
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
         script->set_result_str(result);
         return COLVARS_OK;
         )

CVSCRIPT(colvar_getatomids,
         "Return the list of atom indices used by this colvar\n"
         "indices : array of ints - Atom indices",
         0, 0,
         "",
         script->set_result_int_vec(this_colvar->atom_ids);
         return COLVARS_OK;
         )

CVSCRIPT(colvar_getconfig,
         "Return the configuration string of this colvar\n"
         "conf : string - Current configuration string",
         0, 0,
         "",
         script->set_result_str(this_colvar->get_config());
         return COLVARS_OK;
         )

CVSCRIPT(colvar_getgradients,
         "Return the atomic gradients of this colvar\n"
         "gradients : array of arrays of floats - Atomic gradients",
         0, 0,
         "",
         script->set_result_rvector_vec(this_colvar->atomic_gradients);
         return COLVARS_OK;
         )

CVSCRIPT(colvar_gettotalforce,
         "Return the sum of internal and external forces to this colvar\n"
         "force : float - Total force; matches the colvar dimensionality",
         0, 0,
         "",
         script->set_result_colvarvalue(this_colvar->total_force());
         return COLVARS_OK;
         )

CVSCRIPT(colvar_getvolmapids,
         "Return the list of volumetric map indices used by this colvar",
         0, 0,
         "",
         script->set_result_int_vec(this_colvar->get_volmap_ids());
         return COLVARS_OK;
         )

CVSCRIPT(colvar_help,
         "Get a help summary or the help string of one colvar subcommand\n"
         "help : string - Help string",
         0, 1,
         "command : string - Get the help string of this specific command",
         unsigned char *const cmdobj =
           script->get_colvar_cmd_arg(0, objc, objv);
         if (this_colvar) {
         }
         if (cmdobj) {
           std::string const cmdstr(script->obj_to_str(cmdobj));
           if (cmdstr.size()) {
             script->set_result_str(script->get_command_cmdline_help(colvarscript::use_colvar,
                                                                     cmdstr));
             return cvm::get_error();
           } else {
             return COLVARSCRIPT_ERROR;
           }
         } else {
           script->set_result_str(script->get_cmdline_help_summary(colvarscript::use_colvar));
           return COLVARS_OK;
         }
         )

CVSCRIPT(colvar_modifycvcs,
         "Modify configuration of individual components by passing string arguments",
         1, 1,
         "confs : sequence of strings - New configurations; empty strings are skipped",
         std::vector<std::string> const confs(script->obj_to_str_vector(script->get_colvar_cmd_arg(0, objc, objv)));
         cvm::increase_depth();
         int res = this_colvar->update_cvc_config(confs);
         cvm::decrease_depth();
         if (res != COLVARS_OK) {
           script->add_error_msg("Error setting CVC flags");
           return COLVARSCRIPT_ERROR;
         }
         script->set_result_str("0");
         return COLVARS_OK;
         )

CVSCRIPT(colvar_run_ave,
         "Get the current running average of the value of this colvar\n"
         "value : float or array - Averaged value; matches the colvar dimensionality",
         0, 0,
         "",
         script->set_result_colvarvalue(this_colvar->run_ave());
         return COLVARS_OK;
         )

CVSCRIPT(colvar_set,
         "Set the given feature of this colvar to a new value",
         2, 2,
         "feature : string - Name of the feature\n"
         "value : string - String representation of the new feature value",
         return script->proc_features(this_colvar, objc, objv);
         )

CVSCRIPT(colvar_state,
         "Print a string representation of the feature state of this colvar\n"
         "state : string - The feature state",
         0, 0,
         "",
         this_colvar->print_state();
         return COLVARS_OK;
         )

CVSCRIPT(colvar_type,
         "Get the type description of this colvar\n"
         "type : string - Type description",
         0, 0,
         "",
         script->set_result_str(this_colvar->value().type_desc(this_colvar->value().value_type));
         return COLVARS_OK;
         )

CVSCRIPT(colvar_update,
         "Recompute this colvar and return its up-to-date value\n"
         "value : float or array - Current value; matches the colvar dimensionality",
         0, 0,
         "",
         this_colvar->calc();
         this_colvar->update_forces_energy();
         script->set_result_colvarvalue(this_colvar->value());
         return COLVARS_OK;
         )

CVSCRIPT(colvar_value,
         "Get the current value of this colvar\n"
         "value : float or array - Current value; matches the colvar dimensionality",
         0, 0,
         "",
         script->set_result_colvarvalue(this_colvar->value());
         return COLVARS_OK;
         )

CVSCRIPT(colvar_width,
         "Get the width of this colvar\n"
         "width : float - Value of the width",
         0, 0,
         "",
         script->set_result_str(cvm::to_str(this_colvar->width, 0,
                                            cvm::cv_prec));
         return COLVARS_OK;
         )
