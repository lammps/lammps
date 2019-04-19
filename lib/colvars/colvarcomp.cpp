// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvar.h"
#include "colvarcomp.h"



colvar::cvc::cvc()
  : sup_coeff(1.0),
    sup_np(1),
    b_periodic(false),
    b_try_scalable(true)
{
  init_cvc_requires();
  sup_coeff = 1.0;
  period = 0.0;
  wrap_center = 0.0;
}


colvar::cvc::cvc(std::string const &conf)
  : sup_coeff(1.0),
    sup_np(1),
    b_periodic(false),
    b_try_scalable(true)
{
  init_cvc_requires();
  sup_coeff = 1.0;
  period = 0.0;
  wrap_center = 0.0;
  init(conf);
}


int colvar::cvc::init(std::string const &conf)
{
  if (cvm::debug())
    cvm::log("Initializing cvc base object.\n");

  std::string const old_name(name);

  if (name.size() > 0) {
    cvm::log("Updating configuration for component \""+name+"\"");
  }

  if (get_keyval(conf, "name", name, name)) {
    if (name.size() > 0) {
      description = "cvc \"" + name + "\" of type " + function_type;
    } else {
      description = "unnamed cvc";
    }
    if ((name != old_name) && (old_name.size() > 0)) {
      cvm::error("Error: cannot rename component \""+old_name+
                 "\" after initialization (new name = \""+name+"\")",
                 INPUT_ERROR);
      name = old_name;
    }
  }

  get_keyval(conf, "componentCoeff", sup_coeff, sup_coeff);
  get_keyval(conf, "componentExp", sup_np, sup_np);

  get_keyval(conf, "period", period, period);
  get_keyval(conf, "wrapAround", wrap_center, wrap_center);

  get_keyval_feature(dynamic_cast<colvarparse *>(this), conf, "debugGradients",
                     f_cvc_debug_gradient, false, parse_silent);

  bool b_no_PBC = !is_enabled(f_cvc_pbc_minimum_image); // Enabled by default
  get_keyval(conf, "forceNoPBC", b_no_PBC, b_no_PBC);
  if (b_no_PBC) {
    disable(f_cvc_pbc_minimum_image);
  } else {
    enable(f_cvc_pbc_minimum_image);
  }

  // Attempt scalable calculations when in parallel? (By default yes, if available)
  get_keyval(conf, "scalable", b_try_scalable, b_try_scalable);

  if (cvm::debug())
    cvm::log("Done initializing cvc base object.\n");

  return cvm::get_error();
}


int colvar::cvc::init_total_force_params(std::string const &conf)
{
  if (cvm::get_error()) return COLVARS_ERROR;

  if (get_keyval_feature(this, conf, "oneSiteSystemForce",
                         f_cvc_one_site_total_force, is_enabled(f_cvc_one_site_total_force))) {
    cvm::log("Warning: keyword \"oneSiteSystemForce\" is deprecated: "
             "please use \"oneSiteTotalForce\" instead.\n");
  }
  if (get_keyval_feature(this, conf, "oneSiteTotalForce",
                         f_cvc_one_site_total_force, is_enabled(f_cvc_one_site_total_force))) {
    cvm::log("Computing total force on group 1 only");
  }

  if (! is_enabled(f_cvc_one_site_total_force)) {
    // check whether any of the other atom groups is dummy
    std::vector<cvm::atom_group *>::iterator agi = atom_groups.begin();
    agi++;
    for ( ; agi != atom_groups.end(); agi++) {
      if ((*agi)->b_dummy) {
        provide(f_cvc_inv_gradient, false);
        provide(f_cvc_Jacobian, false);
      }
    }
  }

  return COLVARS_OK;
}


cvm::atom_group *colvar::cvc::parse_group(std::string const &conf,
                                          char const *group_key,
                                          bool optional)
{
  cvm::atom_group *group = NULL;
  std::string group_conf;

  if (key_lookup(conf, group_key, &group_conf)) {
    group = new cvm::atom_group(group_key);

    if (b_try_scalable) {
      if (is_available(f_cvc_scalable_com)
          && is_enabled(f_cvc_com_based)
          && !is_enabled(f_cvc_debug_gradient)) {
        enable(f_cvc_scalable_com);
        enable(f_cvc_scalable);
        // The CVC makes the feature available;
        // the atom group will enable it unless it needs to compute a rotational fit
        group->provide(f_ag_scalable_com);
      }

      // TODO check for other types of parallelism here
    }

    if (group_conf.size() == 0) {
      cvm::error("Error: atom group \""+group->key+
                 "\" is set, but has no definition.\n",
                 INPUT_ERROR);
      return group;
    }

    cvm::increase_depth();
    if (group->parse(group_conf) == COLVARS_OK) {
      register_atom_group(group);
    }
    group->check_keywords(group_conf, group_key);
    if (cvm::get_error()) {
      cvm::error("Error parsing definition for atom group \""+
                 std::string(group_key)+"\"\n.", INPUT_ERROR);
    }
    cvm::decrease_depth();

  } else {
    if (! optional) {
      cvm::error("Error: definition for atom group \""+
                 std::string(group_key)+"\" not found.\n");
    }
  }

  return group;
}


int colvar::cvc::setup()
{
  description = "cvc " + name;
  return COLVARS_OK;
}


colvar::cvc::~cvc()
{
  free_children_deps();
  remove_all_children();
  for (size_t i = 0; i < atom_groups.size(); i++) {
    if (atom_groups[i] != NULL) delete atom_groups[i];
  }
}

void colvar::cvc::read_data()
{
  size_t ig;
  for (ig = 0; ig < atom_groups.size(); ig++) {
    cvm::atom_group &atoms = *(atom_groups[ig]);
    atoms.reset_atoms_data();
    atoms.read_positions();
    atoms.calc_required_properties();
    // each atom group will take care of its own fitting_group, if defined
  }

////  Don't try to get atom velocities, as no back-end currently implements it
//   if (tasks[task_output_velocity] && !tasks[task_fdiff_velocity]) {
//     for (i = 0; i < cvcs.size(); i++) {
//       for (ig = 0; ig < cvcs[i]->atom_groups.size(); ig++) {
//         cvcs[i]->atom_groups[ig]->read_velocities();
//       }
//     }
//   }
}


void colvar::cvc::calc_force_invgrads()
{
  cvm::error("Error: calculation of inverse gradients is not implemented "
             "for colvar components of type \""+function_type+"\".\n",
             COLVARS_NOT_IMPLEMENTED);
}


void colvar::cvc::calc_Jacobian_derivative()
{
  cvm::error("Error: calculation of inverse gradients is not implemented "
             "for colvar components of type \""+function_type+"\".\n",
             COLVARS_NOT_IMPLEMENTED);
}


void colvar::cvc::calc_fit_gradients()
{
  for (size_t ig = 0; ig < atom_groups.size(); ig++) {
    atom_groups[ig]->calc_fit_gradients();
  }
}


void colvar::cvc::debug_gradients()
{
  // this function should work for any scalar cvc:
  // the only difference will be the name of the atom group (here, "group")
  // NOTE: this assumes that groups for this cvc are non-overlapping,
  // since atom coordinates are modified only within the current group

  cvm::log("Debugging gradients for " + description);

  for (size_t ig = 0; ig < atom_groups.size(); ig++) {
    cvm::atom_group *group = atom_groups[ig];
    if (group->b_dummy) continue;

    cvm::rotation const rot_0 = group->rot;
    cvm::rotation const rot_inv = group->rot.inverse();

    cvm::real x_0 = x.real_value;
    if ((x.type() == colvarvalue::type_vector) && (x.size() == 1)) x_0 = x[0];

    // cvm::log("gradients     = "+cvm::to_str (gradients)+"\n");

    cvm::atom_group *group_for_fit = group->fitting_group ? group->fitting_group : group;
    cvm::atom_pos fit_gradient_sum, gradient_sum;

    // print the values of the fit gradients
    if (group->b_rotate || group->b_center) {
      if (group->is_enabled(f_ag_fit_gradients)) {
        size_t j;

        // fit_gradients are in the simulation frame: we should print them in the rotated frame
        cvm::log("Fit gradients:\n");
        for (j = 0; j < group_for_fit->fit_gradients.size(); j++) {
          cvm::log((group->fitting_group ? std::string("refPosGroup") : group->key) +
                  "[" + cvm::to_str(j) + "] = " +
                  (group->b_rotate ?
                    cvm::to_str(rot_0.rotate(group_for_fit->fit_gradients[j])) :
                    cvm::to_str(group_for_fit->fit_gradients[j])));
        }
      }
    }

    // debug the gradients
    for (size_t ia = 0; ia < group->size(); ia++) {

      // tests are best conducted in the unrotated (simulation) frame
      cvm::rvector const atom_grad = (group->b_rotate ?
                                      rot_inv.rotate((*group)[ia].grad) :
                                      (*group)[ia].grad);
      gradient_sum += atom_grad;

      for (size_t id = 0; id < 3; id++) {
        // (re)read original positions
        group->read_positions();
        // change one coordinate
        (*group)[ia].pos[id] += cvm::debug_gradients_step_size;
        group->calc_required_properties();
        calc_value();
        cvm::real x_1 = x.real_value;
        if ((x.type() == colvarvalue::type_vector) && (x.size() == 1)) x_1 = x[0];
        cvm::log("Atom "+cvm::to_str(ia)+", component "+cvm::to_str(id)+":\n");
        cvm::log("dx(actual) = "+cvm::to_str(x_1 - x_0,
                              21, 14)+"\n");
        cvm::real const dx_pred = (group->fit_gradients.size()) ?
          (cvm::debug_gradients_step_size * (atom_grad[id] + group->fit_gradients[ia][id])) :
          (cvm::debug_gradients_step_size * atom_grad[id]);
        cvm::log("dx(interp) = "+cvm::to_str(dx_pred,
                              21, 14)+"\n");
        cvm::log("|dx(actual) - dx(interp)|/|dx(actual)| = "+
                  cvm::to_str(std::fabs(x_1 - x_0 - dx_pred) /
                              std::fabs(x_1 - x_0), 12, 5)+"\n");
      }
    }

    if ((group->is_enabled(f_ag_fit_gradients)) && (group->fitting_group != NULL)) {
      cvm::atom_group *ref_group = group->fitting_group;
      group->read_positions();
      group->calc_required_properties();

      for (size_t ia = 0; ia < ref_group->size(); ia++) {

        // fit gradients are in the unrotated (simulation) frame
        cvm::rvector const atom_grad = ref_group->fit_gradients[ia];
        fit_gradient_sum += atom_grad;

        for (size_t id = 0; id < 3; id++) {
          // (re)read original positions
          group->read_positions();
          ref_group->read_positions();
          // change one coordinate
          (*ref_group)[ia].pos[id] += cvm::debug_gradients_step_size;
          group->calc_required_properties();
          calc_value();

          cvm::real const x_1 = x.real_value;
          cvm::log("refPosGroup atom "+cvm::to_str(ia)+", component "+cvm::to_str (id)+":\n");
          cvm::log("dx(actual) = "+cvm::to_str (x_1 - x_0,
                                21, 14)+"\n");

          cvm::real const dx_pred = cvm::debug_gradients_step_size * atom_grad[id];

          cvm::log("dx(interp) = "+cvm::to_str (dx_pred,
                                21, 14)+"\n");
          cvm::log ("|dx(actual) - dx(interp)|/|dx(actual)| = "+
                    cvm::to_str(std::fabs (x_1 - x_0 - dx_pred) /
                                std::fabs (x_1 - x_0),
                                12, 5)+
                    ".\n");
        }
      }
    }

    cvm::log("Gradient sum: " +  cvm::to_str(gradient_sum) +
          "  Fit gradient sum: " + cvm::to_str(fit_gradient_sum) +
          "  Total " + cvm::to_str(gradient_sum + fit_gradient_sum));
  }
  return;
}


cvm::real colvar::cvc::dist2(colvarvalue const &x1,
                             colvarvalue const &x2) const
{
  return x1.dist2(x2);
}


colvarvalue colvar::cvc::dist2_lgrad(colvarvalue const &x1,
                                     colvarvalue const &x2) const
{
  return x1.dist2_grad(x2);
}


colvarvalue colvar::cvc::dist2_rgrad(colvarvalue const &x1,
                                     colvarvalue const &x2) const
{
  return x2.dist2_grad(x1);
}


void colvar::cvc::wrap(colvarvalue &x) const
{
  return;
}


// Static members

std::vector<colvardeps::feature *> colvar::cvc::cvc_features;
