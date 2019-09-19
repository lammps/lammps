// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include <algorithm>

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
  description = "uninitialized colvar component";
  init_dependencies();
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
  description = "uninitialized colvar component";
  init_dependencies();
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


int colvar::cvc::init_dependencies() {
  size_t i;
  // Initialize static array once and for all
  if (features().size() == 0) {
    for (i = 0; i < colvardeps::f_cvc_ntot; i++) {
      modify_features().push_back(new feature);
    }

    init_feature(f_cvc_active, "active", f_type_dynamic);
//     The dependency below may become useful if we use dynamic atom groups
//     require_feature_children(f_cvc_active, f_ag_active);

    init_feature(f_cvc_scalar, "scalar", f_type_static);

    init_feature(f_cvc_gradient, "gradient", f_type_dynamic);

    init_feature(f_cvc_explicit_gradient, "explicit gradient", f_type_static);
    require_feature_children(f_cvc_explicit_gradient, f_ag_explicit_gradient);

    init_feature(f_cvc_inv_gradient, "inverse gradient", f_type_dynamic);
    require_feature_self(f_cvc_inv_gradient, f_cvc_gradient);

    init_feature(f_cvc_debug_gradient, "debug gradient", f_type_user);
    require_feature_self(f_cvc_debug_gradient, f_cvc_gradient);
    require_feature_self(f_cvc_debug_gradient, f_cvc_explicit_gradient);

    init_feature(f_cvc_Jacobian, "Jacobian derivative", f_type_dynamic);
    require_feature_self(f_cvc_Jacobian, f_cvc_inv_gradient);

    init_feature(f_cvc_com_based, "depends on group centers of mass", f_type_static);

    init_feature(f_cvc_pbc_minimum_image, "use minimum-image distances with PBCs", f_type_user);

    // Compute total force on first site only to avoid unwanted
    // coupling to other colvars (see e.g. Ciccotti et al., 2005)
    init_feature(f_cvc_one_site_total_force, "compute total force from one group", f_type_user);
    require_feature_self(f_cvc_one_site_total_force, f_cvc_com_based);

    init_feature(f_cvc_scalable, "scalable calculation", f_type_static);
    require_feature_self(f_cvc_scalable, f_cvc_scalable_com);

    init_feature(f_cvc_scalable_com, "scalable calculation of centers of mass", f_type_static);
    require_feature_self(f_cvc_scalable_com, f_cvc_com_based);


    // TODO only enable this when f_ag_scalable can be turned on for a pre-initialized group
    // require_feature_children(f_cvc_scalable, f_ag_scalable);
    // require_feature_children(f_cvc_scalable_com, f_ag_scalable_com);

    // check that everything is initialized
    for (i = 0; i < colvardeps::f_cvc_ntot; i++) {
      if (is_not_set(i)) {
        cvm::error("Uninitialized feature " + cvm::to_str(i) + " in " + description);
      }
    }
  }

  // Initialize feature_states for each instance
  // default as available, not enabled
  // except dynamic features which default as unavailable
  feature_states.reserve(f_cvc_ntot);
  for (i = 0; i < colvardeps::f_cvc_ntot; i++) {
    bool avail = is_dynamic(i) ? false : true;
    feature_states.push_back(feature_state(avail, false));
  }

  // Features that are implemented by all cvcs by default
  // Each cvc specifies what other features are available
  feature_states[f_cvc_active].available = true;
  feature_states[f_cvc_gradient].available = true;

  // CVCs are enabled from the start - get disabled based on flags
  enable(f_cvc_active);
  // feature_states[f_cvc_active].enabled = true;

  // Explicit gradients are implemented in mosts CVCs. Exceptions must be specified explicitly.
  // feature_states[f_cvc_explicit_gradient].enabled = true;
  enable(f_cvc_explicit_gradient);

  // Use minimum-image distances by default
  // feature_states[f_cvc_pbc_minimum_image].enabled = true;
  enable(f_cvc_pbc_minimum_image);

  // Features that are implemented by default if their requirements are
  feature_states[f_cvc_one_site_total_force].available = true;

  // Features That are implemented only for certain simulation engine configurations
  feature_states[f_cvc_scalable_com].available = (cvm::proxy->scalable_group_coms() == COLVARS_OK);
  feature_states[f_cvc_scalable].available = feature_states[f_cvc_scalable_com].available;

  return COLVARS_OK;
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


std::vector<std::vector<int> > colvar::cvc::get_atom_lists()
{
  std::vector<std::vector<int> > lists;

  std::vector<cvm::atom_group *>::iterator agi = atom_groups.begin();
  for ( ; agi != atom_groups.end(); ++agi) {
    (*agi)->create_sorted_ids();
    lists.push_back((*agi)->sorted_ids());
    if ((*agi)->is_enabled(f_ag_fitting_group) && (*agi)->is_enabled(f_ag_fit_gradients)) {
      cvm::atom_group &fg = *((*agi)->fitting_group);
      fg.create_sorted_ids();
      lists.push_back(fg.sorted_ids());
    }
  }
  return lists;
}


void colvar::cvc::collect_gradients(std::vector<int> const &atom_ids, std::vector<cvm::rvector> &atomic_gradients)
{
  // Coefficient: d(a * x^n) = a * n * x^(n-1) * dx
  cvm::real coeff = sup_coeff * cvm::real(sup_np) *
    cvm::integer_power(value().real_value, sup_np-1);

  for (size_t j = 0; j < atom_groups.size(); j++) {

    cvm::atom_group &ag = *(atom_groups[j]);

    // If necessary, apply inverse rotation to get atomic
    // gradient in the laboratory frame
    if (ag.b_rotate) {
      cvm::rotation const rot_inv = ag.rot.inverse();

      for (size_t k = 0; k < ag.size(); k++) {
        size_t a = std::lower_bound(atom_ids.begin(), atom_ids.end(),
                                    ag[k].id) - atom_ids.begin();
        atomic_gradients[a] += coeff * rot_inv.rotate(ag[k].grad);
      }

    } else {

      for (size_t k = 0; k < ag.size(); k++) {
        size_t a = std::lower_bound(atom_ids.begin(), atom_ids.end(),
                                    ag[k].id) - atom_ids.begin();
        atomic_gradients[a] += coeff * ag[k].grad;
      }
    }
    if (ag.is_enabled(f_ag_fitting_group) && ag.is_enabled(f_ag_fit_gradients)) {
      cvm::atom_group const &fg = *(ag.fitting_group);
      for (size_t k = 0; k < fg.size(); k++) {
        size_t a = std::lower_bound(atom_ids.begin(), atom_ids.end(),
                                    fg[k].id) - atom_ids.begin();
        // fit gradients are in the unrotated (simulation) frame
        atomic_gradients[a] += coeff * fg.fit_gradients[k];
      }
    }
  }
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
                  cvm::to_str(cvm::fabs(x_1 - x_0 - dx_pred) /
                              cvm::fabs(x_1 - x_0), 12, 5)+"\n");
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
                    cvm::to_str(cvm::fabs (x_1 - x_0 - dx_pred) /
                                cvm::fabs (x_1 - x_0),
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


void colvar::cvc::wrap(colvarvalue &x_unwrapped) const
{
  return;
}


// Static members

std::vector<colvardeps::feature *> colvar::cvc::cvc_features;
