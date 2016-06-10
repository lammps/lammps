// -*- c++ -*-

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
}


colvar::cvc::cvc(std::string const &conf)
  : sup_coeff(1.0),
    sup_np(1),
    b_periodic(false),
    b_try_scalable(true)
{
  if (cvm::debug())
    cvm::log("Initializing cvc base object.\n");

  init_cvc_requires();

  if (get_keyval(conf, "name", this->name, std::string(""), parse_silent)) {
    // Temporary description until child object is initialized
    description = "cvc " + name;
  } else {
    description = "uninitialized cvc";
  }

  get_keyval(conf, "componentCoeff", sup_coeff, 1.0);
  get_keyval(conf, "componentExp", sup_np, 1);

  get_keyval(conf, "period", period, 0.0);
  get_keyval(conf, "wrapAround", wrap_center, 0.0);

  // All cvcs implement this
  provide(f_cvc_debug_gradient);
  {
    bool b_debug_gradient;
    get_keyval(conf, "debugGradients", b_debug_gradient, false, parse_silent);
    if (b_debug_gradient) enable(f_cvc_debug_gradient);
  }

  // Attempt scalable calculations when in parallel? (By default yes, if available)
  get_keyval(conf, "scalable", b_try_scalable, true);

  if (cvm::debug())
    cvm::log("Done initializing cvc base object.\n");
}


cvm::atom_group *colvar::cvc::parse_group(std::string const &conf,
                                          char const *group_key,
                                          bool optional)
{
  cvm::atom_group *group = NULL;

  if (key_lookup(conf, group_key)) {
    group = new cvm::atom_group;
    group->key = group_key;

    if (b_try_scalable) {
      if (is_available(f_cvc_scalable_com) && is_available(f_cvc_com_based)) {
        enable(f_cvc_scalable_com);
        enable(f_cvc_scalable);
        group->enable(f_ag_scalable_com);
        group->enable(f_ag_scalable);
      }

      // TODO check for other types of parallelism here

      if (is_enabled(f_cvc_scalable)) {
        cvm::log("Will enable scalable calculation for group \""+group->key+"\".\n");
      } else {
        cvm::log("Scalable calculation is not available for group \""+group->key+"\" with the current configuration.\n");
      }
    }

    if (group->parse(conf) == COLVARS_OK) {
      atom_groups.push_back(group);
    } else {
      cvm::error("Error parsing definition for atom group \""+
                         std::string(group_key)+"\".\n");
    }
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
  size_t i;
  description = "cvc " + name;

  for (i = 0; i < atom_groups.size(); i++) {
    add_child((cvm::deps *) atom_groups[i]);
  }

  if (b_try_scalable && is_available(f_cvc_scalable)) {
    enable(f_cvc_scalable);
  }

  return COLVARS_OK;
}


colvar::cvc::~cvc()
{
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
  cvm::fatal_error("Error: calculation of inverse gradients is not implemented "
                    "for colvar components of type \""+function_type+"\".\n");
}


void colvar::cvc::calc_Jacobian_derivative()
{
  cvm::fatal_error("Error: calculation of inverse gradients is not implemented "
                    "for colvar components of type \""+function_type+"\".\n");
}


void colvar::cvc::debug_gradients(cvm::atom_group *group)
{
  // this function should work for any scalar variable:
  // the only difference will be the name of the atom group (here, "group")
  // NOTE: this assumes that groups for this cvc are non-overlapping,
  // since atom coordinates are modified only within the current group

  if (group->b_dummy) return;

  cvm::rotation const rot_0 = group->rot;
  cvm::rotation const rot_inv = group->rot.inverse();

  cvm::real x_0 = x.real_value;
  if ((x.type() == colvarvalue::type_vector) && (x.size() == 1)) x_0 = x[0];

  // cvm::log("gradients     = "+cvm::to_str (gradients)+"\n");

  cvm::atom_group *group_for_fit = group->fitting_group ? group->fitting_group : group;
  cvm::atom_pos fit_gradient_sum, gradient_sum;

  // print the values of the fit gradients
  if (group->b_rotate || group->b_center) {
    if (group->b_fit_gradients) {
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

  if ((group->b_fit_gradients) && (group->fitting_group != NULL)) {
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

  return;
}


// Static members

std::vector<cvm::deps::feature *> colvar::cvc::cvc_features;
