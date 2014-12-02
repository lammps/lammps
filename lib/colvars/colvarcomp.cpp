/// -*- c++ -*-

#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvar.h"
#include "colvarcomp.h"



colvar::cvc::cvc()
  : sup_coeff(1.0), sup_np(1),
    b_periodic(false),
    b_inverse_gradients(false),
    b_Jacobian_derivative(false),
    b_debug_gradients(false)
{}


colvar::cvc::cvc(std::string const &conf)
  : sup_coeff(1.0), sup_np(1),
    b_periodic(false),
    b_inverse_gradients(false),
    b_Jacobian_derivative(false),
    b_debug_gradients(false)
{
  if (cvm::debug())
    cvm::log("Initializing cvc base object.\n");

  get_keyval(conf, "name", this->name, std::string(""), parse_silent);

  get_keyval(conf, "componentCoeff", sup_coeff, 1.0);
  get_keyval(conf, "componentExp", sup_np, 1);

  get_keyval(conf, "period", period, 0.0);
  get_keyval(conf, "wrapAround", wrap_center, 0.0);

  get_keyval(conf, "debugGradients", b_debug_gradients, false, parse_silent);

  if (cvm::debug())
    cvm::log("Done initializing cvc base object.\n");
}


void colvar::cvc::parse_group(std::string const &conf,
                               char const *group_key,
                               cvm::atom_group &group,
                               bool optional)
{
  if (key_lookup(conf, group_key)) {
    if (group.parse(conf, group_key) != COLVARS_OK) {
      cvm::error("Error parsing definition for atom group \""+
                         std::string(group_key)+"\".\n");
      return;
    }
  } else {
    if (! optional) {
      cvm::error("Error: definition for atom group \""+
                      std::string(group_key)+"\" not found.\n");
      return;
    }
  }
}


colvar::cvc::~cvc()
{}


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


void colvar::cvc::debug_gradients(cvm::atom_group &group)
{
  // this function should work for any scalar variable:
  // the only difference will be the name of the atom group (here, "group")

  if (group.b_dummy) return;

  cvm::rotation const rot_0 = group.rot;
  cvm::rotation const rot_inv = group.rot.inverse();

  cvm::real x_0 = x.real_value;
  if ((x.type() == colvarvalue::type_vector) && (x.size() == 1)) x_0 = x[0];

  // cvm::log("gradients     = "+cvm::to_str (gradients)+"\n");

  // it only makes sense to debug the fit gradients
  // when the fitting group is the same as this group
  if (group.b_rotate || group.b_center)
    if (group.b_fit_gradients && (group.ref_pos_group == NULL)) {
      group.calc_fit_gradients();
      if (group.b_rotate) {
        // fit_gradients are in the original frame, we should print them in the rotated frame
        for (size_t j = 0; j < group.fit_gradients.size(); j++) {
          group.fit_gradients[j] = rot_0.rotate(group.fit_gradients[j]);
        }
      }
      cvm::log("fit_gradients = "+cvm::to_str(group.fit_gradients)+"\n");
      if (group.b_rotate) {
        for (size_t j = 0; j < group.fit_gradients.size(); j++) {
          group.fit_gradients[j] = rot_inv.rotate(group.fit_gradients[j]);
        }
      }
    }

  for (size_t ia = 0; ia < group.size(); ia++) {

    // tests are best conducted in the unrotated (simulation) frame
    cvm::rvector const atom_grad = group.b_rotate ?
      rot_inv.rotate(group[ia].grad) :
      group[ia].grad;

    for (size_t id = 0; id < 3; id++) {
      // (re)read original positions
      group.read_positions();
      // change one coordinate
      group[ia].pos[id] += cvm::debug_gradients_step_size;
      // (re)do the fit (if defined)
      if (group.b_center || group.b_rotate) {
        group.calc_apply_roto_translation();
      }
      calc_value();
      cvm::real x_1 = x.real_value;
      if ((x.type() == colvarvalue::type_vector) && (x.size() == 1)) x_1 = x[0];
      cvm::log("Atom "+cvm::to_str(ia)+", component "+cvm::to_str(id)+":\n");
      cvm::log("dx(actual) = "+cvm::to_str(x_1 - x_0,
                             21, 14)+"\n");
      //cvm::real const dx_pred = (group.fit_gradients.size() && (group.ref_pos_group == NULL)) ?
      cvm::real const dx_pred = (group.fit_gradients.size()) ?
        (cvm::debug_gradients_step_size * (atom_grad[id] + group.fit_gradients[ia][id])) :
        (cvm::debug_gradients_step_size * atom_grad[id]);
      cvm::log("dx(interp) = "+cvm::to_str(dx_pred,
                             21, 14)+"\n");
      cvm::log("|dx(actual) - dx(interp)|/|dx(actual)| = "+
                cvm::to_str(std::fabs(x_1 - x_0 - dx_pred) /
                             std::fabs(x_1 - x_0), 12, 5)+"\n");
    }
  }

/*
 * The code below is WIP
 */
//   if (group.ref_pos_group != NULL) {
//     cvm::atom_group &ref = *group.ref_pos_group;
//     group.calc_fit_gradients();
//
//     for (size_t ia = 0; ia < ref.size(); ia++) {
//
//       for (size_t id = 0; id < 3; id++) {
//         // (re)read original positions
//         group.read_positions();
//         ref.read_positions();
//         // change one coordinate
//         ref[ia].pos[id] += cvm::debug_gradients_step_size;
//         group.calc_apply_roto_translation();
//         calc_value();
//         cvm::real const x_1 = x.real_value;
//         cvm::log("refPosGroup atom "+cvm::to_str(ia)+", component "+cvm::to_str (id)+":\n");
//         cvm::log("dx(actual) = "+cvm::to_str (x_1 - x_0,
//                                21, 14)+"\n");
//         //cvm::real const dx_pred = (group.fit_gradients.size() && (group.ref_pos_group == NULL)) ?
//         // cvm::real const dx_pred = (group.fit_gradients.size()) ?
//         //  (cvm::debug_gradients_step_size * (atom_grad[id] + group.fit_gradients[ia][id])) :
//         //  (cvm::debug_gradients_step_size * atom_grad[id]);
//         cvm::real const dx_pred = cvm::debug_gradients_step_size * ref.fit_gradients[ia][id];
//         cvm::log("dx(interp) = "+cvm::to_str (dx_pred,
//                                21, 14)+"\n");
//         cvm::log ("|dx(actual) - dx(interp)|/|dx(actual)| = "+
//                   cvm::to_str(std::fabs (x_1 - x_0 - dx_pred) /
//                                std::fabs (x_1 - x_0),
//                                12, 5)+
//                   ".\n");
//       }
//     }
//   }

  return;
}
