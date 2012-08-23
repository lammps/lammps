#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvar.h"
#include "colvarcomp.h"



colvar::cvc::cvc()
  : sup_coeff (1.0), sup_np (1),
    b_periodic (false),
    b_debug_gradients (false),
    b_inverse_gradients (false),
    b_Jacobian_derivative (false)
{}


colvar::cvc::cvc (std::string const &conf)
  : sup_coeff (1.0), sup_np (1),
    b_periodic (false),
    b_debug_gradients (false),
    b_inverse_gradients (false),
    b_Jacobian_derivative (false)
{
  if (cvm::debug())
    cvm::log ("Initializing cvc base object.\n");

  get_keyval (conf, "name", this->name, std::string (""), parse_silent);

  get_keyval (conf, "componentCoeff", sup_coeff, 1.0);
  get_keyval (conf, "componentExp", sup_np, 1);

  get_keyval (conf, "period", period, 0.0);
  get_keyval (conf, "wrapAround", wrap_center, 0.0);

  get_keyval (conf, "debugGradients", b_debug_gradients, false, parse_silent);

  if (cvm::debug())
    cvm::log ("Done initializing cvc base object.\n");
}


void colvar::cvc::parse_group (std::string const &conf, 
                               char const *group_key,
                               cvm::atom_group &group,
                               bool optional)
{
  if (key_lookup (conf, group_key)) {
    group.parse (conf, group_key);
  } else {
    if (! optional) {
      cvm::fatal_error ("Error: definition for atom group \""+
                      std::string (group_key)+"\" not found.\n");
    }
  }
}


colvar::cvc::~cvc()
{}


void colvar::cvc::calc_force_invgrads()
{
  cvm::fatal_error ("Error: calculation of inverse gradients is not implemented "
                    "for colvar components of type \""+function_type+"\".\n");
}


void colvar::cvc::calc_Jacobian_derivative()
{
  cvm::fatal_error ("Error: calculation of inverse gradients is not implemented "
                    "for colvar components of type \""+function_type+"\".\n");
}
