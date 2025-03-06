// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include "colvarcomp.h"


colvar::linearCombination::linearCombination()
{
    set_function_type("linearCombination");
}


int colvar::linearCombination::init(std::string const &conf)
{
    int error_code = cvc::init(conf);
    if (error_code != COLVARS_OK) return error_code;

    // Lookup all available sub-cvcs
    for (auto it_cv_map = colvar::get_global_cvc_map().begin(); it_cv_map != colvar::get_global_cvc_map().end(); ++it_cv_map) {
        if (key_lookup(conf, it_cv_map->first.c_str())) {
            std::vector<std::string> sub_cvc_confs;
            get_key_string_multi_value(conf, it_cv_map->first.c_str(), sub_cvc_confs);
            for (auto it_sub_cvc_conf = sub_cvc_confs.begin(); it_sub_cvc_conf != sub_cvc_confs.end(); ++it_sub_cvc_conf) {
                cv.push_back((it_cv_map->second)());
                cv.back()->init(*(it_sub_cvc_conf));
            }
        }
    }
    // Sort all sub CVs by their names
    std::sort(cv.begin(), cv.end(), colvar::compare_cvc);
    for (auto it_sub_cv = cv.begin(); it_sub_cv != cv.end(); ++it_sub_cv) {
        for (auto it_atom_group = (*it_sub_cv)->atom_groups.begin(); it_atom_group != (*it_sub_cv)->atom_groups.end(); ++it_atom_group) {
            register_atom_group(*it_atom_group);
        }
    }
    // Show useful error messages and prevent crashes if no sub CVC is found
    if (cv.size() == 0) {
       return cvm::error("Error: the CV " + name + " expects one or more nesting components.\n",
                       COLVARS_INPUT_ERROR);
    } else {
        x.type(cv[0]->value());
        x.reset();
    }
    use_explicit_gradients = true;
    for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
        if (!cv[i_cv]->is_enabled(f_cvc_explicit_gradient)) {
            use_explicit_gradients = false;
        }
    }
    if (!use_explicit_gradients) {
        disable(f_cvc_explicit_gradient);
    }
    return error_code;
}

cvm::real colvar::linearCombination::getPolynomialFactorOfCVGradient(size_t i_cv) const {
    cvm::real factor_polynomial = 1.0;
    if (cv[i_cv]->value().type() == colvarvalue::type_scalar) {
        factor_polynomial = cv[i_cv]->sup_coeff * cv[i_cv]->sup_np * cvm::pow(cv[i_cv]->value().real_value, cv[i_cv]->sup_np - 1);
    } else {
        factor_polynomial = cv[i_cv]->sup_coeff;
    }
    return factor_polynomial;
}

colvar::linearCombination::~linearCombination() {
    // Recall the steps we initialize the sub-CVCs:
    // 1. Lookup all sub-CVCs and then register the atom groups for sub-CVCs
    //    in their constructors;
    // 2. Iterate over all sub-CVCs, get the pointers of their atom groups
    //    groups, and register again in the parent (current) CVC.
    // That being said, the atom groups become children of the sub-CVCs at
    // first, and then become children of the parent CVC.
    // So, to destruct this class (parent CVC class), we need to remove the
    // dependencies of the atom groups to the parent CVC at first.
    remove_all_children();
    // Then we remove the dependencies of the atom groups to the sub-CVCs
    // in their destructors.
    for (auto it = cv.begin(); it != cv.end(); ++it) {
        delete (*it);
    }
    // The last step is cleaning up the list of atom groups.
    atom_groups.clear();
}

void colvar::linearCombination::calc_value() {
    x.reset();
    for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
        cv[i_cv]->calc_value();
        colvarvalue current_cv_value(cv[i_cv]->value());
        // polynomial combination allowed
        if (current_cv_value.type() == colvarvalue::type_scalar) {
            x += cv[i_cv]->sup_coeff * (cvm::pow(current_cv_value.real_value, cv[i_cv]->sup_np));
        } else {
            x += cv[i_cv]->sup_coeff * current_cv_value;
        }
    }
}

void colvar::linearCombination::calc_gradients() {
    for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
        cv[i_cv]->calc_gradients();
        if (cv[i_cv]->is_enabled(f_cvc_explicit_gradient)) {
            cvm::real factor_polynomial = getPolynomialFactorOfCVGradient(i_cv);
            for (size_t j_elem = 0; j_elem < cv[i_cv]->value().size(); ++j_elem) {
                for (size_t k_ag = 0 ; k_ag < cv[i_cv]->atom_groups.size(); ++k_ag) {
                    for (size_t l_atom = 0; l_atom < (cv[i_cv]->atom_groups)[k_ag]->size(); ++l_atom) {
                        (*(cv[i_cv]->atom_groups)[k_ag])[l_atom].grad = factor_polynomial * (*(cv[i_cv]->atom_groups)[k_ag])[l_atom].grad;
                    }
                }
            }
        }
    }
}

void colvar::linearCombination::apply_force(colvarvalue const &force) {
    for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
        // If this CV us explicit gradients, then atomic gradients is already calculated
        // We can apply the force to atom groups directly
        if (cv[i_cv]->is_enabled(f_cvc_explicit_gradient)) {
            for (size_t k_ag = 0 ; k_ag < cv[i_cv]->atom_groups.size(); ++k_ag) {
                (cv[i_cv]->atom_groups)[k_ag]->apply_colvar_force(force.real_value);
            }
        } else {
            // Compute factors for polynomial combinations
            cvm::real factor_polynomial = getPolynomialFactorOfCVGradient(i_cv);
            colvarvalue cv_force = force.real_value * factor_polynomial;
            cv[i_cv]->apply_force(cv_force);
        }
    }
}


cvm::real colvar::linearCombination::dist2(colvarvalue const &x1, colvarvalue const &x2) const
{
  return x1.dist2(x2);
}


colvarvalue colvar::linearCombination::dist2_lgrad(colvarvalue const &x1,
                                                   colvarvalue const &x2) const
{
  return x1.dist2_grad(x2);
}


colvarvalue colvar::linearCombination::dist2_rgrad(colvarvalue const &x1,
                                                   colvarvalue const &x2) const
{
  return x2.dist2_grad(x1);
}


void colvar::linearCombination::wrap(colvarvalue & /* x_unwrapped */) const {}



colvar::customColvar::customColvar()
{
    set_function_type("customColvar");
}


int colvar::customColvar::init(std::string const &conf)
{
    int error_code = linearCombination::init(conf);
    if (error_code != COLVARS_OK) return error_code;

    // code swipe from colvar::init_custom_function
    std::string expr_in, expr;
    size_t pos = 0; // current position in config string
#ifdef LEPTON
    std::vector<Lepton::ParsedExpression> pexprs;
    Lepton::ParsedExpression pexpr;
    double *ref;
#endif
    if (key_lookup(conf, "customFunction", &expr_in, &pos)) {
#ifdef LEPTON
        use_custom_function = true;
        cvm::log("This colvar uses a custom function.\n");
        do {
            expr = expr_in;
            if (cvm::debug())
                cvm::log("Parsing expression \"" + expr + "\".\n");
            try {
                pexpr = Lepton::Parser::parse(expr);
                pexprs.push_back(pexpr);
            } catch (...) {
                return cvm::error("Error parsing expression \"" + expr + "\".\n", COLVARS_INPUT_ERROR);
            }
            try {
                value_evaluators.push_back(new Lepton::CompiledExpression(pexpr.createCompiledExpression()));
                // Define variables for cvc values
                for (size_t i = 0; i < cv.size(); ++i) {
                    for (size_t j = 0; j < cv[i]->value().size(); ++j) {
                        std::string vn = cv[i]->name + (cv[i]->value().size() > 1 ? cvm::to_str(j+1) : "");
                        try {
                            ref = &value_evaluators.back()->getVariableReference(vn);
                        } catch (...) {
                            ref = &dev_null;
                            cvm::log("Warning: Variable " + vn + " is absent from expression \"" + expr + "\".\n");
                        }
                        value_eval_var_refs.push_back(ref);
                    }
                }
            } catch (...) {
                return cvm::error("Error compiling expression \"" + expr + "\".\n", COLVARS_INPUT_ERROR);
            }
        } while (key_lookup(conf, "customFunction", &expr_in, &pos));
        // Now define derivative with respect to each scalar sub-component
        for (size_t i = 0; i < cv.size(); ++i) {
            for (size_t j = 0; j < cv[i]->value().size(); ++j) {
                std::string vn = cv[i]->name + (cv[i]->value().size() > 1 ? cvm::to_str(j+1) : "");
                for (size_t c = 0; c < pexprs.size(); ++c) {
                    gradient_evaluators.push_back(new Lepton::CompiledExpression(pexprs[c].differentiate(vn).createCompiledExpression()));
                    for (size_t k = 0; k < cv.size(); ++k) {
                        for (size_t l = 0; l < cv[k]->value().size(); l++) {
                            std::string vvn = cv[k]->name + (cv[k]->value().size() > 1 ? cvm::to_str(l+1) : "");
                            try {
                                ref = &gradient_evaluators.back()->getVariableReference(vvn);
                            } catch (...) {
                                cvm::log("Warning: Variable " + vvn + " is absent from derivative of \"" + expr + "\" wrt " + vn + ".\n");
                                ref = &dev_null;
                            }
                            grad_eval_var_refs.push_back(ref);
                        }
                    }
                }
            }
        }
        if (value_evaluators.size() == 0) {
            return cvm::error("Error: no custom function defined.\n", COLVARS_INPUT_ERROR);
        }
        if (value_evaluators.size() != 1) {
            x.type(colvarvalue::type_vector);
        } else {
            x.type(colvarvalue::type_scalar);
        }
#else
      return cvm::error(
          "customFunction requires the Lepton library, but it is not enabled during compilation.\n"
          "Please refer to the Compilation Notes section of the Colvars manual for more "
          "information.\n",
          COLVARS_NOT_IMPLEMENTED);
#endif
    } else {
        cvm::log("Warning: no customFunction specified.\n");
        cvm::log("Warning: use linear combination instead.\n");
    }
    return error_code;
}

colvar::customColvar::~customColvar() {
#ifdef LEPTON
    for (size_t i = 0; i < value_evaluators.size(); ++i) {
        if (value_evaluators[i] != nullptr) delete value_evaluators[i];
    }
    for (size_t i = 0; i < gradient_evaluators.size(); ++i) {
        if (gradient_evaluators[i] != nullptr) delete gradient_evaluators[i];
    }
#endif
}

void colvar::customColvar::calc_value() {
    if (!use_custom_function) {
        colvar::linearCombination::calc_value();
    } else {
#ifdef LEPTON
        for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
            cv[i_cv]->calc_value();
        }
        x.reset();
        size_t l = 0;
        for (size_t i = 0; i < x.size(); ++i) {
            for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
                const colvarvalue& current_cv_value = cv[i_cv]->value();
                for (size_t j_elem = 0; j_elem < current_cv_value.size(); ++j_elem) {
                    if (current_cv_value.type() == colvarvalue::type_scalar) {
                        *(value_eval_var_refs[l++]) = cv[i_cv]->sup_coeff * (cvm::pow(current_cv_value.real_value, cv[i_cv]->sup_np));
                    } else {
                        *(value_eval_var_refs[l++]) = cv[i_cv]->sup_coeff * current_cv_value[j_elem];
                    }
                }
            }
            x[i] = value_evaluators[i]->evaluate();
        }
#else
        cvm::error("customFunction requires the Lepton library, but it is not enabled during compilation.\n"
                   "Please refer to the Compilation Notes section of the Colvars manual for more information.\n",
                    COLVARS_INPUT_ERROR);
#endif
    }
}

void colvar::customColvar::calc_gradients() {
    if (!use_custom_function) {
        colvar::linearCombination::calc_gradients();
    } else {
#ifdef LEPTON
        size_t r = 0; // index in the vector of variable references
        size_t e = 0; // index of the gradient evaluator
        for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) { // for each CV
            cv[i_cv]->calc_gradients();
            if (cv[i_cv]->is_enabled(f_cvc_explicit_gradient)) {
                const colvarvalue& current_cv_value = cv[i_cv]->value();
                const cvm::real factor_polynomial = getPolynomialFactorOfCVGradient(i_cv);
                for (size_t j_elem = 0; j_elem < current_cv_value.size(); ++j_elem) { // for each element in this CV
                    for (size_t c = 0; c < x.size(); ++c) { // for each custom function expression
                        for (size_t k = 0; k < cv.size(); ++k) { // this is required since we need to feed all CV values to this expression
                            const cvm::real factor_polynomial_k = getPolynomialFactorOfCVGradient(k);
                            for (size_t l = 0; l < cv[k]->value().size(); ++l) {
                                *(grad_eval_var_refs[r++]) = factor_polynomial_k * cv[k]->value()[l];
                            }
                        }
                        const double expr_grad = gradient_evaluators[e++]->evaluate();
                        for (size_t k_ag = 0 ; k_ag < cv[i_cv]->atom_groups.size(); ++k_ag) {
                            for (size_t l_atom = 0; l_atom < (cv[i_cv]->atom_groups)[k_ag]->size(); ++l_atom) {
                                (*(cv[i_cv]->atom_groups)[k_ag])[l_atom].grad = expr_grad * factor_polynomial * (*(cv[i_cv]->atom_groups)[k_ag])[l_atom].grad;
                            }
                        }
                    }
                }
            }
        }
#else
        cvm::error("customFunction requires the Lepton library, but it is not enabled during compilation.\n"
                   "Please refer to the Compilation Notes section of the Colvars manual for more information.\n",
                    COLVARS_INPUT_ERROR);
#endif
    }
}

void colvar::customColvar::apply_force(colvarvalue const &force) {
    if (!use_custom_function) {
        colvar::linearCombination::apply_force(force);
    } else {
#ifdef LEPTON
        size_t r = 0; // index in the vector of variable references
        size_t e = 0; // index of the gradient evaluator
        for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
            // If this CV us explicit gradients, then atomic gradients is already calculated
            // We can apply the force to atom groups directly
            if (cv[i_cv]->is_enabled(f_cvc_explicit_gradient)) {
                for (size_t k_ag = 0 ; k_ag < cv[i_cv]->atom_groups.size(); ++k_ag) {
                    (cv[i_cv]->atom_groups)[k_ag]->apply_colvar_force(force.real_value);
                }
            } else {
                const colvarvalue& current_cv_value = cv[i_cv]->value();
                colvarvalue cv_force(current_cv_value);
                cv_force.reset();
                const cvm::real factor_polynomial = getPolynomialFactorOfCVGradient(i_cv);
                for (size_t j_elem = 0; j_elem < current_cv_value.size(); ++j_elem) {
                    for (size_t c = 0; c < x.size(); ++c) {
                        for (size_t k = 0; k < cv.size(); ++k) {
                            const cvm::real factor_polynomial_k = getPolynomialFactorOfCVGradient(k);
                            for (size_t l = 0; l < cv[k]->value().size(); ++l) {
                                *(grad_eval_var_refs[r++]) = factor_polynomial_k * cv[k]->value()[l];
                            }
                        }
                        cv_force[j_elem] += factor_polynomial * gradient_evaluators[e++]->evaluate() * force.real_value;
                    }
                }
                cv[i_cv]->apply_force(cv_force);
            }
        }
#else
        cvm::error("customFunction requires the Lepton library, but it is not enabled during compilation.\n"
                   "Please refer to the Compilation Notes section of the Colvars manual for more information.\n",
                    COLVARS_INPUT_ERROR);
#endif
    }
}
